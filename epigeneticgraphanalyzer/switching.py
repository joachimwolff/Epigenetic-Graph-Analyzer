"""egaSwitching: enhancer switching and promoter competition between two conditions.

Reads the union graph of egaDifferential (<prefix>_union.graphml: edge attributes status gained / lost / shared, r_a,
r_b, log2FC; node attributes bait, genes, rna_log2FC and one signal per chromatin mark and condition) and, when
present, the element classes of the fragments (<prefix>_elements.tsv).

Contact changes. A contact E-P is lost in B (or weakened) if its status is 'lost' or its log2FC (B over A) is at most
-`--contactThreshold`; it is gained (or strengthened) if its status is 'gained' or its log2FC is at least
+`--contactThreshold`; otherwise it is kept.

Switch event. A fragment E that loses a contact with the baited promoter A and gains a contact with the baited
promoter B (A and B distinct, no gene in common). E is a non-promoter element (main set) or any fragment including
promoters (second set). Every event is reported with the position, element class and chromatin changes of E, the
genes, expression changes, distances and contact strengths of A and B. An event is concordant if the expression of
A goes down and that of B goes up (signs, and in a strict variant both by at least `--rnaThreshold`).

Hub switch (one step). A non-promoter element E leaves the hubs (promoter and its direct partners) of the promoters
it reaches within two steps in A but not in B, and joins those it reaches within two steps in B but not in A. One
row per element with both sets non-empty; concordant if the mean expression change of the left promoters is negative
and that of the joined promoters positive.

Statistics. The concordant fraction is compared with two null models of the expression changes, each with
`--permutations` draws: (i) a permutation of the expression changes among the promoters of the events within strata
of expression level (quintiles) and contact distance class (median distance to their elements: < 100 kb, 100-500 kb,
>= 500 kb); (ii) a rotation of the expression changes of all promoters with expression along each chromosome (gene
order kept, local correlation of expression preserved); (iii) a swap of the roles of the two promoters (lost /
gained), at random and jointly for all events of an element, which keeps the expression of both promoters and their
co-expression and tests the excess of concordant over discordant events. Promoters contacting the same element are
often co-expressed (nearby genes), which lowers both the concordant and the discordant fraction below the
expectation of (i); (ii) and (iii) keep that structure. One-sided permutation P values, (1 + #null >= observed) /
(1 + permutations). Controls: the concordant fraction of ordered promoter pairs that keep their contact with the same
element (no switch; by symmetry half of the pairs with opposite expression changes), and the fraction of lost
promoters going down (gained promoters going up) in switch events against all lost (gained) promoter contacts of
elements without a switch.

Competition. Pairs of baited promoters that share at least one non-promoter element. For each pair the net change of
the shared contacts (gained minus lost) of each promoter; pairs with a different net change are oriented so that
promoter 1 gains more of the shared elements. Competitive if promoter 1 goes up and promoter 2 down; tested against
the same three null models (the role swap per pair).

Ranking. Events are ordered by an evidence score, the sum of the contact change ((min(4, -log2FC lost) + min(4,
log2FC gained)) / 2), the expression change (geometric mean of min(4, -log2FC of A) and min(4, log2FC of B), 0 unless
A goes down and B up) and the chromatin change at E (0.5 x the largest absolute log2 change of a mark, at most 3).

The software reports and indicates; whether an event reflects regulation is for the scientist to decide.
"""

import json
import logging
import os
import re
import time
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from scipy import sparse, stats

from ._cli import base_parser, ensure_parent, setup_logging, write_table

log = logging.getLogger(__name__)

DISTANCE_BREAKS = (100_000, 500_000)
DISTANCE_LABELS = ('<100 kb', '100-500 kb', '>=500 kb')
GRAPHML = '{http://graphml.graphdrawing.org/xmlns}'


# ----------------------------------------------------------------------------------------------------------- input

def _convert(value, kind):
    if value is None:
        return None
    if kind in ('double', 'float'):
        try:
            return float(value)
        except ValueError:
            return np.nan
    if kind in ('long', 'int'):
        try:
            return int(value)
        except ValueError:
            return int(float(value))
    if kind == 'boolean':
        return value.lower() in ('true', '1')
    return value


def read_union_graphml(path):
    """Nodes and edges of a union graph as two DataFrames (streamed, without networkx).

    nodes: index node id, columns chrom, start, end, bait, genes, rna_log2FC and the mark signals; edges: u, v and
    the edge attributes (r_a, r_b, log2FC, status, ...)."""
    keys, node_rows, edge_rows = {}, [], []
    for _, element in ET.iterparse(path, events=('end',)):
        tag = element.tag.replace(GRAPHML, '')
        if tag == 'key':
            keys[element.get('id')] = (element.get('attr.name'), element.get('attr.type'))
        elif tag == 'node':
            row = {'node': element.get('id')}
            for data in element.findall(GRAPHML + 'data'):
                name, kind = keys[data.get('key')]
                row[name] = _convert(data.text if data.text is not None else '', kind) if kind != 'string' \
                    else (data.text or '')
            node_rows.append(row)
            element.clear()
        elif tag == 'edge':
            row = {'u': element.get('source'), 'v': element.get('target')}
            for data in element.findall(GRAPHML + 'data'):
                name, kind = keys[data.get('key')]
                row[name] = _convert(data.text, kind) if kind != 'string' else (data.text or '')
            edge_rows.append(row)
            element.clear()
    nodes = pd.DataFrame(node_rows).set_index('node')
    if 'genes' not in nodes:
        nodes['genes'] = ''
    nodes['genes'] = nodes['genes'].fillna('').astype(str)
    nodes['bait'] = nodes.get('bait', pd.Series(0, index=nodes.index)).fillna(0).astype(int).astype(bool)
    if 'rna_log2FC' not in nodes:
        nodes['rna_log2FC'] = np.nan
    nodes['rna_log2FC'] = pd.to_numeric(nodes['rna_log2FC'], errors='coerce')
    if 'chrom' not in nodes or nodes['chrom'].isna().any():
        parsed = nodes.index.to_series().str.extract(r'^(.+):(\d+)-(\d+)$')
        nodes['chrom'] = nodes.get('chrom', parsed[0]).fillna(parsed[0])
        nodes['start'] = nodes.get('start', parsed[1]).fillna(parsed[1])
        nodes['end'] = nodes.get('end', parsed[2]).fillna(parsed[2])
    nodes['start'] = nodes['start'].astype(np.int64)
    nodes['end'] = nodes['end'].astype(np.int64)
    edges = pd.DataFrame(edge_rows)
    for column in ('r_a', 'r_b', 'log2FC'):
        edges[column] = pd.to_numeric(edges.get(column, np.nan), errors='coerce')
    edges['status'] = edges.get('status', 'shared').fillna('shared')
    return nodes, edges


def condition_names(prefix, nodes=None):
    """Condition names from the prefix (<A>_vs_<B>), else from the mark columns of the nodes."""
    base = os.path.basename(prefix)
    if '_vs_' in base:
        a, b = base.split('_vs_', 1)
        return a, b
    names = []
    for column in nodes.columns if nodes is not None else []:
        if '_' in column and column not in ('rna_log2FC', 'gene_ids'):
            name = column.rsplit('_', 1)[1]
            if name not in names:
                names.append(name)
    return tuple(names[:2]) if len(names) >= 2 else ('A', 'B')


def mark_names(nodes, names):
    a, b = names
    return sorted({c[:-len(a) - 1] for c in nodes.columns if c.endswith('_' + a)}
                  & {c[:-len(b) - 1] for c in nodes.columns if c.endswith('_' + b)})


def read_elements(path, names):
    """Element class of every fragment in A and B from <prefix>_elements.tsv (declared roles, fragment mean)."""
    if not path or not os.path.exists(path):
        return None
    table = pd.read_csv(path, sep='\t')
    column = next((c for c in ('class_declared_mean', 'class_declared_max1kb') if c in table), None)
    if column is None:
        column = next((c for c in table.columns if c.startswith('class_')), None)
    if column is None or 'condition' not in table:
        return None
    out = pd.DataFrame(index=pd.Index(table['node'].unique(), name='node'))
    for name in names:
        part = table[table['condition'] == name].set_index('node')
        out['class_' + name] = part[column]
        if 'element_type' in part and 'element_type' not in out:
            out['element_type'] = part['element_type']
    return out


def expression_levels(nodes, path, gene_column='gene_name', level_columns=None):
    """Expression level of every baited promoter: mean over its genes of log2(baseMean + 1), or of the mean of
    `level_columns` (already log scale) when given or when the table has columns mean_log2_*."""
    if not path:
        return pd.Series(np.nan, index=nodes.index)
    table = pd.read_csv(path, sep='\t')
    if gene_column not in table:
        gene_column = next(c for c in ('gene_name', 'gene', 'genes', 'symbol') if c in table)
    if level_columns:
        level = table[list(level_columns)].astype(float).mean(axis=1)
    elif 'baseMean' in table:
        level = np.log2(table['baseMean'].astype(float) + 1)
    else:
        columns = [c for c in table.columns if c.startswith('mean_log2')]
        level = table[columns].astype(float).mean(axis=1)
    by_gene = pd.Series(level.to_numpy(), index=table[gene_column].astype(str)).groupby(level=0).mean()
    values = {}
    for node, genes in nodes['genes'].items():
        found = [by_gene[g] for g in genes.split(',') if g and g in by_gene.index]
        values[node] = float(np.mean(found)) if found else np.nan
    return pd.Series(values)


# ----------------------------------------------------------------------------------------------- contact changes

def classify_edges(edges, threshold=1.0):
    """'loss', 'gain' or 'kept' for every edge (module doc)."""
    log2 = edges['log2FC'].to_numpy(dtype=float)
    status = edges['status'].to_numpy()
    loss = (status == 'lost') | (log2 <= -threshold)
    gain = (status == 'gained') | (log2 >= threshold)
    change = np.where(loss & ~gain, 'loss', np.where(gain & ~loss, 'gain', 'kept'))
    return pd.Series(change, index=edges.index)


def promoter_mask(nodes, elements=None):
    """Baited promoters (bait with at least one gene)."""
    return nodes['bait'] & nodes['genes'].str.len().gt(0)


def element_is_promoter(nodes, elements=None):
    """Fragments counted as promoters for E: baited promoters and, with an element table, every fragment of type
    'promoter' (TSS within the fragment)."""
    mask = promoter_mask(nodes)
    if elements is not None and 'element_type' in elements:
        typed = elements['element_type'].reindex(nodes.index).eq('promoter')
        mask = mask | typed.fillna(False).to_numpy()
    return mask


def directed_contacts(nodes, edges, threshold=1.0):
    """One row per (fragment E, baited promoter P) contact of the union graph with its change, strengths and
    distance."""
    change = classify_edges(edges, threshold)
    is_promoter = promoter_mask(nodes)
    mid = ((nodes['start'] + nodes['end']) // 2)
    frames = []
    for e_col, p_col in (('u', 'v'), ('v', 'u')):
        keep = is_promoter.reindex(edges[p_col]).to_numpy()
        part = pd.DataFrame({'element': edges[e_col].to_numpy()[keep], 'promoter': edges[p_col].to_numpy()[keep],
                             'status': edges['status'].to_numpy()[keep], 'change': change.to_numpy()[keep],
                             'r_a': edges['r_a'].to_numpy()[keep], 'r_b': edges['r_b'].to_numpy()[keep],
                             'log2FC': edges['log2FC'].to_numpy()[keep]})
        frames.append(part)
    contacts = pd.concat(frames, ignore_index=True)
    contacts = contacts[contacts['element'] != contacts['promoter']]
    same_chrom = nodes['chrom'].reindex(contacts['element']).to_numpy() == \
        nodes['chrom'].reindex(contacts['promoter']).to_numpy()
    contacts = contacts[same_chrom].copy()
    contacts['distance'] = np.abs(mid.reindex(contacts['element']).to_numpy() -
                                  mid.reindex(contacts['promoter']).to_numpy())
    return contacts.reset_index(drop=True)


def _gene_sets(nodes):
    return {n: frozenset(g for g in s.split(',') if g) for n, s in nodes['genes'].items()}


def mark_changes(nodes, marks, names, pseudo=1.0):
    """log2((B + pseudo) / (A + pseudo)) of every mark on every fragment."""
    out = pd.DataFrame(index=nodes.index)
    for mark in marks:
        a = pd.to_numeric(nodes['{}_{}'.format(mark, names[0])], errors='coerce')
        b = pd.to_numeric(nodes['{}_{}'.format(mark, names[1])], errors='coerce')
        out[mark] = np.log2((b + pseudo) / (a + pseudo))
    return out


def switch_events(nodes, edges, names, marks=None, elements=None, contact_threshold=1.0, rna_threshold=1.0,
                  chromatin_threshold=1.0, pseudo=1.0):
    """All switch events (module doc). One row per (E, lost promoter, gained promoter)."""
    marks = mark_names(nodes, names) if marks is None else marks
    contacts = directed_contacts(nodes, edges, contact_threshold)
    lost = contacts[contacts['change'] == 'loss']
    gained = contacts[contacts['change'] == 'gain']
    events = lost.merge(gained, on='element', suffixes=('_lost', '_gained'))
    events = events[events['promoter_lost'] != events['promoter_gained']]
    genes = _gene_sets(nodes)
    same = np.array([bool(genes[a] & genes[b]) for a, b in zip(events['promoter_lost'], events['promoter_gained'])],
                    dtype=bool)
    events = events[~same].copy()
    rename = {'status_lost': 'status_lost', 'change_lost': None, 'change_gained': None}
    events = events.drop(columns=[k for k, v in rename.items() if v is None])
    na, nb = names
    events = events.rename(columns={'r_a_lost': 'r_{}_lost'.format(na), 'r_b_lost': 'r_{}_lost'.format(nb),
                                    'r_a_gained': 'r_{}_gained'.format(na), 'r_b_gained': 'r_{}_gained'.format(nb)})
    e = events['element']
    is_prom = element_is_promoter(nodes, elements)
    events.insert(1, 'element_chrom', nodes['chrom'].reindex(e).to_numpy())
    events.insert(2, 'element_start', nodes['start'].reindex(e).to_numpy())
    events.insert(3, 'element_end', nodes['end'].reindex(e).to_numpy())
    events.insert(4, 'element_type', np.where(is_prom.reindex(e).to_numpy(), 'promoter', 'non-promoter element'))
    events.insert(5, 'element_genes', nodes['genes'].reindex(e).to_numpy())
    if elements is not None:
        for name in names:
            if 'class_' + name in elements:
                events['element_class_' + name] = elements['class_' + name].reindex(e).to_numpy()
    changes = mark_changes(nodes, marks, names, pseudo)
    changed_text, changed_any, max_abs = [], [], []
    sub = changes.reindex(e)
    for mark in marks:
        events['element_{}_{}'.format(mark, na)] = pd.to_numeric(nodes['{}_{}'.format(mark, na)],
                                                                errors='coerce').reindex(e).to_numpy()
        events['element_{}_{}'.format(mark, nb)] = pd.to_numeric(nodes['{}_{}'.format(mark, nb)],
                                                                errors='coerce').reindex(e).to_numpy()
        events['element_{}_log2'.format(mark)] = sub[mark].to_numpy()
    values = sub.to_numpy(dtype=float) if marks else np.zeros((len(events), 0))
    for row in values:
        items = ['{} {:+.1f}'.format(m, v) for m, v in zip(marks, row) if np.isfinite(v) and
                 abs(v) >= chromatin_threshold]
        changed_text.append(', '.join(items))
        changed_any.append(bool(items))
        finite = row[np.isfinite(row)]
        max_abs.append(float(np.abs(finite).max()) if len(finite) else 0.0)
    events['element_chromatin_changed'] = changed_any
    events['element_chromatin_changes'] = changed_text
    if elements is not None and all('class_' + n in elements for n in names):
        events['element_class_changed'] = (events['element_class_' + na].astype(str) !=
                                           events['element_class_' + nb].astype(str)).to_numpy()
    for side in ('lost', 'gained'):
        p = events['promoter_' + side]
        events['genes_' + side] = nodes['genes'].reindex(p).to_numpy()
        events['rna_log2FC_' + side] = nodes['rna_log2FC'].reindex(p).to_numpy()
    rl, rg = events['rna_log2FC_lost'].to_numpy(float), events['rna_log2FC_gained'].to_numpy(float)
    has = np.isfinite(rl) & np.isfinite(rg)
    events['expression_available'] = has
    events['concordant'] = has & (rl < 0) & (rg > 0)
    events['discordant'] = has & (rl > 0) & (rg < 0)
    events['concordant_strict'] = has & (rl <= -rna_threshold) & (rg >= rna_threshold)
    events['strict_contacts'] = (events['log2FC_lost'] <= -contact_threshold) & \
        (events['log2FC_gained'] >= contact_threshold)
    contact = (np.clip(-events['log2FC_lost'].to_numpy(float), 0, 4) +
               np.clip(events['log2FC_gained'].to_numpy(float), 0, 4)) / 2
    # expression evidence: geometric mean of the expected changes (A down, B up), 0 unless both move as expected
    expression = np.sqrt(np.clip(np.nan_to_num(-rl), 0, 4) * np.clip(np.nan_to_num(rg), 0, 4))
    chromatin = 0.5 * np.clip(np.array(max_abs), 0, 3)
    events['evidence_contact'] = contact
    events['evidence_expression'] = expression
    events['evidence_chromatin'] = chromatin
    events['evidence'] = contact + expression + chromatin
    events = events.sort_values(['evidence', 'element'], ascending=[False, True]).reset_index(drop=True)
    events['rank'] = np.arange(1, len(events) + 1)
    order = ['rank', 'element', 'element_chrom', 'element_start', 'element_end', 'element_type', 'element_genes']
    order += [c for c in events.columns if c.startswith('element_class')]
    order += [c for c in events.columns if c.startswith('element_') and c not in order]
    for side in ('lost', 'gained'):
        order += ['promoter_' + side, 'genes_' + side, 'rna_log2FC_' + side, 'distance_' + side, 'status_' + side,
                  'r_{}_{}'.format(na, side), 'r_{}_{}'.format(nb, side), 'log2FC_' + side]
    order += [c for c in events.columns if c not in order]
    return events[order]


# --------------------------------------------------------------------------------------------------- hub switches

def _adjacency(nodes, edges, present):
    index = pd.Series(np.arange(len(nodes)), index=nodes.index)
    u = index.reindex(edges['u'][present]).to_numpy()
    v = index.reindex(edges['v'][present]).to_numpy()
    n = len(nodes)
    m = sparse.coo_matrix((np.ones(len(u) * 2), (np.r_[u, v], np.r_[v, u])), shape=(n, n)).tocsr()
    m.data[:] = 1
    return m


def hub_switches(nodes, edges, names, elements=None):
    """Elements that leave the hubs of some promoters and join the hubs of others (module doc). Returns the table
    and the sparse indicator matrices (elements x promoters) of the left and joined promoters."""
    status = edges['status']
    in_a = _adjacency(nodes, edges, status.isin(['lost', 'shared']).to_numpy())
    in_b = _adjacency(nodes, edges, status.isin(['gained', 'shared']).to_numpy())
    promoter = promoter_mask(nodes).to_numpy()
    element = ~element_is_promoter(nodes, elements).to_numpy()
    p_index = np.flatnonzero(promoter)
    e_index = np.flatnonzero(element)

    def reach(m):
        rows = m[e_index]
        to_promoters = m[:, p_index]
        two = rows[:, p_index] + rows @ to_promoters
        two.data[:] = 1
        two.eliminate_zeros()
        return two.tocsr()

    reach_a, reach_b = reach(in_a), reach(in_b)
    left = (reach_a - reach_a.multiply(reach_b)).tocsr()
    joined = (reach_b - reach_a.multiply(reach_b)).tocsr()
    left.eliminate_zeros()
    joined.eliminate_zeros()
    n_left = np.asarray(left.sum(axis=1)).ravel()
    n_joined = np.asarray(joined.sum(axis=1)).ravel()
    keep = np.flatnonzero((n_left > 0) & (n_joined > 0))
    left, joined = left[keep], joined[keep]
    e_nodes = nodes.index[e_index[keep]]
    p_nodes = nodes.index[p_index]
    rna = nodes['rna_log2FC'].to_numpy(float)[p_index]
    rows = []
    for k, e in enumerate(e_nodes):
        lp = left.indices[left.indptr[k]:left.indptr[k + 1]]
        jp = joined.indices[joined.indptr[k]:joined.indptr[k + 1]]
        lr, jr = rna[lp], rna[jp]
        lm = float(np.nanmean(lr)) if np.isfinite(lr).any() else np.nan
        jm = float(np.nanmean(jr)) if np.isfinite(jr).any() else np.nan
        rows.append({'element': e, 'element_chrom': nodes.at[e, 'chrom'], 'element_start': nodes.at[e, 'start'],
                     'element_end': nodes.at[e, 'end'], 'n_left': len(lp), 'n_joined': len(jp),
                     'genes_left': _join_genes(nodes, p_nodes[lp]), 'genes_joined': _join_genes(nodes, p_nodes[jp]),
                     'mean_rna_log2FC_left': lm, 'mean_rna_log2FC_joined': jm,
                     'expression_available': bool(np.isfinite(lm) and np.isfinite(jm)),
                     'concordant': bool(np.isfinite(lm) and np.isfinite(jm) and lm < 0 and jm > 0),
                     'discordant': bool(np.isfinite(lm) and np.isfinite(jm) and lm > 0 and jm < 0)})
    table = pd.DataFrame(rows, columns=['element', 'element_chrom', 'element_start', 'element_end', 'n_left',
                                        'n_joined', 'genes_left', 'genes_joined', 'mean_rna_log2FC_left',
                                        'mean_rna_log2FC_joined', 'expression_available', 'concordant',
                                        'discordant'])
    return table, left, joined, p_nodes


def _join_genes(nodes, promoters, limit=30):
    genes = []
    for p in promoters:
        for g in nodes.at[p, 'genes'].split(','):
            if g and g not in genes:
                genes.append(g)
    return ','.join(genes[:limit]) + (',...' if len(genes) > limit else '')


# ----------------------------------------------------------------------------------------------------- competition

def competition_pairs(nodes, edges, contact_threshold=1.0, elements=None, max_shared_list=20):
    """Pairs of baited promoters sharing at least one non-promoter element, with the net change of their shared
    contacts; oriented so that promoter_1 gains more (pairs without a difference are kept, delta 0)."""
    contacts = directed_contacts(nodes, edges, contact_threshold)
    element = ~element_is_promoter(nodes, elements)
    contacts = contacts[element.reindex(contacts['element']).to_numpy()]
    contacts = contacts.assign(c=contacts['change'].map({'gain': 1, 'loss': -1, 'kept': 0}).astype(int))
    part = contacts[['element', 'promoter', 'c']]
    pairs = part.merge(part, on='element', suffixes=('_1', '_2'))
    pairs = pairs[pairs['promoter_1'] < pairs['promoter_2']]
    genes = _gene_sets(nodes)
    same = np.array([bool(genes[a] & genes[b]) for a, b in zip(pairs['promoter_1'], pairs['promoter_2'])],
                    dtype=bool)
    pairs = pairs[~same]
    pairs = pairs.assign(g1=(pairs['c_1'] > 0).astype(int), l1=(pairs['c_1'] < 0).astype(int),
                         g2=(pairs['c_2'] > 0).astype(int), l2=(pairs['c_2'] < 0).astype(int))
    grouped = pairs.groupby(['promoter_1', 'promoter_2'], sort=False)
    table = grouped.agg(n_shared=('element', 'size'), gained_1=('g1', 'sum'), lost_1=('l1', 'sum'),
                        gained_2=('g2', 'sum'), lost_2=('l2', 'sum'),
                        shared_elements=('element', lambda s: ';'.join(list(s)[:max_shared_list]) +
                                         (';...' if len(s) > max_shared_list else ''))).reset_index()
    table['net_1'] = table['gained_1'] - table['lost_1']
    table['net_2'] = table['gained_2'] - table['lost_2']
    flip = table['net_2'] > table['net_1']
    for x, y in (('promoter_1', 'promoter_2'), ('gained_1', 'gained_2'), ('lost_1', 'lost_2'), ('net_1', 'net_2')):
        a, b = table.loc[flip, x].copy(), table.loc[flip, y].copy()
        table.loc[flip, x], table.loc[flip, y] = b.to_numpy(), a.to_numpy()
    table['delta'] = table['net_1'] - table['net_2']
    mid = (nodes['start'] + nodes['end']) // 2
    table['promoter_distance'] = np.abs(mid.reindex(table['promoter_1']).to_numpy() -
                                        mid.reindex(table['promoter_2']).to_numpy())
    for k in ('1', '2'):
        p = table['promoter_' + k]
        table['genes_' + k] = nodes['genes'].reindex(p).to_numpy()
        table['rna_log2FC_' + k] = nodes['rna_log2FC'].reindex(p).to_numpy()
    r1, r2 = table['rna_log2FC_1'].to_numpy(float), table['rna_log2FC_2'].to_numpy(float)
    has = np.isfinite(r1) & np.isfinite(r2) & (table['delta'].to_numpy() > 0)
    table['tested'] = has
    table['competitive'] = has & (r1 > 0) & (r2 < 0)
    table['reverse'] = has & (r1 < 0) & (r2 > 0)
    order = ['promoter_1', 'genes_1', 'rna_log2FC_1', 'promoter_2', 'genes_2', 'rna_log2FC_2', 'promoter_distance',
             'n_shared', 'gained_1', 'lost_1', 'net_1', 'gained_2', 'lost_2', 'net_2', 'delta', 'tested',
             'competitive', 'reverse', 'shared_elements']
    table = table[order]
    table = table.sort_values(['delta', 'n_shared'], ascending=False).reset_index(drop=True)
    return table, contacts


# ------------------------------------------------------------------------------------------------------ null models

def role_swap_test(concordant, discordant, clusters, n_perm=10000, seed=0, chunk=500):
    """Asymmetry of concordant over discordant events: the roles of the two promoters (lost / gained, or winner /
    loser) are swapped at random, jointly for all events of a cluster (an element, or a pair), which turns concordant
    into discordant events and back. The null keeps the expression of both promoters and their co-expression.
    Returns (null mean of the concordant fraction, one-sided P, null distribution)."""
    concordant = np.asarray(concordant, dtype=bool)
    discordant = np.asarray(discordant, dtype=bool)
    if not len(concordant):
        return np.nan, np.nan, np.array([])
    codes = pd.factorize(pd.Series(clusters))[0]
    k = codes.max() + 1
    rng = np.random.default_rng(seed)
    c = np.bincount(codes, weights=concordant, minlength=k)
    d = np.bincount(codes, weights=discordant, minlength=k)
    null = []
    done = 0
    while done < n_perm:
        n = min(chunk, n_perm - done)
        flip = rng.random((n, k)) < 0.5
        null.append((np.where(flip, d[None, :], c[None, :]).sum(axis=1)) / len(concordant))
        done += n
    null = np.concatenate(null)
    return float(null.mean()), permutation_p(float(concordant.mean()), null), null


def distance_class(distances):
    return np.digitize(np.asarray(distances, dtype=float), DISTANCE_BREAKS)


def strata_of(pool, levels, distances, n_level_bins=5):
    """Stratum of every promoter of `pool`: expression level quintile (within the pool) x distance class; promoters
    without a level form their own level bin."""
    lv = levels.reindex(pool).to_numpy(float)
    finite = np.isfinite(lv)
    bins = np.full(len(pool), n_level_bins, dtype=int)
    if finite.sum() >= n_level_bins:
        edges = np.unique(np.quantile(lv[finite], np.linspace(0, 1, n_level_bins + 1)[1:-1]))
        bins[finite] = np.digitize(lv[finite], edges)
    dc = distance_class(distances.reindex(pool).to_numpy(float))
    return bins * 10 + dc


def stratified_draws(values, strata, n, rng):
    """`n` permutations of `values` within strata (n x len(values))."""
    values = np.asarray(values, dtype=float)
    order = np.argsort(strata, kind='stable')
    sorted_strata = strata[order]
    keys = rng.random((n, len(values))) + sorted_strata[None, :] * 2.0
    perm = np.argsort(keys, axis=1)
    out = np.empty((n, len(values)))
    out[:, order] = values[order][perm]
    return out


def rotation_draws(values, chroms, positions, n, rng):
    """`n` rotations of `values` along each chromosome (n x len(values))."""
    values = np.asarray(values, dtype=float)
    out = np.empty((n, len(values)))
    for chrom in np.unique(chroms):
        idx = np.flatnonzero(chroms == chrom)
        idx = idx[np.argsort(positions[idx], kind='stable')]
        m = len(idx)
        shifts = rng.integers(0, m, size=n) if m > 1 else np.zeros(n, dtype=int)
        source = (np.arange(m)[None, :] - shifts[:, None]) % m
        out[:, idx] = values[idx][source]
    return out


def permutation_p(observed, null):
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if not np.isfinite(observed) or len(null) == 0:
        return np.nan
    return float((1 + np.sum(null >= observed - 1e-12)) / (1 + len(null)))


def pair_statistics(x_first, x_second, threshold=1.0):
    """Concordant fraction (first < 0, second > 0), strict concordant fraction and mean difference second - first,
    along the last axis."""
    return {'concordant': np.mean((x_first < 0) & (x_second > 0), axis=-1),
            'concordant_strict': np.mean((x_first <= -threshold) & (x_second >= threshold), axis=-1),
            'mean_difference': np.mean(x_second - x_first, axis=-1)}


def null_test(first, second, pool, levels, distances, all_promoters, n_perm=10000, seed=0, threshold=1.0,
              chunk=200):
    """Observed pair statistics of the promoter pairs (first[i], second[i]) and their null distributions under the
    stratified permutation and the chromosome rotation. Returns (observed, null_permutation, null_rotation)."""
    rng = np.random.default_rng(seed)
    rna_all = all_promoters['rna_log2FC']
    pool = pd.Index(pool)
    values = rna_all.reindex(pool).to_numpy(float)
    strata = strata_of(pool, levels, distances)
    pos = pd.Series(np.arange(len(pool)), index=pool)
    i1, i2 = pos.reindex(first).to_numpy(), pos.reindex(second).to_numpy()
    observed = {k: float(v) for k, v in pair_statistics(values[i1], values[i2], threshold).items()}
    perm = {k: [] for k in observed}
    rot = {k: [] for k in observed}
    expressed = all_promoters[np.isfinite(all_promoters['rna_log2FC'].to_numpy(float))]
    rpos = pd.Series(np.arange(len(expressed)), index=expressed.index)
    j1, j2 = rpos.reindex(first).to_numpy(), rpos.reindex(second).to_numpy()
    chroms = expressed['chrom'].to_numpy()
    mids = ((expressed['start'] + expressed['end']) // 2).to_numpy()
    rvalues = expressed['rna_log2FC'].to_numpy(float)
    done = 0
    while done < n_perm:
        n = min(chunk, n_perm - done)
        draws = stratified_draws(values, strata, n, rng)
        for k, v in pair_statistics(draws[:, i1], draws[:, i2], threshold).items():
            perm[k].append(v)
        draws = rotation_draws(rvalues, chroms, mids, n, rng)
        for k, v in pair_statistics(draws[:, j1], draws[:, j2], threshold).items():
            rot[k].append(v)
        done += n
    return observed, {k: np.concatenate(v) for k, v in perm.items()}, {k: np.concatenate(v) for k, v in rot.items()}


def hub_null_test(table, left, joined, p_nodes, nodes, levels, n_perm=10000, seed=0, chunk=100):
    """Concordant fraction of hub switches (mean change of the left promoters < 0 and of the joined > 0) against
    the rotation null (all promoters with expression) and a stratified permutation among the promoters of the hub
    switches (level quintile; the distance class is the two-step reach, not used)."""
    rng = np.random.default_rng(seed)
    rna = nodes['rna_log2FC'].reindex(p_nodes).to_numpy(float)
    finite = np.isfinite(rna)
    L = left[:, finite].tocsr()
    J = joined[:, finite].tocsr()
    keep = (np.asarray(L.sum(axis=1)).ravel() > 0) & (np.asarray(J.sum(axis=1)).ravel() > 0)
    L, J = L[keep], J[keep]
    nl = np.asarray(L.sum(axis=1)).ravel()
    nj = np.asarray(J.sum(axis=1)).ravel()
    values = rna[finite]

    def fraction(X):
        ml = (L @ X) / nl[:, None]
        mj = (J @ X) / nj[:, None]
        return np.mean((ml < 0) & (mj > 0), axis=0)

    observed = float(fraction(values[:, None])[0]) if keep.any() else np.nan
    used = np.flatnonzero((np.asarray(L.sum(axis=0)).ravel() + np.asarray(J.sum(axis=0)).ravel()) > 0)
    sub_nodes = p_nodes[finite]
    strata = strata_of(sub_nodes[used], levels, pd.Series(0.0, index=sub_nodes[used]))
    chroms = nodes['chrom'].reindex(sub_nodes).to_numpy()
    mids = ((nodes['start'] + nodes['end']) // 2).reindex(sub_nodes).to_numpy()
    perm, rot = [], []
    done = 0
    while done < n_perm and keep.any():
        n = min(chunk, n_perm - done)
        X = np.repeat(values[:, None], n, axis=1)
        X[used, :] = stratified_draws(values[used], strata, n, rng).T
        perm.append(fraction(X))
        rot.append(fraction(rotation_draws(values, chroms, mids, n, rng).T))
        done += n
    perm = np.concatenate(perm) if perm else np.array([])
    rot = np.concatenate(rot) if rot else np.array([])
    return int(keep.sum()), observed, perm, rot


# ---------------------------------------------------------------------------------------------------- controls

def stable_pair_control(contacts, rna):
    """Ordered pairs of promoters that keep their contact with the same element: number of ordered pairs with
    expression and the number that are concordant (first down, second up) = unordered pairs of opposite sign."""
    kept = contacts[contacts['change'] == 'kept'].copy()
    kept['rna'] = rna.reindex(kept['promoter']).to_numpy()
    kept = kept[np.isfinite(kept['rna'].to_numpy(float))]
    grouped = kept.groupby('element')['rna']
    down = grouped.apply(lambda s: int((s < 0).sum()))
    up = grouped.apply(lambda s: int((s > 0).sum()))
    k = grouped.size()
    ordered = int((k * (k - 1)).sum())
    concordant = int((down * up).sum())
    return ordered, concordant


def one_sided_control(events, contacts, rna):
    """Fraction of lost promoters going down (gained promoters going up) in switch events against lost (gained)
    contacts of elements without any switch event."""
    switch_elements = set(events['element'])
    rows = {}
    for change, side, sign in (('loss', 'lost', -1), ('gain', 'gained', 1)):
        ev = events[np.isfinite(events['rna_log2FC_' + side].to_numpy(float))]
        ev = ev.drop_duplicates(['element', 'promoter_' + side])
        in_switch = np.sign(ev['rna_log2FC_' + side].to_numpy(float)) == sign
        other = contacts[(contacts['change'] == change) & ~contacts['element'].isin(switch_elements)]
        r = rna.reindex(other['promoter']).to_numpy(float)
        r = r[np.isfinite(r)]
        outside = np.sign(r) == sign
        a, b = int(in_switch.sum()), int(len(in_switch) - in_switch.sum())
        c, d = int(outside.sum()), int(len(outside) - outside.sum())
        p = stats.fisher_exact([[a, b], [c, d]], alternative='greater')[1] if min(a + b, c + d) > 0 else np.nan
        rows[side] = {'switch_n': a + b, 'switch_fraction': a / (a + b) if a + b else np.nan,
                      'other_n': c + d, 'other_fraction': c / (c + d) if c + d else np.nan, 'fisher_p': p}
    return rows


# ------------------------------------------------------------------------------------------------------- biology

def biology_checks(events, competition, gene_sets):
    """Events and competition pairs that involve the genes of named sets ({name: (genes, expected direction)}),
    reported per set: how often a gene of the set is the gaining / losing promoter and how its expression moved."""
    rows = []
    for name, (genes, direction) in gene_sets.items():
        genes = set(genes)

        def hits(column):
            return events[column].fillna('').str.split(',').apply(lambda gs: bool(genes & set(gs))).to_numpy() \
                if len(events) else np.zeros(0, dtype=bool)
        as_gained, as_lost = hits('genes_gained'), hits('genes_lost')
        sub_g, sub_l = events[as_gained], events[as_lost]
        named = sorted({g for gs in pd.concat([sub_g['genes_gained'], sub_l['genes_lost']]).fillna('')
                        for g in gs.split(',') if g in genes})
        if len(competition):
            c1 = competition['genes_1'].fillna('').str.split(',').apply(lambda gs: bool(genes & set(gs))).to_numpy()
            c2 = competition['genes_2'].fillna('').str.split(',').apply(lambda gs: bool(genes & set(gs))).to_numpy()
        else:
            c1 = c2 = np.zeros(0, dtype=bool)
        tested = competition['tested'].to_numpy() if len(competition) else np.zeros(0, dtype=bool)
        rows.append({
            'gene_set': name, 'expected_direction': direction, 'genes_in_events': ','.join(named),
            'events_as_gained': int(as_gained.sum()), 'events_as_lost': int(as_lost.sum()),
            'gained_up': int((sub_g['rna_log2FC_gained'] > 0).sum()),
            'lost_down': int((sub_l['rna_log2FC_lost'] < 0).sum()),
            'concordant_events': int(events.loc[as_gained | as_lost, 'concordant'].sum()),
            'events_expression': int(events.loc[as_gained | as_lost, 'expression_available'].sum()),
            'competition_pairs_winner': int((c1 & tested).sum()), 'competition_pairs_loser': int((c2 & tested).sum()),
            'competitive_winner': int((c1 & tested & competition['competitive'].to_numpy()).sum())
            if len(competition) else 0,
            'competitive_loser': int((c2 & tested & competition['competitive'].to_numpy()).sum())
            if len(competition) else 0})
    return pd.DataFrame(rows)


def parse_gene_sets(items):
    """NAME=GENE1,GENE2[:up|down] or NAME=@file[:up|down] (one gene per line)."""
    sets = {}
    for item in items or []:
        name, spec = item.split('=', 1)
        direction = 'any'
        match = re.match(r'^(.*):(up|down|any)$', spec)
        if match:
            spec, direction = match.group(1), match.group(2)
        if spec.startswith('@'):
            genes = [line.strip() for line in open(spec[1:]) if line.strip()]
        else:
            genes = [g for g in spec.split(',') if g]
        sets[name] = (genes, direction)
    return sets


# -------------------------------------------------------------------------------------------------------- figures

STATUS_COLORS = {'gained': '#b2182b', 'lost': '#2166ac', 'shared': '#9a9994'}
INK = '#1f1f1e'
ELEMENT_COLOR = '#d95f02'
LOST_PROMOTER_COLOR = '#2166ac'
GAINED_PROMOTER_COLOR = '#b2182b'


def _width(strength, thin=0.4, thick=3.0, top=3.0):
    if not np.isfinite(strength) or strength <= 0:
        return thin
    return thin + (thick - thin) * float(np.clip(np.log(strength) / np.log(top), 0, 1)) if strength > 1 else thin


def plot_event_arcs(nodes, edges, event, names, condition, path, pad_fraction=0.15, min_pad=20000):
    """Arc plot of the region of a switch event in one condition (0 = A, 1 = B): the contacts of E, A and B present
    in that condition, width by strength, lost contacts (present in A only) blue dashed, gained (B only) red, shared
    grey; baited promoters as squares, other fragments as circles; E, A and B highlighted."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Arc

    e, a, b = event['element'], event['promoter_lost'], event['promoter_gained']
    mid = (nodes['start'] + nodes['end']) / 2
    key = [mid[e], mid[a], mid[b]]
    span = max(key) - min(key)
    pad = max(min_pad, pad_fraction * span)
    lo, hi = min(key) - pad, max(key) + pad
    chrom = nodes.at[e, 'chrom']
    focus = {e, a, b}
    status_present = ('lost', 'shared') if condition == 0 else ('gained', 'shared')
    sel = edges[(edges['u'].isin(focus) | edges['v'].isin(focus)) & edges['status'].isin(status_present)]
    mu, mv = mid.reindex(sel['u']).to_numpy(), mid.reindex(sel['v']).to_numpy()
    inside = (nodes['chrom'].reindex(sel['u']).to_numpy() == chrom) & \
        (nodes['chrom'].reindex(sel['v']).to_numpy() == chrom) & \
        (np.minimum(mu, mv) >= lo) & (np.maximum(mu, mv) <= hi)
    sel = sel[inside]
    strength_col = 'r_a' if condition == 0 else 'r_b'
    fig, ax = plt.subplots(figsize=(6.2, 2.6))
    longest = max(hi - lo, 1.0)
    drawn = set()
    for _, row in sel.iterrows():
        x0, x1 = sorted((mid[row['u']], mid[row['v']]))
        status = row['status']
        pair = {row['u'], row['v']}
        is_key = pair in ({e, a}, {e, b})
        height = 2 * np.sqrt(max(x1 - x0, 1) / longest)
        ax.add_patch(Arc(((x0 + x1) / 2, 0), x1 - x0, height, theta1=0, theta2=180, color=STATUS_COLORS[status],
                         lw=_width(float(row[strength_col])) * (1.4 if is_key else 1.0),
                         ls=(0, (3, 1.5)) if status == 'lost' else '-', alpha=1.0 if is_key else 0.45,
                         zorder=3 if is_key else 2))
        drawn |= pair
    drawn |= focus
    for n in drawn:
        x = mid[n]
        if not (lo <= x <= hi):
            continue
        bait = bool(nodes.at[n, 'bait'])
        if n == e:
            ax.plot([x], [0], marker='s' if bait else 'o', ms=8, mfc=ELEMENT_COLOR, mec=INK, mew=0.8, zorder=6)
        elif n == a:
            ax.plot([x], [0], marker='s', ms=8, mfc=LOST_PROMOTER_COLOR, mec=INK, mew=0.8, zorder=6)
        elif n == b:
            ax.plot([x], [0], marker='s', ms=8, mfc=GAINED_PROMOTER_COLOR, mec=INK, mew=0.8, zorder=6)
        elif bait:
            ax.plot([x], [0], marker='s', ms=4, color=INK, zorder=5)
        else:
            ax.plot([x], [0], marker='o', ms=3.5, mfc='white', mec='#5f5e5a', mew=0.7, zorder=5)

    def gene(n):
        genes = [g for g in nodes.at[n, 'genes'].split(',') if g]
        named = [g for g in genes if '-AS' not in g and not g.startswith(('Gm', 'LOC', 'RP11', 'AC0', 'AL0'))]
        return (named or genes or [n])[0]

    labels = [(e, 'E', ELEMENT_COLOR, 'normal'), (a, gene(a), LOST_PROMOTER_COLOR, 'italic'),
              (b, gene(b), GAINED_PROMOTER_COLOR, 'italic')]
    for n, text, color, style in labels:
        ax.text(mid[n], -0.22, text, rotation=90, ha='center', va='top', fontsize=7, style=style, color=color,
                fontweight='bold')
    ax.axhline(0, color=INK, lw=0.6, zorder=1)
    ax.set_xlim(lo, hi)
    ax.set_ylim(-1.15, 2.15)
    ax.set_yticks([])
    for side in ('left', 'right', 'top'):
        ax.spines[side].set_visible(False)
    ax.set_xlabel('{} (Mb)'.format(chrom), fontsize=7)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: '{:.2f}'.format(x / 1e6)))
    ax.tick_params(labelsize=6)
    ax.set_title('{}: E {} with A {} ({:+.1f}) and B {} ({:+.1f})'.format(
        names[condition], e, gene(a), float(event['rna_log2FC_lost']), gene(b), float(event['rna_log2FC_gained'])),
        fontsize=7.5, loc='left')
    handles = [plt.Line2D([0], [0], color=STATUS_COLORS['shared'], lw=1.2, label='shared'),
               plt.Line2D([0], [0], color=STATUS_COLORS['lost'], lw=1.2, ls=(0, (3, 1.5)),
                          label='lost in {}'.format(names[1])),
               plt.Line2D([0], [0], color=STATUS_COLORS['gained'], lw=1.2, label='gained in {}'.format(names[1])),
               plt.Line2D([0], [0], marker='s', ls='', color=INK, ms=4, label='baited promoter'),
               plt.Line2D([0], [0], marker='o', ls='', mfc='white', mec='#5f5e5a', ms=4, label='other end'),
               plt.Line2D([0], [0], marker='o', ls='', mfc=ELEMENT_COLOR, mec=INK, ms=6, label='E (switching)'),
               plt.Line2D([0], [0], marker='s', ls='', mfc=LOST_PROMOTER_COLOR, mec=INK, ms=6,
                          label='A (contact lost)'),
               plt.Line2D([0], [0], marker='s', ls='', mfc=GAINED_PROMOTER_COLOR, mec=INK, ms=6,
                          label='B (contact gained)')]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.32), ncol=4, fontsize=5.8,
              frameon=False, handlelength=2.2)
    ax.text(1.0, 1.02, 'line width: interaction strength in {}'.format(names[condition]), transform=ax.transAxes,
            fontsize=5.8, ha='right', va='bottom', color='#5f5e5a')
    fig.subplots_adjust(left=0.03, right=0.98, top=0.88, bottom=0.42)
    fig.savefig(ensure_parent(path))
    plt.close(fig)


def plot_null(null, observed, path, title, xlabel):
    """Histogram of a null distribution with the observed value."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(3.2, 2.4))
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    ax.hist(null, bins=50, color='#9a9994', edgecolor='none')
    ax.axvline(observed, color='#b2182b', lw=1.4)
    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel('permutations', fontsize=7)
    ax.tick_params(labelsize=6)
    ax.set_title(title, fontsize=7.5, loc='left')
    for side in ('right', 'top'):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(ensure_parent(path))
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------- driver

def summarise(label, first, second, pool, levels, distances, promoters, n_perm, seed, threshold):
    observed, perm, rot = null_test(first, second, pool, levels, distances, promoters, n_perm=n_perm, seed=seed,
                                    threshold=threshold)
    row = {'set': label, 'n_tested': len(first)}
    for k, v in observed.items():
        row[k] = v
        row[k + '_null_permutation'] = float(np.mean(perm[k]))
        row[k + '_p_permutation'] = permutation_p(v, perm[k])
        row[k + '_null_rotation'] = float(np.mean(rot[k]))
        row[k + '_p_rotation'] = permutation_p(v, rot[k])
    return row, perm, rot


def run(prefix, out_prefix, names=None, expression=None, gene_column='gene_name', level_columns=None,
        elements_path=None, contact_threshold=1.0, rna_threshold=1.0, chromatin_threshold=1.0, n_perm=10000,
        seed=0, gene_sets=None, top=3, figures=True):
    started = time.time()
    nodes, edges = read_union_graphml(prefix + '_union.graphml')
    names = tuple(names) if names else condition_names(prefix, nodes)
    marks = mark_names(nodes, names)
    log.info('%d nodes, %d edges, conditions %s, marks %s: %.0f s', len(nodes), len(edges), names, marks,
             time.time() - started)
    elements = read_elements(elements_path if elements_path else prefix + '_elements.tsv', names)
    log.info('element classes: %s', 'yes' if elements is not None else 'not available')
    levels = expression_levels(nodes, expression, gene_column, level_columns)
    promoters = nodes[promoter_mask(nodes)]
    rna = promoters['rna_log2FC']

    events = switch_events(nodes, edges, names, marks, elements, contact_threshold, rna_threshold,
                           chromatin_threshold)
    write_table(events, out_prefix + '_switch_events.tsv')
    log.info('%d switch events (%d with a non-promoter element): %.0f s', len(events),
             int((events['element_type'] == 'non-promoter element').sum()), time.time() - started)

    contacts = directed_contacts(nodes, edges, contact_threshold)
    rows, nulls = [], {}
    is_element = ~element_is_promoter(nodes, elements)
    for label, mask in (('non-promoter element', events['element_type'] == 'non-promoter element'),
                        ('any fragment', np.ones(len(events), dtype=bool))):
        sub = events[mask & events['expression_available']]
        pool = pd.Index(pd.unique(pd.concat([sub['promoter_lost'], sub['promoter_gained']])))
        dist = pd.concat([sub[['promoter_lost', 'distance_lost']].set_axis(['p', 'd'], axis=1),
                          sub[['promoter_gained', 'distance_gained']].set_axis(['p', 'd'], axis=1)])
        distances = dist.groupby('p')['d'].median()
        if len(sub):
            row, perm, rot = summarise(label, sub['promoter_lost'].to_numpy(), sub['promoter_gained'].to_numpy(),
                                       pool, levels, distances, promoters, n_perm, seed, rna_threshold)
        else:
            row, perm, rot = {'set': label, 'n_tested': 0}, {}, {}
        all_events = events[mask]
        row.update({'n_events': int(len(all_events)), 'n_elements': int(all_events['element'].nunique()),
                    'n_promoters_lost': int(all_events['promoter_lost'].nunique()),
                    'n_promoters_gained': int(all_events['promoter_gained'].nunique()),
                    'discordant': float(sub['discordant'].mean()) if len(sub) else np.nan,
                    'n_concordant': int(sub['concordant'].sum()), 'n_discordant': int(sub['discordant'].sum()),
                    'element_chromatin_changed': float(all_events['element_chromatin_changed'].mean())
                    if len(all_events) else np.nan,
                    'concordant_chromatin_changed': float(sub.loc[sub['element_chromatin_changed'],
                                                                  'concordant'].mean())
                    if sub['element_chromatin_changed'].any() else np.nan,
                    'concordant_chromatin_kept': float(sub.loc[~sub['element_chromatin_changed'],
                                                               'concordant'].mean())
                    if (~sub['element_chromatin_changed']).any() else np.nan})
        swap_mean, swap_p, swap = role_swap_test(sub['concordant'], sub['discordant'], sub['element'], n_perm, seed)
        row.update({'concordant_null_swap': swap_mean, 'concordant_p_swap': swap_p})
        ctrl_contacts = contacts if label == 'any fragment' else \
            contacts[is_element.reindex(contacts['element']).to_numpy()]
        ordered, concordant = stable_pair_control(ctrl_contacts, rna)
        n_tested = int(row['n_tested'])
        row.update({'stable_pairs_ordered': ordered,
                    'stable_pairs_concordant': concordant / ordered if ordered else np.nan})
        if ordered and n_tested:
            row['stable_pairs_fisher_p'] = stats.fisher_exact(
                [[row['n_concordant'], n_tested - row['n_concordant']],
                 [concordant, ordered - concordant]], alternative='greater')[1]
        one = one_sided_control(all_events, ctrl_contacts, rna)
        for side, values in one.items():
            for k, v in values.items():
                row['{}_{}'.format(side, k)] = v
        rows.append(row)
        nulls[label] = (row.get('concordant', np.nan), {'permutation': perm.get('concordant'),
                                                         'rotation': rot.get('concordant'), 'swap': swap})

    hubs, left, joined, p_nodes = hub_switches(nodes, edges, names, elements)
    write_table(hubs, out_prefix + '_switch_hubs.tsv')
    n_hub, hub_obs, hub_perm, hub_rot = hub_null_test(hubs, left, joined, p_nodes, nodes, levels, n_perm=n_perm,
                                                      seed=seed)
    tested_hubs = hubs[hubs['expression_available']]
    hub_swap_mean, hub_swap_p, hub_swap = role_swap_test(tested_hubs['concordant'], tested_hubs['discordant'],
                                                         tested_hubs['element'], n_perm, seed)
    rows.append({'set': 'hub switch (two steps)', 'n_events': int(len(hubs)), 'n_elements': int(len(hubs)),
                 'n_tested': n_hub, 'concordant': hub_obs,
                 'n_concordant': int(hubs['concordant'].sum()), 'n_discordant': int(hubs['discordant'].sum()),
                 'discordant': float(hubs.loc[hubs['expression_available'], 'discordant'].mean())
                 if hubs['expression_available'].any() else np.nan,
                 'concordant_null_permutation': float(np.mean(hub_perm)) if len(hub_perm) else np.nan,
                 'concordant_p_permutation': permutation_p(hub_obs, hub_perm),
                 'concordant_null_rotation': float(np.mean(hub_rot)) if len(hub_rot) else np.nan,
                 'concordant_p_rotation': permutation_p(hub_obs, hub_rot),
                 'concordant_null_swap': hub_swap_mean, 'concordant_p_swap': hub_swap_p})
    nulls['hub switch'] = (hub_obs, {'permutation': hub_perm, 'rotation': hub_rot, 'swap': hub_swap})
    log.info('%d hub switches: %.0f s', len(hubs), time.time() - started)

    competition, _ = competition_pairs(nodes, edges, contact_threshold, elements)
    write_table(competition, out_prefix + '_competition.tsv')
    tested = competition[competition['tested']]
    comp_row = {'set': 'competition', 'n_events': int(len(competition)), 'n_tested': int(len(tested)),
                'n_concordant': int(tested['competitive'].sum()), 'n_discordant': int(tested['reverse'].sum())}
    if len(tested):
        dist = pd.concat([tested[['promoter_1', 'promoter_distance']].set_axis(['p', 'd'], axis=1),
                          tested[['promoter_2', 'promoter_distance']].set_axis(['p', 'd'], axis=1)])
        pool = pd.Index(pd.unique(dist['p']))
        # competitive: promoter 1 (gains more) up, promoter 2 down = pair statistic with first = 2, second = 1
        row, perm, rot = summarise('competition', tested['promoter_2'].to_numpy(), tested['promoter_1'].to_numpy(),
                                   pool, levels, dist.groupby('p')['d'].median(), promoters, n_perm, seed,
                                   rna_threshold)
        comp_row.update({k: v for k, v in row.items() if k not in ('set', 'n_tested')})
        comp_row['discordant'] = float(tested['reverse'].mean())
        rho, p = stats.spearmanr(tested['delta'], tested['rna_log2FC_1'] - tested['rna_log2FC_2'])
        comp_row.update({'spearman_delta_vs_difference': rho, 'spearman_p': p})
        swap_mean, swap_p, swap = role_swap_test(tested['competitive'], tested['reverse'], np.arange(len(tested)),
                                                 n_perm, seed)
        comp_row.update({'concordant_null_swap': swap_mean, 'concordant_p_swap': swap_p})
        nulls['competition'] = (row['concordant'], {'permutation': perm['concordant'], 'rotation': rot['concordant'],
                                                    'swap': swap})
    rows.append(comp_row)
    summary = pd.DataFrame(rows)
    summary.insert(0, 'comparison', os.path.basename(out_prefix))
    write_table(summary, out_prefix + '_switch_summary.tsv')
    log.info('%d competition pairs (%d tested): %.0f s', len(competition), len(tested), time.time() - started)

    if gene_sets:
        write_table(biology_checks(events, competition, gene_sets), out_prefix + '_switch_biology.tsv')

    if figures:
        base = os.path.join(os.path.dirname(out_prefix), 'figures', os.path.basename(out_prefix))
        main = events[(events['element_type'] == 'non-promoter element') & events['expression_available']]
        chosen = main.drop_duplicates('element').head(top)
        for k, (_, event) in enumerate(chosen.iterrows(), start=1):
            for c in (0, 1):
                plot_event_arcs(nodes, edges, event, names, c,
                                '{}_switch_top{}_{}.pdf'.format(base, k, names[c]))
        for label, key in (('non-promoter element', 'switch_null'), ('hub switch', 'hub_switch_null'),
                           ('competition', 'competition_null')):
            if label not in nulls:
                continue
            obs, draws = nulls[label]
            for kind, null in draws.items():
                if null is not None and len(null):
                    plot_null(null, obs, '{}_{}_{}.pdf'.format(base, key, kind),
                              '{} ({})'.format(label, kind),
                              'fraction concordant' if label != 'competition' else 'fraction competitive')
        np.savez_compressed(out_prefix + '_switch_nulls.npz',
                            **{'{}__{}'.format(label.replace(' ', '_'), kind): np.asarray(null if null is not None
                                                                                          else [])
                               for label, (_, draws) in nulls.items() for kind, null in draws.items()})
    with open(out_prefix + '_switch_run.json', 'w') as handle:
        json.dump({'prefix': prefix, 'names': names, 'marks': marks, 'elements': elements is not None,
                   'contact_threshold': contact_threshold, 'rna_threshold': rna_threshold,
                   'chromatin_threshold': chromatin_threshold, 'permutations': n_perm, 'seed': seed,
                   'expression': expression, 'seconds': round(time.time() - started)}, handle, indent=1)
    log.info('done: %.0f s', time.time() - started)
    return events, summary, competition


def parse_arguments(args=None):
    parser = base_parser('Enhancer switching and promoter competition between the two conditions of an '
                         'egaDifferential run: elements that lose a contact with one promoter and gain one with '
                         'another, the expression changes of both promoters against null models, and promoters that '
                         'share elements.')
    parser.add_argument('--prefix', '-p', required=True, help='Output prefix of egaDifferential (<prefix>_union.graphml)')
    parser.add_argument('--outPrefix', '-o', required=True, help='Output prefix of the tables (figures in '
                                                                  '<dirname>/figures/)')
    parser.add_argument('--names', nargs=2, default=None, help='Condition names A B (default: from the prefix)')
    parser.add_argument('--expression', help='Expression table for the expression level of the promoters '
                                             '(baseMean or mean_log2_* columns); without it, the strata of the '
                                             'permutation use the distance class only')
    parser.add_argument('--geneColumn', default='gene_name')
    parser.add_argument('--levelColumns', nargs='+', default=None,
                        help='Columns of --expression with the (log-scale) expression level, averaged')
    parser.add_argument('--elements', default=None, help='Element classes (default: <prefix>_elements.tsv if present)')
    parser.add_argument('--contactThreshold', type=float, default=1.0,
                        help='|log2 (B / A)| of a contact counted as lost / gained besides the status')
    parser.add_argument('--rnaThreshold', type=float, default=1.0, help='|log2FC| of the strict concordance')
    parser.add_argument('--chromatinThreshold', type=float, default=1.0,
                        help='|log2 change| of a mark at E counted as a chromatin change')
    parser.add_argument('--permutations', type=int, default=10000)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--geneSets', nargs='*', default=None,
                        help='Gene sets to check: NAME=GENE,GENE[:up|down] or NAME=@file[:up|down]')
    parser.add_argument('--top', type=int, default=3, help='Arc plots of the top events (distinct elements)')
    parser.add_argument('--noFigures', action='store_true')
    return parser.parse_args(args)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    run(args.prefix, args.outPrefix, names=args.names, expression=args.expression, gene_column=args.geneColumn,
        level_columns=args.levelColumns, elements_path=args.elements, contact_threshold=args.contactThreshold,
        rna_threshold=args.rnaThreshold, chromatin_threshold=args.chromatinThreshold, n_perm=args.permutations,
        seed=args.seed, gene_sets=parse_gene_sets(args.geneSets), top=args.top, figures=not args.noFigures)


if __name__ == '__main__':
    main()
