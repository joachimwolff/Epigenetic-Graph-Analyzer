"""Structural variants and copy number as a level of the explanation cascade (egaDifferential --structuralVariants,
--svContacts, --copyNumber).

Inputs
    Breakpoints: BEDPE (chrom1 start1 end1 chrom2 start2 end2 [name score strand1 strand2 ...]) with the conditions
    in which a join is present, given by a column 'conditions' (comma-separated condition names), by columns
    present_<condition> (0/1) and / or pairs_<condition> / obs_exp_<condition> (read support, observed over expected),
    or for a plain BEDPE by an 11th column with the condition names; 'condition=file.bedpe' marks every join of the
    file as present in that condition. The output of analysis/capture_sv_calls.py (<prefix>_sv.bedpe) has these
    columns. The inter-chromosomal hotspot table of analysis/capture_sv_cnv.py (chrom1 start1 chrom2 start2 pairs
    expected obs_exp, 1 Mb bin pairs) is read as 'condition=<prefix>_trans_hotspots.tsv', one file per condition:
    bin pairs with at least `min_pairs` pairs and an observed / expected of at least `min_obs_exp` are breakpoints,
    bin pairs within 2 Mb of each other on both sides form one join.
    Junction contacts (optional, <prefix>_sv_contacts.tsv of capture_sv_calls.py): read pairs between a baited
    fragment near a breakend and 100 kb bins near the partner breakend, per condition.
    Copy number (optional): bedGraph of the relative copy number per condition ('condition=file.bedGraph').

Rules (per differentially expressed gene; the window is --svWindow)
    structural variant: the promoter fragment lies within the window of a breakend. The join is gained in B when it
    is present in B only or its observed / expected (pairs when no expected is given) rises by at least
    `sv_threshold` (log2) from A to B, lost in B in the opposite case, present in both otherwise; a join called in
    one condition only whose support changes by less than `sv_threshold` has a 'similar support in both' (the call
    differs by the calling thresholds only) and counts as unchanged. The level explains
    an increase in expression when the join is gained in B and the partner region (the window around the other
    breakend) holds at least one node of the graph with an active element class in B, a decrease when the join is
    lost in B and the partner region holds an active element in A (the enhancers placed next to the promoter by the
    rearrangement are gained or lost). With junction contacts of the promoter fragment, its own contacts across the
    junction must change in the direction of the expression change as well. A join present in both conditions with
    a similar support explains nothing ('present in both') and is reported.
    copy number: the relative copy number at the promoter changes by at least `cn_threshold` (log2 B/A) in the
    direction of the expression change.

The level reports, it does not decide: the breakpoint, the partner region with its element classes, marks and genes,
the read support per condition and the contacts of the promoter across the junction are written for every gene near
a breakpoint (<prefix>_sv_genes.tsv), whether or not it is differential; genes with the promoter on one side and an
active element on the other side are enhancer-hijacking candidates (hijacking_candidate). Orientation of the
rearranged segments (which derivative chromosome joins which side) is not resolved; both directions within the
window are taken.
"""

import logging
import re

import numpy as np
import pandas as pd

from .io import open_text

log = logging.getLogger(__name__)

SV_LEVELS = ('structural_variant', 'copy_number')
ACTIVE_CLASSES = ('active promoter', 'active enhancer')
BEDPE_COLUMNS = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']


def _chrom_number(chrom):
    name = str(chrom)
    return name[3:] if name.lower().startswith('chr') else name


def _chrom_key(chrom):
    name = _chrom_number(chrom)
    return (0, int(name), '') if name.isdigit() else (1, 0, name)


def sv_label(chrom1, chrom2, start1=None, start2=None):
    """'t(8;21)' for an inter-chromosomal join (chromosomes in karyotype order), 'chr8 join 12.3-45.6 Mb' for an
    intra-chromosomal one."""
    if chrom1 != chrom2:
        a, b = sorted([chrom1, chrom2], key=_chrom_key)
        return 't({};{})'.format(_chrom_number(a), _chrom_number(b))
    if start1 is None:
        return '{} join'.format(chrom1)
    return '{} join {:.1f}-{:.1f} Mb'.format(chrom1, min(start1, start2) / 1e6, max(start1, start2) / 1e6)


def _split_item(item):
    """('condition', 'path') of 'condition=path'; (None, path) otherwise (a path may contain '=' only after a '/')."""
    head, sep, tail = item.partition('=')
    if sep and '/' not in head and head:
        return head, tail
    return None, item


def _read_header(path):
    with open_text(path) as handle:
        for line in handle:
            if line.strip():
                return line.rstrip('\n').lstrip('#').split('\t')
    return []


def _hotspot_joins(tables, names, min_pairs, min_obs_exp, gap=2000000):
    """Joins from capture_sv_cnv.py hotspot tables {condition: table}: candidate 1 Mb bin pairs of every condition
    clustered within `gap` on both sides; per condition the strongest bin pair of the cluster."""
    candidates = []
    for condition, table in tables.items():
        hit = table[(table['pairs'] >= min_pairs) & (table['obs_exp'] >= min_obs_exp)]
        candidates += [(r.chrom1, int(r.start1), r.chrom2, int(r.start2)) for r in hit.itertuples(index=False)]
    candidates = sorted(set(candidates))
    groups = []
    for key in candidates:
        for group in groups:
            if any(k[0] == key[0] and k[2] == key[2] and abs(k[1] - key[1]) <= gap and abs(k[3] - key[3]) <= gap
                   for k in group):
                group.append(key)
                break
        else:
            groups.append([key])
    rows = []
    indexed = {c: t.set_index(['chrom1', 'start1', 'chrom2', 'start2']) for c, t in tables.items()}
    for group in groups:
        c1, c2 = group[0][0], group[0][2]
        row = dict(chrom1=c1, chrom2=c2)
        best_pairs, best_key = -1, group[0]
        for condition in names:
            table = indexed.get(condition)
            part = table.reindex(group) if table is not None else pd.DataFrame(columns=['pairs', 'expected',
                                                                                        'obs_exp'])
            pairs = part['pairs'].fillna(0).to_numpy(dtype=float)
            expected = part['expected'].fillna(0).to_numpy(dtype=float)
            ok = (part['pairs'].fillna(0) >= min_pairs) & (part['obs_exp'].fillna(0) >= min_obs_exp)
            row['present_' + condition] = bool(ok.any()) if table is not None else np.nan
            row['pairs_' + condition] = float(pairs.sum()) if table is not None else np.nan
            row['expected_' + condition] = float(expected.sum()) if table is not None else np.nan
            row['obs_exp_' + condition] = float(pairs.sum() / expected.sum()) \
                if table is not None and expected.sum() > 0 else np.nan
            if len(pairs) and pairs.max() > best_pairs:
                best_pairs, best_key = pairs.max(), group[int(np.argmax(pairs))]
        row.update(start1=best_key[1], end1=best_key[1] + 1000000, start2=best_key[3], end2=best_key[3] + 1000000)
        rows.append(row)
    return pd.DataFrame(rows)


def read_structural_variants(items, names, min_pairs=50, min_obs_exp=20.0):
    """Joins of the inputs `items` (paths or 'condition=path', see the module) for the conditions `names` (A, B).

    Returns one row per join: sv_id, name, type (inter / intra), chrom1, start1, end1, chrom2, start2, end2, pos1, pos2
    (breakend midpoints), per condition present_<c> (bool), pairs_<c>, obs_exp_<c> (NaN when not given) and source.
    """
    names = list(names)
    frames, hotspots = [], {}
    for item in items:
        condition, path = _split_item(item)
        if condition is not None and condition not in names:
            raise ValueError('--structuralVariants: condition {} is not one of {}'.format(condition, names))
        header = _read_header(path)
        if {'chrom1', 'start1', 'chrom2', 'start2', 'pairs', 'obs_exp'} <= set(header) and 'end1' not in header:
            if condition is None:
                raise ValueError('{}: a hotspot table of capture_sv_cnv.py needs condition=path'.format(path))
            hotspots[condition] = pd.read_csv(path, sep='\t')
            continue
        has_header = bool(header) and header[0].lower() in ('chrom1', 'chr1', 'chromosome1')
        table = pd.read_csv(path, sep='\t', header=None, comment=None, skiprows=1 if has_header else 0,
                            names=header if has_header else None, dtype={0: str, 3: str} if not has_header else None)
        if not has_header:
            columns = BEDPE_COLUMNS + ['name', 'score', 'strand1', 'strand2', 'conditions']
            table.columns = columns[:table.shape[1]] + ['extra{}'.format(k) for k in
                                                        range(table.shape[1] - len(columns))]
        table = table[~table['chrom1'].astype(str).str.startswith('#')].copy()
        for column in ('start1', 'end1', 'start2', 'end2'):
            table[column] = table[column].astype(np.int64)
        for c in names:
            present = pd.Series(np.nan, index=table.index, dtype=object)
            if condition is not None:
                present[:] = c == condition
            elif 'present_' + c in table:
                present = table['present_' + c].astype(float) > 0
            elif 'conditions' in table:
                present = table['conditions'].fillna('').astype(str).map(
                    lambda v: c in [x.strip() for x in v.split(',')])
            table['present_' + c] = present
            for key in ('pairs', 'obs_exp'):
                if '{}_{}'.format(key, c) not in table:
                    table['{}_{}'.format(key, c)] = np.nan
                table['{}_{}'.format(key, c)] = table['{}_{}'.format(key, c)].astype(float)
        table['source'] = path
        frames.append(table)
    if hotspots:
        joins = _hotspot_joins(hotspots, names, min_pairs, min_obs_exp)
        joins['source'] = ','.join(sorted(_split_item(i)[1] for i in items if _split_item(i)[0] in hotspots))
        frames.append(joins)
    keep = BEDPE_COLUMNS + ['name'] + ['{}_{}'.format(k, c) for c in names for k in ('present', 'pairs', 'obs_exp')] \
        + ['source']
    if not frames:
        return pd.DataFrame(columns=['sv_id', 'type', 'pos1', 'pos2'] + keep)
    svs = pd.concat([f.reindex(columns=keep) for f in frames], ignore_index=True)
    # several files may describe the same join for different conditions: rows with the same breakend bins merged
    svs['chrom1'] = svs['chrom1'].astype(str)
    svs['chrom2'] = svs['chrom2'].astype(str)
    merged = []
    for _, group in svs.groupby(BEDPE_COLUMNS, sort=False):
        row = group.iloc[0].copy()
        for c in names:
            values = group['present_' + c].dropna()
            row['present_' + c] = bool(values.astype(bool).any()) if len(values) else np.nan
            for key in ('pairs', 'obs_exp'):
                values = group['{}_{}'.format(key, c)].dropna()
                row['{}_{}'.format(key, c)] = float(values.max()) if len(values) else np.nan
        row['source'] = ','.join(sorted(set(group['source'].astype(str))))
        merged.append(row)
    svs = pd.DataFrame(merged).reset_index(drop=True)
    svs['type'] = np.where(svs['chrom1'] == svs['chrom2'], 'intra', 'inter')
    svs['pos1'] = ((svs['start1'] + svs['end1']) // 2).astype(np.int64)
    svs['pos2'] = ((svs['start2'] + svs['end2']) // 2).astype(np.int64)
    labels = [sv_label(r.chrom1, r.chrom2, r.pos1, r.pos2) for r in svs.itertuples(index=False)]
    svs['name'] = [str(n) if isinstance(n, str) and n and n != '.' else label
                   for n, label in zip(svs['name'], labels)]
    svs = svs[[c for c in svs if c != 'sv_id']]
    svs.insert(0, 'sv_id', ['sv{}'.format(k + 1) for k in range(len(svs))])
    for c in names:
        svs['present_' + c] = svs['present_' + c].map(lambda v: bool(v) if not pd.isna(v) else False)
    return svs


def sv_status(svs, names, threshold=1.0):
    """Change of every join from A to B: ('gained', +1), ('lost', -1), ('stronger in B', +1), ('weaker in B', -1),
    ('present in both', 0), ('similar support in both', 0: called in one condition only with a support change below
    `threshold`) or ('absent', 0), with the log2 change of the support (obs/exp + 1, else pairs + 1; NaN without
    values in both)."""
    a, b = names
    out = []
    for _, r in svs.iterrows():
        pa, pb = bool(r['present_' + a]), bool(r['present_' + b])
        change = np.nan
        for key in ('obs_exp', 'pairs'):
            va, vb = r.get('{}_{}'.format(key, a), np.nan), r.get('{}_{}'.format(key, b), np.nan)
            if not pd.isna(va) and not pd.isna(vb):
                change = float(np.log2((vb + 1.0) / (va + 1.0)))
                break
        if pa != pb and not pd.isna(change) and abs(change) < threshold:
            # called in one condition only, by the calling thresholds, with a similar support in both
            out.append(('similar support in both', 0, change))
        elif pa and not pb:
            out.append(('lost', -1, change))
        elif pb and not pa:
            out.append(('gained', 1, change))
        elif pa and pb:
            if not pd.isna(change) and abs(change) >= threshold:
                out.append(('stronger in {}'.format(b) if change > 0 else 'weaker in {}'.format(b),
                            int(np.sign(change)), change))
            else:
                out.append(('present in both', 0, change))
        else:
            out.append(('absent', 0, change))
    return pd.DataFrame(out, columns=['status', 'direction', 'support_log2'], index=svs.index)


def read_sv_contacts(path, names):
    """Junction contacts of capture_sv_calls.py: join, bait_chrom, bait_start, bait_end, partner_chrom,
    partner_start, partner_end, pairs_<c>, inter_total_<c> (per condition)."""
    table = pd.read_csv(path, sep='\t')
    for c in names:
        if 'pairs_' + c not in table:
            raise ValueError('{}: no column pairs_{}'.format(path, c))
        if 'inter_total_' + c not in table:
            table['inter_total_' + c] = np.nan
    return table


def read_copy_number(items, names):
    """{condition: DataFrame chrom, start, end, value} of 'condition=file.bedGraph' items."""
    out = {}
    for item in items:
        condition, path = _split_item(item)
        if condition not in names:
            raise ValueError('--copyNumber needs condition=file.bedGraph with a condition of {}'.format(list(names)))
        table = pd.read_csv(path, sep='\t', header=None, comment='#', usecols=[0, 1, 2, 3],
                            names=['chrom', 'start', 'end', 'value'])
        table = table[~table['chrom'].astype(str).str.startswith(('track', 'browser'))]
        table = table.astype({'start': np.int64, 'end': np.int64, 'value': float})
        out[condition] = table.sort_values(['chrom', 'start']).reset_index(drop=True)
    return out


def copy_number_at(table, chrom, position):
    if table is None:
        return np.nan
    part = table[(table['chrom'] == chrom) & (table['start'] <= position) & (table['end'] > position)]
    return float(part['value'].iloc[0]) if len(part) else np.nan


def _middle(data):
    return (int(data['start']) + int(data['end'])) // 2


def flag_sv_edges(graph, svs, window):
    """Mark the edges of `graph` that join the two sides of a breakpoint (one end within `window` of each breakend;
    for an intra-chromosomal join the edge must span the breakpoint): sv=True, sv_name. Returns their number."""
    flagged = 0
    for u, v, data in graph.edges(data=True):
        du, dv = graph.nodes[u], graph.nodes[v]
        for r in svs.itertuples(index=False):
            hit = False
            for (cu, pu), (cv, pv) in (((r.chrom1, r.pos1), (r.chrom2, r.pos2)),
                                       ((r.chrom2, r.pos2), (r.chrom1, r.pos1))):
                if du.get('chrom') == cu and dv.get('chrom') == cv and abs(_middle(du) - pu) <= window and \
                        abs(_middle(dv) - pv) <= window:
                    if cu == cv:
                        lo, hi = sorted((_middle(du), _middle(dv)))
                        if not (lo < min(pu, pv) and hi > max(pu, pv)):
                            continue
                    hit = True
            if hit:
                data['sv'] = True
                data['sv_name'] = r.name
                flagged += 1
                break
    return flagged


def _nodes_near(nodes, chrom, position, window):
    part = nodes[(nodes['chrom'] == chrom) & ((nodes['middle'] - position).abs() <= window)]
    return part.assign(distance=(part['middle'] - position).abs()).sort_values('distance')


def _contacts_of(contacts, node_data, join, partner_chrom, names):
    """Pairs of a promoter fragment across the junction `join` per condition (baits overlapping the fragment)."""
    if contacts is None:
        return None
    part = contacts[(contacts['join'] == join) & (contacts['bait_chrom'] == node_data['chrom']) &
                    (contacts['bait_end'] > int(node_data['start'])) &
                    (contacts['bait_start'] < int(node_data['end'])) & (contacts['partner_chrom'] == partner_chrom)]
    out = {}
    for c in names:
        total = contacts['inter_total_' + c].dropna()
        out[c] = (float(part['pairs_' + c].sum()), float(total.iloc[0]) if len(total) else np.nan)
    return out


def _fmt_mb(bases):
    return '{:.1f} Mb'.format(bases / 1e6) if bases >= 100000 else '{:.0f} kb'.format(bases / 1e3)


def _fmt_position(breakend):
    """'chr8:92.05 Mb' of 'chr8:92050000'."""
    chrom, position = str(breakend).rsplit(':', 1)
    return '{}:{:.2f} Mb'.format(chrom, int(position) / 1e6)


def _support_text(r, names):
    parts = []
    for c in names:
        pairs, oe = r.get('pairs_' + c, np.nan), r.get('obs_exp_' + c, np.nan)
        if pd.isna(pairs) and pd.isna(oe):
            continue
        text = '{} {}'.format(c, '{:,.0f} pairs'.format(pairs) if not pd.isna(pairs) else '')
        if not pd.isna(oe):
            text += ' obs/exp {:.0f}'.format(oe) if oe >= 10 else ' obs/exp {:.1f}'.format(oe)
        parts.append(text.strip())
    return ', '.join(parts)


def structural_genes(graph, change, svs, names, window=1000000, elements=None, marks=None, contacts=None,
                     copy_number=None, differential=None, sv_threshold=1.0, cn_threshold=0.5, element_column=None,
                     pseudo_count=1.0):
    """Genes whose promoter fragment lies within `window` of a breakend, and the copy number at every promoter.

    `change`: node -> RNA log2 B/A of the bait nodes; `svs` of read_structural_variants; `elements` the element
    table of node_elements (class column `element_column`, default class_declared_mean); `marks` mark -> (node ->
    signal A, node -> signal B); `contacts` read_sv_contacts; `copy_number` read_copy_number; `differential` the set
    of differential nodes (None: every node with |log2FC| >= 1 is not assumed, the column is left empty).

    Returns (genes, levels): genes is one row per gene and nearby join (every bait with an expression value, also
    non-differential: sv_genes.tsv); levels is indexed by node with sv_explains, sv_result, sv_detail, sv_by,
    cn_explains, cn_result, cn_detail, cn_by and the fields of the nearest join (for the cascade and the all-level
    table)."""
    a, b = names
    svs = svs.reset_index(drop=True)
    status = sv_status(svs, names, sv_threshold) if len(svs) else pd.DataFrame(columns=['status', 'direction',
                                                                                        'support_log2'])
    nodes = pd.DataFrame([(n, d['chrom'], (int(d['start']) + int(d['end'])) // 2, str(d.get('genes', '') or ''))
                          for n, d in graph.nodes(data=True)], columns=['node', 'chrom', 'middle', 'genes'])
    element_column = element_column or 'class_declared_mean'
    classes = {}
    if elements is not None and element_column in elements:
        for c in names:
            part = elements[elements['condition'] == c]
            classes[c] = dict(zip(part['node'], part[element_column]))
    rows = []
    for node, rna in change.items():
        data = graph.nodes[node]
        middle = _middle(data)
        for i, r in svs.iterrows():
            for side, (chrom, pos, partner_chrom, partner_pos) in enumerate(
                    ((r['chrom1'], r['pos1'], r['chrom2'], r['pos2']), (r['chrom2'], r['pos2'], r['chrom1'], r['pos1']))):
                if data['chrom'] != chrom or abs(middle - pos) > window:
                    continue
                partner = _nodes_near(nodes, partner_chrom, partner_pos, window)
                partner = partner[partner['node'] != node]
                row = dict(node=node, genes=data.get('genes', ''), rna_log2FC=float(rna),
                           differential=(node in differential) if differential is not None else np.nan,
                           sv_id=r['sv_id'], sv_name=r['name'], sv_type=r['type'],
                           breakend='{}:{}'.format(chrom, int(pos)), distance=int(abs(middle - pos)),
                           partner_breakend='{}:{}'.format(partner_chrom, int(partner_pos)),
                           partner_chrom=partner_chrom,
                           sv_status=status.loc[i, 'status'], sv_direction=int(status.loc[i, 'direction']),
                           support_log2=status.loc[i, 'support_log2'])
                for c in names:
                    row['present_' + c] = bool(r['present_' + c])
                    row['pairs_' + c] = r.get('pairs_' + c, np.nan)
                    row['obs_exp_' + c] = r.get('obs_exp_' + c, np.nan)
                row['partner_nodes'] = len(partner)
                genes = [g for g in partner['genes'] if g]
                names_near = []
                for g in genes:
                    for x in g.split(','):
                        if x and x not in names_near:
                            names_near.append(x)
                row['partner_genes'] = ','.join(names_near[:5])
                for c in names:
                    cls = [classes.get(c, {}).get(n, '') for n in partner['node']]
                    counts = pd.Series([x for x in cls if x]).value_counts()
                    row['partner_classes_' + c] = ';'.join('{}:{}'.format(k, v) for k, v in counts.items())
                    row['partner_active_' + c] = int(sum(x in ACTIVE_CLASSES for x in cls))
                for mark, (va, vb) in (marks or {}).items():
                    for c, values in ((a, va), (b, vb)):
                        row['partner_{}_{}'.format(mark, c)] = float(np.mean([values.get(n, 0.0)
                                                                              for n in partner['node']])) \
                            if len(partner) else np.nan
                own = _contacts_of(contacts, data, r['name'], partner_chrom, names)
                if own is not None:
                    for c in names:
                        row['junction_pairs_' + c] = own[c][0]
                    pa, ta = own[a]
                    pb, tb = own[b]
                    row['junction_log2'] = float(np.log2(((pb + pseudo_count) / tb) / ((pa + pseudo_count) / ta))) \
                        if ta and tb and not (pd.isna(ta) or pd.isna(tb)) else \
                        float(np.log2((pb + pseudo_count) / (pa + pseudo_count)))
                present_any = row['present_' + a] or row['present_' + b]
                active_where = [c for c in names if row['present_' + c] and row['partner_active_' + c] > 0]
                own_contacts = own is not None and (own[a][0] + own[b][0]) > 0
                row['hijacking_candidate'] = bool(present_any and active_where and
                                                  (own is None or own_contacts))
                rows.append(row)
    genes = pd.DataFrame(rows)
    # copy number at every promoter
    cn = {}
    if copy_number:
        for node in change:
            data = graph.nodes[node]
            values = {c: copy_number_at(copy_number.get(c), data['chrom'], _middle(data)) for c in names}
            cn[node] = values
    levels = {}
    direction_of = {n: np.sign(v) for n, v in change.items()}
    for node in change:
        d = direction_of[node]
        entry = dict(sv_explains=False, sv_result='no breakpoint nearby',
                     sv_detail='no breakpoint within {}'.format(_fmt_mb(window)), sv_by='', sv_name='',
                     sv_status='', sv_distance=np.nan, sv_partner='', sv_partner_genes='', sv_partner_classes='',
                     sv_support='', sv_junction_log2=np.nan, sv_breakends=0)
        part = genes[genes['node'] == node] if len(genes) else genes
        if len(part):
            texts, results = [], []
            for _, g in part.sort_values('distance').iterrows():
                status_text = g['sv_status']
                direction = g['sv_direction']
                where = names[1] if direction > 0 else names[0]
                active = g['partner_active_' + where] if direction != 0 else 0
                junction = g.get('junction_log2', np.nan)
                concordant_contacts = pd.isna(junction) or np.sign(junction) == d
                if direction == 0:
                    result = 'absent' if status_text == 'absent' else 'present in both'
                    explains = False
                elif direction != d:
                    result, explains = 'discordant', False
                elif active == 0:
                    result, explains = 'no active partner element', False
                elif not concordant_contacts:
                    result, explains = 'junction contacts discordant', False
                else:
                    result, explains = 'concordant', True
                partner_word = 'chromosome {}'.format(_chrom_number(g['partner_chrom'])) if g['sv_type'] == 'inter' \
                    else 'the far side of the join ({})'.format(_fmt_position(g['partner_breakend']))
                near = ' near {}'.format(', '.join(str(g['partner_genes']).split(',')[:2])) \
                    if g['partner_genes'] else ''
                text = 'promoter {} from the {} breakpoint ({} / {})'.format(
                    _fmt_mb(g['distance']), g['sv_name'], _fmt_position(g['breakend']),
                    _fmt_position(g['partner_breakend']))
                present = [c for c in names if g['present_' + c]]
                if status_text == 'present in both':
                    text += ', present in {} and {}'.format(*names)
                elif status_text in ('gained', 'lost'):
                    text += ', present in {} only'.format(present[0] if present else '?')
                else:
                    text += ', {}'.format(status_text)
                support = _support_text(g, names)
                if support:
                    text += ' ({})'.format(support)
                if not pd.isna(junction):
                    pa, pb = g.get('junction_pairs_' + names[0], 0), g.get('junction_pairs_' + names[1], 0)
                    verb = 'gained' if junction > 0 else 'lost' if junction < 0 else 'unchanged'
                    text += '; contacts to {}{} {} ({:.0f} vs {:.0f} read pairs, log2 {:+.1f})'.format(
                        partner_word, near, verb, pa, pb, junction)
                else:
                    text += '; partner {}{}'.format(partner_word, near)
                active_text = ', '.join('{} {}'.format(c, int(g['partner_active_' + c])) for c in names)
                text += '; active elements on the partner side: {} of {} nodes ({})'.format(
                    active_text, int(g['partner_nodes']), 'concordant' if explains else result)
                texts.append(text)
                results.append((explains, result, g))
            explaining = [x for x in results if x[0]]
            chosen = explaining[0] if explaining else results[0]
            g = chosen[2]
            entry.update(sv_explains=bool(explaining), sv_result=chosen[1], sv_detail='; '.join(texts),
                         sv_by='{} {}'.format(g['sv_name'], g['sv_status']) if explaining else '',
                         sv_name=g['sv_name'], sv_status=g['sv_status'], sv_distance=float(g['distance']),
                         sv_partner=g['partner_breakend'], sv_partner_genes=g['partner_genes'],
                         sv_partner_classes=' | '.join('{}: {}'.format(c, g['partner_classes_' + c])
                                                       for c in names),
                         sv_support=_support_text(g, names), sv_junction_log2=g.get('junction_log2', np.nan),
                         sv_breakends=len(part))
        if copy_number:
            va, vb = cn[node][a], cn[node][b]
            entry.update(cn_A=va, cn_B=vb)
            if pd.isna(va) or pd.isna(vb) or va <= 0 or vb <= 0:
                entry.update(cn_explains=False, cn_result='no value', cn_log2=np.nan,
                             cn_detail='no copy number value at the promoter', cn_by='')
            else:
                value = float(np.log2(vb / va))
                explains = abs(value) >= cn_threshold and np.sign(value) == d
                result = 'concordant' if explains else 'discordant' if abs(value) >= cn_threshold else \
                    'below threshold'
                entry.update(cn_explains=bool(explains), cn_result=result, cn_log2=value,
                             cn_by='copy {} (log2 {:+.2f})'.format('gain' if value > 0 else 'loss', value)
                             if explains else '',
                             cn_detail='relative copy number at the promoter {:.2f} in {}, {:.2f} in {} (log2 {:+.2f}'
                                       ', {})'.format(va, a, vb, b, value, result))
        levels[node] = entry
    levels = pd.DataFrame.from_dict(levels, orient='index')
    levels.index.name = 'node'
    return genes, levels


def sv_edge_table(graph, contacts, svs, names):
    """Junction contacts as edges: the promoter fragment (graph node overlapping the bait) to the partner bin, with
    the pairs per condition and a status (gained: pairs in B only, lost: in A only, shared), sv=True."""
    if contacts is None or not len(contacts):
        return pd.DataFrame()
    a, b = names
    baits = {}
    for n, d in graph.nodes(data=True):
        if d.get('bait'):
            baits.setdefault(d['chrom'], []).append((int(d['start']), int(d['end']), n))
    rows = []
    for r in contacts.itertuples(index=False):
        r = r._asdict()
        node = next((n for s, e, n in baits.get(r['bait_chrom'], []) if s < r['bait_end'] and e > r['bait_start']),
                    None)
        if node is None:
            continue
        pa, pb = float(r['pairs_' + a]), float(r['pairs_' + b])
        rows.append(dict(node_a=node, node_b='{}:{}-{}'.format(r['partner_chrom'], r['partner_start'],
                                                               r['partner_end']),
                         genes=graph.nodes[node].get('genes', ''), sv_name=r['join'], pairs_A=pa, pairs_B=pb,
                         status='gained' if pa == 0 and pb > 0 else 'lost' if pb == 0 and pa > 0 else 'shared',
                         sv=True))
    table = pd.DataFrame(rows)
    if len(table):
        table = table.groupby(['node_a', 'node_b', 'genes', 'sv_name', 'status', 'sv'], as_index=False)[
            ['pairs_A', 'pairs_B']].sum()
        table = table.rename(columns={'pairs_A': 'pairs_' + a, 'pairs_B': 'pairs_' + b})
    return table


def add_sv_edges(graph, edges):
    """A copy of `graph` with the junction contacts of sv_edge_table as edges sv=True (partner bins that are no
    node of the graph are added with sv_partner=True)."""
    out = graph.copy()
    for r in edges.itertuples(index=False):
        r = r._asdict()
        if r['node_b'] not in out:
            chrom, rest = r['node_b'].split(':')
            start, end = rest.split('-')
            out.add_node(r['node_b'], chrom=chrom, start=int(start), end=int(end), bait=False, sv_partner=True,
                         genes='')
        attributes = {k: v for k, v in r.items() if k not in ('node_a', 'node_b', 'genes')}
        attributes['sv'] = True
        out.add_edge(r['node_a'], r['node_b'], **attributes)
    return out


def level_text(level, entry):
    """Detail text of the structural variant / copy number level of one gene (entry: a row of the levels table)."""
    if level == 'structural_variant':
        return str(entry.get('sv_detail', ''))
    return str(entry.get('cn_detail', ''))


__all__ = ['SV_LEVELS', 'ACTIVE_CLASSES', 'sv_label', 'read_structural_variants', 'sv_status', 'read_sv_contacts',
           'read_copy_number', 'flag_sv_edges', 'structural_genes', 'sv_edge_table', 'add_sv_edges', 'level_text']


def _unused():
    return re
