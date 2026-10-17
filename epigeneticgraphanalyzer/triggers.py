"""Candidate trigger transcription factors of explained and poised genes.

The explanation cascade of egaDifferential names, for every differentially expressed gene, the level that explains
its change (promoter marks, own contacts, promoter state, marks or states of the contacting fragments, marks at the
nodes k graph steps away) or reports a poised (bivalent) promoter whose trigger is not among the measured layers.
This module adds the sequence layer: it scans the explaining nodes (and the promoter fragment) with position
frequency matrices, tests which motifs are enriched at the explaining nodes of a group of genes against matched
background nodes, and lists per gene the transcription factors whose motif occurs there and whose own expression
changes. The output reports and indicates candidates; it does not decide which factor triggers a gene.

Sequence scanned per node: the window of `window` bp (default 1 kb) with the highest active-mark signal within the
fragment (sum over the active / enhancer mark tracks of both conditions, each track divided by its mean over all
fragments), the whole fragment when it is shorter. Restricting the scan to this window avoids counting motifs in
several kb of mostly unmarked restriction fragment sequence; fragments without signal get their central window.

Explaining nodes of a gene at step k: the nodes at exactly k graph steps (path weight in either condition, as in
egaDifferential) whose own log2 change of an explaining mark (concordant_marks of the level) is concordant with the
expression change and at least the chromatin threshold, counted only when the level explains the gene. At step 1
the state fragment of an explaining step-1 state transition is added. The promoter fragment is always part of the
gene's scanned set; for genes explained by the promoter (marks or state) or poised there it is the explaining node.
"""

import logging
import math
import multiprocessing as mp
import re
from collections import defaultdict

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

log = logging.getLogger(__name__)

BASES = 'ACGT'
# gene name prefixes not used as transcription factor families (large zinc finger groups without a shared motif)
NO_FAMILY = ('ZNF', 'ZFP', 'ZSCAN', 'ZKSCAN', 'ZBTB', 'ZBED')


# --------------------------------------------------------------------------------------------------------------------
# motifs
# --------------------------------------------------------------------------------------------------------------------

def read_jaspar(path):
    """Position frequency matrices of a JASPAR format file: list of dicts id, name, counts (4 x L, rows A C G T)."""
    motifs = []
    current = None
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            parts = line[1:].split(None, 1)
            current = {'id': parts[0], 'name': parts[1].strip() if len(parts) > 1 else parts[0], 'rows': []}
            motifs.append(current)
        else:
            values = re.sub(r'^[ACGT]\s*', '', line).replace('[', ' ').replace(']', ' ').split()
            current['rows'].append([float(v) for v in values])
    for motif in motifs:
        motif['counts'] = np.array(motif.pop('rows'), dtype=float)
        if motif['counts'].shape[0] != 4:
            raise ValueError('motif {} has {} rows'.format(motif['id'], motif['counts'].shape[0]))
    return motifs


def motif_genes(name):
    """Gene symbols of a JASPAR motif name, upper case: dimers 'Pax3::Foxo1' -> [PAX3, FOXO1]; variants
    'NR2F1(var.2)' -> [NR2F1]."""
    name = re.sub(r'\(var\.\d+\)', '', name)
    return [part.strip().upper() for part in name.split('::') if part.strip()]


def tf_family(symbol):
    """Family key of a gene symbol: the letters before the first digit (HOXA2 -> HOX, SOX10 -> SOX, POU5F1 -> POU,
    LHX2 -> LHX, RUNX1 -> RUNX), for symbols without a digit the symbol without its last letter (CEBPA -> CEBP, JUNB
    -> JUN); Hox clusters are merged (HOXA, HOXB, ... -> HOX). None for short keys and for large zinc finger
    groups. A heuristic for reporting related factors, not a structural classification."""
    symbol = symbol.upper()
    match = re.match(r'([A-Z]+)\d', symbol)
    if match:
        key = match.group(1)
    elif len(symbol) > 3 and symbol[-1].isalpha():
        key = symbol[:-1]
    else:
        key = symbol
    if re.fullmatch(r'HOX[A-D]', key):
        key = 'HOX'
    if len(key) < 3 or key.startswith(NO_FAMILY):
        return None
    return key


def nucleotide_background(sequences, pseudo=1.0):
    """A C G T frequencies of a set of sequences (both strands, so A = T and C = G)."""
    counts = np.full(4, pseudo)
    for seq in sequences:
        seq = seq.upper()
        for i, base in enumerate(BASES):
            counts[i] += seq.count(base)
    counts = counts + counts[::-1]
    return list(counts / counts.sum())


def gc_content(seq):
    seq = seq.upper()
    acgt = sum(seq.count(b) for b in BASES)
    return (seq.count('G') + seq.count('C')) / acgt if acgt else np.nan


def _log_odds(motifs, background, pseudocount):
    import MOODS.tools
    forward = [MOODS.tools.log_odds(m['counts'].tolist(), background, pseudocount) for m in motifs]
    return forward, [MOODS.tools.reverse_complement(m) for m in forward]


def motif_thresholds(motifs, background, pvalue=1e-4, pseudocount=1.0):
    """Log-odds score threshold of every motif at P < pvalue under the background (MOODS) and its maximum score."""
    import MOODS.tools
    forward, _ = _log_odds(motifs, background, pseudocount)
    thresholds = np.array([MOODS.tools.threshold_from_p(m, background, pvalue) for m in forward])
    maxima = np.array([MOODS.tools.max_score(m) for m in forward])
    return thresholds, maxima


_SCANNER = {}


def _init_scanner(motifs, background, pvalue, pseudocount):
    import MOODS.scan
    import MOODS.tools
    forward, reverse = _log_odds(motifs, background, pseudocount)
    matrices = forward + reverse
    thresholds = [MOODS.tools.threshold_from_p(m, background, pvalue) for m in matrices]
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, background, thresholds)
    _SCANNER['scanner'] = scanner
    _SCANNER['n'] = len(forward)


def _scan_chunk(chunk):
    scanner, n = _SCANNER['scanner'], _SCANNER['n']
    counts = np.zeros((len(chunk), n), dtype=np.uint16)
    best = np.full((len(chunk), n), np.nan, dtype=np.float32)
    for row, seq in enumerate(chunk):
        if not seq:
            continue
        results = scanner.scan(seq.upper())
        for j in range(n):
            hits = results[j] + results[j + n]
            if hits:
                counts[row, j] = min(len(hits), 65535)
                best[row, j] = max(h.score for h in hits)
    return counts, best


def scan_sequences(sequences, motifs, background, pvalue=1e-4, pseudocount=1.0, threads=1, chunk=500):
    """Scan sequences with every motif on both strands (MOODS, log-odds threshold at P < pvalue under the given
    background). Returns hit counts (uint16, sequences x motifs) and the best hit score (NaN without a hit)."""
    sequences = list(sequences)
    chunks = [sequences[i:i + chunk] for i in range(0, len(sequences), chunk)] or [[]]
    if threads <= 1 or len(chunks) == 1:
        _init_scanner(motifs, background, pvalue, pseudocount)
        parts = [_scan_chunk(c) for c in chunks]
    else:
        ctx = mp.get_context('fork')
        with ctx.Pool(threads, initializer=_init_scanner, initargs=(motifs, background, pvalue, pseudocount)) as pool:
            parts = pool.map(_scan_chunk, chunks)
    if not sequences:
        return np.zeros((0, len(motifs)), dtype=np.uint16), np.zeros((0, len(motifs)), dtype=np.float32)
    return np.vstack([p[0] for p in parts]), np.vstack([p[1] for p in parts])


# --------------------------------------------------------------------------------------------------------------------
# windows and sequences
# --------------------------------------------------------------------------------------------------------------------

def best_window(profile, bin_size, window):
    """Start offset (bp) of the window of `window` bp with the largest summed binned signal; ties and empty
    profiles give the central window."""
    n_bins = len(profile)
    w = max(1, int(round(window / bin_size)))
    if n_bins <= w:
        return 0
    signal = np.nan_to_num(np.asarray(profile, dtype=float))
    sums = np.convolve(signal, np.ones(w), mode='valid')
    centre = (len(sums) - 1) / 2.0
    top = sums.max()
    if top <= 0:
        j = int(round(centre))
    else:
        candidates = np.flatnonzero(sums >= top - 1e-12)
        j = int(candidates[np.argmin(np.abs(candidates - centre))])
    return j * bin_size


def _profile_chunk(task):
    import pyBigWig
    tracks, rows, bin_size = task
    handles = [pyBigWig.open(t) for t in tracks]
    chroms = [set(h.chroms()) for h in handles]
    out = []
    for chrom, start, end in rows:
        n_bins = max(1, int(math.ceil((end - start) / bin_size)))
        profile = np.zeros((len(handles), n_bins))
        for i, h in enumerate(handles):
            if chrom not in chroms[i]:
                continue
            stop = min(end, h.chroms()[chrom])
            if stop <= start:
                continue
            bins = max(1, int(math.ceil((stop - start) / bin_size)))
            values = h.stats(chrom, start, stop, type='mean', nBins=bins)
            profile[i, :bins] = np.array([v if v is not None else 0.0 for v in values])
        out.append(profile)
    for h in handles:
        h.close()
    return out


def node_windows(nodes, tracks, window=1000, bin_size=50, threads=1, chunk=2000):
    """Per node (DataFrame with node, chrom, start, end): the `window` bp sub-window with the highest active-mark
    signal (columns win_start, win_end, win_signal). Every track is divided by its mean over all node bins before
    summing, so that tracks of different scale count alike. Without tracks: the central window."""
    nodes = nodes.reset_index(drop=True)
    rows = list(zip(nodes['chrom'], nodes['start'].astype(int), nodes['end'].astype(int)))
    if tracks:
        tasks = [(list(tracks), rows[i:i + chunk], bin_size) for i in range(0, len(rows), chunk)]
        if threads > 1 and len(tasks) > 1:
            with mp.get_context('fork').Pool(threads) as pool:
                parts = pool.map(_profile_chunk, tasks)
        else:
            parts = [_profile_chunk(t) for t in tasks]
        profiles = [p for part in parts for p in part]
        totals = np.zeros(len(tracks))
        n_bins = 0
        for p in profiles:
            totals += p.sum(axis=1)
            n_bins += p.shape[1]
        scale = np.where(totals > 0, totals / max(n_bins, 1), 1.0)
        summed = [(p / scale[:, None]).sum(axis=0) for p in profiles]
    else:
        summed = [np.zeros(max(1, int(math.ceil((e - s) / bin_size)))) for _, s, e in rows]
    win_start, win_end, win_signal = [], [], []
    for (chrom, start, end), profile in zip(rows, summed):
        length = end - start
        if length <= window:
            ws, we = start, end
            sig = float(profile.sum())
        elif not np.any(profile > 0):
            ws = start + (length - window) // 2
            we = ws + window
            sig = 0.0
        else:
            offset = best_window(profile, bin_size, window)
            ws = min(start + offset, end - window)
            we = ws + window
            w = max(1, int(round(window / bin_size)))
            sig = float(profile[offset // bin_size: offset // bin_size + w].sum())
        win_start.append(int(ws))
        win_end.append(int(we))
        win_signal.append(sig / max(1, int(round(min(length, window) / bin_size))))
    out = nodes.copy()
    out['win_start'] = win_start
    out['win_end'] = win_end
    out['win_signal'] = win_signal
    return out


def fetch_sequences(fasta, windows):
    """Sequences of the windows (pyfaidx Fasta or a dict chrom -> str); empty string for unknown chromosomes."""
    seqs = []
    for chrom, start, end in zip(windows['chrom'], windows['win_start'], windows['win_end']):
        if chrom not in fasta:
            seqs.append('')
            continue
        record = fasta[chrom]
        seq = record[int(start):int(end)]
        seqs.append(str(seq if isinstance(seq, str) else seq.seq).upper())
    return seqs


# --------------------------------------------------------------------------------------------------------------------
# explaining nodes and gene groups
# --------------------------------------------------------------------------------------------------------------------

def graph_marks(graph, names):
    """Marks with a signal of both conditions as node attributes <mark>_<A>, <mark>_<B>."""
    keys = set()
    for _, data in graph.nodes(data=True):
        keys |= set(data)
    a, b = names
    return sorted({k[:-len(a) - 1] for k in keys if k.endswith('_' + a)} &
                  {k[:-len(b) - 1] for k in keys if k.endswith('_' + b)})


def _parse_signs(text):
    signs = {}
    for item in str(text).split(';'):
        if ':' in item:
            mark, value = item.split(':')
            signs[mark] = int(value)
    return signs


def explaining_nodes(graph, levels, names, max_hops=3, chromatin_threshold=1.0, roles='declared', presence='mean',
                     pseudo_count=1.0):
    """Every node of every gene's levels with a flag whether it explains the gene.

    `levels` is <prefix>_explanations_all_levels.tsv of egaDifferential. Rows: gene_node, genes, rna_log2FC, level
    (promoter, step1 .. step<max_hops>), node, explaining (bool), marks (the concordant explaining marks of the node).
    A promoter row is explaining when promoter_marks or promoter_state explains the gene, or when the promoter is
    poised in both conditions. A step-k node is explaining when the level explains and the node's own change of an
    explaining mark is concordant and >= `chromatin_threshold`; non-explaining nodes at the same step are kept
    (explaining False) as the within-gene background.
    """
    from .differential import ActivityOperators
    marks = graph_marks(graph, names)
    a, b = names
    genes = list(dict.fromkeys(levels['node']))
    operators = ActivityOperators(graph, genes, max_hops=max_hops)
    own = {}
    for m in marks:
        va = np.array([float(graph.nodes[n].get('{}_{}'.format(m, a), 0.0) or 0.0) for n in operators.nodes])
        vb = np.array([float(graph.nodes[n].get('{}_{}'.format(m, b), 0.0) or 0.0) for n in operators.nodes])
        own[m] = np.log2((vb + pseudo_count) / (va + pseudo_count))
    row_of = {n: i for i, n in enumerate(operators.sources)}
    explains = 'explains_{}_{}'.format(roles, presence)
    concordant_col = 'concordant_marks_' + roles
    signs_col = 'signs_' + roles
    state_col = 'state_fragment_{}_{}'.format(roles, presence)
    result_col = 'result_{}_{}'.format(roles, presence)
    by_gene = {node: part.set_index('level') for node, part in levels.groupby('node', sort=False)}
    records = []
    for node in genes:
        part = by_gene[node]
        first = part.iloc[0]
        d = np.sign(float(first['rna_log2FC']))
        base = (node, first['genes'], float(first['rna_log2FC']))
        promoter_explains = any(bool(part.at[lv, explains]) for lv in ('promoter_marks', 'promoter_state')
                                if lv in part.index)
        poised_both = 'promoter_state' in part.index and part.at['promoter_state', result_col] == 'poised in both'
        promoter_marks = ''
        if 'promoter_marks' in part.index and bool(part.at['promoter_marks', explains]):
            promoter_marks = str(part.at['promoter_marks', concordant_col])
        records.append(base + ('promoter', node, bool(promoter_explains or poised_both), promoter_marks))
        if node not in row_of:
            continue
        k_row = row_of[node]
        for k in range(1, max_hops + 1):
            level = 'step1_marks' if k == 1 else 'step{}'.format(k)
            H_a, H_b = operators.hop_a[k], operators.hop_b[k]
            columns = np.union1d(H_a.getrow(k_row).indices, H_b.getrow(k_row).indices)
            if not len(columns):
                continue
            hit = np.zeros(len(columns), dtype=bool)
            hit_marks = [[] for _ in columns]
            if level in part.index and bool(part.at[level, explains]):
                signs = _parse_signs(part.at[level, signs_col])
                used = [m for m in str(part.at[level, concordant_col]).split(',') if m in own]
                for m in used:
                    change = own[m][columns]
                    ok = (np.sign(change) * signs.get(m, 1) == d) & (np.abs(change) >= chromatin_threshold)
                    hit |= ok
                    for j in np.flatnonzero(ok):
                        hit_marks[j].append(m)
            extra = None
            if k == 1 and 'step1_state' in part.index and bool(part.at['step1_state', explains]):
                extra = part.at['step1_state', state_col]
            for j, c in enumerate(columns):
                n = operators.nodes[c]
                flag = bool(hit[j]) or n == extra
                mk = ','.join(hit_marks[j]) + (('state' if not hit_marks[j] else ',state') if n == extra else '')
                records.append(base + ('step{}'.format(k), n, flag, mk))
    return pd.DataFrame.from_records(records, columns=['gene_node', 'genes', 'rna_log2FC', 'level', 'node',
                                                       'explaining', 'marks'])


def gene_groups(classes, levels, nodes, max_hops=3, roles='declared', presence='mean', min_genes=5):
    """Groups of genes with their foreground (explaining) nodes and within-gene background nodes.

    Per direction (up / down): promoter_<dir> (promoter marks explain; promoter fragments), step<k>_<dir> (marks at
    step k explain; explaining nodes at step k against the other nodes at step k of the same genes), step1_state_<dir>
    (the poised contacting fragment of an explaining step-1 state transition), poised_both_<dir> (promoter bivalent
    in both conditions: promoter fragments; poised_both_<dir>_steps: their explaining nodes at steps 1..k against
    the other nodes at those steps) and one group per promoter state transition involving a poised state
    (promoter_<from>_to_<to>, e.g. promoter_bivalent_to_active). Promoter groups use the promoter fragments of the
    other differential genes as within-set background. Returns a dict name -> dict(genes, fg, bg_same, kind,
    description); groups with fewer than `min_genes` genes are dropped.
    """
    explains = 'explains_{}_{}'.format(roles, presence)
    direction = dict(zip(classes['node'], np.where(classes['rna_log2FC'] > 0, 'up', 'down')))
    all_promoters = set(classes['node'])
    flag = {}
    for level, part in levels.groupby('level'):
        flag[level] = dict(zip(part['node'], part[explains].astype(bool)))
    by_level = {lv: part for lv, part in nodes.groupby('level')}
    groups = {}

    def add(name, genes, kind, description, steps=None):
        genes = sorted(set(genes))
        if len(genes) < min_genes:
            return
        if kind == 'promoter':
            fg = set(genes)
            bg = all_promoters - fg
        else:
            part = pd.concat([by_level[s] for s in steps if s in by_level])
            part = part[part['gene_node'].isin(genes)]
            fg = set(part.loc[part['explaining'], 'node'])
            bg = set(part.loc[~part['explaining'], 'node']) - fg - set(genes)
        if not fg:
            return
        groups[name] = dict(genes=genes, fg=fg, bg_same=bg, kind=kind, description=description)

    for d in ('up', 'down'):
        in_dir = [g for g in classes['node'] if direction[g] == d]
        add('promoter_' + d, [g for g in in_dir if flag.get('promoter_marks', {}).get(g)], 'promoter',
            'genes {} explained by promoter marks; promoter fragments'.format(d))
        for k in range(1, max_hops + 1):
            level = 'step1_marks' if k == 1 else 'step{}'.format(k)
            add('step{}_{}'.format(k, d), [g for g in in_dir if flag.get(level, {}).get(g)], 'step',
                'genes {} explained by marks at step {}; explaining nodes at step {}'.format(d, k, k),
                steps=['step{}'.format(k)])
        state_genes = [g for g in in_dir if flag.get('step1_state', {}).get(g)]
        if state_genes:
            part = by_level.get('step1')
            if part is not None:
                part = part[part['gene_node'].isin(state_genes) & part['marks'].str.contains('state', na=False)]
                fg = set(part['node'])
                if len(set(part['gene_node'])) >= min_genes and fg:
                    bg = set(by_level['step1'].loc[by_level['step1']['gene_node'].isin(state_genes), 'node']) - fg
                    groups['step1_state_' + d] = dict(
                        genes=sorted(set(part['gene_node'])), fg=fg, bg_same=bg, kind='step',
                        description='genes {} explained by a state transition of a contacting fragment; that '
                                    'fragment'.format(d))
        if 'promoter_transition' in classes:
            poised = [g for g, t in zip(classes['node'], classes['promoter_transition'])
                      if direction[g] == d and t == 'bivalent -> bivalent']
            add('poised_both_' + d, poised, 'promoter',
                'genes {} with a promoter bivalent in both conditions; promoter fragments'.format(d))
            add('poised_both_{}_steps'.format(d), poised, 'step',
                'genes {} with a promoter bivalent in both conditions; their explaining nodes at steps 1-{}'.format(
                    d, max_hops), steps=['step{}'.format(k) for k in range(1, max_hops + 1)])
    if 'promoter_transition' in classes:
        for transition, part in classes.groupby('promoter_transition'):
            match = re.match(r'([\w ]+?) -> ([\w ]+?)(:|$)', str(transition))
            if not match:
                continue
            src, dst = match.group(1).strip(), match.group(2).strip()
            if src == dst or not ({src, dst} & {'bivalent', 'poised enhancer'}):
                continue
            name = 'promoter_{}_to_{}'.format(src.replace(' ', '_'), dst.replace(' ', '_'))
            add(name, list(part['node']), 'promoter',
                'promoter {} -> {} ({}); promoter fragments'.format(src, dst, transition))
    return groups


# --------------------------------------------------------------------------------------------------------------------
# enrichment
# --------------------------------------------------------------------------------------------------------------------

def strata(gc, length, gc_width=0.05, length_breaks=(500, 1000)):
    """Stratum label per window: GC bin of `gc_width` and length bin."""
    gc = np.nan_to_num(np.asarray(gc, dtype=float), nan=-1.0)
    gc_bin = np.floor(gc / gc_width).astype(int)
    length_bin = np.searchsorted(np.asarray(length_breaks), np.asarray(length), side='right')
    return gc_bin * 10 + length_bin


def matched_background(fg, candidates, stratum, ratio=3, rng=None):
    """Indices from `candidates` matched to the stratum distribution of `fg` (GC and length): per stratum up to
    `ratio` background windows per foreground window, sampled without replacement (all when fewer exist)."""
    rng = rng or np.random.default_rng(0)
    fg = np.asarray(fg, dtype=int)
    candidates = np.setdiff1d(np.asarray(candidates, dtype=int), fg)
    need = pd.Series(stratum[fg]).value_counts()
    pool = defaultdict(list)
    for i in candidates:
        pool[stratum[i]].append(i)
    chosen = []
    for s, n in need.items():
        available = pool.get(s, [])
        k = min(len(available), int(ratio * n))
        if k:
            chosen.extend(rng.choice(available, size=k, replace=False))
    return np.array(sorted(chosen), dtype=int)


def benjamini_hochberg(p):
    """Benjamini-Hochberg q values; NaN entries stay NaN and are not counted."""
    p = np.asarray(p, dtype=float)
    q = np.full(len(p), np.nan)
    ok = np.flatnonzero(~np.isnan(p))
    if len(ok):
        values = p[ok]
        order = np.argsort(values)
        ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        out = np.empty(len(values))
        out[order] = np.minimum(ranked, 1.0)
        q[ok] = out
    return q


def motif_enrichment(present, fg, bg):
    """Per motif: windows with at least one hit among foreground and background, odds ratio, one-sided Fisher P
    (enrichment) and BH q across motifs. `present` is windows x motifs boolean."""
    present = np.asarray(present, dtype=bool)
    a = present[fg].sum(axis=0)
    c = present[bg].sum(axis=0)
    n_fg, n_bg = len(fg), len(bg)
    p = np.ones(present.shape[1])
    odds = np.full(present.shape[1], np.nan)
    for j in range(present.shape[1]):
        table = [[a[j], n_fg - a[j]], [c[j], n_bg - c[j]]]
        if n_fg and n_bg:
            p[j] = stats.fisher_exact(table, alternative='greater')[1]
        odds[j] = ((a[j] + 0.5) * (n_bg - c[j] + 0.5)) / ((n_fg - a[j] + 0.5) * (c[j] + 0.5))
    return pd.DataFrame({'fg_hits': a, 'fg_windows': n_fg, 'bg_hits': c, 'bg_windows': n_bg,
                         'fg_fraction': a / max(n_fg, 1), 'bg_fraction': c / max(n_bg, 1),
                         'log2_odds': np.log2(odds), 'p': p, 'q': benjamini_hochberg(p)})


# --------------------------------------------------------------------------------------------------------------------
# candidate triggers per gene
# --------------------------------------------------------------------------------------------------------------------

def expression_lookup(table, gene_column='gene_name', log2fc_column='log2FoldChange', padj_column='padj',
                      base_column=None):
    """Upper-case gene symbol -> (log2FC, padj, base expression); duplicated symbols keep the smallest padj."""
    table = table.dropna(subset=[gene_column]).copy()
    table['_key'] = table[gene_column].astype(str).str.upper()
    table['_p'] = pd.to_numeric(table[padj_column], errors='coerce').fillna(1.0)
    table = table.sort_values('_p').drop_duplicates('_key')
    base = table[base_column] if base_column and base_column in table else pd.Series(np.nan, index=table.index)
    return {k: (float(l), float(p), float(b)) for k, l, p, b in
            zip(table['_key'], pd.to_numeric(table[log2fc_column], errors='coerce'), table['_p'], base)}


def gene_symbol_nodes(graph):
    """Upper-case gene symbol -> set of graph nodes annotated with it."""
    out = defaultdict(set)
    for node, data in graph.nodes(data=True):
        for g in re.split(r'[,;]', str(data.get('genes', '') or '')):
            if g:
                out[g.strip().upper()].add(node)
    return out


def neighbourhood_distance(graph, source, targets, max_hops=3):
    """Smallest number of graph steps from `source` to any node in `targets` (None when farther or absent)."""
    if source not in graph or not targets:
        return None
    if source in targets:
        return 0
    lengths = nx.single_source_shortest_path_length(graph, source, cutoff=max_hops)
    found = [lengths[t] for t in targets if t in lengths]
    return min(found) if found else None


def relation(gene_lfc, tf_lfc):
    """'activator-consistent' when the factor changes with the gene, 'repressor-consistent' when against it."""
    if np.sign(gene_lfc) == np.sign(tf_lfc):
        return 'activator-consistent'
    return 'repressor-consistent'


def rank_candidates(df):
    """Rank candidates of one gene by the mean of the ranks of: number of the gene's scanned nodes carrying the
    motif (more first), best relative motif score, absolute expression change of the factor, graph distance from
    the gene to the factor's own promoter (closer first, absent last) and motif enrichment q of the gene's groups
    (smaller first); candidates of the named factor before family members."""
    if df.empty:
        return df
    hops = df['tf_hops'].astype(float).fillna(99)
    q = df['group_q'].astype(float).fillna(1.0)
    ranks = pd.concat([df['n_nodes_with_motif'].rank(ascending=False), df['best_rel_score'].rank(ascending=False),
                       df['tf_log2FC'].abs().rank(ascending=False), hops.rank(), q.rank()], axis=1)
    df = df.assign(rank_score=ranks.mean(axis=1), _family=(df['evidence'] != 'named factor').astype(int))
    df = df.sort_values(['_family', 'rank_score', 'tf_padj']).drop(columns='_family')
    df['rank'] = np.arange(1, len(df) + 1)
    return df


def gene_triggers(classes, nodes, window_index, counts, best, motifs, maxima, expression, graph, alpha=0.05,
                  group_q=None, symbol_nodes=None, max_hops=3, top=25, family=True):
    """Candidate trigger factors per explained or poised gene.

    Scanned set of a gene: its promoter fragment and its explaining nodes at every step. For every motif with a hit
    there and every factor named by the motif (both partners of a dimer): the factor's expression change; kept when
    the factor changes (padj < alpha; evidence 'named factor'), and with `family`, the changed members of the
    factor's family (evidence 'family of <factor>'). Columns report the motif evidence (nodes with the motif,
    promoter hit, steps with a hit, total hits, best score relative to the motif maximum), the factor (log2FC,
    padj, base expression, relation to the gene: activator- or repressor-consistent), the graph distance from the
    gene's promoter to a node carrying the factor's gene and the smallest enrichment q of the motif in a group of
    the gene. At most `top` candidates per gene (n_candidates gives the total).
    """
    if symbol_nodes is None:
        symbol_nodes = gene_symbol_nodes(graph)
    group_q = group_q or {}
    names = [m['name'] for m in motifs]
    ids = [m['id'] for m in motifs]
    motif_tfs = [motif_genes(n) for n in names]
    families = defaultdict(list)
    if family:
        for symbol, (lfc, padj, _) in expression.items():
            key = tf_family(symbol)
            if key and padj < alpha and not np.isnan(lfc):
                families[key].append(symbol)
    explained = classes.copy()
    keep = nodes[nodes['explaining'] | (nodes['level'] == 'promoter')]
    by_gene = {g: part for g, part in keep.groupby('gene_node', sort=False)}
    poised_col = 'promoter_transition' if 'promoter_transition' in explained else None
    rows = []
    for _, gene in explained.iterrows():
        g = gene['node']
        part = by_gene.get(g)
        if part is None:
            continue
        is_explaining = bool(part['explaining'].any())
        poised = poised_col is not None and str(gene[poised_col]) == 'bivalent -> bivalent'
        if not (is_explaining or poised):
            continue
        part = part[part['explaining'] | (part['level'] == 'promoter')]
        node_list = list(dict.fromkeys(part['node']))
        idx = [window_index[n] for n in node_list if n in window_index]
        if not idx:
            continue
        levels_of = defaultdict(set)
        for n, lv in zip(part['node'], part['level']):
            levels_of[n].add(lv)
        idx_nodes = [n for n in node_list if n in window_index]
        sub = counts[idx]
        sub_best = best[idx]
        promoter_pos = idx_nodes.index(g) if g in idx_nodes else None
        gene_lfc = float(gene['rna_log2FC'])
        hops_cache = {}
        lengths = nx.single_source_shortest_path_length(graph, g, cutoff=max_hops) if g in graph else {}
        candidates = []
        for j in np.flatnonzero(sub.sum(axis=0) > 0):
            with_motif = np.flatnonzero(sub[:, j] > 0)
            steps = sorted({lv for i in with_motif for lv in levels_of[idx_nodes[i]]},
                           key=lambda s: (s != 'promoter', s))
            best_score = float(np.nanmax(sub_best[with_motif, j]))
            evidence = dict(motif_id=ids[j], motif=names[j], n_nodes_scanned=len(idx), n_nodes_with_motif=len(with_motif),
                            promoter_has_motif=bool(promoter_pos is not None and sub[promoter_pos, j] > 0),
                            steps_with_motif=';'.join(steps), hits=int(sub[:, j].sum()), best_score=best_score,
                            best_rel_score=best_score / maxima[j] if maxima[j] > 0 else np.nan,
                            group_q=group_q.get((g, j), np.nan))
            seen = set()
            for tf in motif_tfs[j]:
                members = []
                if tf in expression and expression[tf][1] < alpha and not np.isnan(expression[tf][0]):
                    members.append((tf, 'named factor'))
                key = tf_family(tf) if family else None
                if key:
                    members += [(s, 'family of ' + tf) for s in families.get(key, []) if s != tf]
                for symbol, ev in members:
                    if symbol in seen:
                        continue
                    seen.add(symbol)
                    lfc, padj, base = expression[symbol]
                    if symbol not in hops_cache:
                        found = [lengths[t] for t in symbol_nodes.get(symbol, ()) if t in lengths]
                        hops_cache[symbol] = min(found) if found else None
                    candidates.append(dict(evidence, tf=symbol, evidence=ev, tf_log2FC=lfc, tf_padj=padj,
                                           tf_expression=base, relation=relation(gene_lfc, lfc),
                                           tf_hops=hops_cache[symbol], tf_in_graph=symbol in symbol_nodes))
        if not candidates:
            continue
        table = rank_candidates(pd.DataFrame(candidates))
        n_total = len(table)
        table = table.head(top)
        head = dict(gene=gene.get('genes', ''), gene_node=g, gene_log2FC=gene_lfc,
                    gene_direction='up' if gene_lfc > 0 else 'down',
                    explained_at=gene.get('explained_at', ''), explained_type=gene.get('explained_type', ''),
                    promoter_transition=gene.get(poised_col, '') if poised_col else '',
                    poised_both=poised, n_candidates=n_total)
        for record in table.to_dict('records'):
            rows.append(dict(head, **record))
    columns = ['gene', 'gene_node', 'gene_log2FC', 'gene_direction', 'explained_at', 'explained_type',
               'promoter_transition', 'poised_both', 'rank', 'n_candidates', 'tf', 'evidence', 'relation',
               'tf_log2FC', 'tf_padj', 'tf_expression', 'tf_in_graph', 'tf_hops', 'motif_id', 'motif',
               'n_nodes_with_motif', 'n_nodes_scanned', 'promoter_has_motif', 'steps_with_motif', 'hits',
               'best_score', 'best_rel_score', 'group_q', 'rank_score']
    return pd.DataFrame(rows, columns=columns)
