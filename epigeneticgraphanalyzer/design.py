"""Experiment design from the explanations of egaDifferential: in-silico perturbation of the explaining nodes.

For every explained gene the nodes of its explaining levels are perturbed one at a time and the level statistics are
recomputed with the formulas of egaDifferential (differential.gene_features, explanation_levels, gene_level_hop_test):

  mark levels (promoter, step 1 .. step k): value_m = log2((H_B^(k) v_B,m + 1) / (H_A^(k) v_A,m + 1)), the
      contact-weighted k-hop change of mark m (H^(k) the path weights of the shortest paths of k steps, the edge
      strength r of a condition = score / reference score); statistic S = max_m sign_m * sign(RNA) * value_m over the
      marks with a sign at that level; the level explains the gene when S >= chromatin threshold.
  own contacts: S = sign_contact * sign(RNA) * log2((sum r_B + 1) / (sum r_A + 1)) over the promoter's partners.
  step-1 state: the contacting fragment with the poised transition; the level is lost when it is the only fragment
      with that transition.

Perturbations of a node j:
  flatten  the mark change of j set to 0 (both conditions carry the mean of its two signals; path weights kept),
  silence  the marks of j set to 0 in both conditions (CRISPRi-like loss of the element's marks; paths through j kept),
  remove   j deleted with its edges (deletion; shortest paths and the nodes at every distance recomputed).
The per-gene permutation P of a step level (the null of gene_level_hop_test: mark signals of all non-source nodes
permuted jointly) is recomputed with a reduced number of permutations; the observed statistic before and after a
perturbation is compared with the same null distribution.

Candidates are ranked per gene by the predicted effect and the node's own evidence (design_score); the ranking is a
report, not a decision. A single-condition mode (single_condition_contributions) ranks the nodes of a gene's
neighbourhood by their contribution to its contact-weighted active-mark activity, for cell types without a
differential run.
"""

import logging
import multiprocessing
import re

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .differential import LevelSigns, level_name

log = logging.getLogger(__name__)

PSEUDO_COUNT = 1.0
PERTURBATIONS = ('flatten', 'silence', 'remove')

# --------------------------------------------------------------------------------------------------------------------
# levels and signs
# --------------------------------------------------------------------------------------------------------------------


def level_hops_of(level):
    """Steps of an explanation level name: promoter_* 0, step1_marks / step1_state 1, step<k> k."""
    if level.startswith('promoter'):
        return 0
    return int(re.match(r'step(\d+)', level).group(1))


def feature_level(level):
    """Feature level of a mark level: promoter_marks 'promoter', step1_marks 'direct', step2 'indirect', step<k>
    'hop<k>'."""
    hops = level_hops_of(level)
    return 'promoter' if hops == 0 else level_name(hops)


def mark_level_of(hops):
    return 'promoter_marks' if hops == 0 else 'step1_marks' if hops == 1 else 'step{}'.format(hops)


def parse_signs(text):
    signs = {}
    for item in str(text).split(';'):
        if ':' in item:
            mark, value = item.split(':')
            signs[mark.strip()] = int(value)
    return signs


def level_signs(levels, roles='declared'):
    """LevelSigns of a role set from the signs_<roles> column of an explanations_all_levels table (the signs are the
    same for every gene of a level); the own contacts are the mark 'contact' at the promoter level."""
    out = LevelSigns()
    column = 'signs_' + roles
    for level, part in levels.groupby('level'):
        values = part[column].dropna()
        if values.empty:
            continue
        signs = parse_signs(values.iloc[0])
        if level == 'promoter_contacts':
            out.setdefault('promoter', {})['contact'] = signs.get('contact', 0)
        elif level == 'promoter_marks' or level == 'step1_marks' or re.fullmatch(r'step\d+', level):
            out.setdefault(feature_level(level), {}).update(signs)
    return out


# --------------------------------------------------------------------------------------------------------------------
# graph arrays and neighbourhoods
# --------------------------------------------------------------------------------------------------------------------


class GraphArrays:
    """Index, symmetric adjacency and edge strength matrices of both conditions and the mark vectors of a graph.

    `weights` names the edge attributes of conditions A and B (union graph: r_a, r_b); a single-condition graph uses
    the same attribute (or score / reference score) for both. `marks` maps mark -> (attribute A, attribute B)."""

    def __init__(self, graph, marks, sources=(), weights=('r_a', 'r_b'), reference_score=None):
        self.nodes = list(graph.nodes)
        self.index = {n: i for i, n in enumerate(self.nodes)}
        n = len(self.nodes)
        rows, cols, w_a, w_b = [], [], [], []
        for u, v, data in graph.edges(data=True):
            i, j = self.index[u], self.index[v]
            if i == j:
                continue
            a = float(data.get(weights[0], 0.0) or 0.0)
            b = float(data.get(weights[1], 0.0) or 0.0)
            if reference_score:
                a, b = a / reference_score, b / reference_score
            rows += [i, j]
            cols += [j, i]
            w_a += [a, a]
            w_b += [b, b]
        self.W_a = sp.csr_matrix((w_a, (rows, cols)), shape=(n, n))
        self.W_b = sp.csr_matrix((w_b, (rows, cols)), shape=(n, n))
        self.A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
        self.A.data[:] = 1.0
        self.marks = list(marks)
        self.vectors = {}
        for m, (key_a, key_b) in marks.items():
            va = np.array([float(graph.nodes[x].get(key_a, 0.0) or 0.0) for x in self.nodes])
            vb = np.array([float(graph.nodes[x].get(key_b, 0.0) or 0.0) for x in self.nodes])
            self.vectors[m] = (np.nan_to_num(va), np.nan_to_num(vb))
        self.chrom = np.array([str(graph.nodes[x].get('chrom', x.split(':')[0])) for x in self.nodes], dtype=object)
        self.start = np.array([int(graph.nodes[x].get('start', 0)) for x in self.nodes], dtype=np.int64)
        self.end = np.array([int(graph.nodes[x].get('end', 0)) for x in self.nodes], dtype=np.int64)
        self.genes = np.array([str(graph.nodes[x].get('genes', '') or '') for x in self.nodes], dtype=object)
        self.is_source = np.zeros(n, dtype=bool)
        for s in sources:
            if s in self.index:
                self.is_source[self.index[s]] = True


def neighbourhood(arrays, source, max_hops, removed=None):
    """Path weights of `source` to the nodes at exactly k steps, k = 1 .. max_hops: {k: (node indices, w_A, w_B)}.

    As differential.ActivityOperators: H^(1) = W[source]; H^(k)[m] = sum over the nodes p at k - 1 steps of
    H^(k-1)[p] W[p, m] for the nodes m at exactly k steps (shortest-path distance in the union graph). `removed`
    (node index or set of indices) deletes those nodes with their edges before the distances are computed."""
    removed = set() if removed is None else ({removed} if np.isscalar(removed) else set(removed))
    indptr, indices = arrays.A.indptr, arrays.A.indices
    visited = {source} | removed
    layer = np.array([source], dtype=np.int64)
    h_a = h_b = None
    out = {}
    for k in range(1, max_hops + 1):
        if not len(layer):
            out[k] = (np.zeros(0, dtype=np.int64), np.zeros(0), np.zeros(0))
            continue
        reach = np.unique(np.concatenate([indices[indptr[p]:indptr[p + 1]] for p in layer]))
        new = np.array([m for m in reach if m not in visited], dtype=np.int64)
        if k == 1:
            w_a = np.asarray(arrays.W_a[source, new].todense()).ravel() if len(new) else np.zeros(0)
            w_b = np.asarray(arrays.W_b[source, new].todense()).ravel() if len(new) else np.zeros(0)
        elif len(new):
            sub_a = arrays.W_a[layer][:, new]
            sub_b = arrays.W_b[layer][:, new]
            w_a = np.asarray(sub_a.T @ h_a).ravel()
            w_b = np.asarray(sub_b.T @ h_b).ravel()
        else:
            w_a = w_b = np.zeros(0)
        out[k] = (new, w_a, w_b)
        visited.update(new.tolist())
        layer, h_a, h_b = new, w_a, w_b
    return out


def _log_ratio(x_b, x_a, pseudo_count=PSEUDO_COUNT):
    return np.log2((np.maximum(x_b, 0.0) + pseudo_count) / (np.maximum(x_a, 0.0) + pseudo_count))


def level_values(arrays, hops, source):
    """log2 B/A change of every mark at the promoter (hops is None) or at one step ((cols, w_a, w_b))."""
    out = {}
    for m, (va, vb) in arrays.vectors.items():
        if hops is None:
            out[m] = float(_log_ratio(vb[source], va[source]))
        else:
            cols, w_a, w_b = hops
            out[m] = float(_log_ratio(np.dot(w_b, vb[cols]), np.dot(w_a, va[cols])))
    return out


def statistic(values, signs, direction):
    """max_m sign_m * direction * value_m over the marks with a non-zero sign (NaN without any)."""
    terms = [signs.get(m, 0) * direction * v for m, v in values.items() if signs.get(m, 0)]
    return float(max(terms)) if terms else np.nan


def perturbed_statistics(arrays, hops, signs, direction, mode, pseudo_count=PSEUDO_COUNT):
    """Statistic of a step level after flattening or silencing each of its nodes in turn (vector over the nodes)."""
    cols, w_a, w_b = hops
    best = np.full(len(cols), -np.inf)
    used = False
    for m, (va, vb) in arrays.vectors.items():
        s = signs.get(m, 0)
        if not s:
            continue
        used = True
        a, b = va[cols], vb[cols]
        x_a, x_b = np.dot(w_a, a), np.dot(w_b, b)
        target = (a + b) / 2.0 if mode == 'flatten' else np.zeros(len(cols))
        new_a = x_a - w_a * a + w_a * target
        new_b = x_b - w_b * b + w_b * target
        best = np.maximum(best, s * direction * _log_ratio(new_b, new_a, pseudo_count))
    return best if used else np.full(len(cols), np.nan)


def node_contributions(arrays, hops, signs, direction, marks):
    """Concordant contribution of every node of a step to the explaining marks, (w_B v_B - w_A v_A) sign direction,
    maximised over `marks` (with the mark of the maximum)."""
    cols, w_a, w_b = hops
    best = np.full(len(cols), -np.inf)
    which = np.array([''] * len(cols), dtype=object)
    for m in marks:
        s = signs.get(m, 0)
        if not s or m not in arrays.vectors:
            continue
        va, vb = arrays.vectors[m]
        c = (w_b * vb[cols] - w_a * va[cols]) * s * direction
        better = c > best
        best[better] = c[better]
        which[better] = m
    return best, which


# --------------------------------------------------------------------------------------------------------------------
# per-gene design
# --------------------------------------------------------------------------------------------------------------------

_SHARED = {}


def _map(func, items, threads, shared):
    global _SHARED
    items = list(items)
    _SHARED = shared
    try:
        if threads <= 1 or len(items) <= 1:
            return [func(item) for item in items]
        with multiprocessing.get_context('fork').Pool(min(threads, len(items))) as pool:
            return pool.map(func, items, chunksize=max(1, len(items) // (threads * 8)))
    finally:
        _SHARED = {}


def _state_count(detail):
    match = re.search(r'\((\d+) of [\d,]+ with this transition\)', str(detail))
    return int(match.group(1)) if match else np.nan


def gene_design(arrays, gene, max_hops, signs, threshold=1.0, n_remove=10, explains_col='explains_declared_mean',
                concordant_col='concordant_marks_declared', detail_col='detail_declared_mean',
                state_col='state_fragment_declared_mean', pseudo_count=PSEUDO_COUNT, max_candidates=200):
    """In-silico perturbation of the explaining nodes of one gene.

    `gene` is the gene's part of the explanations_all_levels table indexed by level. Returns (candidates, levels,
    hops): candidate rows (one per node and level), the recomputed statistic of every level and the path weights of
    the explaining step levels (for the permutation null)."""
    first = gene.iloc[0]
    node = first['node']
    if node not in arrays.index:
        return [], [], {}
    source = arrays.index[node]
    d = float(np.sign(first['rna_log2FC']))
    hops = neighbourhood(arrays, source, max_hops)
    explains = {lv: bool(gene.at[lv, explains_col]) if lv in gene.index and not pd.isna(gene.at[lv, explains_col])
                else False for lv in gene.index}
    base = dict(gene_node=node, genes=first['genes'], rna_log2FC=float(first['rna_log2FC']))
    level_rows, candidates, null_hops = [], [], {}

    def stats_of(h):
        """Statistic of every mark level (promoter, steps) of neighbourhood `h`."""
        out = {'promoter_marks': statistic(level_values(arrays, None, source), signs.get('promoter', {}), d)}
        for k in range(1, max_hops + 1):
            out[mark_level_of(k)] = statistic(level_values(arrays, h[k], source), signs.get(level_name(k), {}), d)
        return out

    before = stats_of(hops)
    contact_sign = signs.get('promoter', {}).get('contact', 0)
    cols1, wa1, wb1 = hops[1]
    contact_before = contact_sign * d * float(_log_ratio(wb1.sum(), wa1.sum())) if contact_sign else np.nan
    for lv, value in list(before.items()) + [('promoter_contacts', contact_before)]:
        level_rows.append(dict(base, level=lv, statistic=value, explains_table=explains.get(lv, False),
                               explains_recomputed=bool(value >= threshold) if not pd.isna(value) else False))
    explaining_levels = [lv for lv, e in explains.items() if e]

    def other_levels_explain(level, after=None):
        """Levels other than `level` that still explain (after a removal: the recomputed mark levels)."""
        out = []
        for lv in explaining_levels:
            if lv == level:
                continue
            if after is not None and lv in after:
                if after[lv] >= threshold:
                    out.append(lv)
            else:
                out.append(lv)
        return out

    # own promoter (promoter marks / state): perturbing it flattens the promoter level
    own_levels = [lv for lv in ('promoter_marks', 'promoter_state') if explains.get(lv)]
    if own_levels:
        candidates.append(dict(base, node=node, node_index=source, level='+'.join(own_levels), hops=0,
                               explaining_mark=str(gene.at['promoter_marks', concordant_col])
                               if 'promoter_marks' in own_levels else 'state',
                               contribution=np.nan, path_weight_A=np.nan, path_weight_B=np.nan,
                               statistic_before=before['promoter_marks'], statistic_flatten=0.0,
                               statistic_silence=0.0, statistic_remove=np.nan,
                               level_lost_flatten=True, level_lost_silence=True, level_lost_remove=np.nan,
                               still_explained_flatten=','.join(lv for lv in explaining_levels
                                                                if lv not in own_levels),
                               own_promoter=True))
    # own contacts
    if explains.get('promoter_contacts') and len(cols1) and contact_sign:
        sa, sb = wa1.sum(), wb1.sum()
        mean = (wa1 + wb1) / 2.0
        flat = contact_sign * d * _log_ratio(sb - wb1 + mean, sa - wa1 + mean)
        silence = contact_sign * d * _log_ratio(sb - wb1, sa - wa1)
        contribution = (wb1 - wa1) * contact_sign * d
        for j in np.flatnonzero(contribution > 0):
            candidates.append(dict(base, node=arrays.nodes[cols1[j]], node_index=int(cols1[j]),
                                   level='promoter_contacts', hops=1, explaining_mark='contact',
                                   contribution=float(contribution[j]), path_weight_A=float(wa1[j]),
                                   path_weight_B=float(wb1[j]), statistic_before=contact_before,
                                   statistic_flatten=float(flat[j]), statistic_silence=float(silence[j]),
                                   statistic_remove=float(silence[j]),
                                   level_lost_flatten=bool(flat[j] < threshold),
                                   level_lost_silence=bool(silence[j] < threshold),
                                   level_lost_remove=bool(silence[j] < threshold),
                                   still_explained_flatten=','.join(other_levels_explain('promoter_contacts')),
                                   own_promoter=False))
    # contacting fragment with the poised transition
    if explains.get('step1_state') and state_col in gene.columns:
        fragment = gene.at['step1_state', state_col]
        if isinstance(fragment, str) and fragment in arrays.index:
            count = _state_count(gene.at['step1_state', detail_col])
            j = arrays.index[fragment]
            where = np.flatnonzero(cols1 == j)
            lost = bool(count == 1) if not pd.isna(count) else np.nan
            candidates.append(dict(base, node=fragment, node_index=j, level='step1_state', hops=1,
                                   explaining_mark='state', contribution=np.nan,
                                   path_weight_A=float(wa1[where[0]]) if len(where) else np.nan,
                                   path_weight_B=float(wb1[where[0]]) if len(where) else np.nan,
                                   statistic_before=np.nan, statistic_flatten=np.nan, statistic_silence=np.nan,
                                   statistic_remove=np.nan, level_lost_flatten=lost, level_lost_silence=lost,
                                   level_lost_remove=lost, state_fragments=count,
                                   still_explained_flatten=','.join(other_levels_explain('step1_state')),
                                   own_promoter=False))
    # step levels
    for k in range(1, max_hops + 1):
        lv = mark_level_of(k)
        if not explains.get(lv):
            continue
        h = hops[k]
        if not len(h[0]):
            continue
        null_hops[k] = h
        level_signs_k = signs.get(level_name(k), {})
        used = [m for m in str(gene.at[lv, concordant_col]).split(',') if m in arrays.vectors] or \
            [m for m in arrays.marks if level_signs_k.get(m, 0)]
        contribution, which = node_contributions(arrays, h, level_signs_k, d, used)
        keep = np.flatnonzero(contribution > 0)
        if not len(keep):
            continue
        flat = perturbed_statistics(arrays, h, level_signs_k, d, 'flatten', pseudo_count)
        silence = perturbed_statistics(arrays, h, level_signs_k, d, 'silence', pseudo_count)
        order = keep[np.argsort(-(before[lv] - flat[keep]))][:max_candidates]
        for j in order:
            candidates.append(dict(base, node=arrays.nodes[h[0][j]], node_index=int(h[0][j]), level=lv, hops=k,
                                   explaining_mark=which[j], contribution=float(contribution[j]),
                                   path_weight_A=float(h[1][j]), path_weight_B=float(h[2][j]),
                                   statistic_before=before[lv], statistic_flatten=float(flat[j]),
                                   statistic_silence=float(silence[j]), statistic_remove=np.nan,
                                   level_lost_flatten=bool(flat[j] < threshold),
                                   level_lost_silence=bool(silence[j] < threshold), level_lost_remove=np.nan,
                                   still_explained_flatten=','.join(other_levels_explain(lv)),
                                   own_promoter=False))
    # removal (node and its edges deleted) of the candidates with the largest flatten / silence effect
    scored = [c for c in candidates if not c['own_promoter'] and c['level'] not in ('step1_state',)
              and c['level'] != 'promoter_contacts']
    scored.sort(key=lambda c: -max(c['statistic_before'] - c['statistic_flatten'],
                                   c['statistic_before'] - c['statistic_silence']))
    done = {}
    for c in scored[:n_remove]:
        j = c['node_index']
        if j not in done:
            after_hops = neighbourhood(arrays, source, max_hops, removed=j)
            done[j] = stats_of(after_hops)
            c1 = after_hops[1]
            done[j]['promoter_contacts'] = contact_sign * d * float(_log_ratio(c1[2].sum(), c1[1].sum())) \
                if contact_sign else np.nan
        after = done[j]
        c['statistic_remove'] = after[c['level']]
        c['level_lost_remove'] = bool(after[c['level']] < threshold) if not pd.isna(after[c['level']]) else np.nan
        c['still_explained_remove'] = ','.join(other_levels_explain(c['level'], after))
    return candidates, level_rows, null_hops


def _gene_task(node):
    s = _SHARED
    return gene_design(s['arrays'], s['by_gene'][node], s['max_hops'], s['signs'], s['threshold'], s['n_remove'],
                       s['explains_col'], s['concordant_col'], s['detail_col'], s['state_col'])


# --------------------------------------------------------------------------------------------------------------------
# permutation null of the step levels
# --------------------------------------------------------------------------------------------------------------------


def _null_block(task):
    s = _SHARED
    block, size = task
    rng = np.random.default_rng(np.random.SeedSequence([s['seed'], block]))
    used, position, others = s['used'], s['position'], s['others']
    in_others = position >= 0
    index = np.repeat(used[:, None], size, axis=1)
    for j in range(size):
        donors = others[rng.permutation(len(others))]
        index[in_others, j] = donors[position[in_others]]
    best = np.full((s['H_a'].shape[0], size), -np.inf)
    for m in s['marks']:
        sign = s['signs'].get(m, 0)
        if not sign:
            continue
        va, vb = s['vectors'][m]
        change = np.log2((s['H_b'] @ vb[index] + s['pseudo']) / (s['H_a'] @ va[index] + s['pseudo']))
        np.maximum(best, change * (sign * s['direction'])[:, None], out=best)
    return best


def level_nulls(arrays, hops_by_gene, directions, signs, n_permutations=1000, seed=0, threads=1,
                pseudo_count=PSEUDO_COUNT, block=50):
    """Null statistics of one step level for several genes: {gene node: sorted null values}.

    `hops_by_gene` maps gene node -> (cols, w_A, w_B) of the level. The mark signals of all non-source nodes are
    permuted jointly (gene_level_hop_test), every gene of the level sees the same permutations."""
    genes = list(hops_by_gene)
    if not genes:
        return {}
    n = len(arrays.nodes)
    rows, cols, va_, vb_ = [], [], [], []
    for i, g in enumerate(genes):
        c, wa, wb = hops_by_gene[g]
        rows.append(np.full(len(c), i))
        cols.append(c)
        va_.append(wa)
        vb_.append(wb)
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    used = np.unique(cols)
    local = np.searchsorted(used, cols)
    H_a = sp.csr_matrix((np.concatenate(va_), (rows, local)), shape=(len(genes), len(used)))
    H_b = sp.csr_matrix((np.concatenate(vb_), (rows, local)), shape=(len(genes), len(used)))
    others = np.flatnonzero(~arrays.is_source)
    position = np.full(n, -1, dtype=np.int64)
    position[others] = np.arange(len(others))
    shared = dict(used=used, position=position[used], others=others, H_a=H_a, H_b=H_b, marks=arrays.marks,
                  vectors=arrays.vectors, signs=signs, direction=np.array([directions[g] for g in genes]),
                  seed=seed, pseudo=pseudo_count)
    tasks = [(b, min(block, n_permutations - b * block)) for b in range((n_permutations + block - 1) // block)]
    parts = _map(_null_block, tasks, threads, shared)
    null = np.concatenate(parts, axis=1)
    return {g: np.sort(null[i]) for i, g in enumerate(genes)}


def null_p(null, value):
    """Empirical P of `value` against sorted null values: (#null >= value + 1) / (n + 1)."""
    if null is None or pd.isna(value):
        return np.nan
    # tolerance: the observed value and the null of an unpermuted neighbourhood (only source nodes) agree up to
    # floating-point summation order
    exceed = len(null) - np.searchsorted(null, value - 1e-9 * max(1.0, abs(value)), side='left')
    return float((exceed + 1) / (len(null) + 1))


# --------------------------------------------------------------------------------------------------------------------
# annotation of the candidates
# --------------------------------------------------------------------------------------------------------------------


def read_gene_strands(gtf):
    """Gene name -> (chrom, TSS, strand) from a GTF (gene lines, else the first transcript of a gene name)."""
    import gzip
    opener = gzip.open if str(gtf).endswith('.gz') else open
    out = {}
    with opener(gtf, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] not in ('gene', 'transcript'):
                continue
            match = re.search(r'gene_name "([^"]+)"', parts[8]) or re.search(r'gene_id "([^"]+)"', parts[8])
            if not match:
                continue
            name = match.group(1)
            if name in out and parts[2] == 'transcript':
                continue
            tss = int(parts[3]) if parts[6] == '+' else int(parts[4])
            if name not in out or parts[2] == 'gene':
                out[name] = (parts[0], tss, parts[6])
    return out


class TssIndex:
    """Nearest annotated TSS fragment (elements.tsv tss_genes, else the graph's gene attribute) of a position."""

    def __init__(self, nodes, chroms, starts, ends, genes):
        table = pd.DataFrame(dict(node=nodes, chrom=chroms, mid=(np.asarray(starts) + np.asarray(ends)) // 2,
                                  genes=genes))
        table = table[table['genes'].astype(str).str.len() > 0]
        self.by_chrom = {c: (g['mid'].to_numpy(), g['genes'].to_numpy(), g['node'].to_numpy())
                         for c, g in table.sort_values('mid').groupby('chrom')}

    def nearest(self, chrom, mid, exclude=()):
        if chrom not in self.by_chrom:
            return '', np.nan
        mids, genes, nodes = self.by_chrom[chrom]
        i = np.searchsorted(mids, mid)
        best = ('', np.nan)
        for j in range(max(0, i - 3), min(len(mids), i + 3)):
            if nodes[j] in exclude:
                continue
            dist = abs(int(mids[j]) - int(mid))
            if pd.isna(best[1]) or dist < best[1]:
                best = (genes[j], dist)
        return best


def plain_class(class_a, class_b, names):
    """Element class in plain words: 'active enhancer in both' or 'poised enhancer in A, active enhancer in B'."""
    a, b = str(class_a or 'not classified'), str(class_b or 'not classified')
    if a == 'nan':
        a = 'not classified'
    if b == 'nan':
        b = 'not classified'
    if a == b:
        return '{} in both'.format(a)
    return '{} in {}, {} in {}'.format(a, names[0], b, names[1])


def _kb(distance):
    if pd.isna(distance):
        return ''
    return '{:.0f} kb'.format(distance / 1000.0) if distance >= 1000 else '{:.0f} bp'.format(distance)


def relative_position(gene_chrom, gene_tss, strand, chrom, mid):
    """'120 kb upstream' (5' of the TSS on the gene's strand), 'downstream', or without a strand 'from'."""
    if chrom != gene_chrom:
        return 'on {} (trans)'.format(chrom)
    distance = abs(int(mid) - int(gene_tss))
    if strand not in ('+', '-'):
        return '{} from'.format(_kb(distance))
    upstream = (mid < gene_tss) if strand == '+' else (mid > gene_tss)
    return '{} {} of'.format(_kb(distance), 'upstream' if upstream else 'downstream')


def expected_direction(sign, change, names, gene):
    """Expected expression change of `gene` when the element is silenced (active mark) or deleted (repressive mark):
    (short label, sentence). `sign` is the role sign of the explaining mark, `change` its own log2 B/A change at the
    node; the element is most active (or repressive) in the condition with the larger signal."""
    where = names[1] if change > 0 else names[0]
    if sign > 0:
        return ('down in {}'.format(where),
                'silencing is predicted to lower {} mainly in {}, where the element carries more of the mark'
                .format(gene, where))
    if sign < 0:
        return ('up in {}'.format(where),
                'repressive element (mark higher in {}): CRISPRi does not test it; deletion is predicted to raise {} '
                'mainly in {}'.format(where, gene, where))
    return ('unknown', 'no role sign for the explaining mark')


# --------------------------------------------------------------------------------------------------------------------
# ranking
# --------------------------------------------------------------------------------------------------------------------

EVIDENCE = ('own_change', 'level_p', 'beyond_mirror', 'support_3d', 'contact_both', 'non_promoter', 'not_pleiotropic')


def _ci_low(text):
    match = re.match(r'([\d.]+)-([\d.]+)', str(text))
    return float(match.group(1)) if match else np.nan


def rank_candidates(table, effect='flatten', strong=0.5):
    """design_score and rank per gene. The score is the fraction of the level statistic removed by the perturbation
    (0..1) + 0.25 if the level no longer explains + 0.05 per evidence point (EVIDENCE); the gene's own promoter is
    ranked after every other node (CRISPRi at the promoter silences the gene directly). One row per gene and node is
    kept (the level with the best score; the others in other_levels)."""
    t = table.copy()
    after = t['statistic_' + effect]
    before = t['statistic_before']
    fraction = ((before - after) / before.where(before > 0)).clip(lower=0, upper=1)
    lost = t['level_lost_' + effect].astype('boolean').fillna(False).astype(bool)
    # state level: all of the level or none
    state = t['level'] == 'step1_state'
    fraction = fraction.where(~state, np.where(lost, 1.0, 0.0))
    t['effect_fraction'] = fraction.fillna(0.0)
    points = sum(t['ev_' + e].fillna(False).astype(bool).astype(int) for e in EVIDENCE if 'ev_' + e in t)
    t['evidence_points'] = points
    t['design_score'] = t['effect_fraction'] + 0.25 * lost + 0.05 * points - 2.0 * t['own_promoter'].astype(int)
    t = t.sort_values(['gene_node', 'design_score'], ascending=[True, False])
    other = t.groupby(['gene_node', 'node'])['level'].transform(lambda s: ','.join(s.iloc[1:]))
    t['other_levels'] = other
    t = t.drop_duplicates(['gene_node', 'node'], keep='first')
    t['rank'] = t.groupby('gene_node').cumcount() + 1
    return t


def pleiotropy(table, effect='flatten', strong=0.5):
    """Per node the genes for which it is a strong candidate (removes >= `strong` of a level statistic or the level):
    {node: set of gene labels}."""
    before = table['statistic_before']
    fraction = ((before - table['statistic_' + effect]) / before.where(before > 0)).fillna(0)
    lost = table['level_lost_' + effect].astype('boolean').fillna(False).astype(bool)
    hit = table[((fraction >= strong) | lost) & ~table['own_promoter'].astype(bool)]
    out = {}
    for node, part in hit.groupby('node'):
        out[node] = list(dict.fromkeys(part['genes'].astype(str)))
    return out


# --------------------------------------------------------------------------------------------------------------------
# single-condition mode
# --------------------------------------------------------------------------------------------------------------------


def single_condition_contributions(arrays, sources, max_hops=3, mark=None, pseudo_count=PSEUDO_COUNT):
    """Per gene node and every node within `max_hops` steps: the contribution of the node to the gene's
    contact-weighted activity of one mark in one condition (condition A of `arrays`).

    contribution = path weight * signal; share = contribution / the gene's activity summed over all steps; rank
    within the gene by contribution (the node's part of the active-mark neighbourhood). level_activity X_k = sum of
    the contributions at that step; silence_drop = log2(X_k + 1) - log2(X_k - contribution + 1), the drop of the
    step's log2 activity when the node's mark is set to 0 (the single-condition analogue of the per-level 'silence'
    perturbation), with its own within-gene rank (level_rank; a node alone at its step drops the whole step)."""
    mark = mark or arrays.marks[0]
    v = arrays.vectors[mark][0]
    records = []
    for s in sources:
        if s not in arrays.index:
            continue
        i = arrays.index[s]
        h = neighbourhood(arrays, i, max_hops)
        for k, (cols, w_a, _) in h.items():
            if not len(cols):
                continue
            c = w_a * v[cols]
            x = c.sum()
            drop = np.log2(x + pseudo_count) - np.log2(np.maximum(x - c, 0) + pseudo_count)
            for j in range(len(cols)):
                records.append((s, arrays.nodes[cols[j]], k, float(w_a[j]), float(v[cols[j]]), float(c[j]),
                                float(x), float(drop[j])))
    table = pd.DataFrame.from_records(records, columns=['gene_node', 'node', 'hops', 'path_weight', 'signal',
                                                        'contribution', 'level_activity', 'silence_drop'])
    if table.empty:
        for column in ('share', 'level_rank', 'rank', 'n_nodes'):
            table[column] = []
        return table
    total = table.groupby('gene_node')['contribution'].transform('sum')
    table['share'] = np.where(total > 0, table['contribution'] / total.where(total > 0, 1.0), 0.0)
    table = table.sort_values(['gene_node', 'silence_drop', 'contribution'], ascending=[True, False, False])
    table['level_rank'] = table.groupby('gene_node').cumcount() + 1
    table = table.sort_values(['gene_node', 'contribution', 'silence_drop'], ascending=[True, False, False])
    table['rank'] = table.groupby('gene_node').cumcount() + 1
    table['n_nodes'] = table.groupby('gene_node')['node'].transform('size')
    return table.reset_index(drop=True)


# --------------------------------------------------------------------------------------------------------------------
# text
# --------------------------------------------------------------------------------------------------------------------


def level_phrase(level):
    if level.startswith('promoter_marks') or level.startswith('promoter_state'):
        return 'promoter explanation'
    if level == 'promoter_contacts':
        return 'explanation by the own contacts'
    if level == 'step1_state':
        return 'step-1 state explanation'
    k = level_hops_of(level)
    return '{}-step explanation'.format(k)


def effect_words(fraction, lost):
    if fraction >= 0.9:
        words = 'remove the'
    elif fraction >= 0.5:
        words = 'remove most of the'
    elif fraction >= 0.2:
        words = 'reduce the'
    else:
        words = 'barely change the'
    return words


def _true(value):
    try:
        return (not pd.isna(value)) and bool(value)
    except (TypeError, ValueError):
        return False


def candidate_sentence(r, effect='flatten'):
    """One sentence per candidate (see the module docstring of egaDesign)."""
    gene = str(r['gene_label'])
    if r['own_promoter']:
        return ('CRISPRi at the promoter of {} silences the gene directly and tests the promoter level only'
                .format(gene))
    verb = 'silencing' if r['mark_sign'] >= 0 else 'deleting'
    holds_tss = bool(r['nearest_gene']) and not pd.isna(r['nearest_gene_distance']) and \
        r['nearest_gene_distance'] == 0
    if holds_tss:
        what = '{} the {} of {} ({})'.format(verb, r['element_short'], r['nearest_gene'], r['position'])
        near = ''
    else:
        what = '{} the {} {}'.format(verb, r['element_short'], r['position'])
        near = ' (near {})'.format(r['nearest_gene']) if r['nearest_gene'] and r['nearest_gene'] != gene else ''
    if effect == 'remove' and pd.isna(r.get('statistic_remove', np.nan)):
        effect = 'flatten'
    lost = _true(r['level_lost_' + effect])
    if r['level'] == 'step1_state':
        if lost:
            core = 'is predicted to remove the step-1 state explanation of {} (the only contacting fragment with ' \
                   'this transition)'.format(gene)
        else:
            core = 'leaves {} further contacting fragment(s) with the same transition'.format(
                int(r['state_fragments']) - 1 if not pd.isna(r.get('state_fragments', np.nan)) else 'other')
        text = '{}{} {}'.format(what, near, core)
    else:
        text = '{}{} is predicted to {} {} of {} (statistic {:.1f} -> {:.1f}{})'.format(
            what, near, effect_words(r['effect_fraction'], lost), level_phrase(r['level']),
            gene, r['statistic_before'], r['statistic_' + effect],
            '; the level no longer explains' if lost else '')
    if isinstance(r.get('also_explains'), str) and r['also_explains']:
        text += '; it also contributes to the explanation of {}'.format(r['also_explains'])
    return text
