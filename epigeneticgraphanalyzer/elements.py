"""Element classes of the fragments of a differential run (egaDifferential), per condition, from the marks present and
their roles; transitions of the class between the conditions with the expression change of the genes in contact;
roles of the marks learned per element type; agreement with an external element annotation.

Promoter status
    A fragment is a promoter when it overlaps an annotated transcription start site (the TSS annotation of the run,
    or --annotation), whether it is baited or not: 'baited promoter', 'non-baited promoter' (an other end that holds
    a TSS), 'baited non-promoter' (a bait without an annotated TSS) and 'non-promoter element'. Element type:
    'promoter' (baited or not) or 'non-promoter element'.

Chromatin class
    A mark is present on a fragment in a condition when its signal reaches the presence threshold of the states
    (fold enrichment, statistic fragment mean or maximum 1 kb window). The class follows from the roles of the
    present marks only; no mark name enters the rules. In order of precedence:

    promoters: an active-role and a repressive mark 'bivalent promoter'; an active-role mark 'active promoter'; an
    enhancer-role and a repressive mark 'poised promoter'; an enhancer-role mark 'primed promoter'; a primary
    repressive mark 'repressed promoter'; another repressive mark only 'heterochromatin-like promoter'; an
    architectural mark only 'architectural promoter'; marks without a role (neutral, unknown) only 'promoter with
    marks of no role'; nothing present 'unmarked promoter'.

    non-promoter elements: an active- or enhancer-role and a repressive mark 'poised enhancer'; an active-role mark
    'active enhancer'; an enhancer-role mark only 'primed enhancer' (when an active-role mark is among the tracks, so
    that its absence is measured; otherwise 'enhancer-marked element'); a primary repressive mark 'repressed
    element'; another repressive mark only 'heterochromatin-like element'; an architectural mark only 'architectural
    element'; marks without a role only 'element with marks of no role'; nothing present 'unmarked element'.

    Primary repressive marks are those given as primary (e.g. the Polycomb-type mark), by default the marks declared
    repressive; a mark that is repressive under a role set without being primary (e.g. learned as repressive without
    a declaration) gives 'heterochromatin-like'. When no primary mark is repressive under a role set, every repressive
    mark is primary. The present marks are reported with every class.

    Classes are given under the declared roles and under the learned roles: a baited promoter takes the roles learned
    at the promoter (step 0), every other fragment the roles learned at step 1 (direct contact), the level at which a
    fragment acts on the genes it contacts.
"""

import logging

import numpy as np
import pandas as pd

from .differential import _pool_map, level_hops, level_name
from .io import normalize_chrom, open_text
from .validation import block_resamples

log = logging.getLogger(__name__)

TSS_STATUS = ('baited promoter', 'non-baited promoter', 'baited non-promoter', 'non-promoter element')
ELEMENT_TYPES = ('promoter', 'non-promoter element')
PROMOTER_CLASSES = ('active promoter', 'bivalent promoter', 'poised promoter', 'primed promoter',
                    'repressed promoter', 'heterochromatin-like promoter', 'architectural promoter',
                    'promoter with marks of no role', 'unmarked promoter')
ELEMENT_CLASSES = ('active enhancer', 'poised enhancer', 'primed enhancer', 'enhancer-marked element',
                   'repressed element', 'heterochromatin-like element', 'architectural element',
                   'element with marks of no role', 'unmarked element')
CLASSES = PROMOTER_CLASSES + ELEMENT_CLASSES


def tss_status(graph, tss_genes=None):
    """node -> (TSS status, element type, genes whose TSS lies on the fragment). `tss_genes` maps node -> gene names
    of an other annotation (graph.annotate_genes); default: the 'genes' attribute of the nodes (the TSS annotation of
    the run)."""
    out = {}
    for node, data in graph.nodes(data=True):
        if tss_genes is not None:
            genes = list(tss_genes.get(node, ([], []))[0])
        else:
            genes = [g for g in str(data.get('genes', '') or '').split(',') if g]
        bait = bool(data.get('bait', False))
        if genes:
            status = 'baited promoter' if bait else 'non-baited promoter'
        else:
            status = 'baited non-promoter' if bait else 'non-promoter element'
        out[node] = (status, 'promoter' if genes else 'non-promoter element', ','.join(genes))
    return out


def primary_repressive(roles, primary=None):
    """Primary repressive marks of a role set: `primary` among the marks repressive here, else all repressive."""
    repressive = [m for m, r in roles.items() if r == 'repressive']
    chosen = [m for m in repressive if primary and m in primary]
    return chosen if chosen else repressive


def element_classes(present, roles, is_promoter, primary=None):
    """Class label of every fragment (module doc) from `present` (mark -> bool array), `roles` (mark -> role),
    `is_promoter` (bool array) and the primary repressive marks `primary` (default: all repressive marks)."""
    is_promoter = np.asarray(is_promoter, dtype=bool)
    n = len(is_promoter)
    marks = [m for m in present if m in roles]
    primary_marks = set(primary_repressive(roles, primary))

    def anyof(selected):
        if not selected:
            return np.zeros(n, dtype=bool)
        return np.any([np.asarray(present[m], dtype=bool) for m in selected], axis=0)

    active = anyof([m for m in marks if roles[m] == 'active'])
    enhancer = anyof([m for m in marks if roles[m] == 'enhancer'])
    rep_primary = anyof([m for m in marks if roles[m] == 'repressive' and m in primary_marks])
    rep_other = anyof([m for m in marks if roles[m] == 'repressive' and m not in primary_marks])
    architectural = anyof([m for m in marks if roles[m] == 'architectural'])
    other = anyof([m for m in marks if roles[m] in ('neutral', 'unknown')])
    repressive = rep_primary | rep_other
    active_measured = any(roles[m] == 'active' for m in marks)
    promoter = np.select(
        [active & repressive, active, enhancer & repressive, enhancer, rep_primary, rep_other, architectural, other],
        ['bivalent promoter', 'active promoter', 'poised promoter', 'primed promoter', 'repressed promoter',
         'heterochromatin-like promoter', 'architectural promoter', 'promoter with marks of no role'],
        'unmarked promoter')
    element = np.select(
        [(active | enhancer) & repressive, active, enhancer, rep_primary, rep_other, architectural, other],
        ['poised enhancer', 'active enhancer', 'primed enhancer' if active_measured else 'enhancer-marked element',
         'repressed element', 'heterochromatin-like element', 'architectural element',
         'element with marks of no role'],
        'unmarked element')
    return np.where(is_promoter, promoter, element).astype(object)


def present_marks(present, marks):
    """'H3K4me1+H3K27me3' of every fragment (marks in the given order), 'none' without one."""
    n = len(next(iter(present.values()))) if present else 0
    labels = np.full(n, '', dtype=object)
    for mark in marks:
        labels = np.where(present[mark], np.where(labels == '', mark, labels + '+' + mark), labels)
    return np.where(labels == '', 'none', labels).astype(object)


def node_elements(graph, state_signals_by, threshold, roles_by, status, names=('A', 'B'), primary=None,
                  bait_roles_by=None):
    """Element table: one row per node and condition with tss_status, element_type, tss_genes and per presence
    statistic S (fragment mean 'mean', maximum 1 kb window 'max1kb'): present_marks_S and, per role set R,
    class_R_S. `state_signals_by` maps S -> mark -> (node -> signal in A, node -> signal in B); `roles_by` maps R
    -> roles of the fragments that are not baited promoters; `bait_roles_by` (default `roles_by`) R -> roles of the
    baited promoters."""
    nodes = list(graph.nodes)
    bait_roles_by = bait_roles_by or roles_by
    info = [status[n] for n in nodes]
    status_column = np.array([s for s, _, _ in info], dtype=object)
    is_promoter = np.array([t == 'promoter' for _, t, _ in info], dtype=bool)
    baited_promoter = status_column == 'baited promoter'
    frames = []
    for c, condition in enumerate(names):
        frame = pd.DataFrame({'node': nodes, 'condition': condition, 'tss_status': status_column,
                              'element_type': np.where(is_promoter, 'promoter', 'non-promoter element'),
                              'tss_genes': [g for _, _, g in info]})
        for presence, signals in state_signals_by.items():
            marks = list(signals)
            present = {m: np.array([signals[m][c].get(n, 0.0) for n in nodes], dtype=float) >= threshold
                       for m in marks}
            frame['present_marks_' + presence] = present_marks(present, marks)
            for roles_name in roles_by:
                other = element_classes(present, roles_by[roles_name], is_promoter, primary)
                bait = element_classes(present, bait_roles_by[roles_name], is_promoter, primary)
                frame['class_{}_{}'.format(roles_name, presence)] = np.where(baited_promoter, bait, other)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def class_columns(elements):
    """(role set, presence statistic) of the class_<R>_<S> columns of an element table."""
    out = []
    for column in elements.columns:
        if column.startswith('class_'):
            roles, presence = column[len('class_'):].rsplit('_', 1)
            out.append((roles, presence))
    return out


def element_class_counts(elements, other_ends=None):
    """Fragments per role set, presence statistic, condition, TSS status and class, with the fraction among the
    fragments of that TSS status; `other_ends` (set of nodes that are an other end of a called interaction) adds
    the column other_end_fragments."""
    rows = []
    for roles, presence in class_columns(elements):
        column = 'class_{}_{}'.format(roles, presence)
        for (condition, status, cls), part in elements.groupby(['condition', 'tss_status', column], sort=True):
            total = int(((elements['condition'] == condition) & (elements['tss_status'] == status)).sum())
            row = dict(roles=roles, presence=presence, condition=condition, tss_status=status,
                       element_type='promoter' if status in ('baited promoter', 'non-baited promoter')
                       else 'non-promoter element', element_class=cls, fragments=len(part),
                       fraction_of_status=len(part) / total if total else np.nan)
            if other_ends is not None:
                row['other_end_fragments'] = int(part['node'].isin(other_ends).sum())
            rows.append(row)
    return pd.DataFrame(rows)


def _chromosomes(nodes):
    return np.array([str(n).split(':')[0] for n in nodes], dtype=object)


class _GeneSummary:
    """Mean, median and 95% chromosome block-bootstrap interval of the expression change of gene sets."""

    def __init__(self, nodes, rna, rna_threshold=1.0, n_boot=1000, seed=0, min_genes=5):
        self.rna = np.asarray(rna, dtype=float)
        self.rna_threshold = rna_threshold
        self.min_genes = min_genes
        chroms = _chromosomes(nodes)
        self.labels, resamples = block_resamples(chroms, n_boot, seed)
        self.chrom_index = np.searchsorted(self.labels, chroms)
        self.weights = np.array([np.bincount(r, minlength=len(self.labels)) for r in resamples], dtype=float) \
            if len(self.labels) else np.zeros((0, 0))

    def __call__(self, index):
        index = np.unique(np.asarray(index, dtype=np.int64))
        values = self.rna[index]
        up = int((values >= self.rna_threshold).sum())
        down = int((values <= -self.rna_threshold).sum())
        out = dict(n_genes=len(index), mean_rna_log2FC=float(values.mean()) if len(index) else np.nan,
                   median_rna_log2FC=float(np.median(values)) if len(index) else np.nan, ci_low=np.nan,
                   ci_high=np.nan, n_up=up, n_down=down)
        if len(index) >= self.min_genes and len(self.weights):
            sums = np.bincount(self.chrom_index[index], weights=values, minlength=len(self.labels))
            sizes = np.bincount(self.chrom_index[index], minlength=len(self.labels)).astype(float)
            with np.errstate(invalid='ignore', divide='ignore'):
                boot = (self.weights @ sums) / (self.weights @ sizes)
            if np.isfinite(boot).any():
                out['ci_low'], out['ci_high'] = float(np.nanpercentile(boot, 2.5)), float(np.nanpercentile(boot, 97.5))
        return out


def step_pattern(operators, level):
    """Sparse 0/1 pattern (sources x nodes) of the nodes at a level with a path weight in either condition."""
    H_a, H_b, _ = operators.level(level)
    pattern = ((abs(H_a) > 0).astype(np.int8) + (abs(H_b) > 0).astype(np.int8)).tocsr()
    pattern.eliminate_zeros()
    pattern.data[:] = 1
    return pattern


def element_transitions(elements, operators, features, rna_threshold=1.0, n_boot=1000, seed=0, min_genes=5,
                        differential=None):
    """Transitions of the element class between the conditions and the expression change of the genes in contact.

    For every role set and presence statistic, every fragment has a transition (class in A -> class in B; unchanged
    classes included). Per step k (0: the promoter fragment of the gene itself; k >= 1: the nodes at exactly k
    steps with a path weight in either condition), the genes of `features` (bait genes with an expression change)
    having at least one fragment with that transition at that step: n_genes, mean and median expression change with
    a 95% chromosome block-bootstrap interval of the mean (groups of at least `min_genes`), genes changing at least
    `rna_threshold` up / down, and baseline_mean, the mean over all genes with a node at that step. `differential`
    (set of nodes) adds the numbers of differential genes up and down. Columns: roles, presence, step, transition,
    class_A, class_B, changed, n_fragments (fragments with the transition), n_genes, ..."""
    names = list(dict.fromkeys(elements['condition']))
    data = features.dropna(subset=['rna_log2FC'])
    data = data[data['node'].isin(operators.index)].reset_index(drop=True)
    source_row = {n: i for i, n in enumerate(operators.sources)}
    data = data[data['node'].isin(source_row)].reset_index(drop=True)
    gene_rows = np.array([source_row[n] for n in data['node']], dtype=np.int64)
    summary = _GeneSummary(data['node'], data['rna_log2FC'], rna_threshold, n_boot, seed, min_genes)
    gene_nodes = np.array([operators.index[n] for n in data['node']], dtype=np.int64)
    is_differential = data['node'].isin(differential).to_numpy() if differential is not None else None
    rna = data['rna_log2FC'].to_numpy(dtype=float)
    patterns = {k: step_pattern(operators, level) for k, level in
                ((level_hops(lv), lv) for lv in operators.levels)}
    a = elements[elements['condition'] == names[0]].set_index('node')
    b = elements[elements['condition'] == names[1]].set_index('node')
    rows = []
    for roles, presence in class_columns(elements):
        column = 'class_{}_{}'.format(roles, presence)
        class_a = a[column].reindex(operators.nodes).fillna('').to_numpy(dtype=object)
        class_b = b[column].reindex(operators.nodes).fillna('').to_numpy(dtype=object)
        pairs = pd.Series(list(zip(class_a, class_b)))
        codes, uniques = pd.factorize(pairs)
        fragments = np.bincount(codes, minlength=len(uniques))
        steps = {0: (np.arange(len(data)), gene_nodes)}
        for k, pattern in patterns.items():
            sub = pattern[gene_rows].tocoo()
            steps[k] = (sub.row, sub.col)
        for k in sorted(steps):
            genes, cols = steps[k]
            if not len(genes):
                continue
            baseline = float(np.mean(rna[np.unique(genes)]))
            key = np.unique(genes.astype(np.int64) * len(uniques) + codes[cols])
            gene_of, code_of = key // len(uniques), key % len(uniques)
            order = np.argsort(code_of, kind='mergesort')
            gene_of, code_of = gene_of[order], code_of[order]
            bounds = np.flatnonzero(np.diff(code_of)) + 1
            for members, code in zip(np.split(gene_of, bounds), code_of[np.r_[0, bounds]] if len(code_of) else []):
                ca, cb = uniques[code]
                row = dict(roles=roles, presence=presence, step=k, transition='{} -> {}'.format(ca, cb),
                           class_A=ca, class_B=cb, changed=ca != cb, n_fragments=int(fragments[code]))
                row.update(summary(members))
                row['baseline_mean'] = baseline
                if is_differential is not None:
                    de = members[is_differential[members]]
                    row['differential_up'] = int((rna[de] > 0).sum())
                    row['differential_down'] = int((rna[de] < 0).sum())
                rows.append(row)
    result = pd.DataFrame(rows)
    if len(result):
        result = result.sort_values(['roles', 'presence', 'step', 'n_genes', 'transition'],
                                    ascending=[True, True, True, False, True]).reset_index(drop=True)
    return result


# --------------------------------------------------------------------------------------------------------------------
# external element annotation

def read_element_bed(path):
    """Element annotation (BED with a name column, e.g. ENCODE cCREs or ChromHMM states): chrom, start, end, name."""
    rows = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(('#', 'track', 'browser')) or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4:
                raise ValueError('{}: the element annotation needs a name column (column 4)'.format(path))
            rows.append((normalize_chrom(fields[0]), int(fields[1]), int(fields[2]), fields[3]))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'name'])


def assign_elements(graph, annotation):
    """node -> name of the annotated element with the largest overlap of the fragment ('none' without overlap)."""
    out = {}
    by_chrom = {c: part.sort_values('start') for c, part in annotation.groupby('chrom')}
    nodes = pd.DataFrame([(n, d['chrom'], int(d['start']), int(d['end'])) for n, d in graph.nodes(data=True)],
                         columns=['node', 'chrom', 'start', 'end'])
    for chrom, frag in nodes.groupby('chrom'):
        part = by_chrom.get(chrom)
        if part is None or not len(part):
            out.update({n: 'none' for n in frag['node']})
            continue
        starts = part['start'].to_numpy()
        ends = part['end'].to_numpy()
        labels = part['name'].to_numpy()
        longest = int((ends - starts).max())
        for node, s, e in zip(frag['node'], frag['start'], frag['end']):
            lo = int(np.searchsorted(starts, s - longest, side='left'))
            hi = int(np.searchsorted(starts, e, side='left'))
            overlap = np.minimum(ends[lo:hi], e) - np.maximum(starts[lo:hi], s)
            if len(overlap) and overlap.max() > 0:
                out[node] = str(labels[lo + int(np.argmax(overlap))])
            else:
                out[node] = 'none'
    return out


def element_agreement(elements, external_by_condition):
    """Confusion table of the EGA classes against an external annotation: roles, presence, condition, ega_class,
    external, fragments, fraction_of_ega_class. `external_by_condition` maps condition -> (node -> name)."""
    rows = []
    for roles, presence in class_columns(elements):
        column = 'class_{}_{}'.format(roles, presence)
        for condition, assigned in external_by_condition.items():
            part = elements[elements['condition'] == condition]
            external = part['node'].map(assigned).fillna('none')
            table = pd.crosstab(part[column], external)
            totals = table.sum(axis=1)
            for ega_class in table.index:
                for name in table.columns:
                    n = int(table.loc[ega_class, name])
                    if n:
                        rows.append(dict(roles=roles, presence=presence, condition=condition, ega_class=ega_class,
                                         external=name, fragments=n,
                                         fraction_of_ega_class=n / totals[ega_class]))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------------------------------------
# roles learned per element type

def stratum_masks(elements, operators, roles='declared', presence='mean'):
    """Node strata of the role learning per element type: bool arrays over operators.nodes for 'promoter' (baited or
    not), 'baited promoter', 'non-baited promoter', 'non-promoter element' and, per chromatin class (under `roles`
    and `presence`), 'class: <label>' (the fragment has that class in either condition)."""
    first = elements.drop_duplicates('node').set_index('node')
    status = first['tss_status'].reindex(operators.nodes).fillna('non-promoter element').to_numpy()
    masks = {'promoter': np.isin(status, ['baited promoter', 'non-baited promoter']),
             'baited promoter': status == 'baited promoter',
             'non-baited promoter': status == 'non-baited promoter',
             'non-promoter element': ~np.isin(status, ['baited promoter', 'non-baited promoter'])}
    column = 'class_{}_{}'.format(roles, presence)
    if column in elements:
        per_condition = [part.set_index('node')[column].reindex(operators.nodes).to_numpy(dtype=object)
                         for _, part in elements.groupby('condition', sort=False)]
        labels = {v for values in per_condition for v in values if isinstance(v, str) and v}
        for label in sorted(labels):
            mask = np.zeros(len(operators.nodes), dtype=bool)
            for values in per_condition:
                mask |= values == label
            masks['class: {}'.format(label)] = mask
    return masks


def stratified_features(features, operators, mark_vectors, mask, pseudo_count=1.0):
    """Feature table of learn_mark_roles restricted to the nodes of `mask`: at the promoter, the genes whose own
    fragment is in the stratum (n_promoter 1 / 0, the promoter change unchanged); at every step the contact-weighted
    mark change summed over the nodes of the stratum at exactly k steps, n_<level>: those nodes with a path weight in
    either condition."""
    source_row = {n: i for i, n in enumerate(operators.sources)}
    rows = np.array([source_row[n] for n in features['node']], dtype=np.int64)
    out = pd.DataFrame({'node': features['node'].to_numpy(), 'rna_log2FC': features['rna_log2FC'].to_numpy()})
    weight = mask.astype(float)
    out['n_promoter'] = mask[operators.source_index[rows]].astype(int)
    for mark in mark_vectors:
        out[mark + '_promoter'] = features[mark + '_promoter'].to_numpy()
    for level in operators.levels:
        H_a, H_b, _ = operators.level(level)
        out['n_' + level] = (step_pattern(operators, level) @ weight)[rows].astype(int)
        for mark, vectors in mark_vectors.items():
            # signals of the step feature (contact term: their mean in both conditions)
            vector_a, vector_b = operators.step_vectors(*vectors)
            a = (H_a @ (vector_a * weight))[rows]
            b = (H_b @ (vector_b * weight))[rows]
            out['{}_{}'.format(mark, level)] = np.log2((b + pseudo_count) / (a + pseudo_count))
    return out


def _stratum_roles(name):
    from .states import learn_mark_roles
    s = _SHARED_ROLES
    features = stratified_features(s['features'], s['operators'], s['mark_vectors'], s['masks'][name])
    table = learn_mark_roles(features, s['marks'], s['declared'], s['chromatin_threshold'], s['rna_threshold'],
                             n_boot=s['n_boot'], seed=s['seed'], levels=s['operators'].levels)
    if name.startswith('class: '):
        table = table[table['n_genes'] >= s['min_genes']]
    else:
        table = table[table['n_genes'] > 0]
    return table.assign(element_type=name)


_SHARED_ROLES = {}


def roles_per_element_type(features, operators, mark_vectors, marks, declared, masks, chromatin_threshold=1.0,
                           rna_threshold=1.0, n_boot=1000, seed=0, min_genes=50, threads=1):
    """learn_mark_roles per node stratum (stratum_masks): the same statistics and chromosome block-bootstrap
    intervals, with the column element_type naming the stratum. Rows of chromatin-class strata with fewer than
    `min_genes` genes are dropped; the levels where a stratum has no node are dropped."""
    global _SHARED_ROLES
    _SHARED_ROLES = dict(features=features, operators=operators, mark_vectors=mark_vectors, masks=masks, marks=marks,
                         declared=declared, chromatin_threshold=chromatin_threshold, rna_threshold=rna_threshold,
                         n_boot=n_boot, seed=seed, min_genes=min_genes)
    try:
        tables = _pool_map(_stratum_roles, list(masks), threads, {})
    finally:
        _SHARED_ROLES = {}
    tables = [t for t in tables if len(t)]
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


def node_types(elements):
    """node -> element type ('promoter' / 'non-promoter element') of an element table."""
    first = elements.drop_duplicates('node')
    return dict(zip(first['node'], first['element_type']))


__all__ = ['TSS_STATUS', 'ELEMENT_TYPES', 'CLASSES', 'tss_status', 'element_classes', 'node_elements',
           'element_class_counts', 'element_transitions', 'read_element_bed', 'assign_elements', 'element_agreement',
           'stratum_masks', 'stratified_features', 'roles_per_element_type', 'node_types', 'level_name']
