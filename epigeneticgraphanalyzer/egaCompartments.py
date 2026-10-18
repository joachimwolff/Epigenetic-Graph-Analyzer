"""egaCompartments: A/B compartment estimates from chromatin tracks and capture interactions, optionally evaluated
against a Hi-C eigenvector, and the domain statistics of a graph against TADs."""

import logging

import numpy as np
import pandas as pd

from ._cli import base_parser, ensure_parent, setup_logging, write_table
from .annotation import parse_track_arguments
from .compartments import (bin_bigwig, boundary_enrichment, chromatin_pc1, community_boundaries,
                           compare_to_reference, contact_eigenvector, contacts_to_matrix, domain_boundaries,
                           edges_crossing, flatten_domains, interval_density, log_enrichment, make_bins,
                           read_domains, role_score, smooth_per_chrom, spans_within_domain, standardise_per_chrom)
from .genome import read_chrom_sizes
from .io import read_interactions

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser('A/B compartment estimates per genomic bin from chromatin tracks (first principal '
                         'component oriented by the active marks; active minus repressive role score) and from '
                         'capture interactions (eigenvector of the correlated observed / expected matrix), and '
                         'community / edge statistics of a graph against TADs.', experimental=True)
    parser.add_argument('--chromSizes', required=True, help='Chromosome sizes (UCSC format)')
    parser.add_argument('--binSize', type=int, default=100000, help='Bin size in bp for the chromatin estimates')
    parser.add_argument('--tracks', '-t', nargs='*', default=[], help='name=file.bw[,replicate.bw]; replicates '
                        'are averaged')
    parser.add_argument('--markRoles', nargs='*', default=[],
                        help='name=role; roles active, enhancer, transcription (positive), repressive, '
                             'heterochromatin (negative), other (PC1 only)')
    parser.add_argument('--smooth', type=int, default=0, help='Running mean over 2 * smooth + 1 bins')
    parser.add_argument('--genes', help='BED of genes (or TSS); the gene density phases the contact eigenvector '
                                         '(first of the 3 leading ones unless another correlates 1.5 times more '
                                         'strongly) and orients the PC1 when no active role is declared')
    parser.add_argument('--interactions', '-i', nargs='*', default=[], help='Capture interaction files')
    parser.add_argument('--format', '-f', default='washu', help='Interaction format (see egaCreateGraph)')
    parser.add_argument('--scoreThreshold', type=float, default=None, help='Keep interactions with score >= this')
    parser.add_argument('--contactBinSize', type=int, default=250000, help='Bin size of the contact eigenvector')
    parser.add_argument('--minContacts', type=int, default=3,
                        help='Bins with fewer non-empty matrix entries get no contact eigenvector')
    parser.add_argument('--reference', help='Hi-C eigenvector (bedGraph, positive = A) to evaluate against')
    parser.add_argument('--nodes', help='Node table of a graph (egaCreateGraph *_nodes.tsv) for the domain test')
    parser.add_argument('--edges', help='Edge table of the same graph (*_edges.tsv)')
    parser.add_argument('--domains', help='TADs (BED or Rao et al. domain list) for the domain test')
    parser.add_argument('--shifts', type=int, default=1000, help='Circular shifts of the domain test')
    parser.add_argument('--window', type=int, default=50000, help='Distance in bp to count a boundary as matched')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--outputPrefix', '-o', required=True,
                        help='Writes <prefix>_bins.tsv, <prefix>_compartments.bedGraph, <prefix>_evaluation.tsv, '
                             '<prefix>_loadings.tsv and <prefix>_domains.tsv as applicable')
    return parser.parse_args(args)


def _read_reference(path, bins):
    table = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'value'],
                        comment='#')
    table = table[~table.chrom.astype(str).str.startswith(('track', 'browser'))]
    table = table.astype({'start': np.int64, 'end': np.int64, 'value': float})
    out = np.full(len(bins), np.nan)
    mids = (bins.start.values + bins.end.values) // 2
    for chrom, idx in bins.groupby('chrom', sort=False).indices.items():
        sub = table[table.chrom == chrom].sort_values('start')
        if sub.empty:
            continue
        k = np.searchsorted(sub.start.values, mids[idx], side='right') - 1
        ok = (k >= 0) & (mids[idx] < sub.end.values[np.clip(k, 0, None)])
        out[idx[ok]] = sub.value.values[k[ok]]
    return out


def _gene_density(path, bins):
    genes = pd.read_csv(path, sep='\t', header=None, comment='#')
    return np.log2(interval_density(bins, genes[0].values, genes[1].values, genes[2].values) + 1)


def chromatin_estimates(args, bins, genes):
    tracks = parse_track_arguments(args.tracks)
    roles = dict(entry.split('=', 1) for entry in args.markRoles)
    names = list(tracks)
    features = []
    for name in names:
        values = np.nanmean(np.column_stack([bin_bigwig(path, bins) for path in tracks[name]]), axis=1)
        features.append(log_enrichment(values))
    features = np.column_stack(features)
    role_list = [roles.get(n, 'other') for n in names]
    active = np.array([r in ('active', 'enhancer', 'transcription') for r in role_list])
    pc1, loadings, explained = chromatin_pc1(features, active=active if active.any() else None, orient_by=genes)
    out = {'chromatin_pc1': smooth_per_chrom(pc1, bins.chrom.values, args.smooth)}
    if active.any():
        out['chromatin_role'] = smooth_per_chrom(role_score(features, role_list), bins.chrom.values, args.smooth)
    loads = pd.DataFrame({'track': names, 'role': role_list, 'loading': loadings, 'explained_pc1': explained})
    return out, loads


def contact_estimate(args, sizes, bins, genes_path):
    table = read_interactions(args.interactions, args.format)
    if args.scoreThreshold is not None:
        table = table[table.score >= args.scoreThreshold]
    res = args.contactBinSize
    cbins = make_bins(sizes, res)
    genes = _gene_density(genes_path, cbins) if genes_path else None
    vec = np.full(len(cbins), np.nan)
    for chrom, idx in cbins.groupby('chrom', sort=False).indices.items():
        sub = table[table.chrom == chrom]
        if sub.empty:
            continue
        a = ((sub.bait_start.values + sub.bait_end.values) // 2) // res
        b = ((sub.other_start.values + sub.other_end.values) // 2) // res
        keep = (a < len(idx)) & (b < len(idx))
        weights = sub.score.fillna(1.0).values[keep]
        m = contacts_to_matrix(len(idx), a[keep], b[keep], weights)
        v, _ = contact_eigenvector(m, orient_by=None if genes is None else genes[idx], min_contacts=args.minContacts,
                                   n_candidates=1 if genes is None else 3)
        vec[idx] = v
    # onto the chromatin bins by midpoint
    mids = (bins.start.values + bins.end.values) // 2
    lookup = pd.Series(vec, index=pd.MultiIndex.from_arrays([cbins.chrom.values, cbins.start.values // res]))
    keys = list(zip(bins.chrom.values, mids // res))
    return lookup.reindex(keys).values, (cbins.assign(contact_eigenvector=vec))


def domain_tests(args, sizes):
    nodes = pd.read_csv(args.nodes, sep='\t')
    domains = flatten_domains(read_domains(args.domains))
    nodes['mid'] = (nodes.start + nodes.end) // 2
    spans = {c: (int(g.mid.min()), int(g.mid.max()) + 1) for c, g in nodes.groupby('chrom')}
    cb = community_boundaries(nodes[['chrom', 'start', 'end', 'community']])
    gaps = {c: (g.end.values[:-1] + g.start.values[1:]) // 2 for c, g in nodes.sort_values(['chrom', 'mid'])
            .groupby('chrom')}
    rows = []
    res = boundary_enrichment({c: g.position.values for c, g in cb.groupby('chrom')}, domain_boundaries(domains),
                              spans, window=args.window, n_shifts=args.shifts, seed=args.seed, null_positions=gaps)
    rows.append(dict(test='community_boundaries_at_tad_boundaries', **res))
    com = nodes.groupby('community').agg(chrom=('chrom', 'first'), n_chrom=('chrom', 'nunique'),
                                         start=('start', 'min'), end=('end', 'max'), n=('node', 'size'))
    com = com[(com.n >= 3) & (com.n_chrom == 1)].reset_index(drop=True)
    rows.append(dict(test='communities_within_one_tad',
                     **spans_within_domain(com[['chrom', 'start', 'end']], domains, spans, args.shifts, args.seed)))
    if args.edges:
        edges = pd.read_csv(args.edges, sep='\t', usecols=['node_a', 'node_b'])
        pos = nodes.set_index('node')
        ea, eb = pos.loc[edges.node_a.values], pos.loc[edges.node_b.values]
        same = ea.chrom.values == eb.chrom.values
        ed = pd.DataFrame({'chrom': ea.chrom.values[same], 'mid_a': ea.mid.values[same], 'mid_b': eb.mid.values[same]})
        rows.append(dict(test='edges_crossing_tad_boundaries', **edges_crossing(ed, domains, spans, args.shifts,
                                                                                args.seed)))
    return pd.DataFrame(rows)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    sizes = read_chrom_sizes(args.chromSizes)
    bins = make_bins(sizes, args.binSize)
    genes = _gene_density(args.genes, bins) if args.genes else None
    estimates = {}
    if args.tracks:
        est, loads = chromatin_estimates(args, bins, genes)
        estimates.update(est)
        write_table(loads, args.outputPrefix + '_loadings.tsv')
    if args.interactions:
        estimates['contact_eigenvector'], coarse = contact_estimate(args, sizes, bins, args.genes)
        write_table(coarse, args.outputPrefix + '_contact_eigenvector.tsv')
        if 'chromatin_pc1' in estimates:
            zc = standardise_per_chrom(estimates['chromatin_pc1'], bins.chrom.values)
            zk = standardise_per_chrom(estimates['contact_eigenvector'], bins.chrom.values)
            estimates['combined'] = np.where(np.isfinite(zk), (zc + zk) / 2, zc)
    if estimates:
        frame = bins.copy()
        for key, values in estimates.items():
            frame[key] = values
        write_table(frame, args.outputPrefix + '_bins.tsv')
        main_key = 'combined' if 'combined' in estimates else next(iter(estimates))
        out = frame[['chrom', 'start', 'end', main_key]].dropna()
        out.to_csv(ensure_parent(args.outputPrefix + '_compartments.bedGraph'), sep='\t', header=False, index=False,
                   float_format='%.5g')
        if args.reference:
            ref = _read_reference(args.reference, bins)
            tables = []
            for key, values in estimates.items():
                t = compare_to_reference(values, ref, bins.chrom.values,
                                         threshold='zero' if key == 'contact_eigenvector' else 'median')
                t.insert(0, 'estimate', key)
                tables.append(t)
            write_table(pd.concat(tables, ignore_index=True), args.outputPrefix + '_evaluation.tsv')
    if args.nodes and args.domains:
        write_table(domain_tests(args, sizes), args.outputPrefix + '_domains.tsv')
    if not estimates and not (args.nodes and args.domains):
        log.warning('nothing to do: give --tracks, --interactions or --nodes with --domains')


if __name__ == '__main__':
    main()
