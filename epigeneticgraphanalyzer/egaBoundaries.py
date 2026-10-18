"""egaBoundaries: domain boundaries from capture insulation, CTCF / cohesin binding and CTCF motif orientation."""

import logging

import numpy as np
import pandas as pd

from . import boundaries as bd
from ._cli import base_parser, ensure_parent, setup_logging, write_table
from .compartments import make_bins
from .genome import read_chrom_sizes

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser('Per-bin boundary evidence: capture insulation dip (graph edges spanning each bin against a '
                         'distance-matched mirror expectation), cohesin-bound CTCF with motif orientation (reverse '
                         'sites on the left, forward on the right), convergent CTCF loops among the edges; boundaries '
                         'by a rule with fixed thresholds (no fitted weights).', experimental=True)
    parser.add_argument('--nodes', required=True, help='Node table of a graph (egaCreateGraph *_nodes.tsv)')
    parser.add_argument('--edges', required=True, help='Edge table of the same graph (*_edges.tsv)')
    parser.add_argument('--chromSizes', required=True, help='Chromosome sizes (UCSC format)')
    parser.add_argument('--chromosomes', nargs='*', default=None, help='Restrict to these chromosomes')
    parser.add_argument('--binSize', type=int, default=25000)
    parser.add_argument('--minExpected', type=float, default=5.0,
                        help='Bins with fewer expected spanning edges are not covered (no insulation)')
    parser.add_argument('--dipProminence', type=float, default=0.3,
                        help='Dip: local minimum of the capture insulation with at least this prominence (log2)')
    parser.add_argument('--dipMax', type=float, default=0.0, help='Dip: capture insulation (log2 o/e) at most this')
    parser.add_argument('--ctcfMin', type=float, default=10.0,
                        help='CTCF evidence: cohesin-bound sites with CTCF enrichment >= this')
    parser.add_argument('--window', type=int, default=1, help='CTCF and loop evidence within +- this many bins')
    parser.add_argument('--rule', choices=list(bd.RULES), default=None,
                        help='Rule for the boundary calls (default: dip_and_ctcf_pattern with CTCF, dip without); '
                             'the calls of every rule are columns of <prefix>_bins.tsv')
    parser.add_argument('--ctcf', nargs='*', default=[], help='CTCF ChIP bigWig(s), averaged')
    parser.add_argument('--cohesin', nargs='*', default=[], help='Cohesin (RAD21 / SMC1A / SMC3) bigWig(s)')
    parser.add_argument('--fasta', help='Genome FASTA (plain or gzipped) for the motif orientation')
    parser.add_argument('--motifs', help='JASPAR-format motif file')
    parser.add_argument('--motifId', default='MA0139', help='CTCF motif ID (prefix match)')
    parser.add_argument('--motifPvalue', type=float, default=1e-4)
    parser.add_argument('--minEnrichment', type=float, default=10.0,
                        help='CTCF summits: minimum signal over the genome-wide mean')
    parser.add_argument('--reference', help='BED of reference (Hi-C) boundaries; boundary = interval midpoint')
    parser.add_argument('--outputPrefix', '-o', required=True,
                        help='Writes <prefix>_bins.tsv, <prefix>_capture_insulation.bedGraph, <prefix>_boundaries.bed, '
                             '<prefix>_domains.bed and, with CTCF, <prefix>_ctcf_sites.bed and <prefix>_loops.tsv; '
                             'with --reference <prefix>_evaluation.tsv')
    return parser.parse_args(args)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    sizes = read_chrom_sizes(args.chromSizes)
    if args.chromosomes:
        sizes = {c: s for c, s in sizes.items() if c in set(args.chromosomes)}
    nodes = pd.read_csv(args.nodes, sep='\t')
    edges = pd.read_csv(args.edges, sep='\t')
    nodes = nodes[nodes.chrom.isin(sizes)]
    sizes = {c: s for c, s in sizes.items() if c in set(nodes.chrom)}
    bins = make_bins(sizes, args.binSize)
    et = bd.edge_table(nodes, edges)
    ins = bd.capture_insulation(bins, et, min_expected=args.minExpected)
    log.info('%d cis edges, %.1f %% of bins covered', len(et), 100 * ins.covered.mean())
    sites = loops = None
    if args.ctcf:
        sites = bd.ctcf_sites(args.ctcf, args.cohesin, sizes, min_enrichment=args.minEnrichment)
        if not args.cohesin:
            sites['cohesin_bound'] = True
        if args.fasta and args.motifs:
            pwm = bd.log_odds(bd.read_jaspar(args.motifs, args.motifId))
            sites = bd.annotate_motifs(sites, args.fasta, pwm, bd.score_threshold(pwm, args.motifPvalue))
        else:
            sites['orientation'] = ''
        nd = nodes.set_index('node')
        left_a = et.a.values <= et.b.values
        la = np.where(left_a, et.node_a.values, et.node_b.values)
        rb = np.where(left_a, et.node_b.values, et.node_a.values)
        loops = pd.DataFrame({'chrom': et.chrom.values, 'left_node': la, 'right_node': rb,
                              'left_start': nd.start.loc[la].values, 'left_end': nd.end.loc[la].values,
                              'right_start': nd.start.loc[rb].values, 'right_end': nd.end.loc[rb].values,
                              'left_mid': np.minimum(et.a.values, et.b.values),
                              'right_mid': np.maximum(et.a.values, et.b.values), 'weight': et.weight.values})
        loops['cls'], loops['strength'] = bd.classify_loops(loops, sites)
        bed = sites.assign(start=sites.summit - 50, end=sites.summit + 50,
                           name=np.where(sites.cohesin_bound, 'CTCF_cohesin', 'CTCF'),
                           strand=sites.orientation.replace({'both': '.', '': '.', '+': '+', '-': '-'}))
        bed[['chrom', 'start', 'end', 'name', 'ctcf', 'strand']].to_csv(
            ensure_parent(args.outputPrefix + '_ctcf_sites.bed'), sep='\t', header=False, index=False,
            float_format='%.3g')
        write_table(loops, args.outputPrefix + '_loops.tsv')
    ev = bd.bin_evidence(bins, sites, ins, loops, dip_prominence=args.dipProminence, dip_max=args.dipMax,
                         ctcf_min=args.ctcfMin, window_bins=args.window)
    raw = ins.drop(columns=['covered', 'insulation']).rename(columns={'insulation_w': 'capture_insulation_w',
                                                                      'insulation_global': 'capture_insulation_global'})
    table = pd.concat([bins, raw, ev], axis=1)
    codes = bd.chrom_codes(bins)
    rules = list(bd.RULES) if sites is not None else ['dip']
    rule = args.rule or ('dip_and_ctcf_pattern' if sites is not None else 'dip')
    if rule not in rules:
        raise SystemExit('rule {} needs CTCF input (--ctcf)'.format(rule))
    for r in rules:
        flags, score = bd.rule_flags(ev, r)
        table['call_' + r] = np.isin(np.arange(len(table)), bd.calls_from_flags(flags, score, codes))
        table['score_' + r] = score
    calls = np.flatnonzero(table['call_' + rule].values)
    table['score'] = table['score_' + rule]
    table['boundary'] = np.isin(np.arange(len(table)), calls)
    write_table(table, args.outputPrefix + '_bins.tsv')
    cov = table[table.covered]
    cov[['chrom', 'start', 'end', 'capture_insulation']].to_csv(
        ensure_parent(args.outputPrefix + '_capture_insulation.bedGraph'), sep='\t', header=False, index=False,
        float_format='%.4g')
    b = table.iloc[calls]
    b[['chrom', 'start', 'end']].assign(name='boundary', score=b.score.values).to_csv(
        args.outputPrefix + '_boundaries.bed', sep='\t', header=False, index=False, float_format='%.4g')
    dom = []
    for chrom, g in b.groupby('chrom', sort=False):
        mids = ((g.start + g.end) // 2).values
        for s, e in zip(mids[:-1], mids[1:]):
            dom.append((chrom, int(s), int(e)))
    pd.DataFrame(dom, columns=['chrom', 'start', 'end']).to_csv(args.outputPrefix + '_domains.bed', sep='\t',
                                                                header=False, index=False)
    log.info('%d boundaries', len(calls))
    if args.reference:
        ref = pd.read_csv(args.reference, sep='\t', header=None, usecols=[0, 1, 2], comment='#')
        ref = ref[ref[0].isin(sizes)]
        pos = {c: ((g[1] + g[2]) // 2).values for c, g in ref.groupby(0)}
        rb = bd.positions_to_bins(bins, pos)
        rows = []
        score = np.nan_to_num(table.score.values, nan=-np.inf)
        for region, mask in (('all', None), ('covered', table.covered.values.astype(bool))):
            sel = np.ones(len(table), dtype=bool) if mask is None else mask
            # AUPRC of the rule's score against the bins holding a reference boundary (does not depend on tolerance)
            auprc = bd.average_precision(score[sel], bd.bin_labels(len(table), rb)[sel])
            for tol in (0, 1, 2):
                m = bd.match_boundaries(calls, rb, codes, tolerance=tol, mask=mask)
                rows.append({'region': region, 'tolerance_bins': tol, **m, 'auprc': auprc})
        write_table(pd.DataFrame(rows), args.outputPrefix + '_evaluation.tsv')


if __name__ == '__main__':
    main()
