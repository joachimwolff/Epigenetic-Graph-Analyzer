"""egaBoundaries: domain boundaries from capture insulation, CTCF / cohesin binding and CTCF motif orientation."""

import json
import logging

import numpy as np
import pandas as pd

from . import boundaries as bd
from ._cli import base_parser, ensure_parent, setup_logging, write_table
from .compartments import make_bins
from .genome import read_chrom_sizes

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser('Capture insulation (graph edges spanning each bin against a distance-matched mirror '
                         'expectation), cohesin-bound CTCF sites with motif orientation, convergent CTCF loops among '
                         'the edges, and domain boundaries: local minima of the capture insulation, or the calls of '
                         'a logistic model (boundaries_evaluate.py writes trained models as JSON).')
    parser.add_argument('--nodes', required=True, help='Node table of a graph (egaCreateGraph *_nodes.tsv)')
    parser.add_argument('--edges', required=True, help='Edge table of the same graph (*_edges.tsv)')
    parser.add_argument('--chromSizes', required=True, help='Chromosome sizes (UCSC format)')
    parser.add_argument('--chromosomes', nargs='*', default=None, help='Restrict to these chromosomes')
    parser.add_argument('--binSize', type=int, default=25000)
    parser.add_argument('--minExpected', type=float, default=5.0,
                        help='Bins with fewer expected spanning edges are not covered (no insulation)')
    parser.add_argument('--prominence', type=float, default=0.3,
                        help='Prominence (log2) of capture insulation minima called as boundaries')
    parser.add_argument('--ctcf', nargs='*', default=[], help='CTCF ChIP bigWig(s), averaged')
    parser.add_argument('--cohesin', nargs='*', default=[], help='Cohesin (RAD21 / SMC1A / SMC3) bigWig(s)')
    parser.add_argument('--fasta', help='Genome FASTA (plain or gzipped) for the motif orientation')
    parser.add_argument('--motifs', help='JASPAR-format motif file')
    parser.add_argument('--motifId', default='MA0139', help='CTCF motif ID (prefix match)')
    parser.add_argument('--motifPvalue', type=float, default=1e-4)
    parser.add_argument('--minEnrichment', type=float, default=10.0,
                        help='CTCF summits: minimum signal over the genome-wide mean')
    parser.add_argument('--model', help='Model JSON (features, coefficients, threshold, radius_bins); the bin size '
                                        'should match the one the model was trained at')
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
    feats = bd.bin_features(bins, sites, ins, loops)
    raw = ins.drop(columns='covered').rename(columns={'insulation': 'capture_insulation',
                                                      'insulation_w': 'capture_insulation_w',
                                                      'insulation_global': 'capture_insulation_global'})
    table = pd.concat([bins, raw, feats], axis=1)     # feats: model features (insulation 0 where not covered)
    table['covered'] = feats.covered.values > 0
    codes = bd.chrom_codes(bins)
    if args.model:
        with open(args.model) as h:
            spec = json.load(h)
        model = bd.LogisticModel.from_dict(spec)
        table['score'] = model.predict(feats)
        calls = bd.peak_calls(table.score.values, codes, spec['threshold'], spec.get('radius_bins', 1))
    else:
        table['score'] = -table.capture_insulation
        calls = bd.local_minima(ins.insulation.values, codes, args.prominence)
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
        for region, mask in (('all', None), ('covered', table.covered.values.astype(bool))):
            for tol in (0, 1, 2):
                m = bd.match_boundaries(calls, rb, codes, tolerance=tol, mask=mask)
                rows.append({'region': region, 'tolerance_bins': tol, **m})
            sel = np.ones(len(table), dtype=bool) if mask is None else mask
            score = np.nan_to_num(table.score.values, nan=-np.inf)
            rows[-1]['auprc'] = bd.average_precision(score[sel], bd.bin_labels(len(table), rb)[sel])
        write_table(pd.DataFrame(rows), args.outputPrefix + '_evaluation.tsv')
