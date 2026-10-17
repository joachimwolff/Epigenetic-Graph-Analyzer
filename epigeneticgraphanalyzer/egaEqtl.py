"""egaEqtl: are eQTL variant-gene pairs connected in the promoter interaction graph at 1, 2 or 3 steps more often
than distance-matched non-eQTL pairs?

The eQTL tables (`--leads`: the tested variant-gene pairs, e.g. the lead variant per gene or fine-mapped variants;
`--significant`: all significant pairs of the data set, used to exclude eQTL pairs from the controls) are
tab-separated with the columns chrom, pos (1-based) and gene_id, in the assembly of the graph.

Outputs <prefix>_mapping.tsv (pairs kept and excluded per step, graph mapping), <prefix>_pairs.tsv.gz (eQTL pairs
and controls with their graph features), <prefix>_fractions.tsv, <prefix>_odds.tsv and <prefix>_models.tsv
(column `graph`: observed, or mirrored with --mirror).
"""

import logging

import numpy as np
import pandas as pd

from ._cli import base_parser, setup_logging, write_table
from .eqtl import (FragmentMap, add_features, conditional_models, fraction_table, gene_promoters,
                   graph_on_fragments, lift_fragments, logistic_models, map_pairs, odds_table, overlay_nodes,
                   read_graph_tables, read_rmap, sample_controls, significant_pairs, tile_bins, variant_fragments)
from .genome import read_chrom_sizes
from .io import read_tss_from_gtf

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser('Test eQTL variant-gene pairs for graph paths of 1, 2 and 3 steps against distance-matched '
                         'non-eQTL pairs.')
    parser.add_argument('--leads', required=True, help='Tested eQTL pairs (chrom, pos, gene_id)')
    parser.add_argument('--significant', required=True, help='All significant pairs of the data set (chrom, pos, '
                                                             'gene_id); also the source of control fragments')
    parser.add_argument('--nodes', required=True, help='<graph>_nodes.tsv of egaCreateGraph')
    parser.add_argument('--edges', required=True, help='<graph>_edges.tsv of egaCreateGraph')
    parser.add_argument('--rmap', help='Restriction fragments (rmap or BED, 1-based closed)')
    parser.add_argument('--liftRmap', help='Chain file to lift the --rmap fragments to the graph assembly')
    parser.add_argument('--binSize', type=int, help='Fixed bins as fragments (virtual capture) instead of --rmap')
    parser.add_argument('--chromSizes', help='Chromosome sizes for --binSize')
    parser.add_argument('--baits', help='Baited fragments (baitmap or BED) in the graph assembly')
    parser.add_argument('--baitGtf', help='Baits are the fragments holding a TSS of --baitGeneTypes in this GTF')
    parser.add_argument('--baitGeneTypes', default='protein_coding')
    parser.add_argument('--gtf', required=True, help='Gene annotation (GTF) of the graph assembly for the promoters')
    parser.add_argument('--geneTypes', default='protein_coding,lincRNA,antisense')
    parser.add_argument('--maxDistance', type=int, default=1000000)
    parser.add_argument('--controls', type=int, default=5, help='Controls per eQTL pair and control type')
    parser.add_argument('--window', type=float, default=0.1, help='Distance matching window, log10 units')
    parser.add_argument('--lengthWindow', type=float, default=0.15,
                        help='Fragment length matching window of control fragments, log10 units (negative: off)')
    parser.add_argument('--noBaitMatching', action='store_true',
                        help='Do not match control fragments on the bait status of the eQTL fragment')
    parser.add_argument('--bootstrap', type=int, default=1000, help='Block-bootstrap resamples for fractions, odds')
    parser.add_argument('--modelBootstrap', type=int, default=200, help='Block-bootstrap resamples for the models')
    parser.add_argument('--mirror', action='store_true', help='Also evaluate the mirrored control graph')
    parser.add_argument('--label', default='', help='Label written into every output table')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--outputPrefix', '-o', required=True)
    return parser.parse_args(args)


def read_pairs(path):
    df = pd.read_csv(path, sep='\t', usecols=['chrom', 'pos', 'gene_id'], dtype={'chrom': str, 'gene_id': str})
    df['pos'] = df['pos'].astype(np.int64)
    return df


def build_fragments(args, nodes):
    if args.binSize:
        fragments = tile_bins(read_chrom_sizes(args.chromSizes), args.binSize)
    else:
        fragments = read_rmap(args.rmap)
        if args.liftRmap:
            fragments = lift_fragments(fragments, args.liftRmap)
    return FragmentMap(overlay_nodes(fragments, nodes))


def bait_mask(args, fmap, nodes, node_index):
    bait = np.zeros(fmap.n, dtype=bool)
    if args.baits:
        design = read_rmap(args.baits)
        index = fmap.locate(design['chrom'], (design['start'] + design['end']) // 2)
        bait[index[index >= 0]] = True
    if args.baitGtf:
        tss = read_tss_from_gtf(args.baitGtf, gene_types=args.baitGeneTypes.split(','))
        index = fmap.locate(tss['chrom'], tss['tss'] + 1)
        bait[index[index >= 0]] = True
    from_graph = node_index[(node_index >= 0) & nodes['bait'].to_numpy()]
    bait[from_graph] = True
    return bait


def evaluate(table, graph, promoters, fmap, bait, args, label):
    features = add_features(table, graph, promoters, fmap, bait)
    fractions = fraction_table(features, args.bootstrap, args.seed).assign(graph=label)
    odds = odds_table(features, args.bootstrap, args.seed).assign(graph=label)
    models = pd.concat([logistic_models(features, args.modelBootstrap, args.seed), conditional_models(features)],
                       ignore_index=True).assign(graph=label)
    return features, fractions, odds, models


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    nodes, edges = read_graph_tables(args.nodes, args.edges)
    fmap = build_fragments(args, nodes)
    graph, node_index, info = graph_on_fragments(fmap, nodes, edges)
    bait = bait_mask(args, fmap, nodes, node_index)
    log.info('fragment map: %d fragments, %d baits; graph %s', fmap.n, int(bait.sum()), info)
    tss = read_tss_from_gtf(args.gtf, gene_types=args.geneTypes.split(','))
    promoters = gene_promoters(tss, fmap, bait)
    leads = read_pairs(args.leads)
    significant = pd.concat([read_pairs(args.significant), leads], ignore_index=True)
    pairs, mapping = map_pairs(leads, fmap, promoters, all_genes=set(tss['gene_id']), max_distance=args.maxDistance)
    sig = significant_pairs(significant, fmap)
    universe = sorted(set(significant['gene_id']) & set(promoters))
    candidates = variant_fragments(significant, fmap)
    table = sample_controls(pairs, promoters, universe, candidates, sig, n_controls=args.controls,
                            window=args.window, seed=args.seed,
                            fragment_length=fmap.length(np.arange(fmap.n)) if args.lengthWindow >= 0 else None,
                            length_window=args.lengthWindow, bait=None if args.noBaitMatching else bait)
    extra = dict(info, fragments=fmap.n, baits=int(bait.sum()), genes_with_baited_promoter=len(promoters),
                 control_gene_universe=len(universe), control_fragment_candidates=len(candidates),
                 eqtl_pairs_with_control_gene=int(table.loc[table['group'] == 'control_gene', 'set_id'].nunique()),
                 eqtl_pairs_with_control_fragment=int(table.loc[table['group'] == 'control_fragment',
                                                                'set_id'].nunique()),
                 control_gene_pairs=int((table['group'] == 'control_gene').sum()),
                 control_fragment_pairs=int((table['group'] == 'control_fragment').sum()))
    mapping = pd.concat([mapping, pd.DataFrame(list(extra.items()), columns=['step', 'n'])], ignore_index=True)
    write_table(mapping.assign(label=args.label), args.outputPrefix + '_mapping.tsv')
    log.info('%d eQTL pairs, %d controls', len(pairs), len(table) - len(pairs))

    graphs = [('observed', graph)]
    if args.mirror:
        graphs.append(('mirrored', graph.mirrored(fmap, bait)))
    outputs = {'fractions': [], 'odds': [], 'models': []}
    feature_tables = []
    for label, g in graphs:
        features, fractions, odds, models = evaluate(table, g, promoters, fmap, bait, args, label)
        feature_tables.append(features.assign(graph=label))
        outputs['fractions'].append(fractions)
        outputs['odds'].append(odds)
        outputs['models'].append(models)
        log.info('%s graph done', label)
    pd.concat(feature_tables, ignore_index=True).assign(label=args.label).to_csv(
        args.outputPrefix + '_pairs.tsv.gz', sep='\t', index=False)
    for name, frames in outputs.items():
        write_table(pd.concat(frames, ignore_index=True).assign(label=args.label),
                    args.outputPrefix + '_{}.tsv'.format(name))


if __name__ == '__main__':
    main()
