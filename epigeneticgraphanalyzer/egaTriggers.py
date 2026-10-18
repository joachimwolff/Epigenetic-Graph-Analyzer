"""egaTriggers: candidate trigger transcription factors of the genes explained by egaDifferential.

For genes whose explanation includes a poised (bivalent) promoter ('poised, trigger not among the measured layers')
or whose explaining nodes lie k graph steps away, egaTriggers adds the sequence layer. Every node of the union graph
is represented by the 1 kb window with the highest active-mark signal within its fragment (--windowTracks; the whole
fragment when shorter) and scanned with all motifs of a JASPAR file (MOODS, log-odds threshold at P < --pvalue under
the nucleotide composition of all scanned windows). Outputs (prefix -o):

  <out>_motif_enrichment.tsv  per gene group (promoter / step k explained genes up or down, poised in both, promoter
                              state transitions involving a poised state) and background: motifs at the explaining
                              nodes against background windows matched for GC content and length, either the
                              non-explaining nodes at the same steps of the same genes (promoter groups: promoters of
                              the other differential genes) ('same_genes') or random nodes with active-mark signal
                              at least the median ('random_marked'; bait nodes for promoter groups). One-sided Fisher
                              test on the windows with at least one hit, BH across motifs within group and
                              background.
  <out>_gene_triggers.tsv     per explained or poised gene the ranked candidate factors: motif present at the gene's
                              promoter fragment or explaining nodes, factor (or a member of its family) changed in
                              expression (padj < --alpha), with the relation to the gene (activator- or
                              repressor-consistent), the motif evidence, the factor's graph distance to the gene and
                              the motif's enrichment q in the gene's groups, and the gene-level motif P (binomial:
                              nodes of the gene with the motif against the motif's frequency in marked windows).
                              Candidates are reported, not decided.
  <out>_explaining_nodes.tsv.gz  every node of every gene's levels with the explaining flag.
  <out>_windows.tsv.gz        scanned window per node with signal and GC content.
  <out>_scan.npz              hit counts and best scores (nodes x motifs), reused when present (--cache).
"""

import json
import logging
import os
import time

import numpy as np
import pandas as pd

from ._cli import base_parser, ensure_parent, setup_logging, write_table
from .graph import read_graph
from .triggers import (explaining_nodes, expression_lookup, fetch_sequences, gc_content, gene_groups,
                       gene_symbol_nodes, gene_triggers, matched_background, motif_enrichment, motif_genes,
                       motif_thresholds, node_windows, nucleotide_background, read_jaspar, scan_sequences, strata)

log = logging.getLogger(__name__)


def parse_arguments(args=None):
    parser = base_parser(__doc__.split('\n\n')[0], experimental=True)
    parser.add_argument('--prefix', '-p', required=True, help='Output prefix of egaDifferential (union graph, '
                        'gene classes and explanations_all_levels)')
    parser.add_argument('--names', nargs=2, default=None,
                        help='Condition names A B (default: from the prefix <A>_vs_<B>)')
    parser.add_argument('--genome', required=True, help='Genome FASTA (uncompressed or bgzip, indexed by pyfaidx)')
    parser.add_argument('--motifs', required=True, help='Motifs in JASPAR format')
    parser.add_argument('--windowTracks', nargs='*', default=[],
                        help='bigWig tracks of the active / enhancer marks (both conditions) used to place the '
                             'scanned window within each fragment')
    parser.add_argument('--window', type=int, default=1000, help='Scanned window per fragment (bp)')
    parser.add_argument('--pvalue', type=float, default=1e-4, help='Motif hit threshold (P of the log-odds score)')
    parser.add_argument('--expression', required=True, help='Expression table of all genes (DESeq2 or similar)')
    parser.add_argument('--geneColumn', default='gene_name')
    parser.add_argument('--log2FCColumn', default='log2FoldChange')
    parser.add_argument('--padjColumn', default='padj')
    parser.add_argument('--baseColumn', default='baseMean', help='Expression level column (reported only)')
    parser.add_argument('--alpha', type=float, default=0.05, help='Factor expression change: padj below')
    parser.add_argument('--maxHops', type=int, default=3)
    parser.add_argument('--chromatinThreshold', type=float, default=1.0)
    parser.add_argument('--roles', default='declared', choices=['declared', 'learned'])
    parser.add_argument('--presence', default='mean', choices=['mean', 'max1kb'])
    parser.add_argument('--minGenes', type=int, default=5, help='Smallest gene group tested')
    parser.add_argument('--minWindows', type=int, default=10, help='Smallest number of foreground windows tested')
    parser.add_argument('--backgroundRatio', type=float, default=3.0,
                        help='Background windows per foreground window and GC / length stratum')
    parser.add_argument('--topPerGene', type=int, default=25)
    parser.add_argument('--noFamily', action='store_true', help='Do not report family members of the factors')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--cache', action='store_true', help='Reuse <out>_windows.tsv.gz and <out>_scan.npz')
    parser.add_argument('--outPrefix', '-o', required=True)
    return parser.parse_args(args)


def open_genome(path):
    import pyfaidx
    return pyfaidx.Fasta(path, as_raw=True, sequence_always_upper=True, rebuild=False)


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    out = args.outPrefix
    ensure_parent(out)
    names = args.names or os.path.basename(args.prefix).split('_vs_')
    graph = read_graph(args.prefix + '_union.graphml')
    levels = pd.read_csv(args.prefix + '_explanations_all_levels.tsv', sep='\t')
    classes = pd.read_csv(args.prefix + '_gene_classes.tsv', sep='\t')
    motifs = read_jaspar(args.motifs)
    log.info('%d nodes, %d genes, %d motifs', graph.number_of_nodes(), len(classes), len(motifs))

    t = time.time()
    nodes = explaining_nodes(graph, levels, names, max_hops=args.maxHops,
                             chromatin_threshold=args.chromatinThreshold, roles=args.roles, presence=args.presence)
    nodes.to_csv(out + '_explaining_nodes.tsv.gz', sep='\t', index=False)
    log.info('explaining nodes: %d rows, %d explaining (%.0f s)', len(nodes), int(nodes['explaining'].sum()),
             time.time() - t)

    # windows and sequences of every graph node
    window_path, scan_path = out + '_windows.tsv.gz', out + '_scan.npz'
    if args.cache and os.path.exists(window_path):
        windows = pd.read_csv(window_path, sep='\t')
    else:
        t = time.time()
        table = pd.DataFrame([(n, d['chrom'], int(d['start']), int(d['end']), bool(d.get('bait', False)))
                              for n, d in graph.nodes(data=True)], columns=['node', 'chrom', 'start', 'end', 'bait'])
        windows = node_windows(table, args.windowTracks, window=args.window, threads=args.threads)
        genome = open_genome(args.genome)
        sequences = fetch_sequences(genome, windows)
        windows['gc'] = [gc_content(s) for s in sequences]
        windows['length'] = windows['win_end'] - windows['win_start']
        windows.to_csv(window_path, sep='\t', index=False)
        log.info('windows of %d nodes (%.0f s)', len(windows), time.time() - t)
    if args.cache and os.path.exists(scan_path):
        data = np.load(scan_path, allow_pickle=False)
        counts, best, background = data['counts'], data['best'], list(data['background'])
    else:
        if 'sequences' not in locals():
            sequences = fetch_sequences(open_genome(args.genome), windows)
        background = nucleotide_background(sequences)
        t = time.time()
        counts, best = scan_sequences(sequences, motifs, background, pvalue=args.pvalue, threads=args.threads)
        np.savez_compressed(scan_path, counts=counts, best=best, background=np.array(background),
                            motif_ids=np.array([m['id'] for m in motifs]))
        log.info('scanned %d windows with %d motifs (%.0f s), background %s', len(windows), len(motifs),
                 time.time() - t, ', '.join('{} {:.3f}'.format(b, f) for b, f in zip('ACGT', background)))
    thresholds, maxima = motif_thresholds(motifs, background, pvalue=args.pvalue)
    window_index = {n: i for i, n in enumerate(windows['node'])}
    present = counts > 0
    stratum = strata(windows['gc'].to_numpy(), windows['length'].to_numpy())
    signal = windows['win_signal'].to_numpy(dtype=float)
    marked = signal >= max(np.median(signal), 1e-12) if np.any(signal > 0) else np.ones(len(windows), dtype=bool)
    bait = windows['bait'].astype(str).isin(['True', '1', 'true']).to_numpy()

    # enrichment per gene group and background
    groups = gene_groups(classes, levels, nodes, max_hops=args.maxHops, roles=args.roles, presence=args.presence,
                         min_genes=args.minGenes)
    rng = np.random.default_rng(args.seed)
    comparison = os.path.basename(args.prefix)
    motif_ids = [m['id'] for m in motifs]
    motif_names = [m['name'] for m in motifs]
    results = []
    summary = {}
    for name, group in groups.items():
        fg = np.array(sorted(window_index[n] for n in group['fg'] if n in window_index), dtype=int)
        fg = fg[~np.isnan(windows['gc'].to_numpy()[fg])] if len(fg) else fg
        if len(fg) < args.minWindows:
            continue
        same = np.array(sorted(window_index[n] for n in group['bg_same'] if n in window_index), dtype=int)
        pool = np.flatnonzero(marked & bait) if group['kind'] == 'promoter' else np.flatnonzero(marked)
        for background_name, candidates in (('same_genes', same), ('random_marked', pool)):
            bg = matched_background(fg, candidates, stratum, ratio=args.backgroundRatio, rng=rng)
            if len(bg) < args.minWindows:
                continue
            table = motif_enrichment(present, fg, bg)
            table.insert(0, 'motif', motif_names)
            table.insert(0, 'motif_id', motif_ids)
            table.insert(0, 'tf_genes', [','.join(motif_genes(n)) for n in motif_names])
            for column, value in (('fg_gc', float(np.nanmean(windows['gc'].to_numpy()[fg]))),
                                  ('bg_gc', float(np.nanmean(windows['gc'].to_numpy()[bg]))),
                                  ('n_genes', len(group['genes'])), ('background', background_name),
                                  ('description', group['description']), ('group', name),
                                  ('comparison', comparison)):
                table.insert(0, column, value)
            results.append(table)
            summary['{}/{}'.format(name, background_name)] = dict(
                genes=len(group['genes']), fg=int(len(fg)), bg=int(len(bg)),
                significant=int((table['q'] < 0.05).sum()))
    enrichment = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
    if not enrichment.empty:
        enrichment = enrichment.sort_values(['group', 'background', 'p']).reset_index(drop=True)
    write_table(enrichment, out + '_motif_enrichment.tsv')
    log.info('enrichment: %s', json.dumps(summary))

    # group q per gene and motif: the larger q of the two backgrounds, smallest over the gene's groups
    group_q = {}
    if not enrichment.empty:
        wide = enrichment.pivot_table(index=['group', 'motif_id'], columns='background', values='q', aggfunc='min')
        wide['q'] = wide.max(axis=1, skipna=False)
        column = {m: j for j, m in enumerate(motif_ids)}
        for (name, motif_id), q in wide['q'].dropna().items():
            for gene in groups[name]['genes']:
                key = (gene, column[motif_id])
                if q < group_q.get(key, np.inf):
                    group_q[key] = q

    expression = expression_lookup(pd.read_csv(args.expression, sep='\t'), args.geneColumn, args.log2FCColumn,
                                   args.padjColumn, args.baseColumn)
    triggers = gene_triggers(classes, nodes, window_index, counts, best, motifs, maxima, expression, graph,
                             alpha=args.alpha, group_q=group_q, symbol_nodes=gene_symbol_nodes(graph),
                             max_hops=args.maxHops, top=args.topPerGene, family=not args.noFamily,
                             motif_background=present[marked].mean(axis=0))
    triggers.insert(0, 'comparison', comparison)
    write_table(triggers, out + '_gene_triggers.tsv')
    log.info('candidate triggers: %d rows for %d genes', len(triggers), triggers['gene_node'].nunique()
             if len(triggers) else 0)


if __name__ == '__main__':
    main()
