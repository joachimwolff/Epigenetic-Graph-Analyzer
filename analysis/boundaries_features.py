"""Per-bin inputs of the boundary rules for every data set (27_boundaries.sh, after boundaries_hic.py).

Data sets (graph, CTCF and cohesin ChIP, genome):
  GM12878  GM_chicago graph (hg19, nodes lifted to GRCh38 by midpoint), ENCODE CTCF + RAD21 (GRCh38)
  K562     K562 graph (GRCh38), ENCODE CTCF + RAD21
  H1       H1_jung graph (hg19, lifted), ENCODE CTCF + RAD21
  rubin_D0/D3/D6  keratinocyte day graphs (hg19), CTCF (one sample per day) and SMC1A (two donors, averaged), hg19
CTCF sites: summits of the CTCF signal (100 bp windows, smoothed over 300 bp, local maximum within +-500 bp,
>= 10 x the genome-wide mean and >= 4 x the mean over +-10 kb), cohesin enrichment = maximum of the cohesin signal
within +-200 bp of the summit over its genome-wide mean (cohesin-bound: >= 3 and >= 2 x its mean over +-10 kb),
CTCF motif = JASPAR 2024 MA0139.2, best hit within +-150 bp on either strand at P < 1e-4 (exact score distribution,
uniform background). Convergent cohesin-CTCF loops among the graph edges (end = node interval +- 1 kb).
Per resolution (10, 25, 50 kb): capture insulation (mirror and global null) and the Hi-C insulation track of the same
cell / day (boundaries_hic.py) as reference.

Usage: python boundaries_features.py <raw dir> <results dir> [threads]
Output: <results>/boundaries/cache/sites_<set>.tsv.gz, tables/sites_qc.tsv, tables/loop_orientation.tsv,
        cache/loops_<set>.tsv.gz, cache/features_<set>_<res>kb.tsv.gz
"""

import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epigeneticgraphanalyzer import boundaries as bd  # noqa: E402
from epigeneticgraphanalyzer.compartments import make_bins  # noqa: E402
from epigeneticgraphanalyzer.genome import LiftOver, read_chrom_sizes  # noqa: E402

CHROMS = ['chr{}'.format(c) for c in list(range(1, 23)) + ['X']]
RESOLUTIONS = (10000, 25000, 50000)
ENC = {'GM12878': os.path.expanduser('~/data/hic-multiome/hg38_gm12878'),
       'K562': '/mnt/processed/wolffjoa/hic-multiome/hg38_k562',
       'H1': '/mnt/processed/wolffjoa/hic-multiome/hg38_h1'}
HG38_FASTA = os.path.expanduser('~/data/hic-multiome/genome/hg38.fa')


def datasets(raw, results):
    rubin = os.path.join(raw, 'rubin2017_keratinocyte', 'GSE84657_RAW')
    tc = os.path.join(results, 'timecourse_rubin2017', 'graphs')
    sets = []
    for cell, graph, asm in (('GM12878', 'GM_chicago', 'hg19'), ('K562', 'K562', 'hg38'), ('H1', 'H1', 'hg19')):
        sets.append({'name': cell, 'assembly': 'hg38', 'graph_assembly': asm,
                     'nodes': os.path.join(results, 'graphs', graph + '_nodes.tsv'),
                     'edges': os.path.join(results, 'graphs', graph + '_edges.tsv'),
                     'ctcf': [os.path.join(ENC[cell], 'CTCF.bigwig')],
                     'cohesin': [os.path.join(ENC[cell], 'RAD21.bigwig')],
                     'hic': 'hic_insulation_{}'.format(cell)})
    gsm = {0: ('47', '50', '53'), 3: ('48', '51', '54'), 6: ('49', '52', '55')}
    for day, (c, s1, s2) in gsm.items():
        sets.append({'name': 'rubin_D{}'.format(day), 'assembly': 'hg19', 'graph_assembly': 'hg19',
                     'nodes': os.path.join(tc, 'D{}_nodes.tsv'.format(day)),
                     'edges': os.path.join(tc, 'D{}_edges.tsv'.format(day)),
                     'ctcf': [os.path.join(rubin, 'GSM22472{}_CTCF_ChIPseq_Day{}.bigWig'.format(c, day))],
                     'cohesin': [os.path.join(rubin, 'GSM22472{}_SMC1A_ChIPseq_Donor1_Day{}.bigWig'.format(s1, day)),
                                 os.path.join(rubin, 'GSM22472{}_SMC1A_ChIPseq_Donor2_Day{}.bigWig'.format(s2, day))],
                     'hic': 'hic_insulation_rubin_D{}'.format(day)})
    return sets


def _read(args):
    path, sizes = args
    return path, bd.read_windows(path, sizes)


def sites_for(sets, sizes, fasta, pwm, threshold, cache, threads):
    todo = [s for s in sets if not os.path.exists(os.path.join(cache, 'sites_{}.tsv.gz'.format(s['name'])))]
    if todo:
        paths = sorted({p for s in todo for p in s['ctcf'] + s['cohesin']})
        with Pool(min(threads, len(paths))) as pool:
            windows = dict(pool.map(_read, [(p, sizes) for p in paths]))
        all_sites = []
        for s in todo:
            sites = bd.ctcf_sites(s['ctcf'], s['cohesin'], sizes, windows=windows)
            sites['set'] = s['name']
            all_sites.append(sites)
            print(s['name'], len(sites), 'sites', int(sites.cohesin_bound.sum()), 'cohesin-bound', flush=True)
        del windows
        sites = bd.annotate_motifs(pd.concat(all_sites, ignore_index=True), fasta, pwm, threshold)
        for name, g in sites.groupby('set'):
            g.drop(columns='set').to_csv(os.path.join(cache, 'sites_{}.tsv.gz'.format(name)), sep='\t', index=False,
                                         float_format='%.4g')
    return {s['name']: pd.read_csv(os.path.join(cache, 'sites_{}.tsv.gz'.format(s['name'])), sep='\t',
                                   keep_default_na=False, na_values=['nan', 'NaN'])
            for s in sets}


def lifted_nodes(nodes, source, target, chain_dir):
    """Node table in the target assembly (midpoint liftover, length kept; unmapped / other chromosome dropped)."""
    if source == target:
        return nodes.copy()
    lift = LiftOver(os.path.join(chain_dir, '{}To{}.over.chain.gz'.format(source, target.capitalize())))
    mid = ((nodes.start + nodes.end) // 2).values
    half = ((nodes.end - nodes.start) // 2).values
    new = np.full(len(nodes), -1, dtype=np.int64)
    for i, (c, m) in enumerate(zip(nodes.chrom.values, mid)):
        hit = lift.position(c, int(m))
        if hit is not None and hit[0] == c:
            new[i] = hit[1]
    out = nodes.assign(start=new - half, end=new + half)
    return out[new >= 0].reset_index(drop=True)


def loops_for(edges, nodes, sites):
    nd = nodes.set_index('node')
    left_is_a = edges.a.values <= edges.b.values
    la = np.where(left_is_a, edges.node_a.values, edges.node_b.values)
    rb = np.where(left_is_a, edges.node_b.values, edges.node_a.values)
    ends = pd.DataFrame({'chrom': edges.chrom.values,
                         'left_start': nd.start.loc[la].values, 'left_end': nd.end.loc[la].values,
                         'right_start': nd.start.loc[rb].values, 'right_end': nd.end.loc[rb].values,
                         'left_mid': np.minimum(edges.a.values, edges.b.values),
                         'right_mid': np.maximum(edges.a.values, edges.b.values),
                         'weight': edges.weight.values})
    cls, strength, ambiguous = bd.classify_loops(ends, sites, return_ambiguous=True)
    ends['cls'] = cls
    ends['strength'] = strength
    ends['ambiguous'] = ambiguous
    return ends


def orientation_qc(name, loops, sites, rng):
    """Orientation classes of edges with an oriented cohesin-bound site at both ends (ends with sites of both
    orientations excluded), against the classes expected from the orientation frequency of single sites."""
    both = loops[loops.cls.isin(['convergent', 'divergent', 'tandem_forward', 'tandem_reverse']) & ~loops.ambiguous]
    obs = both.cls.value_counts()
    s = sites[sites.cohesin_bound & sites.orientation.isin(['+', '-'])]
    p_f = float((s.orientation == '+').mean()) if len(s) else np.nan
    rows = []
    n = len(both)
    exp = {'convergent': p_f * (1 - p_f), 'divergent': (1 - p_f) * p_f, 'tandem_forward': p_f * p_f,
           'tandem_reverse': (1 - p_f) * (1 - p_f)}
    for c in ('convergent', 'divergent', 'tandem_forward', 'tandem_reverse'):
        o = int(obs.get(c, 0))
        rows.append({'set': name, 'class': c, 'n': o, 'fraction': o / n if n else np.nan,
                     'expected_fraction_single_sites': exp[c], 'n_edges_both_ends': n,
                     'n_edges': len(loops)})
    return rows


def main():
    raw, results = sys.argv[1], sys.argv[2]
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    out = os.path.join(results, 'boundaries')
    cache = os.path.join(out, 'cache')
    tables = os.path.join(out, 'tables')
    hic_dir = os.path.join(out, 'hic')
    for d in (cache, tables):
        os.makedirs(d, exist_ok=True)
    pfm = bd.read_jaspar(os.path.join(raw, 'motifs', 'JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt'),
                         'MA0139')
    pwm = bd.log_odds(pfm)
    threshold = bd.score_threshold(pwm, 1e-4)
    print('MA0139 log-odds threshold (P < 1e-4):', threshold, flush=True)
    sets = datasets(raw, results)
    sizes = {asm: {c: s for c, s in read_chrom_sizes(os.path.join(raw, 'annotation', asm + '.chrom.sizes')).items()
                   if c in CHROMS} for asm in ('hg38', 'hg19')}
    fasta = {'hg38': HG38_FASTA, 'hg19': os.path.join(raw, 'annotation', 'hg19.fa.gz')}
    sites = {}
    for asm in ('hg38', 'hg19'):
        sub = [s for s in sets if s['assembly'] == asm]
        sites.update(sites_for(sub, sizes[asm], fasta[asm], pwm, threshold, cache, threads))
    qc, orient = [], []
    rng = np.random.default_rng(0)
    for s in sets:
        st = sites[s['name']]
        cb = st[st.cohesin_bound]
        qc.append({'set': s['name'], 'n_sites': len(st), 'n_cohesin_bound': len(cb),
                   'motif_fraction': float((st.orientation != '').mean()),
                   'motif_fraction_cohesin_bound': float((cb.orientation != '').mean()),
                   'forward': int((cb.orientation == '+').sum()), 'reverse': int((cb.orientation == '-').sum()),
                   'both': int((cb.orientation == 'both').sum()),
                   'median_ctcf_enrichment': float(st.ctcf.median())})
        nodes = pd.read_csv(s['nodes'], sep='\t')
        edges = pd.read_csv(s['edges'], sep='\t', usecols=['node_a', 'node_b', 'weight'])
        nodes = lifted_nodes(nodes, s['graph_assembly'], s['assembly'], os.path.join(raw, 'annotation'))
        nodes = nodes[nodes.chrom.isin(CHROMS)]
        et = bd.edge_table(nodes, edges)
        loops = loops_for(et, nodes, st)
        loops.to_csv(os.path.join(cache, 'loops_{}.tsv.gz'.format(s['name'])), sep='\t', index=False)
        orient += orientation_qc(s['name'], loops, st, rng)
        qc[-1].update(n_edges_cis=len(et), n_convergent_loops=int((loops.cls == 'convergent').sum()))
        for res in RESOLUTIONS:
            target = os.path.join(cache, 'features_{}_{}kb.tsv.gz'.format(s['name'], res // 1000))
            if os.path.exists(target):
                continue
            bins = make_bins(sizes[s['assembly']], res)
            ins = bd.capture_insulation(bins, et)
            table = pd.concat([bins, ins], axis=1)
            hic_path = os.path.join(hic_dir, '{}_{}kb.tsv.gz'.format(s['hic'], res // 1000))
            if os.path.exists(hic_path):
                hic = pd.read_csv(hic_path, sep='\t')
                hic = hic[['chrom', 'start', 'insulation']].rename(columns={'insulation': 'hic_insulation'})
                table = table.merge(hic,
                                    on=['chrom', 'start'], how='left')
            else:
                table['hic_insulation'] = np.nan
            table.to_csv(target, sep='\t', index=False, float_format='%.5g')
            print(s['name'], res, 'covered fraction', round(float(ins.covered.mean()), 3), flush=True)
    pd.DataFrame(qc).to_csv(os.path.join(tables, 'sites_qc.tsv'), sep='\t', index=False, float_format='%.4g')
    pd.DataFrame(orient).to_csv(os.path.join(tables, 'loop_orientation.tsv'), sep='\t', index=False,
                                float_format='%.4g')
    print(pd.DataFrame(qc).to_string())
    print(pd.DataFrame(orient).to_string())


if __name__ == '__main__':
    main()
