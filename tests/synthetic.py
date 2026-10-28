"""Synthetic data sets with planted structure for the test suite."""

import os

import numpy as np
import pandas as pd

CHROM_SIZES = {'chr1': 4000000, 'chr2': 3000000}
FRAGMENT = 5000


def fragments(chrom):
    size = CHROM_SIZES[chrom]
    return [(chrom, s, min(size, s + FRAGMENT) - 1) for s in range(0, size, FRAGMENT)]


def frag_id(fragment):
    return '{}.{}.{}'.format(*fragment)


def node(fragment):
    return '{}:{}-{}'.format(*fragment)


def write_chrom_sizes(path):
    with open(path, 'w') as handle:
        for chrom, size in CHROM_SIZES.items():
            handle.write('{}\t{}\n'.format(chrom, size))
    return path


def simulate_capture(directory, seed=0, bait_every=40, planted=None, n_noise=1, name='sample'):
    """CHiCAGO ibed file of called interactions: three planted partners per bait (scores 7 to 15), one weak
    noise call per bait (score 5 to 6) and promoter-promoter calls between neighbouring baits.

    Returns the ibed path, the planted (chrom, bait index, other index) pairs and the baits.
    """
    rng = np.random.default_rng(seed)
    baits = []
    for chrom in CHROM_SIZES:
        frags = fragments(chrom)
        for i in range(10, len(frags) - 10, bait_every):
            baits.append((chrom, i))
    bait_set = set(baits)
    planted_pairs = set() if planted is None else set(planted)
    if planted is None:
        for chrom, i in baits:
            n_frags = len(fragments(chrom))
            for _ in range(3):
                j = i + int(rng.choice([-1, 1]) * rng.integers(6, 60))
                if 0 <= j < n_frags and j != i and (chrom, j) not in bait_set:
                    planted_pairs.add((chrom, i, j))
    calls = {(c, i, j): float(rng.uniform(7, 15)) for c, i, j in planted_pairs}
    for chrom, i in baits:
        n_frags = len(fragments(chrom))
        for _ in range(n_noise):
            j = i + int(rng.choice([-1, 1]) * rng.integers(3, 120))
            if 0 <= j < n_frags and (chrom, j) not in bait_set:
                calls.setdefault((chrom, i, j), float(rng.uniform(5, 6)))
    for (chrom, i), (chrom2, k) in zip(baits[:-1], baits[1:]):
        if chrom == chrom2 and rng.random() < 0.5:
            calls[(chrom, i, k)] = float(rng.uniform(6, 10))
    rows = []
    for (chrom, i, j), score in sorted(calls.items()):
        bait, other = fragments(chrom)[i], fragments(chrom)[j]
        rows.append((chrom, bait[1], bait[2], 'GENE_{}'.format(i), chrom, other[1], other[2],
                     'GENE_{}'.format(j) if (chrom, j) in bait_set else '.', int(rng.integers(10, 80)),
                     round(score, 2)))
    columns = ['bait_chr', 'bait_start', 'bait_end', 'bait_name', 'otherEnd_chr', 'otherEnd_start', 'otherEnd_end',
               'otherEnd_name', 'N_reads', 'score']
    path = os.path.join(directory, name + '.ibed')
    pd.DataFrame(rows, columns=columns).to_csv(path, sep='\t', index=False)
    return path, planted_pairs, baits


def write_cool(path, interactions_ibed):
    """Fragment-level cool file (h5py, cooler layout) holding the ibed read counts as pixels."""
    import h5py
    table = pd.read_csv(interactions_ibed, sep='\t')
    bins = [f for chrom in CHROM_SIZES for f in fragments(chrom)]
    index = {(c, s): k for k, (c, s, _) in enumerate(bins)}
    a = np.array([index[(c, s)] for c, s in zip(table['bait_chr'], table['bait_start'])])
    b = np.array([index[(c, s)] for c, s in zip(table['otherEnd_chr'], table['otherEnd_start'])])
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    order = np.lexsort((hi, lo))
    with h5py.File(path, 'w') as handle:
        handle['chroms/name'] = np.array(list(CHROM_SIZES), dtype='S')
        handle['chroms/length'] = np.array(list(CHROM_SIZES.values()))
        handle['bins/chrom'] = np.array([list(CHROM_SIZES).index(c) for c, _, _ in bins])
        handle['bins/start'] = np.array([s for _, s, _ in bins])
        handle['bins/end'] = np.array([e + 1 for _, _, e in bins])
        handle['pixels/bin1_id'] = lo[order]
        handle['pixels/bin2_id'] = hi[order]
        handle['pixels/count'] = table['N_reads'].to_numpy()[order]
    return path


def write_tss_bed(path, baits):
    rows = []
    for k, (chrom, i) in enumerate(baits):
        start = i * FRAGMENT + 1000
        rows.append((chrom, start, start + 1, 'GENE{}'.format(k), 0, '+'))
    pd.DataFrame(rows).to_csv(path, sep='\t', header=False, index=False)
    return path


def write_bigwig(path, values_by_fragment, default=1.0):
    """bigWig with one value per 5 kb fragment; `values_by_fragment` maps (chrom, index) -> value."""
    import pyBigWig
    handle = pyBigWig.open(path, 'w')
    handle.addHeader(list(CHROM_SIZES.items()))
    for chrom in CHROM_SIZES:
        frags = fragments(chrom)
        starts = [f[1] for f in frags]
        ends = [min(CHROM_SIZES[chrom], f[1] + FRAGMENT) for f in frags]
        values = [float(values_by_fragment.get((chrom, i), default)) for i in range(len(frags))]
        handle.addEntries([chrom] * len(frags), starts, ends=ends, values=values)
    handle.close()
    return path


def simulate_sprite(path, hubs, bin_size=25000, n_background=40000, n_hub_clusters=30, seed=0):
    """SPRITE cluster file: distance-decaying background clusters plus clusters holding complete hubs.

    `hubs` is a list of (chrom, [positions]) that co-occur as a whole.
    """
    rng = np.random.default_rng(seed)
    lines = []
    chroms = list(CHROM_SIZES)
    for k in range(n_background):
        chrom = chroms[rng.integers(0, len(chroms))]
        size = CHROM_SIZES[chrom]
        centre = int(rng.integers(0, size))
        n = int(rng.integers(2, 7))
        positions = np.clip(centre + (rng.geometric(1 / 8.0, n) * rng.choice([-1, 1], n) * bin_size), 0, size - 1)
        lines.append('C{}\t'.format(k) + '\t'.join('{}:{}'.format(chrom, int(p)) for p in positions))
    for h, (chrom, positions) in enumerate(hubs):
        for k in range(n_hub_clusters):
            jitter = rng.integers(0, 2000, len(positions))
            extra = [int(rng.integers(0, CHROM_SIZES[chrom]))]
            lines.append('H{}_{}\t'.format(h, k) + '\t'.join(
                '{}:{}'.format(chrom, int(p + j)) for p, j in zip(list(positions) + extra, list(jitter) + [0])))
    rng.shuffle(lines)
    with open(path, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')
    return path


def simulate_gam(path, hubs, window=25000, n_tubes=300, seed=0, slice_bins=12, hub_probability=0.25):
    """GAM segregation table: every tube holds a few contiguous slices; hubs co-segregate as a whole."""
    rng = np.random.default_rng(seed)
    windows = []
    for chrom, size in CHROM_SIZES.items():
        for start in range(0, size, window):
            windows.append((chrom, start, min(size, start + window)))
    index = {(c, s): k for k, (c, s, _) in enumerate(windows)}
    table = np.zeros((len(windows), n_tubes), dtype=np.int8)
    for tube in range(n_tubes):
        for _ in range(6):
            first = int(rng.integers(0, len(windows) - slice_bins))
            hits = first + np.nonzero(rng.random(slice_bins) < 0.5)[0]
            table[hits, tube] = 1
        for chrom, positions in hubs:
            if rng.random() < hub_probability:
                for p in positions:
                    table[index[(chrom, (p // window) * window)], tube] = 1
    df = pd.DataFrame(table, columns=['tube{}'.format(t) for t in range(n_tubes)])
    df.insert(0, 'stop', [w[2] for w in windows])
    df.insert(0, 'start', [w[1] for w in windows])
    df.insert(0, 'chrom', [w[0] for w in windows])
    df.to_csv(path, sep='\t', index=False)
    return path
