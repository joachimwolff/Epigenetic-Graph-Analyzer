"""Readers for called promoter capture Hi-C interactions and gene annotation.

EGA starts from interactions called by CHiCAGO or another caller (ibed, WashU,
BEDPE, published tables) or read out of a cool matrix at the baited bins. All
interaction readers return a pandas DataFrame in a common long format (one row
per bait-other pair, cis only):

    chrom, bait_start, bait_end, other_start, other_end, other_is_bait,
    distance, count, expected, score

`count` and `expected` are NaN when the source does not provide them; `score`
is the source-specific interaction strength (for example a CHiCAGO score or the
log(observed/expected) of Mifsud et al.).
"""

import logging

import numpy as np
import pandas as pd

from .genome import normalize_chrom, open_text

log = logging.getLogger(__name__)

INTERACTION_COLUMNS = ['chrom', 'bait_start', 'bait_end', 'other_start', 'other_end',
                       'other_is_bait', 'distance', 'count', 'expected', 'score']

INTERACTION_FORMATS = ('ibed', 'washu', 'bedpe', 'mifsud', 'homer', 'cool')


def _finalize(df):
    df = df.copy()
    df['chrom'] = df['chrom'].map(normalize_chrom)
    for column in ('bait_start', 'bait_end', 'other_start', 'other_end'):
        df[column] = df[column].astype(np.int64)
    if 'distance' not in df.columns or df['distance'].isna().any():
        df['distance'] = ((df['other_start'] + df['other_end']) // 2 -
                          (df['bait_start'] + df['bait_end']) // 2).abs()
    for column in ('count', 'expected', 'score'):
        if column not in df.columns:
            df[column] = np.nan
    if 'other_is_bait' not in df.columns:
        df['other_is_bait'] = False
    df['other_is_bait'] = df['other_is_bait'].astype(bool)
    return df[INTERACTION_COLUMNS].reset_index(drop=True)


def read_mifsud(path):
    """Read Mifsud et al. 2015 TS5 significant promoter-other interactions."""
    df = pd.read_csv(path, sep='\t')
    df.columns = ['bait_chrom', 'bait_start', 'bait_end', 'symbol', 'ensembl', 'expression_quartile',
                  'other_chrom', 'other_start', 'other_end', 'count', 'log_obs_exp']
    df = df[df['bait_chrom'] == df['other_chrom']]
    out = pd.DataFrame({
        'chrom': df['bait_chrom'], 'bait_start': df['bait_start'], 'bait_end': df['bait_end'],
        'other_start': df['other_start'], 'other_end': df['other_end'], 'other_is_bait': False,
        'count': df['count'].astype(float), 'score': df['log_obs_exp'].astype(float),
    })
    return _finalize(out)


def read_ibed(path):
    """Read a CHiCAGO ibed file (bait/other coordinates, N_reads, score)."""
    df = pd.read_csv(path, sep='\t', compression='infer')
    df = df[df['bait_chr'].astype(str) == df['otherEnd_chr'].astype(str)]
    other_is_bait = df['otherEnd_name'].astype(str) != '.'
    out = pd.DataFrame({
        'chrom': df['bait_chr'].astype(str), 'bait_start': df['bait_start'], 'bait_end': df['bait_end'],
        'other_start': df['otherEnd_start'], 'other_end': df['otherEnd_end'], 'other_is_bait': other_is_bait,
        'count': df['N_reads'].astype(float), 'score': df['score'].astype(float),
    })
    return _finalize(out)


def read_washu(path):
    """Read a CHiCAGO WashU export, bait first, in either layout:

    * WashU text: 'chrom,start,end<TAB>chrom,start,end<TAB>score'
    * WashU browser (long-range) track: 'chrom<TAB>start<TAB>end<TAB>chrom:start-end,score[<TAB>...]'

    Baits are the fragments listed first; an other end that is itself listed
    as a bait anywhere in the file is flagged as baited.
    """
    raw = pd.read_csv(path, sep='\t', header=None, compression='infer', dtype=str)
    if raw.shape[1] >= 4 and raw[3].str.contains(':').all():
        other_fields = raw[3].str.extract(r'^([^:]+):(\d+)-(\d+),([^,]+)$')
        df = pd.DataFrame({'bait': raw[0] + ',' + raw[1] + ',' + raw[2],
                           'other': other_fields[0] + ',' + other_fields[1] + ',' + other_fields[2],
                           'score': other_fields[3]})
    else:
        df = raw.iloc[:, :3].copy()
        df.columns = ['bait', 'other', 'score']
    bait = df['bait'].str.split(',', expand=True)
    other = df['other'].str.split(',', expand=True)
    cis = bait[0] == other[0]
    bait_keys = set(df['bait'])
    out = pd.DataFrame({
        'chrom': bait[0][cis], 'bait_start': bait[1][cis].astype(np.int64), 'bait_end': bait[2][cis].astype(np.int64),
        'other_start': other[1][cis].astype(np.int64), 'other_end': other[2][cis].astype(np.int64),
        'other_is_bait': df['other'][cis].isin(bait_keys), 'score': df['score'][cis].astype(float),
    })
    return _finalize(out)


def read_bedpe(path):
    """Read a generic BEDPE-like table: bait chrom/start/end, other chrom/start/end, score[, count]."""
    df = pd.read_csv(path, sep='\t', header=None, comment='#', compression='infer')
    if df.shape[1] < 7:
        raise ValueError('BEDPE input needs at least 7 columns, got {}'.format(df.shape[1]))
    df = df[df[0].astype(str) == df[3].astype(str)]
    out = pd.DataFrame({
        'chrom': df[0].astype(str), 'bait_start': df[1], 'bait_end': df[2],
        'other_start': df[4], 'other_end': df[5], 'other_is_bait': False,
        'score': df[6].astype(float),
    })
    if df.shape[1] >= 8:
        out['count'] = pd.to_numeric(df[7], errors='coerce')
    return _finalize(out)


def read_homer(path):
    """Read a HOMER analyzeHiC interaction table (as deposited for capture Hi-C by Assi et al. 2019).

    HOMER centres its windows on the reads of each library, so window starts differ between samples; every window
    is snapped to the fixed grid given by its peak id ('chrN-<grid start>') with the window size of the table, so
    that the windows of different samples are the same nodes. The bait is not marked in these tables; pass the
    capture design as `baits` to read_interactions. Score is -log10 of the HOMER p-value, count the interaction
    reads; the HOMER FDR is kept in column `fdr` (see --maxFdr of egaCreateGraph).
    """
    columns = {1: 'peak1', 2: 'chrom', 3: 'start1', 4: 'end1', 7: 'peak2', 8: 'chrom2', 14: 'count', 15: 'expected',
               17: 'logp', 18: 'fdr'}
    df = pd.read_csv(path, sep='\t', usecols=list(columns), compression='infer')
    df.columns = [columns[i] for i in sorted(columns)]
    df = df[df['chrom'].astype(str) == df['chrom2'].astype(str)]
    size = int((df['end1'] - df['start1']).median())
    grid1 = df['peak1'].str.rsplit('-', n=1).str[1].astype(np.int64)
    grid2 = df['peak2'].str.rsplit('-', n=1).str[1].astype(np.int64)
    out = pd.DataFrame({'chrom': df['chrom'].to_numpy(), 'bait_start': grid1.to_numpy(),
                        'bait_end': grid1.to_numpy() + size - 1, 'other_start': grid2.to_numpy(),
                        'other_end': grid2.to_numpy() + size - 1, 'count': df['count'].to_numpy(),
                        'expected': df['expected'].to_numpy(), 'score': (-df['logp'] / np.log(10)).to_numpy()})
    out = _finalize(out)
    out['fdr'] = df['fdr'].to_numpy()
    return out


def read_baits(path):
    """Baited fragments: BED (chrom, start, end[, ...]) or CHiCAGO .baitmap; returns chrom, start, end."""
    df = pd.read_csv(path, sep='\t', header=None, comment='#', compression='infer', usecols=[0, 1, 2])
    df.columns = ['chrom', 'start', 'end']
    df['chrom'] = df['chrom'].map(normalize_chrom)
    return df


def _overlaps_bait(chrom, start, end, baits):
    """Boolean array: interval i overlaps a bait of `baits` (DataFrame of read_baits)."""
    chrom = np.asarray(chrom)
    start = np.asarray(start, dtype=np.int64)
    end = np.asarray(end, dtype=np.int64)
    result = np.zeros(len(chrom), dtype=bool)
    for name, sites in baits.groupby('chrom'):
        sites = sites.sort_values('start')
        bait_start = sites['start'].to_numpy()
        bait_end = np.maximum.accumulate(sites['end'].to_numpy())
        mask = chrom == name
        # the last bait starting before the interval end must reach into the interval
        index = np.searchsorted(bait_start, end[mask], side='left') - 1
        ok = index >= 0
        hit = np.zeros(mask.sum(), dtype=bool)
        hit[ok] = bait_end[index[ok]] > start[mask][ok]
        result[mask] = hit
    return result


def apply_baits(table, baits):
    """Set bait status from a capture design: pairs are oriented bait first, pairs without a baited end are
    dropped and pairs listed in both orientations are kept once (largest score)."""
    table = table.copy()
    first = _overlaps_bait(table['chrom'], table['bait_start'], table['bait_end'], baits)
    second = _overlaps_bait(table['chrom'], table['other_start'], table['other_end'], baits)
    swap = ~first & second
    for a, b in (('bait_start', 'other_start'), ('bait_end', 'other_end')):
        table.loc[swap, [a, b]] = table.loc[swap, [b, a]].to_numpy()
    table['other_is_bait'] = first & second
    keep = first | second
    if (~keep).any():
        log.info('%d pairs without a baited fragment dropped', int((~keep).sum()))
    table = table[keep]
    # bait-bait pairs: upstream fragment first, so both orientations collapse
    both = table['other_is_bait'] & (table['other_start'] < table['bait_start'])
    for a, b in (('bait_start', 'other_start'), ('bait_end', 'other_end')):
        table.loc[both, [a, b]] = table.loc[both, [b, a]].to_numpy()
    table = table.sort_values('score', ascending=False).drop_duplicates(
        ['chrom', 'bait_start', 'bait_end', 'other_start', 'other_end'])
    return table.sort_values(['chrom', 'bait_start', 'other_start']).reset_index(drop=True)


def _cool_group(handle, path):
    """The cooler group of a .cool file or of 'file.mcool::resolutions/<n>'."""
    if '::' in path:
        return handle[path.split('::', 1)[1]]
    return handle


def read_cool(path, baits, field='count', min_value=None):
    """Bait-other pairs of a cool matrix (e.g. a fragment-level capture Hi-C or CHiCAGO score matrix).

    Every cis pixel touching a bin that overlaps a bait in `baits` (DataFrame of read_baits) becomes one pair
    with the pixel value of `field` as score (and as count for the 'count' field); with `min_value` only
    pixels at or above it are kept. Pixels between two baited bins are bait-bait pairs.
    """
    import h5py
    file_path = path.split('::', 1)[0]
    with h5py.File(file_path, 'r') as handle:
        group = _cool_group(handle, path)
        chrom_names = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in group['chroms/name'][:]])
        bin_chrom = chrom_names[group['bins/chrom'][:]]
        bin_start = group['bins/start'][:].astype(np.int64)
        bin_end = group['bins/end'][:].astype(np.int64)
        bin1 = group['pixels/bin1_id'][:]
        bin2 = group['pixels/bin2_id'][:]
        value = group['pixels/' + field][:].astype(float)
    bin_chrom = np.array([normalize_chrom(c) for c in bin_chrom])
    is_bait = np.zeros(len(bin_start), dtype=bool)
    for chrom, sites in baits.groupby('chrom'):
        on_chrom = np.nonzero(bin_chrom == chrom)[0]
        if len(on_chrom) == 0:
            continue
        starts = bin_start[on_chrom]
        first = np.searchsorted(starts, sites['start'].to_numpy(), side='right') - 1
        last = np.searchsorted(starts, sites['end'].to_numpy(), side='left')
        for lo, hi in zip(np.maximum(first, 0), last):
            is_bait[on_chrom[lo:max(hi, lo + 1)]] = True
    keep = (bin_chrom[bin1] == bin_chrom[bin2]) & (bin1 != bin2) & (is_bait[bin1] | is_bait[bin2])
    if min_value is not None:
        keep &= value >= min_value
    bin1, bin2, value = bin1[keep], bin2[keep], value[keep]
    # the bait is listed first; for bait-bait pixels the upstream bin
    swap = ~is_bait[bin1]
    bait_bin = np.where(swap, bin2, bin1)
    other_bin = np.where(swap, bin1, bin2)
    out = pd.DataFrame({
        'chrom': bin_chrom[bait_bin], 'bait_start': bin_start[bait_bin], 'bait_end': bin_end[bait_bin],
        'other_start': bin_start[other_bin], 'other_end': bin_end[other_bin], 'other_is_bait': is_bait[other_bin],
        'score': value, 'count': value if field == 'count' else np.nan,
    })
    log.info('%s: %d bait-other pixels of %d baited bins', path, len(out), int(is_bait.sum()))
    return _finalize(out)


def read_interactions(paths, file_format, baits=None, cool_field='count', min_value=None):
    """Dispatch to the reader of `file_format`. `paths` is a list of input files.

    `baits` (DataFrame of read_baits) is required for 'cool'; for the other formats it overrides the bait
    status given by the file (see apply_baits).
    """
    if file_format == 'cool':
        if baits is None:
            raise ValueError('cool input needs the baited fragments (--baits)')
        frames = [read_cool(path, baits, field=cool_field, min_value=min_value) for path in paths]
        return pd.concat(frames, ignore_index=True)
    readers = {'mifsud': read_mifsud, 'ibed': read_ibed, 'washu': read_washu, 'bedpe': read_bedpe,
               'homer': read_homer}
    if file_format not in readers:
        raise ValueError('Unknown interaction format {}'.format(file_format))
    frames = [readers[file_format](path) for path in paths]
    table = pd.concat(frames, ignore_index=True)
    return apply_baits(table, baits) if baits is not None else table


def read_chicdiff(path):
    """Chicdiff results exported as a table (data.table::fwrite of <prefix>_results.Rds).

    Returns bait and other-end coordinates with log2FoldChange, pvalue and the weighted adjusted p-value
    (weighted_padj, or padj when absent). Other ends are Chicdiff regions and may span several fragments.
    """
    df = pd.read_csv(path, sep=None, engine='python', compression='infer')
    missing = [c for c in ('baitchr', 'baitstart', 'baitend', 'OEchr', 'OEstart', 'OEend', 'log2FoldChange')
               if c not in df.columns]
    if missing:
        raise ValueError('{} misses Chicdiff columns {}'.format(path, missing))
    padj = df['weighted_padj'] if 'weighted_padj' in df.columns else df.get('padj', pd.Series(np.nan, index=df.index))
    pvalue = df['weighted_pvalue'] if 'weighted_pvalue' in df.columns else df.get('pvalue', pd.Series(np.nan, index=df.index))
    df = df[df['baitchr'].astype(str) == df['OEchr'].astype(str)]
    return pd.DataFrame({
        'chrom': df['baitchr'].astype(str).map(normalize_chrom).to_numpy(),
        'bait_start': df['baitstart'].astype(np.int64).to_numpy(), 'bait_end': df['baitend'].astype(np.int64).to_numpy(),
        'other_start': df['OEstart'].astype(np.int64).to_numpy(), 'other_end': df['OEend'].astype(np.int64).to_numpy(),
        'log2FC': df['log2FoldChange'].astype(float).to_numpy(),
        'pvalue': pvalue.loc[df.index].astype(float).to_numpy(), 'padj': padj.loc[df.index].astype(float).to_numpy(),
    })


def read_tss_from_gtf(path, gene_types=('protein_coding', 'lincRNA', 'antisense')):
    """Transcription start sites of all transcripts of the selected gene types in a GTF file.

    Returns a DataFrame with chrom, tss, strand, gene_id (without version), gene_name, gene_type.
    """
    rows = []
    gene_types = set(gene_types) if gene_types else None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'transcript':
                continue
            attributes = {}
            for item in fields[8].split(';'):
                item = item.strip()
                if not item or ' ' not in item:
                    continue
                key, value = item.split(' ', 1)
                attributes.setdefault(key, value.strip('"'))
            gene_type = attributes.get('gene_type', attributes.get('gene_biotype'))
            if gene_types is not None and gene_type not in gene_types:
                continue
            strand = fields[6]
            tss = int(fields[3]) - 1 if strand == '+' else int(fields[4]) - 1
            rows.append((normalize_chrom(fields[0]), tss, strand,
                         attributes.get('gene_id', '').split('.')[0],
                         attributes.get('gene_name', ''), gene_type))
    tss = pd.DataFrame(rows, columns=['chrom', 'tss', 'strand', 'gene_id', 'gene_name', 'gene_type'])
    return tss.drop_duplicates(['chrom', 'tss', 'gene_id']).reset_index(drop=True)


def read_tss_bed(path):
    """TSS or promoter BED file (chrom, start, end, name[, score, strand]); the TSS is the midpoint."""
    df = pd.read_csv(path, sep='\t', header=None, comment='#', compression='infer')
    name = df[3].astype(str) if df.shape[1] > 3 else pd.Series([''] * len(df))
    strand = df[5].astype(str) if df.shape[1] > 5 else pd.Series(['.'] * len(df))
    return pd.DataFrame({
        'chrom': df[0].map(normalize_chrom), 'tss': ((df[1] + df[2]) // 2).astype(np.int64),
        'strand': strand, 'gene_id': name, 'gene_name': name, 'gene_type': '',
    })


def read_expression(path, gene_column=None, value_columns=None):
    """Read a gene expression or differential expression table.

    The gene column is detected among common names unless given. Returns a
    DataFrame indexed by gene identifier (version suffix removed).
    """
    # index_col=False keeps the first column as data when rows carry a trailing delimiter
    df = pd.read_csv(path, sep=None, engine='python', compression='infer', index_col=False)
    if gene_column is None:
        for candidate in ('gene_name', 'Gene', 'gene', 'symbol', 'gene_id', 'GeneName', 'Unnamed: 0'):
            if candidate in df.columns:
                gene_column = candidate
                break
        else:
            gene_column = df.columns[0]
    df[gene_column] = df[gene_column].astype(str).str.replace(r'^(ENS[A-Z]*G\d+)\.\d+$', r'\1', regex=True)
    df = df.drop_duplicates(gene_column).set_index(gene_column)
    if value_columns is not None:
        df = df[list(value_columns)]
    return df
