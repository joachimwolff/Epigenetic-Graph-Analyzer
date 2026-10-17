"""Convert Pore-C / HiPore-C concatemers into the multi-way cluster format of egaBuildMultiway (SPRITE format).

Every read (concatemer) becomes one group: `read_id TAB chrom:pos TAB chrom:pos ...`, with pos the midpoint of each
aligned fragment, optionally lifted to another assembly with a UCSC chain (fragments that do not lift, or lift to
another chromosome, are dropped). Reads with fewer than two fragments are skipped. The input is a table with one row
per aligned fragment (parquet, or tab/comma separated, gzipped or not), or a pairs file whose rows carry the read id
(both ends of each pair are added to the read). Rows of one read need not be adjacent: reads are grouped in memory per
chunk of the input sorted by read id, so the input should be sorted or grouped by read (pore-c outputs are).

Usage: python porec_to_clusters.py <input> <output.txt.gz> --read COL --chrom COL --start COL [--end COL]
                                   [--chrom2 COL --start2 COL [--end2 COL]] [--chain hg38ToHg19.over.chain.gz]
                                   [--minMapq N --mapq COL] [--keep COL=VALUE] [--sep TAB]

--keep keeps only rows whose column has the given value, e.g. --keep pass_filter=True for pore-c-py tables.
"""

import argparse
import gzip
import sys

import pandas as pd


def read_chunks(path, columns, sep, chunksize=2000000):
    if path.endswith('.parquet') or path.endswith('.pq'):
        import pyarrow.parquet as pq
        handle = pq.ParquetFile(path)
        for batch in handle.iter_batches(batch_size=chunksize, columns=columns):
            yield batch.to_pandas()
    else:
        header = 'infer'
        comment = None
        if path.endswith(('.pairs', '.pairs.gz')):
            comment = '#'
            header = None
        for chunk in pd.read_csv(path, sep=sep, chunksize=chunksize, comment=comment, header=header,
                                 usecols=None if header is None else columns, low_memory=False):
            if header is None:
                chunk.columns = [str(c) for c in chunk.columns]
                chunk = chunk[columns]
            yield chunk


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--read', required=True)
    parser.add_argument('--chrom', required=True)
    parser.add_argument('--start', required=True)
    parser.add_argument('--end')
    parser.add_argument('--chrom2')
    parser.add_argument('--start2')
    parser.add_argument('--end2')
    parser.add_argument('--mapq')
    parser.add_argument('--minMapq', type=int, default=0)
    parser.add_argument('--chain')
    parser.add_argument('--keep', help='COL=VALUE: keep only rows with this value in COL')
    parser.add_argument('--sep', default='\t')
    args = parser.parse_args()
    lift = None
    if args.chain:
        from pyliftover import LiftOver
        lift = LiftOver(args.chain)
    ends = [(args.chrom, args.start, args.end)]
    if args.chrom2:
        ends.append((args.chrom2, args.start2, args.end2))
    keep_col, keep_value = args.keep.split('=', 1) if args.keep else (None, None)
    columns = [args.read] + [c for end in ends for c in end if c] + ([args.mapq] if args.mapq else []) \
        + ([keep_col] if keep_col else [])
    pending = {}
    reads = written = dropped = 0

    def flush(keep_last=None):
        nonlocal reads
        for read_id in [r for r in pending if r != keep_last]:
            positions = sorted(set(pending.pop(read_id)))
            if len(positions) >= 2:
                out.write('{}\t{}\n'.format(read_id, '\t'.join('{}:{}'.format(c, p) for c, p in positions)))
                reads += 1

    with gzip.open(args.output, 'wt') as out:
        for chunk in read_chunks(args.input, columns, args.sep):
            if args.mapq:
                # non-numeric values are repeated header lines inside concatenated CSV exports
                chunk = chunk[pd.to_numeric(chunk[args.mapq], errors='coerce') >= args.minMapq]
            if keep_col:
                chunk = chunk[chunk[keep_col].astype(str) == keep_value]
            for chrom_col, start_col, end_col in ends:
                starts = chunk[start_col].astype('int64')
                middle = (starts + chunk[end_col].astype('int64')) // 2 if end_col else starts
                for read_id, chrom, pos in zip(chunk[args.read].astype(str), chunk[chrom_col].astype(str), middle):
                    if not chrom.startswith('chr'):
                        chrom = 'chr' + chrom
                    if lift is not None:
                        hit = lift.convert_coordinate(chrom, int(pos))
                        if not hit or hit[0][0] != chrom:
                            dropped += 1
                            continue
                        pos = hit[0][1]
                    pending.setdefault(read_id, []).append((chrom, int(pos)))
                    written += 1
            last = str(chunk[args.read].iloc[-1]) if len(chunk) else None
            flush(keep_last=last)
        flush()
    print('{:,} reads with >= 2 fragments, {:,} fragments kept, {:,} not lifted'.format(reads, written, dropped),
          file=sys.stderr)


if __name__ == '__main__':
    main()
