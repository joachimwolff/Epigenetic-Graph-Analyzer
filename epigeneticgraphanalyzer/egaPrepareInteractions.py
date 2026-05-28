import argparse
import os

import pandas as pd


EXPECTED_COLS = [
    'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'value'
]


def parse_arguments():
    parser = argparse.ArgumentParser(description='Prepare and merge interaction tables into a standard 7-column format')
    parser.add_argument('--inputs', '-i', type=str, nargs='+', required=True, help='Input interaction TSV/CSV files')
    parser.add_argument('--output', '-o', type=str, required=True, help='Output TSV path')
    parser.add_argument('--sep', type=str, default='\t', help='Output separator (default: tab)')
    return parser.parse_args()


def _read_any_table(path):
    return pd.read_csv(path, sep=None, engine='python', header=None)


def main():
    args = parse_arguments()
    frames = []
    for path in args.inputs:
        df = _read_any_table(path)
        if df.shape[1] < 7:
            raise ValueError(f'File {path} has only {df.shape[1]} columns; expected at least 7.')
        df = df.iloc[:, :7].copy()
        df.columns = EXPECTED_COLS
        frames.append(df)

    merged = pd.concat(frames, ignore_index=True)
    merged = merged.dropna(subset=['chrom1', 'start1', 'chrom2', 'start2'])

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    merged.to_csv(args.output, sep=args.sep, index=False)


if __name__ == '__main__':
    main()
