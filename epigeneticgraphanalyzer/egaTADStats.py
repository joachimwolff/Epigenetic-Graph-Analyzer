import argparse
import os

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description='Compute summary statistics for TAD intervals')
    parser.add_argument('--tadBed', '-t', required=True, help='TAD BED file (chrom, start, end)')
    parser.add_argument('--output', '-o', required=True, help='Output TSV path')
    return parser.parse_args()


def main():
    args = parse_arguments()
    df = pd.read_csv(args.tadBed, sep='\t', header=None, comment='#')
    if df.shape[1] < 3:
        raise ValueError('Expected at least 3 columns in TAD BED file.')

    df = df.iloc[:, :3].copy()
    df.columns = ['chrom', 'start', 'end']
    df['length'] = df['end'].astype(int) - df['start'].astype(int)

    summary = pd.DataFrame([
        {
            'tad_count': int(len(df)),
            'mean_length': float(df['length'].mean()),
            'median_length': float(df['length'].median()),
            'min_length': int(df['length'].min()),
            'max_length': int(df['length'].max()),
        }
    ])

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    summary.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
