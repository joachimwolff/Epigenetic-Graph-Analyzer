import argparse
import hashlib
import os

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description='Verify file checksums from a manifest')
    parser.add_argument('--manifest', '-m', required=True, help='TSV/CSV manifest with filepath and md5 columns')
    parser.add_argument('--output', '-o', required=True, help='Output report TSV')
    parser.add_argument('--pathCol', default='filepath', help='Manifest path column name')
    parser.add_argument('--md5Col', default='md5', help='Manifest md5 column name')
    return parser.parse_args()


def compute_md5(path):
    md5 = hashlib.md5()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            md5.update(chunk)
    return md5.hexdigest()


def main():
    args = parse_arguments()
    manifest = pd.read_csv(args.manifest, sep=None, engine='python')

    if args.pathCol not in manifest.columns or args.md5Col not in manifest.columns:
        raise ValueError(f'Manifest must include columns {args.pathCol} and {args.md5Col}.')

    rows = []
    for _, row in manifest.iterrows():
        path = str(row[args.pathCol])
        expected = str(row[args.md5Col]).lower()

        exists = os.path.exists(path)
        actual = compute_md5(path) if exists else ''
        rows.append({
            'filepath': path,
            'exists': int(exists),
            'expected_md5': expected,
            'actual_md5': actual,
            'match': int(exists and expected == actual),
        })

    report = pd.DataFrame(rows)
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    report.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
