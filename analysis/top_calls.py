"""The n highest-scoring calls of an ibed file (CHiCAGO score), to compare call sets of equal size.

Only cis calls within 2 Mb are ranked, the pairs that the benchmark compares and that loop callers report.

Usage: python top_calls.py <calls.ibed> <n or reference.ibed> <output.ibed>
With a reference ibed, n is its number of calls.
"""

import os
import sys

import pandas as pd


def main():
    calls_path, n, output = sys.argv[1:4]
    calls = pd.read_csv(calls_path, sep='\t')
    calls = calls[(calls['bait_chr'] == calls['otherEnd_chr'])
                  & ((calls['otherEnd_start'] - calls['bait_start']).abs() <= 2000000)]
    n = len(pd.read_csv(n, sep='\t')) if os.path.exists(n) else int(n)
    top = calls.sort_values('score', ascending=False, kind='mergesort').head(n)
    top.to_csv(output, sep='\t', index=False)
    print('{} of {} calls, score >= {}'.format(len(top), len(calls), top['score'].min()))


if __name__ == '__main__':
    main()
