import hashlib
import os
import tempfile
import unittest
from unittest import mock

import pandas as pd

from epigeneticgraphanalyzer import (
    egaCheckDownloads,
    egaFetchPromoters,
    egaPrepareInteractions,
    egaReconstructCutSites,
    egaTADStats,
)


class TestCliUtilities(unittest.TestCase):
    def test_prepare_interactions_merges_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            in1 = os.path.join(tmpdir, 'a.tsv')
            in2 = os.path.join(tmpdir, 'b.tsv')
            out = os.path.join(tmpdir, 'merged.tsv')

            pd.DataFrame([
                ['chr1', 0, 10, 'chr1', 20, 30, 5.0],
            ]).to_csv(in1, sep='\t', header=False, index=False)
            pd.DataFrame([
                ['chr2', 0, 10, 'chr2', 20, 30, 2.0],
            ]).to_csv(in2, sep='\t', header=False, index=False)

            argv = ['egaPrepareInteractions', '--inputs', in1, in2, '--output', out]
            with mock.patch('sys.argv', argv):
                egaPrepareInteractions.main()

            merged = pd.read_csv(out, sep='\t')
            self.assertEqual(len(merged), 2)
            self.assertEqual(list(merged.columns), egaPrepareInteractions.EXPECTED_COLS)

    def test_fetch_promoters_creates_bed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            annotation = os.path.join(tmpdir, 'annotation.tsv')
            out = os.path.join(tmpdir, 'promoters.bed')

            pd.DataFrame([
                {'Chromosome': 'chr1', 'Start': 1000, 'End': 2000, 'Strand': '+', 'Name': 'GeneA'},
                {'Chromosome': 'chr1', 'Start': 4000, 'End': 5000, 'Strand': '-', 'Name': 'GeneB'},
            ]).to_csv(annotation, sep='\t', index=False)

            argv = ['egaFetchPromoters', '--annotation', annotation, '--output', out, '--upstream', '100', '--downstream', '50']
            with mock.patch('sys.argv', argv):
                egaFetchPromoters.main()

            bed = pd.read_csv(out, sep='\t', header=None)
            self.assertEqual(len(bed), 2)
            self.assertEqual(bed.iloc[0, 3], 'GeneA')
            self.assertEqual(bed.iloc[1, 3], 'GeneB')

    def test_reconstruct_cut_sites(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            chrom = os.path.join(tmpdir, 'chrom.sizes')
            out = os.path.join(tmpdir, 'cuts.bed')
            with open(chrom, 'w', encoding='utf-8') as handle:
                handle.write('chr1\t10\n')

            argv = ['egaReconstructCutSites', '--chromSizes', chrom, '--resolution', '4', '--output', out]
            with mock.patch('sys.argv', argv):
                egaReconstructCutSites.main()

            with open(out, 'r', encoding='utf-8') as handle:
                lines = handle.read().strip().splitlines()
            self.assertEqual(len(lines), 3)
            self.assertEqual(lines[0], 'chr1\t0\t4')

    def test_tad_stats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tad = os.path.join(tmpdir, 'tad.bed')
            out = os.path.join(tmpdir, 'stats.tsv')
            with open(tad, 'w', encoding='utf-8') as handle:
                handle.write('chr1\t0\t10\nchr1\t10\t30\n')

            argv = ['egaTADStats', '--tadBed', tad, '--output', out]
            with mock.patch('sys.argv', argv):
                egaTADStats.main()

            summary = pd.read_csv(out, sep='\t')
            self.assertEqual(int(summary.loc[0, 'tad_count']), 2)
            self.assertEqual(int(summary.loc[0, 'min_length']), 10)
            self.assertEqual(int(summary.loc[0, 'max_length']), 20)

    def test_check_downloads(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            data = os.path.join(tmpdir, 'data.txt')
            manifest = os.path.join(tmpdir, 'manifest.tsv')
            out = os.path.join(tmpdir, 'report.tsv')

            with open(data, 'w', encoding='utf-8') as handle:
                handle.write('abc')

            md5 = hashlib.md5(b'abc').hexdigest()
            pd.DataFrame([{'filepath': data, 'md5': md5}]).to_csv(manifest, sep='\t', index=False)

            argv = ['egaCheckDownloads', '--manifest', manifest, '--output', out]
            with mock.patch('sys.argv', argv):
                egaCheckDownloads.main()

            report = pd.read_csv(out, sep='\t')
            self.assertEqual(int(report.loc[0, 'exists']), 1)
            self.assertEqual(int(report.loc[0, 'match']), 1)


if __name__ == '__main__':
    unittest.main()
