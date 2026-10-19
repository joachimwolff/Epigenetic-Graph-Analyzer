"""The minimal example (example/run_example.sh: Md / PA2 at the Hoxa locus) reproduces example/expected/.

About one minute; marked slow (deselect with -m 'not slow').
"""

import os
import shutil
import subprocess
import sys

import pandas as pd
import pytest

EXAMPLE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'example')
FILES = ('graphs.tsv', 'cascade.tsv', 'hoxa3_summary.tsv', 'reports.tsv')


@pytest.mark.slow
@pytest.mark.skipif(shutil.which('bash') is None or not os.path.isdir(os.path.join(EXAMPLE, 'data')),
                    reason='bash or example data missing')
def test_example_matches_expected(tmp_path):
    env = dict(os.environ, PYTHON=sys.executable)
    run = subprocess.run(['bash', os.path.join(EXAMPLE, 'run_example.sh'), str(tmp_path / 'output'), '2'],
                         env=env, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    assert 'Hoxa3: levels with an observation: promoter_state, step1_state, step2' in run.stdout
    for name in FILES:
        observed = pd.read_csv(tmp_path / 'output' / 'summary' / name, sep='\t')
        expected = pd.read_csv(os.path.join(EXAMPLE, 'expected', name), sep='\t')
        pd.testing.assert_frame_equal(observed, expected, check_exact=False, rtol=1e-5, obj=name)
    hoxa3 = pd.read_csv(tmp_path / 'output' / 'summary' / 'hoxa3_summary.tsv', sep='\t')
    step2 = hoxa3[hoxa3['level'] == 'step2']
    # Hoxa3 is explained two steps away, by the H3K27ac gain at the Hoxa2 promoter (as in the genome-wide run)
    assert step2['explains_declared_mean'].all()
    assert step2['strongest_declared'].str.startswith('Hoxa2,Hoxaas2 promoter H3K27ac').all()
