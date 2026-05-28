import os
import subprocess
import sys
import unittest


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
BIN_DIR = os.path.join(REPO_ROOT, 'bin')


class TestBinWrappers(unittest.TestCase):
    def test_all_bin_wrappers_help(self):
        wrappers = sorted(
            os.path.join(BIN_DIR, name)
            for name in os.listdir(BIN_DIR)
            if os.path.isfile(os.path.join(BIN_DIR, name))
        )

        self.assertGreater(len(wrappers), 0)

        env = os.environ.copy()
        env['PYTHONPATH'] = REPO_ROOT + os.pathsep + env.get('PYTHONPATH', '')

        for wrapper in wrappers:
            with self.subTest(wrapper=os.path.basename(wrapper)):
                result = subprocess.run(
                    [sys.executable, wrapper, '--help'],
                    cwd=REPO_ROOT,
                    env=env,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=60,
                )
                self.assertEqual(
                    result.returncode,
                    0,
                    msg=f'Wrapper failed: {wrapper}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}',
                )
                self.assertIn('usage', result.stdout.lower())


if __name__ == '__main__':
    unittest.main()
