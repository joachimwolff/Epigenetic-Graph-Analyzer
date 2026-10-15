"""Shared command line helpers."""

import argparse
import logging
import os

from ._version import __version__


def base_parser(description):
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--verbose', action='store_true', help='Print progress information')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    return parser


def setup_logging(verbose):
    logging.basicConfig(level=logging.INFO if verbose else logging.WARNING,
                        format='%(asctime)s %(name)s %(levelname)s: %(message)s')


def ensure_parent(path):
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    return path


def write_table(df, path):
    ensure_parent(path)
    df.to_csv(path, sep='\t', index=False)
