"""Annotation of graph nodes with chromatin signal and gene expression."""

import logging

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def bigwig_node_signal(graph, bigwig_path, flank=0):
    """Mean bigWig signal over every node fragment (optionally extended by `flank` bp). Missing data -> 0."""
    import pyBigWig
    values = {}
    handle = pyBigWig.open(bigwig_path)
    try:
        chrom_sizes = handle.chroms()
        for node, data in graph.nodes(data=True):
            chrom = data['chrom']
            if chrom not in chrom_sizes:
                values[node] = 0.0
                continue
            start = max(0, int(data['start']) - flank)
            end = min(chrom_sizes[chrom], int(data['end']) + flank)
            if end <= start:
                values[node] = 0.0
                continue
            stat = handle.stats(chrom, start, end, type='mean', exact=True)[0]
            values[node] = float(stat) if stat is not None else 0.0
    finally:
        handle.close()
    return values


def bigwig_interval_signal(bigwig_path, chroms, starts, ends):
    """Mean bigWig signal over arbitrary intervals (0 outside the chromosome or without data)."""
    import pyBigWig
    values = np.zeros(len(starts))
    handle = pyBigWig.open(bigwig_path)
    try:
        sizes = handle.chroms()
        for i, (chrom, start, end) in enumerate(zip(chroms, starts, ends)):
            if chrom not in sizes:
                continue
            start = max(0, int(start))
            end = min(sizes[chrom], int(end))
            if end <= start:
                continue
            stat = handle.stats(chrom, start, end, type='mean', exact=True)[0]
            values[i] = float(stat) if stat is not None else 0.0
    finally:
        handle.close()
    return values


def bigwig_interval_max_window(bigwig_path, chroms, starts, ends, window=1000):
    """Maximum over windows of the mean bigWig signal within every interval: the interval is divided into
    ceil(length / `window`) windows of equal width (at most `window` bp each) and the largest window mean is
    returned (the interval mean for intervals up to `window` bp; 0 outside the chromosome or without data)."""
    import pyBigWig
    values = np.zeros(len(starts))
    handle = pyBigWig.open(bigwig_path)
    try:
        sizes = handle.chroms()
        for i, (chrom, start, end) in enumerate(zip(chroms, starts, ends)):
            if chrom not in sizes:
                continue
            start = max(0, int(start))
            end = min(sizes[chrom], int(end))
            if end <= start:
                continue
            bins = max(1, min(end - start, int(np.ceil((end - start) / float(window)))))
            stat = handle.stats(chrom, start, end, type='mean', nBins=bins, exact=True)
            found = [float(s) for s in stat if s is not None]
            values[i] = max(found) if found else 0.0
    finally:
        handle.close()
    return values


def annotate_signals(graph, tracks, flank=0):
    """Attach one attribute per track. `tracks` maps attribute name -> list of bigWig replicates."""
    for name, paths in tracks.items():
        per_replicate = [bigwig_node_signal(graph, path, flank) for path in paths]
        for node in graph.nodes:
            graph.nodes[node][name] = float(np.mean([rep[node] for rep in per_replicate]))
        log.info('annotated %s from %d file(s)', name, len(paths))


def parse_track_arguments(entries):
    """Parse 'name=file1,file2' command line entries into an ordered dict."""
    tracks = {}
    for entry in entries or []:
        if '=' not in entry:
            raise ValueError('Track must be given as name=file[,file...]: {}'.format(entry))
        name, files = entry.split('=', 1)
        tracks[name] = [f for f in files.split(',') if f]
    return tracks


def gene_values_to_nodes(graph, gene_values, key='gene_ids', aggregate='max_abs'):
    """Map a gene -> value Series onto bait nodes by their annotated genes.

    `key` selects the node attribute holding comma separated identifiers
    ('gene_ids' or 'genes'). Several genes on one fragment are aggregated by the
    value with the largest magnitude (default) or the mean.
    """
    result = {}
    lookup = gene_values.to_dict()
    for node, data in graph.nodes(data=True):
        identifiers = [g for g in str(data.get(key, '')).split(',') if g]
        values = [lookup[g] for g in identifiers if g in lookup and not pd.isna(lookup[g])]
        if not values:
            continue
        if aggregate == 'mean':
            result[node] = float(np.mean(values))
        else:
            result[node] = float(values[int(np.argmax(np.abs(values)))])
    return result


def log2_ratio(a, b, pseudo_count):
    return np.log2((np.asarray(b, dtype=float) + pseudo_count) / (np.asarray(a, dtype=float) + pseudo_count))
