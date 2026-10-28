"""Process-parallel map over chunks of work that read large shared data (e.g. the multi-way incidence matrix).

Workers are forked after the shared data is stored in this module (`share`), so they read the parent's arrays
through copy-on-write pages instead of receiving pickled copies. Chunks are defined by the data (fixed batch sizes,
chromosomes, items), never by the number of workers, and results come back in chunk order, so every output is
independent of the thread count. Random draws stay in the calling process.
"""

import multiprocessing
import os

MAX_THREADS = 24
_SHARED = {}


def default_threads():
    """CPUs available to this process, at most MAX_THREADS."""
    try:
        available = len(os.sched_getaffinity(0))
    except AttributeError:  # pragma: no cover - platforms without sched_getaffinity
        available = os.cpu_count() or 1
    return max(1, min(MAX_THREADS, available))


def resolve_threads(threads):
    """None or <= 0: default_threads(); else the given number."""
    return default_threads() if threads is None or threads <= 0 else int(threads)


def share(**objects):
    """Store read-only objects for the workers of the next `map_chunks` call."""
    _SHARED.update(objects)


def shared(name):
    return _SHARED[name]


def map_chunks(function, chunks, threads=1):
    """[function(chunk) for chunk in chunks], in forked worker processes when threads > 1.

    `function` must be a module-level function; it reads large inputs with `shared(name)`.
    """
    chunks = list(chunks)
    threads = min(int(threads or 1), len(chunks))
    # workers do not start pools of their own (daemonic processes cannot have children)
    if threads <= 1 or multiprocessing.current_process().daemon or \
            'fork' not in multiprocessing.get_all_start_methods():
        return [function(chunk) for chunk in chunks]
    context = multiprocessing.get_context('fork')
    with context.Pool(threads) as pool:
        return pool.map(function, chunks, chunksize=1)
