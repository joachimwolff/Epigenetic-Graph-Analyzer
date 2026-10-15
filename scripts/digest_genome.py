"""Restriction fragments of a genome as BED (1-based closed coordinates of the CHiCAGO .rmap convention).

Usage: python digest_genome.py <genome.fa(.gz)> <site, e.g. AAGCTT> <cut offset, e.g. 1 for A^AGCTT> <output.bed>
       [chromosome regex, default ^chr([0-9]+|X|Y)$]

Each fragment runs from the base after one cut to the base before the next, the first fragment of a
chromosome starts at 1 and the last ends at the chromosome length.
"""

import gzip
import re
import sys


def chromosomes(path):
    opener = gzip.open if path.endswith('.gz') else open
    name, parts = None, []
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip().upper())
    if name is not None:
        yield name, ''.join(parts)


def main():
    fasta, site, offset, output = sys.argv[1], sys.argv[2].upper(), int(sys.argv[3]), sys.argv[4]
    pattern = re.compile(sys.argv[5] if len(sys.argv) > 5 else r'^chr([0-9]+|X|Y)$')
    with open(output, 'w') as out:
        for name, sequence in chromosomes(fasta):
            if not pattern.match(name):
                continue
            cuts = [m.start() + offset for m in re.finditer('(?={})'.format(site), sequence)]
            start = 1
            for cut in cuts:
                # the fragment ends with the base before the cut (1-based position cut)
                out.write('{}\t{}\t{}\n'.format(name, start, cut))
                start = cut + 1
            out.write('{}\t{}\t{}\n'.format(name, start, len(sequence)))


if __name__ == '__main__':
    main()
