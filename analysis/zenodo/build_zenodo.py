"""Zenodo records of the EGA manuscript: a file list, a manifest with sizes and sha256 and a README per record, one
.tar.gz per record written directly from the file list (no copies; the archive holds the record directory with the
README and MANIFEST next to the data), and the files to upload.

Records (each below the 50 GB Zenodo default):
  01_software                 source archive of the release tag, environment files, REPRODUCE.md
  02_results                  results of all analyses (run_all.sh output and the revision analyses)
  03_processed_inputs         inputs processed by us from public raw data (CHiCAGO calls, capture design, salmon)
  04_multiway_intermediates   multi-way incidence matrices and Pore-C / HiPore-C clusters (hours of compute)

Left out on purpose (public or rebuilt in minutes by the scripts): genome FASTA files, Chicdiff .Rds intermediates,
CHiCAGO .Rds and .chinput files, all-pair score tables, rescaled time-course bigWig files, BAM files, the salmon index.

Usage: python build_zenodo.py <out_dir> [--tag v1.0.0] [--repo <EGA checkout>] [--threads 8] [--no-tar]
                              [--only RECORD ...]   (rebuild these records, keep the others of records.json)
<out_dir>/records/<record>/  MANIFEST.tsv, README.md
<out_dir>/upload/            <record>.tar.gz, <record>.README.md, SHA256SUMS, README.md, records.json
"""

import argparse
import fnmatch
import hashlib
import json
import os
import subprocess
import tarfile
from concurrent.futures import ThreadPoolExecutor

EGA = '/mnt/processed/wolffjoa/ega'
HOME_DATA = os.path.expanduser('~/data')

EXCLUDE = ['*.fa', '*.fa.fai', '*.Rds', '*.rds', '*.chinput', '*.bam', '*.bai', '*_scores.txt', 'hicx_scores.txt',
           '*.scaled.bw', '*/reviewer_controls_killed_run/*', '*/porec/*_clusters.txt.gz', '*/multiway/*.npz',
           '*.pyc', '*/__pycache__/*', '*/salmon_index*/*']

RECORDS = {
    '02_results': {
        'title': 'Results of all analyses',
        'sources': [(EGA + '/results_final', 'results_final'),
                    (EGA + '/results/revision_bioinf', 'revision_analyses'),
                    (EGA + '/results/absolute_rates', 'revision_analyses/absolute_rates')],
        'description': 'Output of analysis/run_all.sh (results_final: graphs, multi-way evaluation, differential '
                       'comparisons, explanation levels, reports, boundary crossing, further modules, supplementary '
                       'tables, figures and number_check.md) and of the revision analyses (revision_analyses: '
                       'decomposed step features, weighting variants, held-out roles, expression baselines, graph '
                       'sensitivity, Hoxa2 Hi-C test, runtime).',
    },
    '03_processed_inputs': {
        'title': 'Inputs processed from public raw data',
        'sources': [(EGA + '/raw/chicago_calls/H1_jung2019_reprocessed', 'H1_jung2019_chicago'),
                    (HOME_DATA + '/ega_aml/kasumi1/design', 'kasumi1/capture_design'),
                    (HOME_DATA + '/ega_aml/kasumi1/chicago', 'kasumi1/chicago'),
                    (HOME_DATA + '/ega_aml/kasumi1/rna_salmon', 'kasumi1/rna_salmon'),
                    (HOME_DATA + '/ega_aml/kasumi1/call_chicago.sh', 'kasumi1/call_chicago.sh'),
                    (HOME_DATA + '/ega_aml/kasumi1/map_runs.sh', 'kasumi1/map_runs.sh'),
                    (HOME_DATA + '/ega_virtual_capture/GM12878', 'virtual_capture/GM12878'),
                    (HOME_DATA + '/ega_virtual_capture/H1', 'virtual_capture/H1'),
                    (HOME_DATA + '/ega_virtual_capture/K562', 'virtual_capture/K562'),
                    (HOME_DATA + '/ega_virtual_capture/run_vc.sh', 'virtual_capture/run_vc.sh')],
        'description': 'Promoter capture Hi-C of H1-hESC (Jung et al. 2019, GSE86189) mapped with HiCUP and called '
                       'with CHiCAGO; Kasumi-1 after RUNX1-ETO knockdown (Ptasinska et al. 2019, GSE117108): capture '
                       'design of HindIII fragments with a transcription start site, CHiCAGO calls per replicate '
                       'down to a score of 1, and salmon quantification of the RNA-seq (GSE54478); virtual promoter '
                       'capture called from the Hi-C and Micro-C of GM12878, K562 and H1-hESC (4DN).',
    },
    '04_multiway_intermediates': {
        'title': 'Multi-way intermediates',
        'sources': [(EGA + '/results_final/multiway', 'multiway'),
                    (EGA + '/results_final/porec', 'porec_clusters')],
        'include_only': ['*.npz', '*_clusters.txt.gz', '*convert.log', '*_build.log'],
        'description': 'Sparse locus-by-group incidence matrices (EGA .npz) of GM12878 SPRITE, H1-hESC GAM, '
                       'GM12878 and K562 HiPore-C and GM12878 Pore-C, and the Pore-C and HiPore-C reads converted '
                       'to clusters of restriction fragments. They are rebuilt from the public data by '
                       'analysis/09_porec.sh and egaBuildMultiway, which takes several hours.',
    },
}


def excluded(path, rules):
    return any(fnmatch.fnmatch(path, r) or fnmatch.fnmatch(os.path.basename(path), r) for r in rules)


def collect(record):
    """(source file, archive path) of a record."""
    items = []
    only = record.get('include_only')
    for source, target in record['sources']:
        if os.path.isfile(source):
            items.append((source, target))
            continue
        for root, _, files in os.walk(source):
            for name in sorted(files):
                path = os.path.join(root, name)
                rel = os.path.relpath(path, source)
                if only:
                    if not any(fnmatch.fnmatch(name, r) for r in only):
                        continue
                elif excluded(path, EXCLUDE):
                    continue
                items.append((path, os.path.join(target, rel)))
    return items


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 24), b''):
            digest.update(block)
    return digest.hexdigest()


def write_tar(items, name, directory, target, threads):
    """<name>/<archive path> for every item plus README.md and MANIFEST.tsv of the record, gzip by pigz."""
    with open(target, 'wb') as out:
        pigz = subprocess.Popen(['pigz', '-p', str(threads)], stdin=subprocess.PIPE, stdout=out)
        with tarfile.open(fileobj=pigz.stdin, mode='w|', format=tarfile.PAX_FORMAT) as tar:
            for extra in ('README.md', 'MANIFEST.tsv'):
                path = os.path.join(directory, extra)
                if os.path.exists(path):
                    tar.add(path, arcname=os.path.join(name, extra))
            for source, arcname in items:
                tar.add(os.path.realpath(source), arcname=os.path.join(name, arcname), recursive=False)
        pigz.stdin.close()
        if pigz.wait():
            raise SystemExit('pigz failed for ' + name)


def write_manifest(items, directory, threads):
    with ThreadPoolExecutor(threads) as pool:
        digests = list(pool.map(sha256, [s for s, _ in items]))
    rows = [(t, os.path.getsize(s), d) for (s, t), d in zip(items, digests)]
    with open(os.path.join(directory, 'MANIFEST.tsv'), 'w') as out:
        out.write('path\tbytes\tsha256\n')
        for t, n, d in rows:
            out.write('{}\t{}\t{}\n'.format(t, n, d))
    return rows


def readme(name, record, rows, tag):
    total = sum(n for _, n, _ in rows)
    top = {}
    for t, n, _ in rows:
        key = t.split('/')[0] if '/' not in t.split('/', 1)[-1] else '/'.join(t.split('/')[:2])
        top[key] = top.get(key, 0) + n
    lines = ['# Epigenetic Graph Analyzer (EGA) {}: {}'.format(tag, record['title']), '',
             record['description'], '',
             '{} files, {:.2f} GB before compression. MANIFEST.tsv lists every file with its size and sha256.'.format(
                 len(rows), total / 1e9), '', '## Contents', '', '| Directory | GB |', '|---|---|']
    lines += ['| {} | {:.3f} |'.format(k, v / 1e9) for k, v in sorted(top.items())]
    lines += ['', 'Software, environment and the mapping of every figure and table to its script: record '
              '01_software (REPRODUCE.md). Public data and accessions: DATA_SOURCES.md in the upload.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('out')
    parser.add_argument('--tag', default='v1.0.0')
    parser.add_argument('--repo', default=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--no-tar', action='store_true')
    parser.add_argument('--only', nargs='+', help='Records to rebuild; the others are kept from records.json')
    args = parser.parse_args()
    records_dir, upload = os.path.join(args.out, 'records'), os.path.join(args.out, 'upload')
    os.makedirs(records_dir, exist_ok=True)
    os.makedirs(upload, exist_ok=True)
    summary = {}
    previous = os.path.join(upload, 'records.json')
    if args.only and os.path.exists(previous):
        summary = json.load(open(previous))
    wanted = (lambda name: not args.only or name in args.only)
    built = []
    # 01_software: source archive of the tag and the reproduction documents
    record_items = {}
    software = os.path.join(records_dir, '01_software')
    if wanted('01_software'):
        os.makedirs(software, exist_ok=True)
        archive = os.path.join(software, 'Epigenetic-Graph-Analyzer-{}.tar.gz'.format(args.tag))
        subprocess.run(['git', '-C', args.repo, 'archive', '--format=tar.gz',
                        '--prefix=Epigenetic-Graph-Analyzer-{}/'.format(args.tag), '-o', archive, args.tag], check=True)
        for name in ('REPRODUCE.md', 'environment.yml', 'CITATION.cff', 'README.md'):
            out = subprocess.run(['git', '-C', args.repo, 'show', '{}:{}'.format(args.tag, name)], check=True,
                                 capture_output=True).stdout
            open(os.path.join(software, name), 'wb').write(out)
        items = [(os.path.join(software, f), f) for f in sorted(os.listdir(software)) if f != 'MANIFEST.tsv']
        rows = write_manifest(items, software, args.threads)
        summary['01_software'] = {'files': len(rows), 'bytes': sum(n for _, n, _ in rows)}
        record_items['01_software'] = items
        built.append('01_software')
    for name, record in RECORDS.items():
        if not wanted(name):
            continue
        built.append(name)
        directory = os.path.join(records_dir, name)
        os.makedirs(directory, exist_ok=True)
        items = collect(record)
        rows = write_manifest(items, directory, args.threads)
        open(os.path.join(directory, 'README.md'), 'w').write(readme(name, record, rows, args.tag))
        summary[name] = {'files': len(rows), 'bytes': sum(n for _, n, _ in rows)}
        record_items[name] = items
        print('{}: {} files, {:.2f} GB'.format(name, len(rows), summary[name]['bytes'] / 1e9), flush=True)
    if not args.no_tar:
        for name in built:
            target = os.path.join(upload, name + '.tar.gz')
            write_tar(record_items[name], name, os.path.join(records_dir, name), target, args.threads)
            readme_path = os.path.join(records_dir, name, 'README.md')
            if os.path.exists(readme_path):
                open(os.path.join(upload, name + '.README.md'), 'w').write(open(readme_path).read())
            summary[name]['packed_bytes'] = os.path.getsize(target)
            summary[name]['sha256'] = sha256(target)
            print('{}.tar.gz: {:.2f} GB'.format(name, summary[name]['packed_bytes'] / 1e9), flush=True)
        with open(os.path.join(upload, 'SHA256SUMS'), 'w') as out:
            for name, info in summary.items():
                out.write('{}  {}.tar.gz\n'.format(info['sha256'], name))
    json.dump(summary, open(os.path.join(upload, 'records.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
