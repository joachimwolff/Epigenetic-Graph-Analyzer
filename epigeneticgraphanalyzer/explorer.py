"""egaExplorer: one self-contained, offline HTML file per comparison to search and browse a differential run.

The explorer reads the outputs of egaDifferential (<prefix>_*: gene classes, all-level explanations, elements,
communities, mark roles, states, the union graph) and, when present, the enhancer switching events
(<results>/enhancer_switching/<name>_switch_events.tsv). It writes one HTML file with inline CSS and JavaScript and the data
embedded as gzip-compressed JSON (base64), decompressed in the browser; the page makes no network request.

The text of every gene ('At a glance', evidence by level, chromatin states, element classes, hub and community) is
written by the generators of egaReport, so that the explorer and the reports read the same. The figures (arc view,
neighbourhood graph, explanation ladder) are drawn in the browser as SVG from the embedded graph.

A gene opens directly with the URL hash '#gene=NAME' (optionally '&k=2' for the steps of the neighbourhood);
'#view=genes', '#view=communities' and '#view=switching' open the browse views.
"""

import argparse
import base64
import datetime
import gzip
import json
import logging
import math
import os
import re
import sys
import time

import numpy as np
import pandas as pd

from ._cli import base_parser, setup_logging
from .egaTadCrossing import add_tad_arguments
from ._version import __version__
from . import egaReport as report
from .egaPlotDifferential import IGNORED_GENE_PREFIXES
from .egaReport import (GeneStory, Places, Report, Run, esc, explaining_nodes, gene_names, is_step, isnan, level_hops,
                        level_label, level_order, level_short, prov_legend, read_strands, sensitive_only, supported,
                        truthy)

log = logging.getLogger(__name__)

TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'templates', 'explorer')
SWITCH_COLUMNS = ('rank', 'element', 'element_type', 'element_genes', 'element_class_changed', 'genes_lost',
                  'rna_log2FC_lost', 'distance_lost', 'promoter_lost', 'genes_gained', 'rna_log2FC_gained',
                  'distance_gained', 'promoter_gained', 'concordant', 'evidence')
STATUS_CODE = {'shared': 0, 'gained': 1, 'lost': 2}
# support class of an edge (egaDifferential attribute support): stable, threshold-sensitive, not testable
SUPPORT_CODE = {'stable': 0, 'threshold-sensitive': 1, 'not testable': 2}


# --------------------------------------------------------------------------------------------------------------------
# command line

def parse_arguments(args=None):
    parser = base_parser('Self-contained offline HTML explorer of a differential run of egaDifferential: gene search '
                         'with the text of egaReport, interactive arc view and neighbourhood graph, and browse views '
                         'of the differential genes, communities and enhancer switching events.')
    parser.add_argument('--prefix', '-p', required=True, help='Output prefix of egaDifferential')
    parser.add_argument('--explainPrefix', default=None,
                        help='Output prefix of egaExplain (default: <prefix dir>/explain/<prefix name> if present)')
    parser.add_argument('--names', nargs=2, default=None, help='Condition names A B (default: from the run)')
    parser.add_argument('--labels', nargs=2, default=None, help='Readable names of the conditions A B')
    parser.add_argument('--annotation', default=None,
                        help='Gene annotation (GTF, GTF.gz or BED with strand) for upstream / downstream')
    parser.add_argument('--rechip', default=None, help='Re-ChIP signal per bait (default as egaReport)')
    parser.add_argument('--switching', default=None,
                        help='Enhancer switching events (default: <results>/enhancer_switching/<name>_switch_events'
                             '.tsv if present; "none" to leave out)')
    parser.add_argument('--maxSwitchEvents', type=int, default=10000,
                        help='Switching events embedded, by rank (the concordant ones first)')
    parser.add_argument('--differentialOnly', action='store_true',
                        help='Search index of the differential genes only (default: every baited gene of the graph)')
    parser.add_argument('--maxNodes', type=int, default=80,
                        help='Default largest neighbourhood drawn (changeable in the page)')
    parser.add_argument('--maxLog2', type=float, default=2.0, help='Colour scale limit of the log2 changes of nodes')
    parser.add_argument('--rnaThreshold', type=float, default=None, help='As egaReport')
    parser.add_argument('--chromatinThreshold', type=float, default=None, help='As egaReport')
    parser.add_argument('--log', nargs='*', default=None, help='Run logs (default: <prefix>.log)')
    parser.add_argument('--name', default=None, help='Name of the comparison (default: the prefix name)')
    add_tad_arguments(parser)
    parser.add_argument('--output', '-o', required=True,
                        help='Output HTML file, or a directory (then <dir>/<name>_explorer.html)')
    return parser.parse_args(args)


# --------------------------------------------------------------------------------------------------------------------
# helpers

GENE_SPAN = re.compile(r'<span class="gene">([^<]*)</span>')


def gene_tags(markup):
    """Gene names in the text of egaReport as <x-g> elements (shorter; the page makes them links)."""
    return GENE_SPAN.sub(r'<x-g>\1</x-g>', str(markup or ''))


def number(value, digits=4):
    """A float for the JSON with `digits` significant digits, None for missing values."""
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        return None
    if value == 0:
        return 0
    rounded = round(value, digits - 1 - int(math.floor(math.log10(abs(value)))))
    return int(rounded) if rounded == int(rounded) and abs(rounded) < 1e15 else rounded


class Capture:
    """Stands in for egaReport.Document: keeps the markup of paragraphs and headings; notes (the method text that
    is the same for every gene) are kept once in `notes` and referenced by number."""

    def __init__(self, notes):
        self.notes = notes
        self.parts = []
        self.refs = []

    def p(self, markup, cls=None):
        if cls == 'note' and self.is_shared_note(markup):
            self.ref(markup)
            return
        self.parts.append('<p{}>{}</p>'.format(' class="{}"'.format(cls) if cls else '', gene_tags(markup)))

    @staticmethod
    def is_shared_note(markup):
        return markup.startswith(('States are called', 'Element classes follow', 'Boundaries: '))

    def ref(self, markup):
        if markup not in self.notes:
            self.notes[markup] = len(self.notes)
        if self.notes[markup] not in self.refs:
            self.refs.append(self.notes[markup])

    def html(self, markup, text=None):
        self.parts.append(gene_tags(markup))

    def heading(self, title, level=3, anchor=None, cls=None):
        self.parts.append('<h4>{}</h4>'.format(gene_tags(title)))

    def section(self, *args, **kwargs):
        pass

    def take(self):
        out = {'h': ''.join(self.parts)}
        if self.refs:
            out['n'] = list(self.refs)
        self.parts, self.refs = [], []
        return out if out['h'] or self.refs else None


def default_path(prefix, folder, suffix, explicit):
    if explicit is not None:
        return None if explicit.lower() == 'none' else explicit
    results = os.path.dirname(os.path.dirname(os.path.abspath(prefix)))
    path = os.path.join(results, folder, os.path.basename(prefix) + suffix)
    return path if os.path.exists(path) else None


# --------------------------------------------------------------------------------------------------------------------
# the embedded data

class Builder:
    def __init__(self, run, places, args):
        self.run = run
        self.places = places
        self.args = args
        self.report = Report(run, argparse.Namespace(maxNodes=args.maxNodes, maxLog2=args.maxLog2), False, places)
        self.graph = run.graph
        self.node_ids = list(self.graph.nodes)
        self.index = {n: k for k, n in enumerate(self.node_ids)}
        self.notes = {}

    # graph ----------------------------------------------------------------------------------------------------------
    def graph_data(self):
        graph, run = self.graph, self.run
        chroms, chrom_index = [], {}
        classes, class_index = [], {}

        def code(table, index, value):
            value = '' if isnan(value) else str(value)
            if value not in index:
                index[value] = len(table)
                table.append(value)
            return index[value]
        elements = self.places.elements
        nodes = dict(c=[], s=[], e=[], b=[], g=[], r=[], ca=[], cb=[], t=[])
        marks = {m: ([], []) for m in run.marks}
        statuses, status_index = [], {}
        for n in self.node_ids:
            d = graph.nodes[n]
            nodes['c'].append(code(chroms, chrom_index, d.get('chrom')))
            start, end = int(d['start']), int(d['end'])
            nodes['s'].append(start)
            nodes['e'].append(end - start)
            nodes['b'].append(1 if d.get('bait') else 0)
            genes = d.get('genes', '')
            nodes['g'].append('' if isnan(genes) else str(genes))
            nodes['r'].append(number(d.get('rna_log2FC'), 3))
            info = elements.get(n)
            nodes['ca'].append(code(classes, class_index, info[0] if info else ''))
            nodes['cb'].append(code(classes, class_index, info[1] if info else ''))
            nodes['t'].append(code(statuses, status_index, info[2] if info else ''))
            for m in run.marks:
                for k, name in enumerate(run.names):
                    marks[m][k].append(number(d.get('{}_{}'.format(m, name)), 3))
        nodes['m'] = {m: list(v) for m, v in marks.items()}
        edges = dict(u=[], v=[], s=[], a=[], b=[])
        support = []
        for u, v, d in graph.edges(data=True):
            edges['u'].append(self.index[u])
            edges['v'].append(self.index[v])
            edges['s'].append(STATUS_CODE.get(d.get('status', 'shared'), 0))
            edges['a'].append(number(d.get('r_a', 0) or 0, 3))
            edges['b'].append(number(d.get('r_b', 0) or 0, 3))
            support.append(SUPPORT_CODE.get(d.get('support'), -1))
        if any(c >= 0 for c in support):
            # support class per edge (0 stable, 1 threshold-sensitive, 2 not testable)
            edges['t'] = [max(c, 0) for c in support]
        return dict(chroms=chroms, classes=classes, tss=statuses, nodes=nodes, edges=edges)

    # genes ----------------------------------------------------------------------------------------------------------
    def gene_index(self):
        """Gene name -> promoter fragments: the differential fragments carrying the name (as egaReport) and, unless
        --differentialOnly, the baited fragments of every other gene of the graph."""
        out = {}
        for node, genes in zip(self.run.classes['node'], self.run.classes['genes']):
            for name in gene_names(genes):
                out.setdefault(name, []).append(node)
        differential = set(out)
        if not self.args.differentialOnly:
            for n, d in self.graph.nodes(data=True):
                if d.get('bait'):
                    for name in gene_names(d.get('genes', '')):
                        if name not in differential:
                            out.setdefault(name, []).append(n)
        return out, differential

    def ladder(self, story):
        if not story.fragments:
            return None
        levels = sorted({lv for _, rows in story.fragments for lv in rows}, key=level_order)
        columns = []
        for node, rows in story.fragments:
            cells = {lv: [1 if r['explains'] else 0, 1 if supported(r) else 0, 1 if sensitive_only(r) else 0]
                     for lv, r in rows.items()}
            columns.append({'l': story.column_label(node), 'c': cells})
        support = {}
        _, rows = story.fragments[0]
        for lv, row in rows.items():
            if is_step(lv) and not isnan(row.get('support')):
                support[lv] = report.sig(row['support'], 2)
        return {'levels': levels, 'labels': [level_label(lv) for lv in levels], 'cols': columns, 'support': support}

    def evidence(self, story):
        if not story.fragments:
            return []
        with_support = any(not isnan(r.get('support')) and level_hops(lv) > 0 for _, rows in story.fragments
                           for lv, r in rows.items())
        out = []
        for node, rows in story.fragments:
            frame, classes = story.evidence_frame(node, rows, with_support)
            out.append({'h': gene_tags(story.fragment_heading(node)) if len(story.fragments) > 1 else '',
                        'node': self.index.get(node), 'cols': list(frame.columns),
                        'rows': [[gene_tags(v) for v in row] for row in frame.astype(str).values.tolist()],
                        'cls': classes, 'learned': [gene_tags(x) for x in story.learned_items(rows)]})
        return out

    def explaining(self, story, nodes):
        """Nodes that carry the explanation per step (egaReport.explaining_nodes): [[node, step, mark, log2]]."""
        seeds = [n for n in nodes if n in self.graph]
        rows = self.run.levels[self.run.levels['node'].isin(nodes)] if len(self.run.levels) else []
        if not seeds or not len(rows):
            return [], None
        direction = np.sign(story.rna) if not isnan(story.rna) else 1.0
        levels = rows.loc[rows['explains'], 'level'].tolist()
        found, which = explaining_nodes(self.run, seeds, levels, direction)
        out = [[self.index[n], int(which[n][0]), which[n][1], number(which[n][2], 3)] for k in sorted(found)
               for n in found[k]]
        mark = None
        for lv in sorted(levels, key=level_order, reverse=True):
            row = rows[rows['level'] == lv].iloc[0]
            if row.get('concordant_marks'):
                mark = row['concordant_marks'][0]
                break
        return out, mark

    def gene_entry(self, gene, nodes, differential):
        run, rep = self.run, self.report
        story = GeneStory(run, self.places, gene, nodes)
        cap = Capture(self.notes)
        entry = {'n': [self.index[n] for n in nodes if n in self.index], 'd': 1 if story.fragments else 0,
                 'r': number(story.rna, 3), 'glance': gene_tags(story.glance())}
        strand = self.places.strands.get(gene)
        if strand in ('+', '-'):
            entry['strand'] = strand
        if story.fragments:
            entry['ladder'] = self.ladder(story)
            entry['evidence'] = self.evidence(story)
            rep.state_block(cap, story)
            entry['states'] = cap.take()
            explaining, mark = self.explaining(story, nodes)
            if explaining:
                entry['x'] = explaining
            if mark:
                entry['mark'] = mark
        if run.tad_analysis is not None:
            rep.tad_block(cap, gene, nodes, story)
            entry['tads'] = cap.take()
        rep.partner_elements(cap, gene, nodes)
        entry['elements'] = cap.take()
        if any(n in self.graph for n in nodes):
            rep.hub_block(cap, gene, nodes, draw=False)
            entry['hub'] = cap.take()
        if run.gene_communities is not None:
            rep.gene_community_text(cap, gene, nodes)
            entry['community'] = cap.take()
        return {k: v for k, v in entry.items() if v is not None}

    def genes(self):
        index, differential = self.gene_index()
        out = {}
        start = time.time()
        for k, (gene, nodes) in enumerate(sorted(index.items())):
            out[gene] = self.gene_entry(gene, nodes, gene in differential)
            if (k + 1) % 2000 == 0:
                log.info('%d of %d genes (%.0f s)', k + 1, len(index), time.time() - start)
        return out

    # browse tables --------------------------------------------------------------------------------------------------
    def gene_table(self):
        """One row per differential promoter fragment: genes, expression change, first explaining level, explaining
        levels, best q of a step, supported (significant beyond the mirror), promoter class in A and B."""
        run = self.run
        rows = []
        levels = run.levels
        grouped = {node: part for node, part in levels.groupby('node')} if len(levels) else {}
        for node, genes, rna in zip(run.classes['node'], run.classes['genes'], run.classes['rna_log2FC']):
            part = grouped.get(node)
            explaining, first, best_q, best_level, beyond = [], '', None, '', 0
            if part is not None:
                part = part.sort_values('order')
                hits = part[part['explains']]
                explaining = [str(lv) for lv in hits['level']]
                first = explaining[0] if explaining else ''
                steps = [r for _, r in hits.iterrows() if is_step(r['level']) and not isnan(r['q'])]
                if steps:
                    top = min(steps, key=lambda r: (float(r['q']), r['order']))
                    best_q, best_level = number(top['q'], 3), str(top['level'])
                supported_levels = [str(r['level']) for _, r in hits.iterrows() if supported(r)]
                beyond = 1 if supported_levels else 0
            info = self.places.elements.get(node)
            names = gene_names(genes)
            rows.append([self.index.get(node), names[0] if names else str(genes), ','.join(names), number(rna, 3),
                         first, explaining, best_q, best_level, beyond,
                         info[0] if info else '', info[1] if info else ''])
        columns = ['node', 'gene', 'genes', 'rna', 'first', 'levels', 'q', 'q_level', 'supported', 'class_a',
                   'class_b']
        return {'cols': columns, 'rows': rows}

    def communities(self):
        run = self.run
        if run.community_transitions is None or run.gene_communities is None:
            return None
        table = run.gene_communities
        members = dict(node=[], a=[], b=[], d=[])
        for row in table.itertuples():
            if row.node not in self.index:
                continue
            members['node'].append(self.index[row.node])
            members['a'].append(None if isnan(row.community_a) else int(row.community_a))
            members['b'].append(None if isnan(row.community_b) else int(row.community_b))
            members['d'].append(1 if truthy(getattr(row, 'differential', False)) else 0)
        # communities named by the differential member with the largest expression change (as egaReport)
        names = {}
        diff = table[table['differential'].map(truthy)] if 'differential' in table else table
        diff = diff.assign(abs_fc=diff['rna_log2FC'].abs()).sort_values('abs_fc', ascending=False)
        for side in ('a', 'b'):
            column = 'community_' + side
            for community, genes in zip(diff[column], diff['genes']):
                if isnan(community) or (side, int(community)) in names:
                    continue
                label = gene_names(genes)
                if label:
                    names[(side, int(community))] = label[0]
        coherence = {}
        if run.community_expression is not None:
            for r in run.community_expression.itertuples():
                coherence[(int(r.community_a), int(r.community_b))] = (int(r.n_genes), number(r.mean_rna_log2FC, 3),
                                                                       number(r.coherence_p, 3),
                                                                       number(r.coherence_q, 3),
                                                                       number(r.rewired_fraction, 3))
        rows = []
        for r in run.community_transitions.itertuples():
            a = None if isnan(r.community_a) else int(r.community_a)
            b = None if isnan(r.community_b) else int(r.community_b)
            c = coherence.get((a, b), (None,) * 5)
            name = names.get(('a', a)) if a is not None else None
            name = name or (names.get(('b', b)) if b is not None else None) or ''
            rows.append([str(r.side), a, None if isnan(r.size_a) else int(r.size_a), str(r.transition), b,
                         None if isnan(r.size_b) else int(r.size_b), number(r.jaccard, 3), name] + list(c))
        columns = ['side', 'a', 'size_a', 'transition', 'b', 'size_b', 'jaccard', 'name', 'n_genes', 'mean_rna',
                   'p', 'q', 'rewired']
        return {'cols': columns, 'rows': rows, 'members': members}

    def switching(self, path):
        if not path or not os.path.exists(path):
            return None
        table = pd.read_csv(path, sep='\t', low_memory=False)
        total = len(table)
        if 'rank' in table:
            table = table.sort_values('rank')
        table = table.head(self.args.maxSwitchEvents)
        names = self.run.names
        columns = [c for c in SWITCH_COLUMNS if c in table] + \
            [c for c in ('element_class_{}'.format(names[0]), 'element_class_{}'.format(names[1])) if c in table]
        rows = []
        for r in table[columns].itertuples(index=False):
            row = []
            for c, v in zip(columns, r):
                if c in ('element', 'promoter_lost', 'promoter_gained'):
                    row.append(self.index.get(v))
                elif isinstance(v, (bool, np.bool_)):
                    row.append(1 if v else 0)
                elif isinstance(v, (float, np.floating)):
                    row.append(number(v, 3))
                elif isinstance(v, (int, np.integer)):
                    row.append(int(v))
                else:
                    row.append('' if isnan(v) else str(v))
            rows.append(row)
        log.info('switching: %d of %d events from %s', len(rows), total, path)
        return {'cols': [c.replace('element_class_' + names[0], 'class_a').replace('element_class_' + names[1],
                                                                                   'class_b') for c in columns],
                'rows': rows, 'total': total, 'file': os.path.basename(path)}

    def tads(self):
        """Boundaries counted as crossed (--tads), per chromosome: [position, strength rank, strong (top decile)]."""
        analysis = self.run.tad_analysis
        if analysis is None:
            return None
        b = analysis.tads.boundaries
        b = b[b['strength_pct'].fillna(1.0) >= analysis.min_strength_pct]
        out = {}
        for chrom, g in b.groupby('chrom', sort=False):
            out[chrom] = [[int(p), None if isnan(r) else int(r), 1 if (not isnan(q) and q >= 0.9) else 0]
                          for p, r, q in zip(g['pos'], g['strength_rank'], g['strength_pct'])]
        return {'chroms': out, 'n': int(len(analysis.tads)), 'file': os.path.basename(str(analysis.tads.source))}

    # meta -----------------------------------------------------------------------------------------------------------
    def meta(self, n_genes):
        run = self.run
        counts = run.edge_counts()
        cascade = {}
        if run.cascade is not None and 'genes' in run.cascade:
            cascade = {str(r.level): int(r.genes) for r in run.cascade.itertuples() if not isnan(r.genes)}
        roles = run.declared_roles()
        return {'name': run.name, 'names': list(run.names), 'labels': list(run.labels), 'kind': run.kind,
                'title': run.condition_title(), 'marks': list(run.marks),
                'roles': {m: roles.get(m, '') for m in run.marks}, 'maxHops': int(run.max_hops),
                'rnaThreshold': run.rna_threshold, 'chromatinThreshold': run.chromatin_threshold,
                'maxLog2': self.args.maxLog2, 'maxNodes': self.args.maxNodes,
                'nDifferential': int(len(run.classes)), 'nGenes': n_genes,
                'nNodes': self.graph.number_of_nodes(), 'nEdges': self.graph.number_of_edges(),
                'edgeCounts': {str(k): int(v) for k, v in counts.items()}, 'cascade': cascade,
                'levelLabels': {lv: level_label(lv) for lv in self.level_names()},
                'levelShort': {lv: level_short(lv) for lv in self.level_names()},
                'levelOrder': {lv: level_order(lv) for lv in self.level_names()},
                'version': __version__, 'date': datetime.date.today().isoformat(),
                'levelsSource': run.levels_source, 'strands': bool(self.places.strands),
                'provenance': gene_tags(prov_legend()),
                'supportDefinition': self.report.support_definition() if run.support_edges is not None else '',
                'ignored': list(IGNORED_GENE_PREFIXES),
                'notes': [gene_tags(n) for n, _ in sorted(self.notes.items(), key=lambda x: x[1])]}

    def level_names(self):
        names = set(self.run.levels['level']) if len(self.run.levels) else set()
        names |= {'promoter_marks', 'contacts', 'promoter_state', 'step1', 'step1_state'}
        names |= {'step{}'.format(k) for k in range(2, self.run.max_hops + 1)}
        return sorted(names, key=level_order)


def encode(data):
    raw = json.dumps(data, separators=(',', ':'), ensure_ascii=False, allow_nan=False).encode('utf-8')
    packed = gzip.compress(raw, compresslevel=9, mtime=0)
    return base64.b64encode(packed).decode('ascii'), len(raw)


def decode(text):
    """The embedded data of an explorer page (for tests and scripts)."""
    match = re.search(r'<script type="application/octet-stream" id="ega-data">([^<]*)</script>', text)
    return json.loads(gzip.decompress(base64.b64decode(match.group(1))).decode('utf-8'))


def search_index(text):
    match = re.search(r'<script type="application/json" id="ega-index">(.*?)</script>', text, flags=re.S)
    return json.loads(match.group(1))


def script_json(data):
    """JSON safe inside a <script> element."""
    return json.dumps(data, separators=(',', ':'), ensure_ascii=False).replace('</', '<\\/')


def render(meta, index, payload):
    with open(os.path.join(TEMPLATE_DIR, 'explorer.html')) as handle:
        template = handle.read()
    with open(os.path.join(TEMPLATE_DIR, 'explorer.css')) as handle:
        css = handle.read()
    with open(os.path.join(TEMPLATE_DIR, 'explorer.js')) as handle:
        js = handle.read()
    replacements = {'title': esc('{}: explorer'.format(meta['title'])), 'css': css, 'js': js,
                    'index': script_json(index), 'data': payload, 'heading': esc(meta['title']),
                    'version': esc(meta['version']), 'date': esc(meta['date'])}
    # one pass, so that nothing inserted is read as a placeholder
    return re.sub(r'\{\{(\w+)\}\}', lambda m: replacements.get(m.group(1), m.group(0)), template)


def build(args):
    run = Run(args.prefix, args.explainPrefix, args.names, args.rnaThreshold, args.chromatinThreshold, args.log,
              args.name)
    run.set_labels(args.labels)
    run.read_rechip(args.rechip)
    run.set_tads(args)
    log.info('%s: conditions %s, marks %s, %d differential fragments', run.name, run.names, run.marks,
             len(run.classes))
    places = Places(run.graph, elements=run.element_lookup(), names=run.names)
    builder = Builder(run, places, args)
    index, _ = builder.gene_index()
    if args.annotation:
        places.strands = read_strands(args.annotation, list(index))
        log.info('strand of %d of %d genes from %s', len(places.strands), len(index), args.annotation)
    genes = builder.genes()
    data = {'graph': builder.graph_data(), 'genes': genes, 'table': builder.gene_table(),
            'communities': builder.communities(),
            'switching': builder.switching(default_path(args.prefix, 'enhancer_switching', '_switch_events.tsv',
                                                        args.switching))}
    tads = builder.tads()
    if tads is not None:
        data['tads'] = tads
    data['meta'] = builder.meta(len(genes))
    search = {'meta': {k: data['meta'][k] for k in ('name', 'names', 'labels', 'title', 'nDifferential', 'nGenes',
                                                    'nNodes', 'nEdges', 'version', 'date')},
              'genes': [[g, e['d'], e.get('r')] for g, e in sorted(genes.items(), key=lambda x: x[0].lower())]}
    payload, raw_size = encode(data)
    log.info('embedded data: %.1f MB JSON, %.1f MB gzip + base64', raw_size / 1e6, len(payload) / 1e6)
    return render(data['meta'], search, payload), run


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    page, run = build(args)
    path = args.output
    if os.path.isdir(path) or path.endswith(os.sep) or not path.endswith('.html'):
        os.makedirs(path, exist_ok=True)
        path = os.path.join(path, '{}_explorer.html'.format(report.slug(run.name)))
    else:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w') as handle:
        handle.write(page)
    log.info('wrote %s (%.1f MB)', path, os.path.getsize(path) / 1e6)
    sys.stdout.write('{}\n'.format(path))
    return path


if __name__ == '__main__':
    main()
