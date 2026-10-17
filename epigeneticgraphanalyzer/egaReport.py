"""egaReport: a self-contained HTML and PDF report of a differential run (egaDifferential, egaExplain).

The report only reads the outputs of egaDifferential (<prefix>_*) and, when present, of egaExplain
(<explain prefix>_*). Every number shown is taken from these files or counted from them; the figures are drawn from
<prefix>_union.graphml with the drawing functions of egaPlotDifferential.

Main report (default): (a) summary of the comparison, (b) explanation cascade (first level and any level, level
combinations, mark roles declared and learned per level, chromatin state transitions), (c) major findings: per
cascade level the genes with the strongest support, each with all explaining levels in rank order, interaction arcs in
both conditions and their change and the neighbourhood graph with the explaining nodes, (d) hubs and communities:
members changing in the same and in the opposite direction, genes that left or joined, drawn in both conditions,
(e) methods, parameters and versions.

Gene report (--gene NAME [NAME ...]): the same structure restricted to each gene, one page per gene and an index.
Community report (--community ID): one page for a community of condition B (or A with --communitySide A).

The per-gene explanation of every level comes from <prefix>_explanations_all_levels.tsv when egaDifferential wrote
it; otherwise it is derived from <prefix>_gene_classes.tsv (mark changes at the promoter, own interactions,
promoter state, step 1, step 1 state, step k), the per-gene permutation tests (<prefix>_indirect_gene_test.tsv,
<prefix>_hop<k>_gene_test.tsv) and the 3D support per step of <prefix>_explanation_cascade.tsv, with the
thresholds --rnaThreshold and --chromatinThreshold (the egaDifferential defaults unless given).

The PDF is printed by Chromium in headless mode (snap or system installation); without it the report is written
with matplotlib (PdfPages, plain text pages and the figures), which the log and the report state.
"""

import datetime
import gzip
import html
import io
import json
import logging
import os
import platform
import re
import shutil
import string
import subprocess
import sys
import tempfile
import textwrap

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import colors as mcolors  # noqa: E402
from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.patches import Arc  # noqa: E402

from ._cli import base_parser, setup_logging  # noqa: E402
from ._version import __version__  # noqa: E402
from .egaPlotDifferential import (IGNORED_GENE_PREFIXES, INK, STATUS_COLORS,  # noqa: E402
                                  draw_arcs, draw_graph, edge_norm, edge_strength, edge_width, mark_change, separate,
                                  place_labels, shape_handles, spread_layout, status_handles, width_handles)
from .egaPlotRegion import ordered_subgraph  # noqa: E402
from .graph import read_graph  # noqa: E402

log = logging.getLogger(__name__)

TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'templates')
FONT_FAMILY = ['DejaVu Sans']
# colour-blind-safe categorical colours (Okabe and Ito) for marks; interactions keep the red-blue convention
MARK_COLORS = ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#56B4E9', '#D55E00', '#000000']
BAR_COLOR = '#3d6f99'
EXPLAINING = '#6a3d9a'
FADED = '#d6d5cf'
# membership changes of hubs and communities: colour and line style differ (colour-blind safe on the red-blue fill)
LEFT_COLOR = '#E69F00'
JOINED_COLOR = '#009E73'
ROLE_SIGN = {'active': 1, 'enhancer': 1, 'repressive': -1}
FRAGMENT_COLUMNS = ('fragment', 'promoter fragment', 'node', 'path', 'genes')
CHROMIUM_CANDIDATES = ('/snap/bin/chromium', 'chromium', 'chromium-browser', 'google-chrome', 'google-chrome-stable')

matplotlib.rcParams.update({'svg.fonttype': 'none', 'font.family': 'sans-serif', 'font.sans-serif': FONT_FAMILY,
                            'svg.hashsalt': 'egaReport', 'pdf.fonttype': 42})


# --------------------------------------------------------------------------------------------------------------------
# command line

def parse_arguments(args=None):
    parser = base_parser('Self-contained HTML and PDF report of a differential run of egaDifferential (and egaExplain): '
                         'summary, explanation cascade, major findings with figures, hubs and communities, methods; '
                         'or reports of single genes (--gene) or of one community (--community).')
    parser.add_argument('--prefix', '-p', required=True, help='Output prefix of egaDifferential')
    parser.add_argument('--explainPrefix', default=None,
                        help='Output prefix of egaExplain (default: <prefix dir>/explain/<prefix name> if present)')
    parser.add_argument('--names', nargs=2, default=None,
                        help='Condition names A B (default: read from <prefix>_state_marks.tsv, the edges or the graph)')
    parser.add_argument('--gene', nargs='+', default=None, help='Gene report(s) for these gene names')
    parser.add_argument('--community', type=int, default=None, help='Community report for this community ID')
    parser.add_argument('--communitySide', choices=('A', 'B'), default='B',
                        help='Condition whose community IDs --community refers to')
    parser.add_argument('--top', type=int, default=3,
                        help='Genes shown per cascade level in the major findings (main report)')
    parser.add_argument('--maxCommunities', type=int, default=3,
                        help='Coordinated communities drawn in the main report')
    parser.add_argument('--maxNodes', type=int, default=80,
                        help='Largest neighbourhood drawn; larger neighbourhoods keep the promoter, its strongest '
                             'partners and the paths to the explaining nodes')
    parser.add_argument('--rnaThreshold', type=float, default=None,
                        help='|log2FC| of a differential gene (default: from the run log, else 1, as egaDifferential)')
    parser.add_argument('--chromatinThreshold', type=float, default=None,
                        help='|log2FC| of an explaining mark change (default: from the run log, else 1)')
    parser.add_argument('--maxLog2', type=float, default=2.0, help='Colour scale limit of the log2 changes of nodes')
    parser.add_argument('--log', nargs='*', default=None,
                        help='Run logs or scripts with the egaDifferential command line (default: <prefix>.log)')
    parser.add_argument('--name', default=None, help='Name of the comparison (default: the prefix name)')
    parser.add_argument('--pdf', choices=('auto', 'chromium', 'matplotlib', 'none'), default='auto',
                        help='PDF engine: Chromium headless (auto: if available), matplotlib, or no PDF')
    parser.add_argument('--chromium', default=None, help='Chromium executable (default: snap or PATH)')
    parser.add_argument('--outputDir', '-o', required=True, help='Output directory')
    return parser.parse_args(args)


# --------------------------------------------------------------------------------------------------------------------
# formatting

def esc(value):
    return html.escape('' if value is None else str(value), quote=True)


def isnan(value):
    try:
        return value is None or (isinstance(value, float) and np.isnan(value)) or pd.isna(value)
    except (TypeError, ValueError):
        return False


def fmt_num(value, digits=2, signed=False):
    if isnan(value):
        return ''
    try:
        value = float(value)
    except (TypeError, ValueError):
        return str(value)
    if signed:
        return '{:+.{d}f}'.format(value, d=digits)
    return '{:.{d}f}'.format(value, d=digits)


def fmt_int(value):
    if isnan(value):
        return ''
    try:
        return '{:,}'.format(int(round(float(value))))
    except (TypeError, ValueError):
        return str(value)


def fmt_p(value):
    if isnan(value):
        return ''
    value = float(value)
    if value == 0:
        return '0'
    if value < 0.001:
        return '{:.1e}'.format(value)
    return '{:.3g}'.format(value)


def fmt_pct(value):
    return '' if isnan(value) else '{:.1f}%'.format(100 * float(value))


def gene_span(name):
    return '<span class="gene">{}</span>'.format(esc(name))


def first_gene(genes):
    names = [g for g in str(genes if not isnan(genes) else '').split(',') if g and g != 'nan']
    informative = [g for g in names if not g.startswith(IGNORED_GENE_PREFIXES)]
    return (informative or names or [''])[0]


def slug(text):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(text)).strip('_') or 'x'


def truthy(value):
    if isnan(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() in ('true', '1', 'yes', 'y')
    return bool(value)


# --------------------------------------------------------------------------------------------------------------------
# levels of the explanation

def level_order(level):
    """Rank order of a level: promoter marks, contacts, promoter state, step 1, step 1 state, step 2, ..."""
    fixed = {'promoter_marks': 0, 'contacts': 1, 'promoter_state': 2, 'step1': 3, 'step1_state': 4}
    if level in fixed:
        return fixed[level]
    match = re.fullmatch(r'step(\d+)', str(level))
    if match:
        return 3 + 2 * (int(match.group(1)) - 1)
    return 1000


def level_label(level):
    labels = {'promoter_marks': 'promoter marks', 'contacts': 'contacts of the promoter',
              'promoter_state': 'promoter state', 'step1': 'step 1 (direct partners)',
              'step1_state': 'step 1 state (contacting fragments)'}
    if level in labels:
        return labels[level]
    match = re.fullmatch(r'step(\d+)', str(level))
    if match:
        return 'step {} (nodes at exactly {} steps)'.format(match.group(1), match.group(1))
    return str(level)


def cascade_level(level):
    """Cascade level (gene, step1, step2, ...) of an explanation level."""
    if level in ('promoter_marks', 'contacts', 'promoter_state'):
        return 'gene'
    if level == 'step1_state':
        return 'step1'
    return level


def level_hops(level):
    match = re.fullmatch(r'step(\d+)(_state)?', str(level))
    return int(match.group(1)) if match else 0


def normalise_level(value):
    """Level names of an all-level table mapped onto the names used here."""
    text = re.sub(r'[\s_-]+', '', str(value).lower())
    aliases = {'promoter': 'promoter_marks', 'promotermarks': 'promoter_marks', 'genemarks': 'promoter_marks',
               'marks': 'promoter_marks', 'step0': 'promoter_marks', 'contacts': 'contacts', 'contact': 'contacts',
               'interactions': 'contacts', 'owninteractions': 'contacts', 'promotercontacts': 'contacts',
               'promoterstate': 'promoter_state', 'step1states': 'step1_state',
               'state': 'promoter_state', 'genestate': 'promoter_state', 'step1': 'step1', 'direct': 'step1',
               'step1marks': 'step1', 'step1state': 'step1_state', 'indirect': 'step2'}
    if text in aliases:
        return aliases[text]
    match = re.fullmatch(r'(?:step|hop)(\d+)(marks)?', text)
    if match:
        return 'step{}'.format(int(match.group(1)))
    return str(value)


def pick(columns, *names):
    lower = {c.lower(): c for c in columns}
    for name in names:
        if name.lower() in lower:
            return lower[name.lower()]
    return None


# --------------------------------------------------------------------------------------------------------------------
# reading a run

def read_tsv(path):
    if path and os.path.exists(path):
        try:
            return pd.read_csv(path, sep='\t', low_memory=False)
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
    return None


def read_json(path):
    if path and os.path.exists(path):
        with open(path) as handle:
            return json.load(handle)
    return None


class Run:
    """The outputs of one egaDifferential run (and its egaExplain report), read lazily."""

    def __init__(self, prefix, explain_prefix=None, names=None, rna_threshold=None, chromatin_threshold=None,
                 logs=None, name=None):
        self.prefix = prefix
        self.name = name or os.path.basename(prefix)
        if explain_prefix is None:
            candidate = os.path.join(os.path.dirname(prefix), 'explain', os.path.basename(prefix))
            explain_prefix = candidate if os.path.exists(candidate + '_genes.tsv') else None
        self.explain_prefix = explain_prefix
        self.used = {}
        self.classes = self.table('_gene_classes.tsv', required=True)
        self.classes['gene'] = self.classes['genes'].map(first_gene)
        self.features = self.table('_gene_features.tsv')
        self.cascade = self.table('_explanation_cascade.tsv')
        self.cascade_learned = self.table('_explanation_cascade_learned.tsv')
        self.types = self.table('_explanation_types.tsv')
        self.roles = self.table('_mark_roles.tsv')
        self.state_transitions = self.table('_state_transitions.tsv')
        self.bivalent = self.table('_bivalent_transitions.tsv')
        self.state_marks = self.table('_state_marks.tsv')
        self.state_thresholds = self.table('_state_thresholds.tsv')
        self.community_expression = self.table('_community_expression.tsv')
        self.community_transitions = self.table('_community_transitions.tsv')
        self.gene_communities = self.table('_gene_community_changes.tsv')
        self.indirect_only = self.table('_indirect_only_genes.tsv')
        self.nested = self.table('_nested_models.tsv')
        self.all_levels_raw = self.table('_explanations_all_levels.tsv')
        self.combinations = self.table('_explanation_combinations.tsv')
        self.community_summary = self.json('_community_summary.json')
        self.permutation = self.json('_permutation.json')
        self.tests = {}
        for path in sorted(os.listdir(os.path.dirname(os.path.abspath(prefix)))):
            match = re.fullmatch(re.escape(os.path.basename(prefix)) + r'_(direct|indirect|hop\d+)_gene_test\.tsv',
                                 path)
            if match:
                self.tests[match.group(1)] = self.table('_{}_gene_test.tsv'.format(match.group(1)))
        self.explain_genes = self.explain_table('_genes.tsv')
        self.explain_communities = self.explain_table('_communities.tsv')
        self.edge_chromatin = self.explain_table('_edge_chromatin.tsv')
        self.log_paths = logs if logs is not None else [p for p in (prefix + '.log',) if os.path.exists(p)]
        self.log_lines = []
        for path in self.log_paths:
            with open(path, errors='replace') as handle:
                self.log_lines += [line.rstrip('\n') for line in handle]
            self.used['log: ' + os.path.basename(path)] = path
        self.command = self.find_command()
        self.rna_threshold = rna_threshold if rna_threshold is not None else self.command_value('--rnaThreshold', 1.0)
        self.chromatin_threshold = chromatin_threshold if chromatin_threshold is not None else \
            self.command_value('--chromatinThreshold', 1.0)
        self.names = list(names) if names else self.condition_names()
        self.marks = self.find_marks()
        self.max_hops = max([int(m.group(1)) for c in self.classes.columns
                             for m in [re.fullmatch(r'n_hop(\d+)', c)] if m] + [2])
        self.signs = self.declared_signs()
        self.support = self.step_support()
        self._graph = None
        self._edge_counts = None
        self.levels_variant = None
        self.level_counts = self.table('_explanation_levels.tsv')
        self.levels, self.levels_source = self.build_levels()

    # files ----------------------------------------------------------------------------------------------------------
    def table(self, suffix, required=False):
        path = self.prefix + suffix
        table = read_tsv(path)
        if table is None and required:
            raise FileNotFoundError('{} not found: is {} the prefix of an egaDifferential run?'.format(path,
                                                                                                     self.prefix))
        if table is not None:
            self.used[os.path.basename(path)] = path
        return table

    def json(self, suffix):
        path = self.prefix + suffix
        data = read_json(path)
        if data is not None:
            self.used[os.path.basename(path)] = path
        return data

    def explain_table(self, suffix):
        if not self.explain_prefix:
            return None
        path = self.explain_prefix + suffix
        table = read_tsv(path)
        if table is not None:
            self.used[os.path.basename(path)] = path
        return table

    @property
    def graph(self):
        if self._graph is None:
            path = self.prefix + '_union.graphml'
            log.info('reading %s', path)
            self._graph = read_graph(path)
            self.used[os.path.basename(path)] = path
        return self._graph

    def edge_counts(self):
        if self._edge_counts is None:
            path = self.prefix + '_edges.tsv.gz'
            if os.path.exists(path):
                status = pd.read_csv(path, sep='\t', usecols=['status'])['status']
                self._edge_counts = status.value_counts().to_dict()
                self.used[os.path.basename(path)] = path
            else:
                counts = {}
                for _, _, s in self.graph.edges(data='status'):
                    counts[s or 'shared'] = counts.get(s or 'shared', 0) + 1
                self._edge_counts = counts
        return self._edge_counts

    def find_command(self):
        for line in self.log_lines:
            if 'egaDifferential' in line and '--' in line:
                return line.strip()
        return ''

    def command_value(self, option, default):
        match = re.search(re.escape(option) + r'[ =]+([-0-9.eE]+)', self.command)
        return float(match.group(1)) if match else default

    def condition_names(self):
        if self.state_marks is not None and 'condition' in self.state_marks:
            names = list(dict.fromkeys(self.state_marks['condition'].astype(str)))
            if len(names) >= 2:
                return names[:2]
        path = self.prefix + '_edges.tsv.gz'
        if os.path.exists(path):
            with gzip.open(path, 'rt') as handle:
                header = handle.readline().rstrip('\n').split('\t')
            names = [c[len('score_'):] for c in header if c.startswith('score_')]
            if len(names) >= 2:
                return names[:2]
        from .egaPlotDifferential import condition_names
        return condition_names(self.graph)[:2]

    def find_marks(self):
        if self.roles is not None and 'mark' in self.roles:
            marks = [m for m in dict.fromkeys(self.roles['mark'].astype(str)) if m != 'contact']
            if marks:
                return marks
        return sorted({c[:-len('_promoter')] for c in self.classes.columns if c.endswith('_promoter')
                       and not c.startswith('promoter')})

    def declared_roles(self):
        roles = {}
        if self.roles is not None and 'declared_role' in self.roles:
            for row in self.roles.itertuples():
                roles.setdefault(row.mark, row.declared_role)
        for mark in self.marks:
            roles.setdefault(mark, 'repressive' if mark in ('H3K27me3', 'H3K9me3') else 'active')
        return roles

    def declared_signs(self):
        return {m: ROLE_SIGN.get(r, 0) for m, r in self.declared_roles().items() if m != 'contact'}

    def learned(self, mark, roles_level):
        if self.roles is None or 'learned_role' not in self.roles:
            return None
        rows = self.roles[(self.roles['mark'] == mark) & (self.roles['level'] == roles_level)]
        return rows.iloc[0] if len(rows) else None

    def step_support(self):
        """3D proximity support per step from the cascade table: {k: (support, interval)}."""
        support = {}
        if self.cascade is not None and 'proximity_support' in self.cascade:
            for row in self.cascade.itertuples():
                match = re.fullmatch(r'step(\d+)', str(row.level))
                if match and not isnan(row.proximity_support):
                    support[int(match.group(1))] = (float(row.proximity_support),
                                                    str(getattr(row, 'proximity_ci', '') or ''))
        return support

    # gene lookup ----------------------------------------------------------------------------------------------------
    def gene_nodes(self, gene):
        """Baited promoter fragments carrying `gene` (differential genes first, then any bait of the graph)."""
        hits = self.classes[self.classes['genes'].astype(str).str.split(',').map(lambda g: gene in g)]
        if len(hits):
            return hits['node'].tolist(), True
        if self.features is not None:
            hits = self.features[self.features['genes'].astype(str).str.split(',').map(lambda g: gene in g)]
            if len(hits):
                return hits['node'].tolist(), False
        nodes = [n for n, d in self.graph.nodes(data=True) if d.get('bait') and gene in str(d.get('genes', '')).split(',')]
        return nodes, False

    # levels ---------------------------------------------------------------------------------------------------------
    def build_levels(self):
        if self.all_levels_raw is not None and len(self.all_levels_raw):
            try:
                return self.levels_from_file(self.all_levels_raw), \
                    os.path.basename(self.prefix) + '_explanations_all_levels.tsv'
            except (KeyError, ValueError) as error:
                log.warning('cannot read the all-level table (%s); deriving the levels from the gene classes', error)
        return self.derive_levels(), os.path.basename(self.prefix) + '_gene_classes.tsv (derived)'

    def levels_from_file(self, raw):
        """Rows of <prefix>_explanations_all_levels.tsv (egaDifferential): one per gene and level, with the result of
        every level under a role set (declared, learned) and presence statistic (mean, max1kb) of the states. The
        report shows the declared roles with the mean statistic and adds the learned roles and max1kb results."""
        columns = list(raw.columns)
        c_node, c_gene, c_level = pick(columns, 'node'), pick(columns, 'gene', 'genes'), pick(columns, 'level')
        if c_level is None or (c_node is None and c_gene is None):
            raise KeyError('level and node / gene columns')
        variants = [(m.group(1), m.group(2)) for c in columns
                    for m in [re.fullmatch(r'explains_(\w+?)_(mean|max1kb)', c)] if m]
        roles = 'declared' if any(r == 'declared' for r, _ in variants) else (variants[0][0] if variants else None)
        presence = 'mean' if (roles, 'mean') in variants else \
            next((p for r, p in variants if r == roles), None)
        key = '{}_{}'.format(roles, presence) if roles else None
        self.levels_variant = key
        learned_key = 'learned_{}'.format(presence) if ('learned', presence) in variants and roles != 'learned' \
            else None
        max_key = '{}_max1kb'.format(roles) if (roles, 'max1kb') in variants and presence != 'max1kb' else None

        def col(generic, specific=None):
            return pick(columns, *(([specific] if specific else []) + list(generic)))

        c = dict(explains=col(('explains', 'explaining', 'explained'), 'explains_' + key if key else None),
                 detail=col(('detail', 'explanation', 'text'), 'detail_' + key if key else None),
                 result=col(('result',), 'result_' + key if key else None),
                 transition=col(('transition',), 'transition_' + key if key else None),
                 p=col(('p', 'pvalue', 'p_value'), 'p_' + roles if roles else None),
                 q=col(('q', 'qvalue', 'q_value', 'fdr'), 'q_' + roles if roles else None),
                 beyond=col(('beyond_mirror',), 'beyond_mirror_' + roles if roles else None),
                 concordant=col(('concordant_marks', 'marks'), 'concordant_marks_' + roles if roles else None),
                 strongest=col(('strongest',), 'strongest_' + roles if roles else None),
                 n_concordant=col(('n_concordant_nodes',), 'n_concordant_nodes_' + roles if roles else None),
                 changes=col(('mark_changes', 'changes')), mirror=col(('mirror_changes',)),
                 nodes=col(('nodes_tested', 'nodes', 'n_nodes')), support=col(('proximity_support', 'support_3d')),
                 support_ci=col(('proximity_support_ci', 'support_3d_ci', 'proximity_ci')),
                 rna=col(('rna_log2FC',)), gained=col(('contacts_gained',)), lost=col(('contacts_lost',)),
                 shared=col(('contacts_shared',)), contact=col(('contact_log2',)))

        def value(r, name):
            return r[c[name]] if c[name] else np.nan

        def parsed(text):
            out = {}
            for item in str(text).split(';'):
                if ':' in item:
                    mark, number = item.rsplit(':', 1)
                    try:
                        out[mark] = float(number)
                    except ValueError:
                        continue
            return out

        used = {v for v in c.values() if v} | {c_node, c_gene, c_level}
        classes = self.classes.set_index('node')
        rows = []
        for _, r in raw.iterrows():
            node = r[c_node] if c_node else None
            genes = r[c_gene] if c_gene else (classes.loc[node, 'genes'] if node in classes.index else '')
            if node is None:
                match = self.classes[self.classes['genes'].astype(str) == str(genes)]
                node = match['node'].iloc[0] if len(match) else str(genes)
            level = normalise_level(r[c_level])
            state = level in ('promoter_state', 'step1_state')
            row = dict(node=node, genes=genes, gene=first_gene(genes), level=level, level_label=level_label(level),
                       cascade_level=cascade_level(level), order=level_order(level), declared='', extra={})
            row['rna_log2FC'] = value(r, 'rna') if not isnan(value(r, 'rna')) else \
                (classes.loc[node, 'rna_log2FC'] if node in classes.index else np.nan)
            row['explains'] = truthy(value(r, 'explains'))
            row['nodes'] = value(r, 'nodes')
            row['p'], row['q'] = value(r, 'p'), value(r, 'q')
            row['beyond_mirror'] = truthy(value(r, 'beyond')) if not isnan(value(r, 'beyond')) else np.nan
            row['mirror'] = '' if isnan(row['beyond_mirror']) else \
                ('beyond the mirror' if row['beyond_mirror'] else 'not beyond the mirror (mirrored loci change the '
                                                                  'same way)')
            row['support'] = np.nan if state else value(r, 'support')
            row['support_ci'] = '' if state or isnan(value(r, 'support_ci')) else str(value(r, 'support_ci'))
            result = value(r, 'result')
            result = '' if isnan(result) else str(result)
            marks = parsed(value(r, 'changes')) if not isnan(value(r, 'changes')) else {}
            mirror = parsed(value(r, 'mirror')) if not isnan(value(r, 'mirror')) else {}
            if level == 'contacts':
                text = '{} gained, {} lost, {} shared; summed strength log2 {}'.format(
                    fmt_int(value(r, 'gained')), fmt_int(value(r, 'lost')), fmt_int(value(r, 'shared')),
                    fmt_num(value(r, 'contact'), signed=True)) if c['gained'] else str(value(r, 'detail'))
                row['changes'] = row['all_changes'] = text
                row['detail'] = '' if row['explains'] else result
            elif state:
                detail = value(r, 'detail')
                row['changes'] = row['all_changes'] = '' if isnan(detail) else str(detail)
                row['detail'] = ''
            else:
                concordant = [m for m in str(value(r, 'concordant') if not isnan(value(r, 'concordant')) else '')
                              .split(',') if m]
                row['changes'] = ', '.join('{} {}'.format(m, fmt_num(marks.get(m), signed=True)) for m in concordant)
                row['all_changes'] = ', '.join('{} {}{}'.format(m, fmt_num(v, signed=True),
                                                                 ' (mirror {})'.format(fmt_num(mirror[m], signed=True))
                                                                 if m in mirror else '') for m, v in marks.items())
                parts = []
                if not isnan(value(r, 'n_concordant')) and level != 'promoter_marks':
                    parts.append('{} concordant node{}'.format(fmt_int(value(r, 'n_concordant')),
                                                               '' if fmt_int(value(r, 'n_concordant')) == '1' else 's'))
                strongest = value(r, 'strongest')
                if row['explains'] and isinstance(strongest, str) and strongest:
                    parts.append('strongest {}'.format(strongest))
                if not row['explains'] and result:
                    parts.append(result)
                row['detail'] = '; '.join(parts)
            # the learned roles and the max1kb presence statistic
            notes = []
            if learned_key and 'explains_' + learned_key in raw:
                learned = truthy(r['explains_' + learned_key])
                learned_result = r.get('result_' + learned_key, '')
                notes.append('{}{}'.format('explains' if learned else 'does not explain',
                                           ' ({})'.format(learned_result) if isinstance(learned_result, str) and
                                           learned_result else ''))
            row['learned'] = '; '.join(notes)
            if max_key and state and 'explains_' + max_key in raw:
                other = r.get('result_' + max_key, '')
                row['extra']['max 1 kb statistic'] = '{}{}'.format(
                    'explains' if truthy(r['explains_' + max_key]) else 'does not explain',
                    ' ({})'.format(other) if isinstance(other, str) and other else '')
            for column in columns:
                if column not in used and not column.endswith(('_learned_mean', '_learned_max1kb', '_learned')) \
                        and not isnan(r[column]) and column not in ('rank', 'hops'):
                    if roles and (column.endswith('_' + roles) or column.endswith('_' + key) or
                                  (max_key and column.endswith('_' + max_key))):
                        continue
                    row['extra'][column] = r[column]
            rows.append(row)
        table = pd.DataFrame(rows)
        return self.assign_ranks(table).reset_index(drop=True)

    @staticmethod
    def assign_ranks(table):
        table = table.sort_values(['node', 'order']).copy()
        ranks = []
        for _, part in table.groupby('node', sort=False):
            k = 0
            for explains in part['explains']:
                if explains:
                    k += 1
                    ranks.append(k)
                else:
                    ranks.append(np.nan)
        table['rank'] = ranks
        return table

    def roles_text(self, roles_level, contact=False):
        declared = self.declared_roles()
        marks = ['contact'] if contact else self.marks
        parts_declared, parts_learned = [], []
        for mark in marks:
            parts_declared.append('{} {}'.format(mark, declared.get(mark, 'active')))
            row = self.learned(mark, roles_level)
            if row is not None and (truthy(row.get('disagrees', False)) or contact):
                parts_learned.append('{} {}'.format(mark, row['learned_role']))
        learned = '; '.join(parts_learned) if parts_learned else ('as declared' if self.roles is not None else '')
        return '; '.join(parts_declared), learned

    def concordant(self, value, sign, direction):
        return not isnan(value) and np.sign(float(value)) * sign == direction and \
            abs(float(value)) >= self.chromatin_threshold

    def derive_levels(self):
        """One row per differential gene and level from the gene classes (no all-level table)."""
        rows = []
        support = self.support
        tests = {}
        for name, test in self.tests.items():
            if test is not None and len(test):
                tests[name] = test.set_index('node')
        hop_names = {1: 'direct', 2: 'indirect'}
        for _, r in self.classes.iterrows():
            node, direction = r['node'], np.sign(float(r['rna_log2FC']))
            base = dict(node=node, genes=r['genes'], gene=r['gene'], rna_log2FC=r['rna_log2FC'], extra={})

            def add(level, explains, changes, all_changes='', nodes=np.nan, p=np.nan, q=np.nan, beyond=np.nan,
                    mirror='', sup=np.nan, sup_ci='', detail='', roles_level=None, contact=False):
                declared, learned = self.roles_text(roles_level, contact) if roles_level else ('', '')
                rows.append(dict(base, level=level, level_label=level_label(level), cascade_level=cascade_level(level),
                                 order=level_order(level), explains=bool(explains), changes=changes,
                                 all_changes=all_changes, nodes=nodes, p=p, q=q, beyond_mirror=beyond, mirror=mirror,
                                 support=sup, support_ci=sup_ci, declared=declared, learned=learned, detail=detail))

            # promoter marks
            values = {m: r.get(m + '_promoter', np.nan) for m in self.marks}
            hit = [m for m in self.marks if self.concordant(values[m], self.signs.get(m, 1), direction)]
            add('promoter_marks', bool(hit), ', '.join('{} {}'.format(m, fmt_num(values[m], signed=True)) for m in hit),
                ', '.join('{} {}'.format(m, fmt_num(v, signed=True)) for m, v in values.items() if not isnan(v)),
                nodes=1, sup=1.0 if support else np.nan, sup_ci='NA' if support else '', roles_level='promoter')
            # own interactions of the promoter
            contact = r.get('contact_direct', np.nan)
            counts = [r.get('interactions_' + s, np.nan) for s in ('gained', 'lost', 'shared')]
            text = '{} gained, {} lost, {} shared'.format(*[fmt_int(c) for c in counts]) \
                if not any(isnan(c) for c in counts) else ''
            text += '{}summed strength log2 {}'.format('; ' if text else '', fmt_num(contact, signed=True)) \
                if not isnan(contact) else ''
            add('contacts', self.concordant(contact, 1, direction), text, text,
                nodes=sum(c for c in counts if not isnan(c)) if not all(isnan(c) for c in counts) else np.nan,
                roles_level='promoter', contact=True)
            # promoter state
            transition = r.get('promoter_transition', np.nan)
            if not isnan(transition) and transition != 'states not callable':
                explains = truthy(r.get('promoter_explains')) or truthy(r.get('promoter_unresolved'))
                kind = 'poised, trigger not measured' if truthy(r.get('promoter_unresolved')) else \
                    ('state transition' if truthy(r.get('promoter_explains')) else '')
                add('promoter_state', explains, '{}{}'.format(transition, ' ({})'.format(kind) if kind else ''),
                    str(r.get('promoter_state_summary', '') if not isnan(r.get('promoter_state_summary')) else ''),
                    nodes=1, roles_level='promoter')
            # steps
            for k in range(1, self.max_hops + 1):
                name = hop_names.get(k, 'hop{}'.format(k))
                values = {m: r.get('{}_{}'.format(m, name), np.nan) for m in self.marks}
                mirrors = {m: r.get('{}_{}_mirror'.format(m, name), np.nan) for m in self.marks}
                level_signs = {m: self.level_sign(m, k) for m in self.marks}
                hit = [m for m in self.marks if self.concordant(values[m], level_signs[m], direction)]
                mirror_hit = [m for m in hit if self.concordant(mirrors[m], level_signs[m], direction)]
                have_mirror = any(not isnan(v) for v in mirrors.values())
                beyond = (len(mirror_hit) < len(hit)) if hit and have_mirror else np.nan
                p = q = np.nan
                if name in tests and node in tests[name].index:
                    p = tests[name].loc[node, name + '_p']
                    q = tests[name].loc[node, name + '_q'] if name + '_q' in tests[name] else np.nan
                    p = p.iloc[0] if isinstance(p, pd.Series) else p
                    q = q.iloc[0] if isinstance(q, pd.Series) else q
                parts = ['{} {}{}'.format(m, fmt_num(values[m], signed=True),
                                          ' (mirror {})'.format(fmt_num(mirrors[m], signed=True))
                                          if not isnan(mirrors[m]) else '') for m in self.marks if not isnan(values[m])]
                detail = ''
                if k == 1 and not isnan(r.get('step1_fragments', np.nan)):
                    detail = 'concordant at {} of {} contacting fragments'.format(
                        fmt_int(r.get('step1_concordant_fragments')), fmt_int(r.get('step1_fragments')))
                    if isinstance(r.get('step1_strongest'), str) and r.get('step1_strongest'):
                        detail += '; strongest {}'.format(r['step1_strongest'])
                sup, sup_ci = support.get(k, (np.nan, ''))
                mirror = '' if isnan(beyond) else ('beyond the mirror' if beyond else
                                                   'not beyond the mirror (mirrored loci change the same way)')
                add('step{}'.format(k), bool(hit),
                    ', '.join('{} {}'.format(m, fmt_num(values[m], signed=True)) for m in hit), '; '.join(parts),
                    nodes=r.get('n_' + name, np.nan), p=p, q=q, beyond=beyond, mirror=mirror, sup=sup, sup_ci=sup_ci,
                    detail=detail, roles_level='step{}'.format(k))
                if k == 1:
                    transition = r.get('step1_transition', np.nan)
                    if not isnan(transition) and transition not in ('states not callable', 'no poised fragment'):
                        explains = truthy(r.get('step1_explains')) or truthy(r.get('step1_unresolved'))
                        kind = 'poised, trigger not measured' if truthy(r.get('step1_unresolved')) else \
                            ('state transition' if truthy(r.get('step1_explains')) else '')
                        fragment = r.get('step1_state_fragment', '')
                        add('step1_state', explains,
                            '{}{}'.format(transition, ' ({})'.format(kind) if kind else ''),
                            str(r.get('step1_state_summary', '') if not isnan(r.get('step1_state_summary')) else ''),
                            nodes=r.get('step1_transition_fragments', np.nan), sup=sup, sup_ci=sup_ci,
                            detail='fragment {}'.format(fragment) if isinstance(fragment, str) and fragment else '',
                            roles_level='step1')
        table = pd.DataFrame(rows)
        if not len(table):
            return pd.DataFrame(columns=['node', 'gene', 'level', 'rank', 'explains', 'order', 'cascade_level'])
        return self.assign_ranks(table).reset_index(drop=True)

    def level_sign(self, mark, k):
        """Declared sign of a mark at step k (enhancer and active +1, repressive -1)."""
        return self.signs.get(mark, 1)

    # derived views --------------------------------------------------------------------------------------------------
    def first_level(self):
        """Node -> first cascade level (explained_at of egaDifferential)."""
        if 'explained_at' in self.classes:
            return dict(zip(self.classes['node'], self.classes['explained_at'].astype(str)))
        first = {}
        for node, part in self.levels[self.levels['explains']].groupby('node'):
            first[node] = cascade_level(part.sort_values('order')['level'].iloc[0])
        return first

    def any_level_counts(self):
        """Genes for which a level explains, at any rank (cascade levels gene, step1, ...)."""
        if not len(self.levels):
            return pd.Series(dtype=int)
        explains = self.levels[self.levels['explains']]
        counts = explains.groupby('cascade_level')['node'].nunique()
        order = sorted(counts.index, key=lambda c: (c != 'gene', level_hops(c)))
        return counts.reindex(order)

    def variant_rows(self, table):
        """Rows of the role set and presence statistic shown in the report (declared, mean), without those columns."""
        if {'roles', 'presence'} <= set(table.columns):
            roles, presence = (self.levels_variant or 'declared_mean').split('_', 1)
            part = table[(table['roles'] == roles) & (table['presence'] == presence)]
            table = (part if len(part) else table).drop(columns=['roles', 'presence'])
        return table.reset_index(drop=True)

    def level_combinations(self):
        if self.combinations is not None and len(self.combinations):
            table = self.variant_rows(self.combinations)
            if 'combination' in table:
                table['combination'] = table['combination'].map(
                    lambda c: ' + '.join(level_short(normalise_level(x)) for x in str(c).split(';')))
            return table, True
        explains = self.levels[self.levels['explains']]
        combos = explains.groupby('node')['cascade_level'].agg(
            lambda c: ' + '.join(sorted(set(c), key=lambda x: (x != 'gene', level_hops(x)))))
        all_nodes = set(self.classes['node'])
        none = len(all_nodes - set(combos.index))
        rna = self.classes.set_index('node')['rna_log2FC']
        table = combos.to_frame('levels').assign(up=lambda d: (rna.reindex(d.index) > 0).to_numpy())
        result = table.groupby('levels').agg(genes=('up', 'size'), genes_up=('up', 'sum')).reset_index()
        result['genes_down'] = result['genes'] - result['genes_up']
        if none:
            result = pd.concat([result, pd.DataFrame([{'levels': 'none', 'genes': none, 'genes_up': np.nan,
                                                        'genes_down': np.nan}])], ignore_index=True)
        result['fraction_of_differential'] = result['genes'] / len(all_nodes)
        return result.sort_values('genes', ascending=False).reset_index(drop=True), False


# --------------------------------------------------------------------------------------------------------------------
# HTML document

class Document:
    """Sections, paragraphs, numbered figures (captions below) and numbered tables (captions below)."""

    def __init__(self, title, keep_figures=False):
        self.title = title
        self.parts = []
        self.toc = []
        self.figure_count = 0
        self.table_count = 0
        self.keep_figures = keep_figures
        self.figures = []      # (Figure objects, caption) for the matplotlib PDF
        self.text_blocks = []  # plain text in order, for the matplotlib PDF

    def section(self, anchor, title, level=2):
        self.toc.append((level, anchor, title))
        self.parts.append('<h{l} id="{a}">{t}</h{l}>'.format(l=level, a=esc(anchor), t=esc(title)))
        self.text_blocks.append(('heading', title, level))

    def heading(self, title, level=4):
        self.parts.append('<h{l}>{t}</h{l}>'.format(l=level, t=title))
        self.text_blocks.append(('heading', re.sub('<[^>]+>', '', title), level))

    def html(self, markup, text=None):
        self.parts.append(markup)
        plain = text if text is not None else re.sub(r'<[^>]+>', ' ', markup)
        plain = html.unescape(re.sub(r'\s+', ' ', plain)).strip()
        if plain:
            self.text_blocks.append(('text', plain, 0))

    def p(self, markup, cls=None):
        self.html('<p{}>{}</p>'.format(' class="{}"'.format(cls) if cls else '', markup))

    def figure(self, figures, caption, cols=1):
        """Panels (matplotlib Figures) as inline SVG in a grid, caption below."""
        self.figure_count += 1
        number = self.figure_count
        panels = []
        for k, fig in enumerate(figures):
            panels.append('<div class="panel">{}</div>'.format(svg_markup(fig, 'f{}p{}'.format(number, k))))
        self.parts.append('<figure id="figure-{n}"><div class="panels cols-{c}">{p}</div><figcaption><span '
                          'class="label">Figure {n}.</span> {cap}</figcaption></figure>'.format(
                              n=number, c=cols, p=''.join(panels), cap=caption))
        plain = 'Figure {}. {}'.format(number, html.unescape(re.sub(r'<[^>]+>', '', caption)))
        if self.keep_figures:
            self.figures.append((figures, plain))
            self.text_blocks.append(('figure', len(self.figures) - 1, 0))
        else:
            for fig in figures:
                fig.clear()
        return number

    def table(self, frame, caption, formats=None, numeric=(), long=(), escape=True, max_rows=None, nowrap=()):
        """A table with its caption below; `formats` maps columns to formatting functions."""
        self.table_count += 1
        number = self.table_count
        formats = formats or {}
        frame = frame if max_rows is None else frame.head(max_rows)
        nowrap = set(nowrap) | {c for c in ('level', 'mark', 'category', 'gene', 'kind', 'side') if c in frame.columns}
        head = ''.join('<th class="{}">{}</th>'.format('num' if c in numeric else '',
                                                      esc(str(c).replace('_lt_', ' < ').replace('_', ' ')))
                       for c in frame.columns)
        body = []
        for _, row in frame.iterrows():
            cells = []
            for c in frame.columns:
                value = row[c]
                text = formats[c](value) if c in formats else ('' if isnan(value) else str(value))
                text = esc(text) if escape else text
                cls = 'num' if c in numeric else ('long' if c in long else ('key' if c in nowrap else
                                                                          ('frag' if c in FRAGMENT_COLUMNS else '')))
                cells.append('<td class="{}">{}</td>'.format(cls, text))
            body.append('<tr>{}</tr>'.format(''.join(cells)))
        self.parts.append('<div class="table-wrap" id="table-{n}"><div class="table-scroll"><table><thead><tr>{h}</tr>'
                          '</thead><tbody>{b}</tbody></table></div><p class="caption"><span class="label">Table {n}.'
                          '</span> {cap}</p></div>'.format(n=number, h=head, b=''.join(body), cap=caption))
        plain_rows = [' | '.join(str(c) for c in frame.columns)]
        for _, row in frame.iterrows():
            plain_rows.append(' | '.join(html.unescape(re.sub('<[^>]+>', '', formats[c](row[c]) if c in formats
                                                              else ('' if isnan(row[c]) else str(row[c]))))
                                         for c in frame.columns))
        self.text_blocks.append(('table', '\n'.join(plain_rows), 0))
        self.text_blocks.append(('text', 'Table {}. {}'.format(number, html.unescape(re.sub('<[^>]+>', '', caption))), 0))
        return number

    def render(self, heading, meta):
        with open(os.path.join(TEMPLATE_DIR, 'report.css')) as handle:
            css = handle.read()
        with open(os.path.join(TEMPLATE_DIR, 'report.html')) as handle:
            template = string.Template(handle.read())
        items = []
        for level, anchor, title in self.toc:
            items.append('<li class="{}"><a href="#{}">{}</a></li>'.format('sub' if level > 2 else '', esc(anchor),
                                                                          esc(title)))
        return template.substitute(version=esc(__version__), title=esc(self.title), css=css,
                                   toc='<ol>{}</ol>'.format(''.join(items)), heading=heading, meta=meta,
                                   body='\n'.join(self.parts))


def new_figure(width, height):
    fig = Figure(figsize=(width, height))
    FigureCanvasAgg(fig)
    return fig


def svg_markup(fig, uid):
    """Inline SVG of a figure: no XML prolog, no metadata, no namespace URLs, ids unique per figure."""
    buffer = io.StringIO()
    fig.savefig(buffer, format='svg', bbox_inches='tight', pad_inches=0.04, metadata={'Date': None})
    text = buffer.getvalue()
    text = re.sub(r'<\?xml[^>]*\?>\s*', '', text)
    text = re.sub(r'<!DOCTYPE[^>]*>\s*', '', text)
    text = re.sub(r'<metadata>.*?</metadata>\s*', '', text, flags=re.S)
    text = re.sub(r'\s+xmlns(:\w+)?="[^"]*"', '', text)
    text = re.sub(r'id="([^"]+)"', r'id="{}-\1"'.format(uid), text)
    text = re.sub(r'url\(#([^)]+)\)', r'url(#{}-\1)'.format(uid), text)
    text = re.sub(r'href="#([^"]+)"', r'href="#{}-\1"'.format(uid), text)
    match = re.search(r'<svg[^>]*\swidth="([0-9.]+)pt"', text)
    width = float(match.group(1)) if match else None
    text = re.sub(r'(<svg[^>]*?)\swidth="[^"]*"', r'\1', text, count=1)
    text = re.sub(r'(<svg[^>]*?)\sheight="[^"]*"', r'\1', text, count=1)
    style = ' style="max-width:{:.0f}pt"'.format(width) if width else ''
    text = text.replace('<svg', '<svg role="img"{}'.format(style), 1)
    return text


# --------------------------------------------------------------------------------------------------------------------
# figures

def node_middle(graph, node):
    data = graph.nodes[node]
    return (int(data['start']) + int(data['end'])) / 2


def distances_from(graph, seeds, cutoff):
    """Shortest-path steps from the nearest seed, up to `cutoff`."""
    distance = {s: 0 for s in seeds if s in graph}
    frontier = list(distance)
    for step in range(1, cutoff + 1):
        nxt = []
        for n in frontier:
            for w in graph.neighbors(n):
                if w not in distance:
                    distance[w] = step
                    nxt.append(w)
        frontier = nxt
    return distance


def explaining_nodes(run, seeds, levels, direction, limit=12):
    """Nodes that carry the explanation at each explaining level: the promoter for promoter levels, fragments at
    exactly k steps with a concordant change of at least the chromatin threshold for step k (mark changes from the
    graph attributes, log2 (B + 1) / (A + 1))."""
    graph = run.graph
    steps = sorted({level_hops(lv) for lv in levels if level_hops(lv) > 0})
    if not steps:
        return {}, {}
    distance = distances_from(graph, seeds, max(steps))
    found, which = {}, {}
    for k in steps:
        candidates = []
        for n, d in distance.items():
            if d != k:
                continue
            for mark in run.marks:
                change = mark_change(graph.nodes[n], mark, run.names)
                if run.concordant(change, run.level_sign(mark, k), direction):
                    candidates.append((abs(change), n, mark, change))
        candidates.sort(key=lambda x: (-x[0], x[1]))
        chosen = []
        for _, n, mark, change in candidates:
            if n not in chosen:
                chosen.append(n)
                which[n] = (k, mark, change)
            if len(chosen) >= limit:
                break
        found[k] = chosen
    return found, which


def neighbourhood_subgraph(run, seeds, hops, explaining, max_nodes):
    """Up to `hops` steps around the seeds; larger neighbourhoods are reduced to the seeds, their strongest direct
    partners and one shortest path to each explaining node."""
    graph = run.graph
    distance = distances_from(graph, seeds, max(hops, 1))
    if len(distance) <= max_nodes:
        nodes = set(distance)
    else:
        nodes = set(seeds)
        partners = [n for n, d in distance.items() if d == 1]
        strength = {n: max(edge_strength(graph.edges[s, n]) for s in seeds if graph.has_edge(s, n)) for n in partners}
        partners.sort(key=lambda n: (-strength[n], n))
        nodes |= set(partners[:max(10, max_nodes // 3)])
        for target in explaining:
            source = min((s for s in seeds if s in graph), key=lambda s: distance.get(target, 99))
            try:
                nodes |= set(nx.shortest_path(graph, source, target))
            except nx.NetworkXNoPath:
                continue
    chrom = {graph.nodes[s]['chrom'] for s in seeds}
    nodes = {n for n in nodes if graph.nodes[n].get('chrom') in chrom}
    sub = ordered_subgraph(graph, nodes)
    return ordered_subgraph(sub, [n for n in sub.nodes if sub.degree[n] > 0 or n in seeds])


def highlight_handles(gene_label='gene promoter', explaining_label='explaining node'):
    return [plt.Line2D([0], [0], marker='s', color='white', markeredgecolor='#000000', markeredgewidth=1.8, ms=7, lw=0,
                       label=gene_label),
            plt.Line2D([0], [0], marker='o', color='white', markeredgecolor=EXPLAINING, markeredgewidth=1.6, ms=6,
                       lw=0, label=explaining_label)]


def overlay_highlights(ax, sub, positions, seeds, explaining):
    for n in explaining:
        if n in positions:
            bait = sub.nodes[n].get('bait')
            ax.scatter(*positions[n], marker='s' if bait else 'o', s=95 if bait else 40, facecolor='none',
                       edgecolor=EXPLAINING, lw=1.5, zorder=4)
    for n in seeds:
        if n in positions:
            ax.scatter(*positions[n], marker='s', s=100, facecolor='none', edgecolor='#000000', lw=1.8, zorder=4)


def arcs_figure(run, sub, seeds, explaining):
    """Interaction arcs of condition A, condition B and their change (gained above, lost below)."""
    names = run.names
    fig = new_figure(7.0, 4.3)
    grid = fig.add_gridspec(3, 1, height_ratios=[1, 1, 2.3], hspace=0.12)
    axes = [fig.add_subplot(grid[k, 0]) for k in range(3)]
    middle = {n: node_middle(sub, n) for n in sub.nodes}
    span = (min(middle.values()), max(middle.values()))
    pad = 0.03 * (span[1] - span[0] or 1.0)
    longest = max((abs(middle[u] - middle[v]) for u, v in sub.edges), default=1.0) or 1.0
    for ax, which in zip(axes[:2], ('A', 'B')):
        count = 0
        for u, v, data in sub.edges(data=True):
            status = data.get('status', 'shared')
            if which == 'A' and status == 'gained' or which == 'B' and status == 'lost':
                continue
            count += 1
            x0, x1 = sorted((middle[u], middle[v]))
            ax.add_patch(Arc(((x0 + x1) / 2, 0), x1 - x0, 2 * np.sqrt((x1 - x0) / longest), theta1=0, theta2=180,
                             color='#7f7f7a', lw=edge_width(edge_strength(data, which)), alpha=0.9))
        ax.axhline(0, color=INK, lw=0.6)
        for n in sub.nodes:
            if sub.nodes[n].get('bait'):
                ax.plot([middle[n]], [0], marker='s', ms=3, color=INK, zorder=3)
        ax.set_xlim(span[0] - pad, span[1] + pad)
        ax.set_ylim(-0.08, 1.15)
        ax.set_yticks([])
        ax.set_xticks([])
        for side in ('left', 'right', 'top', 'bottom'):
            ax.spines[side].set_visible(False)
        ax.text(0.0, 0.98, '{}: {} interactions'.format(names[0 if which == 'A' else 1], count),
                transform=ax.transAxes, fontsize=7, va='top', color=INK)
    draw_arcs(axes[2], sub, names, priority=set(seeds))
    for ax in axes:
        for n in explaining:
            if n in middle:
                ax.plot([middle[n]], [0], marker='o', ms=4.5, mfc='white', mec=EXPLAINING, mew=1.2, zorder=4)
        for n in seeds:
            if n in middle:
                ax.plot([middle[n]], [0], marker='s', ms=5, mfc=INK, mec='#000000', mew=1.0, zorder=5)
    handles = [plt.Line2D([0], [0], color='#7f7f7a', lw=2, label='interaction of the condition'),
               plt.Line2D([0], [0], color=STATUS_COLORS['gained'], lw=2, label='gained (only in {})'.format(names[1])),
               plt.Line2D([0], [0], color=STATUS_COLORS['lost'], lw=2, label='lost (only in {})'.format(names[0])),
               plt.Line2D([0], [0], color=STATUS_COLORS['shared'], lw=2, label='shared')] + width_handles() + \
        [plt.Line2D([0], [0], marker='s', ls='', color=INK, ms=4, label='baited promoter'),
         plt.Line2D([0], [0], marker='o', ls='', mfc='white', mec=EXPLAINING, mew=1.2, ms=5, label='explaining node')]
    axes[2].legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.42), fontsize=6, frameon=False,
                   ncol=3, handlelength=2.2)
    return fig


def graph_panel_legend(fig, ax, mode, names, mark, max_log2, highlight=True):
    """Legend and colour bar below one graph panel, with the highlight outlines."""
    extra = highlight_handles() if highlight else []
    if mode in ('A', 'B'):
        handles = width_handles() + shape_handles() + extra
        offset = -0.04
    elif mode == 'edges':
        handles = status_handles(names) + width_handles()[::2] + shape_handles() + extra
        cax = ax.inset_axes([0.05, -0.1, 0.9, 0.03])
        bar = fig.colorbar(plt.cm.ScalarMappable(norm=edge_norm, cmap='RdBu_r'), cax=cax, orientation='horizontal')
        bar.set_label('log2 change {} / {} of shared interactions'.format(names[1], names[0]), fontsize=6.5)
        bar.ax.tick_params(labelsize=6)
        offset = -0.32
    else:
        handles = shape_handles(neutral=False) + extra
        norm = mcolors.TwoSlopeNorm(vmin=-max_log2, vcenter=0, vmax=max_log2)
        cax = ax.inset_axes([0.05, -0.1, 0.9, 0.03])
        bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='RdBu_r'), cax=cax, orientation='horizontal')
        bar.set_label('log2 change {} / {}: expression (squares){}'.format(
            names[1], names[0], ', {} (circles)'.format(mark) if mark else ''), fontsize=6.5)
        bar.ax.tick_params(labelsize=6)
        offset = -0.32
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, offset), fontsize=6, frameon=False, ncol=2,
              handlelength=2.2)


def neighbourhood_figures(run, sub, seeds, explaining, mark, max_log2):
    """Neighbourhood graph in condition A, condition B, the change of the interactions and the change of expression
    and of the explaining mark; one panel each, on one layout."""
    names = run.names
    positions = spread_layout(sub, aspect=1.2)
    status = {}
    for _, _, d in sub.edges(data=True):
        status[d.get('status', 'shared')] = status.get(d.get('status', 'shared'), 0) + 1
    titles = {'A': '{}: {} interactions'.format(names[0], status.get('shared', 0) + status.get('lost', 0)),
              'B': '{}: {} interactions'.format(names[1], status.get('shared', 0) + status.get('gained', 0)),
              'edges': 'change: {} gained, {} lost, {} shared'.format(status.get('gained', 0), status.get('lost', 0),
                                                                      status.get('shared', 0)),
              'nodes': 'expression{} change'.format(' and {}'.format(mark) if mark else '')}
    figures = []
    for mode in ('A', 'B', 'edges', 'nodes'):
        fig = new_figure(3.5, 3.3)
        ax = fig.add_subplot(1, 1, 1)
        draw_graph(ax, sub, positions, names, mark, max_log2, mode=mode,
                   candidates=set(seeds) if mode == 'nodes' else (),
                   drivers=set(explaining) if mode == 'nodes' else (), title=titles[mode])
        if mode != 'nodes':
            overlay_highlights(ax, sub, positions, seeds, explaining)
        graph_panel_legend(fig, ax, mode, names, mark, max_log2)
        figures.append(fig)
    return figures


def hub_figures(run, sub, positions, members_a, members_b, focus, title_a, title_b, max_log2, members_label='hub',
                max_labels=18):
    """A hub or community in condition A and condition B side by side on one layout: interactions of each condition
    (widths by strength), promoters coloured by the expression change, members that left (orange dotted outline) or
    joined (green outline) between the conditions, members of the other condition only faded."""
    names = run.names
    left = members_a - members_b
    joined = members_b - members_a
    cmap = plt.get_cmap('RdBu_r')
    norm = mcolors.TwoSlopeNorm(vmin=-max_log2, vcenter=0, vmax=max_log2)
    figures = []
    for which, members, title in (('A', members_a, title_a), ('B', members_b, title_b)):
        fig = new_figure(3.5, 3.4)
        ax = fig.add_subplot(1, 1, 1)
        for u, v, data in sub.edges(data=True):
            status = data.get('status', 'shared')
            if which == 'A' and status == 'gained' or which == 'B' and status == 'lost':
                continue
            ax.plot([positions[u][0], positions[v][0]], [positions[u][1], positions[v][1]], color='#8c8b85',
                    lw=edge_width(edge_strength(data, which)), alpha=0.85, zorder=1)
        labelled, seen = [], set()
        for n, data in sub.nodes(data=True):
            bait = bool(data.get('bait'))
            inside = n in members or not bait
            value = data.get('rna_log2FC')
            if bait and value is not None and not isnan(value):
                face = cmap(norm(np.clip(float(value), -max_log2, max_log2)))
            else:
                face = '#ffffff'
            edge, width, style = INK, 0.5, '-'
            if n in left:
                edge, width, style = LEFT_COLOR, 1.9, (0, (1.6, 1.0))
            elif n in joined:
                edge, width = JOINED_COLOR, 1.9
            if n in focus:
                edge, width, style = '#000000', 2.0, '-'
            ax.scatter(*positions[n], marker='s' if bait else 'o', s=48 if bait else 12, color=face, edgecolor=edge,
                       lw=width, linestyle=style, alpha=1.0 if inside else 0.35, zorder=3 if bait else 2)
            if bait:
                label = first_gene(data.get('genes', ''))
                if label and label not in seen and (n in members_a | members_b):
                    seen.add(label)
                    labelled.append((n, label))
        ax.set_axis_off()
        xy = np.array(list(positions.values()))
        ax.set_xlim(xy[:, 0].min() - 0.12, xy[:, 0].max() + 0.12)
        ax.set_ylim(xy[:, 1].min() - 0.08, xy[:, 1].max() + 0.15)
        if len(labelled) > max_labels:
            # the gene of the report, then the promoters with the largest expression change, then left / joined
            def priority(item):
                n = item[0]
                value = sub.nodes[n].get('rna_log2FC')
                change = abs(float(value)) if value is not None and not isnan(value) else -1.0
                return (n not in focus, -change, n not in left and n not in joined, n)
            labelled = sorted(labelled, key=priority)[:max_labels]
        place_labels(ax, labelled, positions)
        ax.set_title(title, loc='left', fontsize=7.5)
        cax = ax.inset_axes([0.05, -0.1, 0.9, 0.03])
        bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='RdBu_r'), cax=cax, orientation='horizontal')
        bar.set_label('expression log2 change {} / {} (promoters)'.format(names[1], names[0]), fontsize=6.5)
        bar.ax.tick_params(labelsize=6)
        handles = [plt.Line2D([0], [0], marker='s', ls='', color='#ffffff', markeredgecolor=INK, ms=6,
                              label='baited promoter'),
                   plt.Line2D([0], [0], marker='o', ls='', color='#ffffff', markeredgecolor=INK, ms=4,
                              label='other end (OE)'),
                   ax.scatter([], [], marker='s', s=40, facecolor='#ffffff', edgecolor=LEFT_COLOR, lw=1.9,
                              linestyle=(0, (1.6, 1.0)), label='left the {} in {}'.format(members_label, names[1])),
                   ax.scatter([], [], marker='s', s=40, facecolor='#ffffff', edgecolor=JOINED_COLOR, lw=1.9,
                              label='joined the {} in {}'.format(members_label, names[1])),
                   plt.Line2D([0], [0], marker='s', ls='', color='#ffffff', markeredgecolor='#000000',
                              markeredgewidth=2.0, ms=6, label='gene of the report'),
                   plt.Line2D([0], [0], marker='s', ls='', color='#ffffff', markeredgecolor=INK, alpha=0.35, ms=6,
                              label='member of the other condition only')] + width_handles()[::2]
        ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.32), fontsize=6, frameon=False, ncol=2,
                  handlelength=2.0)
        figures.append(fig)
    return figures


def cascade_figure(run):
    """Genes per cascade level: first explaining level (egaDifferential cascade) and any explaining level."""
    cascade = run.cascade
    levels = ['gene'] + ['step{}'.format(k) for k in range(1, run.max_hops + 1)]
    labels = {'gene': 'gene (promoter, contacts, state)', 'unexplained': 'unexplained'}
    first = cascade.set_index('level')['genes'] if cascade is not None else pd.Series(dtype=float)
    total = int(first.reindex(levels + ['unexplained']).fillna(0).sum()) if len(first) else len(run.classes)
    any_level = run.any_level_counts()
    rows = levels + ['unexplained']
    names = [labels.get(r, 'step {}'.format(level_hops(r))) for r in rows]
    fig = new_figure(7.0, 0.42 * len(rows) + 1.0)
    axes = [fig.add_subplot(1, 2, k + 1) for k in range(2)]
    for ax, values, title in ((axes[0], first.reindex(rows).fillna(0), 'first explaining level'),
                              (axes[1], any_level.reindex(rows).fillna(0), 'any explaining level')):
        if title.startswith('any'):
            values['unexplained'] = first.get('unexplained', np.nan) if len(first) else np.nan
        y = np.arange(len(rows))[::-1]
        ax.barh(y, values.to_numpy(dtype=float), height=0.62, color=BAR_COLOR, lw=0)
        top = float(np.nanmax(values.to_numpy(dtype=float))) if len(values) else 1.0
        for yy, v in zip(y, values.to_numpy(dtype=float)):
            if not np.isnan(v):
                ax.text(v + top * 0.02, yy, '{:,} ({:.0f}%)'.format(int(v), 100 * v / total if total else 0),
                        va='center', fontsize=6.5, color=INK)
        ax.set_yticks(y)
        ax.set_yticklabels(names if ax is axes[0] else [''] * len(names), fontsize=7)
        ax.set_xlim(0, top * 1.45 if top > 0 else 1)
        ax.tick_params(axis='x', labelsize=6.5)
        ax.set_xlabel('differential genes', fontsize=7)
        ax.set_title(title, loc='left', fontsize=7.5)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        ax.grid(axis='x', color='#e6e5e0', lw=0.6)
        ax.set_axisbelow(True)
    fig.subplots_adjust(wspace=0.08)
    return fig, total


def roles_figure(run):
    """Spearman correlation of every mark change with the expression change per level, with the bootstrap interval;
    filled markers: learned role equals the declared role, open: it differs."""
    roles = run.roles[run.roles['versus'] == 'expression'] if 'versus' in run.roles else run.roles
    roles = roles[roles['mark'] != 'contact']
    levels = list(dict.fromkeys(roles['level']))
    marks = list(dict.fromkeys(roles['mark']))
    fig = new_figure(6.6, 2.7)
    ax = fig.add_subplot(1, 1, 1)
    width = 0.7 / max(len(marks), 1)
    for k, mark in enumerate(marks):
        color = MARK_COLORS[k % len(MARK_COLORS)]
        part = roles[roles['mark'] == mark].set_index('level').reindex(levels)
        x = np.arange(len(levels)) - 0.35 + width * (k + 0.5)
        rho = part['rho'].to_numpy(dtype=float)
        low = part['rho_ci_low'].to_numpy(dtype=float)
        high = part['rho_ci_high'].to_numpy(dtype=float)
        ax.vlines(x, low, high, color=color, lw=1.4)
        differs = part['disagrees'].map(truthy).to_numpy() if 'disagrees' in part else np.zeros(len(x), bool)
        ax.scatter(x[~differs], rho[~differs], color=color, s=26, zorder=3, label=mark)
        ax.scatter(x[differs], rho[differs], facecolor='white', edgecolor=color, lw=1.4, s=26, zorder=3)
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_xticks(np.arange(len(levels)))
    ax.set_xticklabels([str(lv).replace('step', 'step ') for lv in levels], fontsize=7)
    ax.set_ylabel('Spearman rho with expression change', fontsize=7)
    ax.tick_params(axis='y', labelsize=6.5)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.grid(axis='y', color='#e6e5e0', lw=0.6)
    ax.set_axisbelow(True)
    handles, _ = ax.get_legend_handles_labels()
    handles += [plt.Line2D([0], [0], marker='o', ls='', color=INK, ms=5, label='learned role = declared role'),
                plt.Line2D([0], [0], marker='o', ls='', mfc='white', mec=INK, ms=5, label='learned role differs'),
                plt.Line2D([0], [0], color=INK, lw=1.4, label='bootstrap interval')]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=4, fontsize=6.5, frameon=False)
    return fig


# --------------------------------------------------------------------------------------------------------------------
# text of the explanation of one gene

def level_sentence(row, names):
    """One item of the ranked explanation: level, what changed, nodes, P / q, mirror, 3D support."""
    parts = []
    changes = row.get('changes') if row.get('explains') else row.get('all_changes')
    if isinstance(changes, str) and changes:
        parts.append(esc(changes))
    if isinstance(row.get('detail'), str) and row.get('detail'):
        parts.append(esc(row['detail']))
    if not isnan(row.get('nodes')) and row.get('level') not in ('promoter_marks', 'promoter_state'):
        parts.append('{} node{}'.format(fmt_int(row['nodes']), '' if fmt_int(row['nodes']) == '1' else 's'))
    if not isnan(row.get('p')):
        parts.append('P = {}'.format(fmt_p(row['p'])) + ('' if isnan(row.get('q')) else
                                                         ', q = {}'.format(fmt_p(row['q']))))
    if isinstance(row.get('mirror'), str) and row.get('mirror'):
        parts.append(esc(row['mirror']))
    if not isnan(row.get('support')) and level_hops(row.get('level')) > 0:
        parts.append('3D support {}{}'.format(fmt_num(row['support']), ' ({})'.format(esc(row['support_ci']))
                                              if isinstance(row.get('support_ci'), str) and row['support_ci']
                                              not in ('', 'NA') else ''))
    return '; '.join(parts)


def ranked_list(rows, names):
    """Ordered list of all levels of a gene: explaining levels in rank order, then the others."""
    rows = rows.sort_values('order')
    items = []
    explaining = rows[rows['explains']].sort_values(['rank', 'order'])
    for _, row in explaining.iterrows():
        items.append('<li><strong>{}</strong> <span class="tag yes">explains, rank {}</span>: {}</li>'.format(
            esc(row['level_label']), fmt_int(row['rank']), level_sentence(row, names) or 'change recorded'))
    for _, row in rows[~rows['explains']].iterrows():
        sentence = level_sentence(row, names)
        items.append('<li class="miss"><strong>{}</strong> <span class="tag">does not explain</span>{}</li>'.format(
            esc(row['level_label']), ': ' + sentence if sentence else ''))
    return '<ol class="levels">{}</ol>'.format(''.join(items))


def level_short(level):
    short = {'promoter_marks': 'promoter marks', 'contacts': 'contacts', 'promoter_state': 'promoter state',
             'step1': 'step 1', 'step1_state': 'step 1 state'}
    return short.get(level, 'step {}'.format(level_hops(level)) if re.fullmatch(r'step\d+', str(level)) else level)


def evidence_table(doc, rows, gene, names, declared=''):
    rows = rows.sort_values(['explains', 'rank', 'order'], ascending=[False, True, True])

    def mirror(row):
        if not isnan(row['beyond_mirror']):
            return 'beyond' if truthy(row['beyond_mirror']) else 'not beyond'
        return row['mirror'] if isinstance(row['mirror'], str) else ''

    def support(row):
        if isnan(row['support']) or level_hops(row['level']) == 0:
            return ''
        ci = row['support_ci'] if isinstance(row['support_ci'], str) and row['support_ci'] not in ('', 'NA') else ''
        return '{}{}'.format(fmt_num(row['support']), ' ({})'.format(ci) if ci else '')

    changes = 'changes (log2 {} / {})'.format(names[1], names[0])
    frame = pd.DataFrame({
        'rank': rows['rank'], 'level': rows['level'].map(level_short),
        'explains': rows['explains'].map({True: 'yes', False: 'no'}), changes: rows['all_changes'],
        'nodes': rows['nodes'], 'P': rows['p'], 'q': rows['q'], 'mirror': rows.apply(mirror, axis=1),
        '3D support': rows.apply(support, axis=1), 'learned roles': rows['learned']})
    frame = frame[[c for c in frame.columns if not frame[c].map(lambda v: isnan(v) or v == '').all()]]
    doc.table(frame, 'Evidence of every level for {}, explaining levels in rank order. Changes are log2 fold changes '
                     '{} / {} (mark signal at the promoter or contact-weighted over the nodes of a step, values at the '
                     'mirrored loci in brackets); P and q of the per-gene permutation test of the step; mirror: '
                     'whether the explaining marks change differently at the mirrored loci; 3D support: probability '
                     'that graph neighbours at this step are closer in the nucleus than distance-matched pairs (0.5 '
                     'without relation), with its interval; learned roles: the result of the level under the roles '
                     'learned from the data (all-level table), or the marks whose learned role differs from the '
                     'declared one (derived levels){}.'.format(
                         gene_span(gene), names[1], names[0],
                         ' (declared roles: {})'.format(esc(declared)) if declared else ''),
              formats={'rank': fmt_int, 'nodes': fmt_int, 'P': fmt_p, 'q': fmt_p},
              numeric=('rank', 'nodes', 'P', 'q'), long=(changes,), nowrap=('level', '3D support', 'mirror'))


# --------------------------------------------------------------------------------------------------------------------
# report builder

class Report:
    def __init__(self, run, args, keep_figures):
        self.run = run
        self.args = args
        self.keep_figures = keep_figures

    # shared blocks ------------------------------------------------------------------------------------------------
    def gene_figures(self, doc, gene, nodes, rows):
        """Arcs and neighbourhood graph of a gene, the explaining nodes highlighted."""
        run = self.run
        seeds = [n for n in nodes if n in run.graph]
        if not seeds:
            doc.p('The promoter of {} is not in the union graph.'.format(gene_span(gene)), cls='note')
            return [], {}
        direction = np.sign(float(rows['rna_log2FC'].iloc[0])) if len(rows) and not isnan(rows['rna_log2FC'].iloc[0]) \
            else 1.0
        levels = rows.loc[rows['explains'], 'level'].tolist() if len(rows) else []
        found, which = explaining_nodes(run, seeds, levels, direction)
        explaining = [n for k in sorted(found) for n in found[k]]
        hops = max([level_hops(lv) for lv in levels] + [1])
        sub = neighbourhood_subgraph(run, seeds, hops, explaining, self.args.maxNodes)
        explaining = [n for n in explaining if n in sub]
        if sub.number_of_edges() == 0:
            doc.p('No interactions of {} in either condition.'.format(gene_span(gene)), cls='note')
            return explaining, which
        steps = sorted(found)
        what = '; '.join('step {}: {}'.format(k, ', '.join('{} ({} {})'.format(
            first_gene(run.graph.nodes[n].get('genes', '')) or n, which[n][1], fmt_num(which[n][2], signed=True))
            for n in found[k] if n in sub)) for k in steps if any(n in sub for n in found[k]))
        doc.figure([arcs_figure(run, sub, seeds, explaining)],
                   'Interactions around {g} in {a} (top), {b} (middle) and their change (bottom; gained in {b} above '
                   'the axis in red, lost below in blue, shared in grey; line width by interaction strength relative '
                   'to the calling threshold). Squares on the axis: baited promoters, the promoter of {g} in black; '
                   'open purple circles: explaining nodes. {n} nodes, {e} interactions{r}.'.format(
                       g=gene_span(gene), a=esc(run.names[0]), b=esc(run.names[1]), n=sub.number_of_nodes(),
                       e=sub.number_of_edges(), r=' (reduced neighbourhood)' if sub.number_of_nodes() >=
                       self.args.maxNodes else ''))
        mark = None
        for lv in sorted(levels, key=level_order, reverse=True):
            row = rows[rows['level'] == lv].iloc[0]
            names = [m for m in run.marks if isinstance(row['changes'], str) and m in row['changes']]
            if names:
                mark = names[0]
                break
        mark = mark or (run.marks[0] if run.marks else None)
        doc.figure(neighbourhood_figures(run, sub, seeds, explaining, mark, self.args.maxLog2),
                   'Neighbourhood graph of {g} up to {h} step{s} on one layout: interactions in {a} and in {b} (grey, '
                   'width by strength), the change of the interactions (shared: log2 change {b} / {a}, gained and lost '
                   'dashed) and the expression change of the promoters (squares) with the change of {m} at the other '
                   'ends (circles). Black outline: promoter of {g}; purple outline: explaining nodes ({w}).'.format(
                       g=gene_span(gene), h=hops, s='' if hops == 1 else 's', a=esc(run.names[0]),
                       b=esc(run.names[1]), m=esc(mark or 'the mark'),
                       w=esc(what) if what else 'none at a step; the gene is explained at its promoter'), cols=2)
        return explaining, which

    def hub_block(self, doc, gene, nodes, draw=True):
        """Hub of a gene: its promoter, the partners and the other promoters contacting these partners, in A and B."""
        run = self.run
        graph = run.graph
        seeds = [n for n in nodes if n in graph]
        if not seeds:
            return
        members = {}
        for which, excluded in (('A', 'gained'), ('B', 'lost')):
            partners = {w for s in seeds for w in graph.neighbors(s)
                        if graph.edges[s, w].get('status', 'shared') != excluded}
            promoters = {p for w in partners for p in graph.neighbors(w)
                         if graph.nodes[p].get('bait') and graph.edges[w, p].get('status', 'shared') != excluded}
            members[which] = (partners, (promoters | {w for w in partners if graph.nodes[w].get('bait')}) - set(seeds))
        promoters_a, promoters_b = members['A'][1], members['B'][1]
        partners_a, partners_b = members['A'][0], members['B'][0]
        all_promoters = promoters_a | promoters_b
        rna = {n: graph.nodes[n].get('rna_log2FC') for n in all_promoters}
        differential = set(run.classes['node'])
        direction = np.sign(float(graph.nodes[seeds[0]].get('rna_log2FC', 0) or 0))
        same = sorted((n for n in all_promoters if n in differential and rna[n] is not None
                       and np.sign(float(rna[n])) == direction), key=lambda n: -abs(float(rna[n])))
        opposite = sorted((n for n in all_promoters if n in differential and rna[n] is not None
                           and np.sign(float(rna[n])) == -direction), key=lambda n: -abs(float(rna[n])))
        left = sorted(promoters_a - promoters_b, key=lambda n: node_middle(graph, n))
        joined = sorted(promoters_b - promoters_a, key=lambda n: node_middle(graph, n))

        def names_of(items, with_value=True):
            if not items:
                return 'none'
            out = []
            for n in items:
                label = first_gene(graph.nodes[n].get('genes', '')) or n
                value = rna.get(n, graph.nodes[n].get('rna_log2FC'))
                if with_value and value is not None and not isnan(value):
                    cls = 'up' if float(value) > 0 else 'down'
                    out.append('{} <span class="{}">{}</span>'.format(gene_span(label), cls,
                                                                       fmt_num(value, signed=True)))
                else:
                    out.append(gene_span(label))
            return ', '.join(out)

        verb = 'up' if direction > 0 else 'down'
        doc.html('<dl class="lists">'
                 '<dt>Hub</dt><dd>{a} partners and {pa} other promoters in {na}; {b} partners and {pb} other promoters '
                 'in {nb} ({g} {v}, log2 {fc}).</dd>'
                 '<dt>Same direction as {g} (differential promoters of the hub)</dt><dd>{same}</dd>'
                 '<dt>Opposite direction</dt><dd>{opp}</dd>'
                 '<dt>Left the hub in {nb}</dt><dd>{left}</dd>'
                 '<dt>Joined the hub in {nb}</dt><dd>{joined}</dd></dl>'.format(
                     a=len(partners_a), pa=len(promoters_a), b=len(partners_b), pb=len(promoters_b),
                     na=esc(run.names[0]), nb=esc(run.names[1]), g=gene_span(gene), v=verb,
                     fc=fmt_num(graph.nodes[seeds[0]].get('rna_log2FC'), signed=True),
                     same=names_of(same[:25]) + (' and {} more'.format(len(same) - 25) if len(same) > 25 else ''),
                     opp=names_of(opposite[:25]) + (' and {} more'.format(len(opposite) - 25)
                                                    if len(opposite) > 25 else ''),
                     left=names_of(left[:30]) + (' and {} more'.format(len(left) - 30) if len(left) > 30 else ''),
                     joined=names_of(joined[:30]) + (' and {} more'.format(len(joined) - 30)
                                                     if len(joined) > 30 else '')))
        if not draw:
            return
        nodes_drawn = set(seeds) | partners_a | partners_b | all_promoters
        if len(nodes_drawn) > self.args.maxNodes:
            keep = set(seeds) | set(same[:10]) | set(opposite[:10]) | set(left[:10]) | set(joined[:10])
            keep |= {w for w in partners_a | partners_b if any(graph.has_edge(w, p) for p in keep)}
            nodes_drawn = set(sorted(keep, key=lambda n: (n not in seeds, n))[:self.args.maxNodes + 20])
        sub = ordered_subgraph(graph, nodes_drawn)
        sub = ordered_subgraph(sub, [n for n in sub.nodes if sub.degree[n] > 0])
        if sub.number_of_edges() == 0:
            return
        positions = separate(spread_layout(sub, aspect=1.2), sub, min_distance=0.06)
        figures = hub_figures(run, sub, positions, (promoters_a | partners_a | set(seeds)) & set(sub.nodes),
                              (promoters_b | partners_b | set(seeds)) & set(sub.nodes), set(seeds),
                              '{}: hub of {}'.format(run.names[0], gene), '{}: hub of {}'.format(run.names[1], gene),
                              self.args.maxLog2)
        doc.figure(figures, 'Hub of {g} in {a} (left) and {b} (right) on one layout: the promoter of {g} (black '
                            'outline), its partners and the other promoters contacting them. Promoters coloured by the '
                            'expression change {b} / {a}; orange dotted outline: left the hub in {b}; green outline: '
                            'joined the hub in {b}; faded: present only in the other condition.{r}'.format(
                                g=gene_span(gene), a=esc(run.names[0]), b=esc(run.names[1]),
                                r=' Reduced to the members named above and their partners.' if
                                sub.number_of_nodes() >= self.args.maxNodes else ''), cols=2)

    def community_rows(self, side, community):
        table = self.run.gene_communities
        if table is None:
            return None
        column = 'community_a' if side == 'A' else 'community_b'
        return table[table[column] == community]

    def community_block(self, doc, community_a, community_b, focus=(), info=None, draw=True):
        """Members of a community in A and its counterpart in B: expression, marks, left and joined, figure."""
        run = self.run
        table = run.gene_communities
        if table is None:
            doc.p('No community table (<code>{}_gene_community_changes.tsv</code>).'.format(esc(run.name)), cls='note')
            return
        members_a = set(table.loc[table['community_a'] == community_a, 'node']) if community_a >= 0 else set()
        members_b = set(table.loc[table['community_b'] == community_b, 'node']) if community_b >= 0 else set()
        rows = table[table['node'].isin(members_a | members_b)].copy()
        rows['in'] = np.where(rows['node'].isin(members_a) & rows['node'].isin(members_b), 'both',
                              np.where(rows['node'].isin(members_a), 'left in {}'.format(run.names[1]),
                                       'joined in {}'.format(run.names[1])))
        rows['gene'] = rows['genes'].map(first_gene)
        mean = rows['rna_log2FC'].mean()
        direction = np.sign(mean) if not isnan(mean) else 1.0
        diff = rows['differential'].map(truthy) if 'differential' in rows else rows['rna_log2FC'].abs() >= \
            run.rna_threshold
        same = rows[diff & (np.sign(rows['rna_log2FC']) == direction)].sort_values('rna_log2FC', key=np.abs,
                                                                                   ascending=False)
        opposite = rows[diff & (np.sign(rows['rna_log2FC']) == -direction)].sort_values('rna_log2FC', key=np.abs,
                                                                                        ascending=False)

        def listing(frame, limit=25):
            if not len(frame):
                return 'none'
            items = ['{} <span class="{}">{}</span>'.format(gene_span(r.gene or r.node), 'up' if r.rna_log2FC > 0
                                                             else 'down', fmt_num(r.rna_log2FC, signed=True))
                     for r in frame.head(limit).itertuples()]
            return ', '.join(items) + (' and {} more'.format(len(frame) - limit) if len(frame) > limit else '')

        left = rows[rows['in'].str.startswith('left')]
        joined = rows[rows['in'].str.startswith('joined')]

        def plain(frame, limit=30):
            names = [g for g in frame['gene'] if g]
            return ', '.join(gene_span(g) for g in names[:limit]) + (' and {} more'.format(len(names) - limit)
                                                                    if len(names) > limit else '') if names else 'none'

        transition = info.get('transition', '') if info else ''
        doc.html('<dl class="lists"><dt>Members</dt><dd>{na} community {ca}: {a} promoters; {nb} community {cb}: '
                 '{b} promoters; {both} in both; transition: {t}; mean expression change {m}.</dd>'
                 '<dt>Same direction as the community mean (differential)</dt><dd>{same}</dd>'
                 '<dt>Opposite direction</dt><dd>{opp}</dd><dt>Left in {nb}</dt><dd>{left}</dd>'
                 '<dt>Joined in {nb}</dt><dd>{joined}</dd></dl>'.format(
                     na=esc(run.names[0]), nb=esc(run.names[1]), ca=community_a if community_a >= 0 else 'none',
                     cb=community_b if community_b >= 0 else 'none', a=len(members_a), b=len(members_b),
                     both=len(members_a & members_b), t=esc(transition or 'not recorded'),
                     m=fmt_num(mean, signed=True), same=listing(same), opp=listing(opposite), left=plain(left),
                     joined=plain(joined)))
        # member table with expression and promoter mark changes
        features = run.features.set_index('node') if run.features is not None else None
        show = rows.sort_values('rna_log2FC', key=np.abs, ascending=False)
        frame = pd.DataFrame({'gene': show['gene'], 'promoter fragment': show['node'], 'membership': show['in'],
                              'expression log2FC': show['rna_log2FC'],
                              'differential': show['differential'].map(lambda v: 'yes' if truthy(v) else 'no')
                              if 'differential' in show else ''})
        if features is not None:
            for mark in run.marks:
                column = mark + '_promoter'
                if column in features:
                    frame['{} promoter'.format(mark)] = features[column].reindex(show['node']).to_numpy()
        numeric = ['expression log2FC'] + ['{} promoter'.format(m) for m in run.marks]
        doc.table(frame, 'Promoters of the community, by absolute expression change (up to 30 shown): membership in {} '
                         'and {}, expression log2 fold change {} / {}, and the log2 change of the marks at the '
                         'promoter.'.format(esc(run.names[0]), esc(run.names[1]), esc(run.names[1]),
                                            esc(run.names[0])),
                  formats={c: (lambda v: fmt_num(v, signed=True)) for c in numeric}, numeric=numeric, max_rows=30)
        if not draw:
            return
        graph = run.graph
        baits = [n for n in members_a | members_b if n in graph]
        if not baits:
            return
        # promoters of the community and the other ends that connect at least two of them (or the promoter's
        # strongest partner if it has no shared one)
        counts = {}
        for n in baits:
            for w in graph.neighbors(n):
                if not graph.nodes[w].get('bait'):
                    counts[w] = counts.get(w, 0) + 1
        connectors = sorted((w for w, c in counts.items() if c >= 2), key=lambda w: (-counts[w], w))
        nodes = set(baits) | set(connectors[:max(0, self.args.maxNodes - len(baits))])
        if len(baits) > self.args.maxNodes:
            keep = set(focus) | set(same['node'].head(12)) | set(opposite['node'].head(12)) | \
                set(left['node'].head(12)) | set(joined['node'].head(12))
            nodes = keep | {w for w in connectors if sum(graph.has_edge(w, p) for p in keep) >= 2}
        sub = ordered_subgraph(graph, nodes)
        sub = ordered_subgraph(sub, [n for n in sub.nodes if sub.degree[n] > 0])
        if sub.number_of_edges() == 0:
            doc.p('The promoters of the community share no interacting fragment in the union graph.', cls='note')
            return
        positions = separate(spread_layout(sub, aspect=1.2), sub, min_distance=0.06)
        figures = hub_figures(run, sub, positions, (members_a | set(n for n in sub if not sub.nodes[n].get('bait')))
                              & set(sub.nodes), (members_b | set(n for n in sub if not sub.nodes[n].get('bait')))
                              & set(sub.nodes), set(focus),
                              '{}: community {}'.format(run.names[0], community_a if community_a >= 0 else 'none'),
                              '{}: community {}'.format(run.names[1], community_b if community_b >= 0 else 'none'),
                              self.args.maxLog2, members_label='community')
        doc.figure(figures, 'Community in {a} (left) and its counterpart in {b} (right) on one layout: member promoters '
                            '(squares, coloured by the expression change {b} / {a}) and the other ends shared by at '
                            'least two of them, with the interactions of each condition. Orange dotted outline: '
                            'promoters that left the community in {b}; green outline: promoters that joined; faded: '
                            'member in the other condition only.{r}'.format(
                                a=esc(run.names[0]), b=esc(run.names[1]),
                                r=' Reduced to the named members and their shared other ends.'
                                if len(baits) > self.args.maxNodes else ''), cols=2)

    def gene_community_text(self, doc, gene, nodes):
        run = self.run
        if run.gene_communities is None:
            return None
        rows = run.gene_communities[run.gene_communities['node'].isin(nodes)]
        if not len(rows):
            doc.p('{} is not assigned to a community.'.format(gene_span(gene)), cls='note')
            return None
        r = rows.iloc[0]

        def genes(value, limit=30):
            names = [g for g in str(value).split(',') if g and g != 'nan']
            if not names:
                return 'none'
            text = ', '.join(gene_span(g) for g in names[:limit])
            return text + (' and {} more'.format(len(names) - limit) if len(names) > limit else '')

        doc.p('Community of {g}: {s} (community {ca} of {na}, {sa} promoters; community {cb} of {nb}, {sb} promoters; '
              'member Jaccard {j}; transition of the {na} community: {t}). Genes that joined its community in {nb}: '
              '{jo}. Genes that left: {le}.'.format(
                  g=gene_span(gene), s=esc(r.get('status', '')), ca=fmt_int(r.get('community_a')),
                  cb=fmt_int(r.get('community_b')), na=esc(run.names[0]), nb=esc(run.names[1]),
                  sa=fmt_int(r.get('size_a')), sb=fmt_int(r.get('size_b')), j=fmt_num(r.get('member_jaccard')),
                  t=esc(r.get('community_a_transition', '') if not isnan(r.get('community_a_transition')) else
                        'none'), jo=genes(r.get('joined_genes')), le=genes(r.get('left_genes'))))
        return r

    def gene_finding(self, doc, gene, nodes, anchor=None, figures=True, hub=True, level_heading=4):
        """Explanation of a gene: ranked levels, evidence, arcs, neighbourhood, hub."""
        run = self.run
        rows = run.levels[run.levels['node'].isin(nodes)]
        for k, node in enumerate(nodes):
            part = rows[rows['node'] == node]
            if not len(part):
                continue
            genes = part['genes'].iloc[0]
            first = run.first_level().get(node, '')
            label = 'promoter fragment <span class="mono">{}</span> ({})'.format(esc(node), esc(genes))
            doc.html('<p>{}: expression log2 fold change {} / {} {}; first level of the cascade: {}. {} level{} '
                     'explain{}.'
                     '</p>'.format(label, esc(run.names[1]), esc(run.names[0]),
                                   fmt_num(part['rna_log2FC'].iloc[0], signed=True), esc(first or 'none'),
                                   int(part['explains'].sum()), '' if int(part['explains'].sum()) == 1 else 's',
                                   's' if int(part['explains'].sum()) == 1 else ''))
            doc.html(ranked_list(part, run.names))
            classes = run.classes.set_index('node')
            if node in classes.index:
                for column, title in (('explanation_detail', 'egaDifferential (declared roles)'),
                                      ('explanation_detail_learned', 'egaDifferential (learned roles)')):
                    value = classes.loc[node, column] if column in classes else np.nan
                    if isinstance(value, str) and value:
                        doc.p('{}: {}'.format(esc(title), esc(value)), cls='note')
        if figures:
            self.gene_figures(doc, gene, nodes, rows)
        if hub:
            self.hub_block(doc, gene, nodes, draw=figures)

    # main report ----------------------------------------------------------------------------------------------------
    def main_report(self):
        run = self.run
        a, b = run.names
        doc = Document('{} vs {}: differential graph report'.format(a, b), self.keep_figures)
        self.summary(doc)
        self.cascade_section(doc)
        findings = self.findings_section(doc)
        self.hubs_section(doc, findings)
        self.methods(doc)
        return doc

    def summary(self, doc):
        run = self.run
        a, b = run.names
        doc.section('summary', 'Summary of the comparison')
        counts = run.edge_counts()
        total_edges = sum(counts.values())
        n_diff = len(run.classes)
        up = int((run.classes['rna_log2FC'] > 0).sum())
        tested = len(run.features) if run.features is not None else np.nan
        transitions = (run.community_summary or {}).get('transitions', {})
        stats = [('{:,}'.format(n_diff), 'differential genes ({:,} up, {:,} down in {})'.format(up, n_diff - up, b)),
                 (fmt_int(tested), 'promoters with an expression change'),
                 ('{:,}'.format(counts.get('gained', 0)), 'interactions gained in {}'.format(b)),
                 ('{:,}'.format(counts.get('lost', 0)), 'interactions lost in {}'.format(b)),
                 ('{:,}'.format(counts.get('shared', 0)), 'interactions shared'),
                 ('{:,}'.format(sum(transitions.values())) if transitions else '', 'community transitions')]
        doc.html('<div class="stats">{}</div>'.format(''.join(
            '<div class="stat"><div class="value">{}</div><div class="label">{}</div></div>'.format(esc(v), esc(t))
            for v, t in stats if v)))
        marks = run.declared_roles()
        doc.p('Comparison of {a} (condition A) with {b} (condition B), output prefix <code>{p}</code>. The union graph '
              'holds {e:,} interactions: {g:,} gained in {b}, {l:,} lost in {b} and {s:,} shared. Chromatin marks with '
              'their declared roles: {m}. Differential genes: |log2 fold change| of at least {r:g}{pa}; mark changes '
              'explain a gene from |log2 fold change| {c:g} on.'.format(
                  a=esc(a), b=esc(b), p=esc(run.prefix), e=total_edges, g=counts.get('gained', 0),
                  l=counts.get('lost', 0), s=counts.get('shared', 0),
                  m=esc(', '.join('{} {}'.format(k, v) for k, v in marks.items() if k != 'contact')),
                  r=run.rna_threshold, pa=' and the adjusted P threshold of the run' if run.command and
                  '--padjColumn' in run.command else '', c=run.chromatin_threshold))
        rows = []
        if run.permutation:
            perm = run.permutation
            rows.append(('genes in the models (promoter, direct and indirect partners)', fmt_int(perm.get('n_genes'))))
            rows.append(('R2 gain of the indirect level (observed)', fmt_num(perm.get('observed_indirect_delta_r2'), 5)))
            rows.append(('R2 gain of the indirect level (permutations, mean)',
                         fmt_num(perm.get('null_indirect_delta_r2_mean'), 5)))
            rows.append(('empirical P of the indirect R2 gain', fmt_p(perm.get('empirical_p_delta_r2'))))
            rows.append(('permutations', fmt_int(perm.get('n_permutations'))))
            mirror = perm.get('mirror_comparison', {})
            for level in ('direct', 'indirect', 'hop3'):
                key = '{}_graph_minus_mirror'.format(level)
                if key in mirror:
                    rows.append(('R2 gain graph minus mirror, {} level (95% interval)'.format(level),
                                 '{} ({} to {})'.format(fmt_num(mirror[key], 5), fmt_num(mirror.get(level + '_ci_low'), 5),
                                                         fmt_num(mirror.get(level + '_ci_high'), 5))))
        summary = run.community_summary or {}
        if transitions:
            rows.append(('community transitions', ', '.join('{} {:,}'.format(k, v) for k, v in transitions.items())))
        if summary.get('gene_status'):
            rows.append(('promoters by community fate', ', '.join('{} {:,}'.format(k, v)
                                                                 for k, v in summary['gene_status'].items())))
        if 'communities_tested' in summary:
            rows.append(('community pairs tested for coordinated expression', fmt_int(summary['communities_tested'])))
            rows.append(('coordinated community pairs (q < 0.1)',
                         fmt_int(summary.get('coherent_community_pairs_q_lt_0.1'))))
        switching = summary.get('switching_vs_differential', {})
        if switching:
            rows.append(('community switching vs differential expression: odds ratio (Fisher P)',
                         '{} ({})'.format(fmt_num(switching.get('odds_ratio')), fmt_p(switching.get('fisher_p')))))
        if rows:
            doc.table(pd.DataFrame(rows, columns=['quantity', 'value']),
                      'Graph-level tests and community changes of the run (<code>{n}_permutation.json</code>, '
                      '<code>{n}_community_summary.json</code>).'.format(n=esc(run.name)), long=('value',))
        if run.nested is not None and len(run.nested):
            doc.table(run.nested, 'Nested regression models of the expression change on the mark changes at the '
                                  'promoter, the direct partners and the indirect partners (<code>{}_nested_models.tsv'
                                  '</code>).'.format(esc(run.name)),
                      formats={'r2': lambda v: fmt_num(v, 4), 'delta_r2': lambda v: fmt_num(v, 5),
                               'f_statistic': fmt_num, 'f_pvalue': fmt_p, 'n_genes': fmt_int},
                      numeric=('n_genes', 'r2', 'n_parameters', 'delta_r2', 'f_statistic', 'f_pvalue'))

    def cascade_section(self, doc):
        run = self.run
        a, b = run.names
        doc.section('cascade', 'Explanation cascade')
        doc.p('Every differential gene is explained at the first level with a concordant change: the promoter (mark '
              'changes, the own interactions, the promoter chromatin state), step 1 (direct partners and their '
              'chromatin states), then the nodes at exactly 2, 3, ... steps. The any-level view counts every level '
              'with a concordant change, whatever its rank. Per-gene levels are read from <code>{}</code>.'.format(
                  esc(run.levels_source)))
        if run.cascade is not None:
            fig, total = cascade_figure(run)
            doc.figure([fig], 'Differential genes ({:,}) by the first level that explains their expression change '
                              '(left, <code>{n}_explanation_cascade.tsv</code>) and by every level that explains it '
                              '(right; a gene counts at each level with a concordant change). Percentages of all '
                              'differential genes.'.format(total, n=esc(run.name)))
            cascade = run.cascade.copy()
            columns = [c for c in ('level', 'genes', 'genes_reaching', 'fraction_of_differential',
                                   'fraction_of_remaining', 'beyond_mirror', 'p_lt_0.05', 'q_lt_0.1', 'q_lt_0.05',
                                   'proximity_support', 'proximity_ci', 'median_nodes_tested') if c in cascade]
            doc.table(cascade[columns], 'Explanation cascade: genes first explained at each level, genes reaching it '
                                        '(not explained before), the fractions, genes beyond the mirror (the explaining '
                                        'marks do not all change the same way at the mirrored loci), genes with a '
                                        'per-gene permutation P below 0.05 and q below 0.1 and 0.05, the 3D proximity '
                                        'support of the step with its interval, and the median number of nodes '
                                        'tested.',
                      formats={'fraction_of_differential': fmt_pct, 'fraction_of_remaining': fmt_pct,
                               'beyond_mirror': fmt_int, 'p_lt_0.05': fmt_int, 'q_lt_0.1': fmt_int,
                               'q_lt_0.05': fmt_int, 'proximity_support': fmt_num, 'median_nodes_tested': fmt_num,
                               'genes': fmt_int, 'genes_reaching': fmt_int},
                      numeric=('genes', 'genes_reaching', 'fraction_of_differential', 'fraction_of_remaining',
                               'beyond_mirror', 'p_lt_0.05', 'q_lt_0.1', 'q_lt_0.05', 'proximity_support',
                               'median_nodes_tested'))
            if run.cascade_learned is not None:
                learned = run.cascade_learned[[c for c in ('level', 'genes', 'fraction_of_differential',
                                                           'beyond_mirror') if c in run.cascade_learned]]
                merged = cascade[['level', 'genes']].merge(learned, on='level', how='outer',
                                                           suffixes=(' (declared roles)', ' (learned roles)'))
                doc.table(merged, 'The cascade under the declared mark roles and under the roles learned at every '
                                  'level (<code>{}_explanation_cascade_learned.tsv</code>).'.format(esc(run.name)),
                          formats={'genes (declared roles)': fmt_int, 'genes (learned roles)': fmt_int,
                                   'fraction_of_differential': fmt_pct, 'beyond_mirror': fmt_int},
                          numeric=('genes (declared roles)', 'genes (learned roles)', 'fraction_of_differential',
                                   'beyond_mirror'))
        if run.level_counts is not None and len(run.level_counts):
            counts = run.variant_rows(run.level_counts)
            columns = [c for c in ('rank', 'level', 'genes_any', 'genes_first', 'fraction_any', 'fraction_first',
                                   'any_transition', 'any_poised_in_both', 'any_beyond_mirror', 'any_p_lt_0.05',
                                   'any_q_lt_0.05', 'first_bh_q_lt_0.05', 'proximity_support', 'median_nodes_tested')
                       if c in counts]
            ints = [c for c in columns if c.startswith(('genes_', 'any_', 'first_'))]
            doc.table(counts[columns], 'Genes explained at every level whatever the earlier levels (any) and genes '
                                       'first explained there (first), with the state transitions and poised states, '
                                       'the genes beyond the mirror and with P &lt; 0.05 and q &lt; 0.05 at any rank, '
                                       'the 3D support and the median nodes tested (<code>{}_explanation_levels.tsv'
                                       '</code>, declared roles, mean presence statistic).'.format(esc(run.name)),
                      formats=dict({c: fmt_int for c in ints}, fraction_any=fmt_pct, fraction_first=fmt_pct,
                                   proximity_support=fmt_num, median_nodes_tested=fmt_num),
                      numeric=tuple(c for c in columns if c not in ('level',)))
        doc.section('combinations', 'Level combinations', level=3)
        combos, from_file = run.level_combinations()
        if from_file:
            doc.table(combos, 'The most frequent combinations of explaining levels per gene '
                              '(<code>{}_explanation_combinations.tsv</code>, declared roles, mean presence statistic), '
                              'first 15 rows.'.format(esc(run.name)), max_rows=15,
                      formats={'genes': fmt_int, 'fraction_of_differential': fmt_pct},
                      numeric=('rank', 'n_levels', 'genes', 'fraction_of_differential'))
        else:
            doc.table(combos, 'Combinations of explaining levels per differential gene (any rank), counted from the '
                              'per-gene levels; genes up and down in {}. First 15 rows.'.format(esc(b)),
                      formats={'genes': fmt_int, 'genes_up': fmt_int, 'genes_down': fmt_int,
                               'fraction_of_differential': fmt_pct},
                      numeric=('genes', 'genes_up', 'genes_down', 'fraction_of_differential'), max_rows=15)
        if run.types is not None and len(run.types):
            doc.table(run.types, 'Gene level and step 1 by explanation type, by mark (genes with the mark changing '
                                 'concordantly, gained and lost) and by chromatin state transition '
                                 '(<code>{}_explanation_types.tsv</code>).'.format(esc(run.name)),
                      formats={c: fmt_int for c in ('genes', 'genes_up', 'genes_down', 'mark_gained', 'mark_lost')},
                      numeric=('genes', 'genes_up', 'genes_down', 'mark_gained', 'mark_lost'))
        if run.roles is not None and len(run.roles):
            doc.section('roles', 'Mark roles per level', level=3)
            doc.figure([roles_figure(run)], 'Association of every mark change with the expression change at the '
                                            'promoter and at each step (Spearman rho with the chromosome block-'
                                            'bootstrap interval, <code>{}_mark_roles.tsv</code>). Filled markers: the '
                                            'role learned at this level equals the declared role; open markers: it '
                                            'differs.'.format(esc(run.name)))
            roles = run.roles
            frame = pd.DataFrame({'mark': roles['mark'], 'level': roles['level'],
                                  'versus': roles.get('versus', ''), 'declared': roles.get('declared_role', ''),
                                  'learned': roles.get('learned_role', ''), 'genes': roles.get('n_genes', np.nan)})

            def interval(value, low, high):
                return '{} ({} to {})'.format(fmt_num(value), fmt_num(low), fmt_num(high)) if not isnan(value) else ''
            for key, label in (('rho', 'rho (interval)'), ('concordance', 'concordance (interval)')):
                if key in roles:
                    frame[label] = [interval(v, lo, hi) for v, lo, hi in
                                    zip(roles[key], roles[key + '_ci_low'], roles[key + '_ci_high'])]
            if 'disagrees' in roles:
                frame['differs'] = roles['disagrees'].map(lambda v: 'yes' if truthy(v) else 'no')
            doc.table(frame, 'Declared and learned role of every mark per level (versus: the quantity the mark change '
                             'is compared with), with the rank correlation and the concordance of strong changes and '
                             'their chromosome block-bootstrap intervals (<code>{}_mark_roles.tsv</code>).'.format(
                                 esc(run.name)),
                      formats={'genes': fmt_int}, numeric=('genes',))
        doc.section('states', 'Chromatin state transitions', level=3)
        if run.state_transitions is not None and len(run.state_transitions):
            table = run.state_transitions
            if 'changed' in table:
                table = table[table['changed'].map(truthy)]
            if 'n_genes' in table:
                table = table.sort_values('n_genes', ascending=False)
            columns = [c for c in ('level', 'transition', 'n_fragments', 'n_genes', 'mean_rna_log2FC', 'ci_low',
                                   'ci_high', 'n_up', 'n_down', 'concordance', 'label') if c in table]
            doc.table(table[columns], 'Changed mark combinations (learned states) at the promoter (step 0) and at the '
                                      'contacting fragments (step 1), by number of genes: fragments, genes, mean '
                                      'expression change of the genes with its interval, genes up and down, and the '
                                      'concordance (<code>{}_state_transitions.tsv</code>). First 20 rows.'.format(
                                          esc(run.name)),
                      formats={'n_fragments': fmt_int, 'n_genes': fmt_int, 'n_up': fmt_int, 'n_down': fmt_int,
                               'mean_rna_log2FC': lambda v: fmt_num(v, signed=True),
                               'ci_low': lambda v: fmt_num(v, signed=True), 'ci_high': lambda v: fmt_num(v, signed=True),
                               'concordance': fmt_num},
                      numeric=('n_fragments', 'n_genes', 'mean_rna_log2FC', 'ci_low', 'ci_high', 'n_up', 'n_down',
                               'concordance'), max_rows=20)
        if run.bivalent is not None and len(run.bivalent):
            threshold = run.command_value('--stateThreshold', 2.0)
            table = run.bivalent
            if 'threshold' in table and (table['threshold'] == threshold).any():
                table = table[table['threshold'] == threshold]
            table = table.sort_values('genes', ascending=False) if 'genes' in table else table
            doc.table(table, 'Declared states (bivalent promoter, poised enhancer) at the presence threshold {:g}: '
                             'transitions at the promoter and at step 1 with the genes up and down and the concordant '
                             'genes (<code>{}_bivalent_transitions.tsv</code>). First 20 rows.'.format(
                                 threshold, esc(run.name)),
                      formats={'threshold': fmt_num, 'genes': fmt_int, 'genes_up': fmt_int, 'genes_down': fmt_int,
                               'concordant': fmt_int, 'supported_concordant': fmt_int,
                               'fraction_concordant': fmt_pct, 'direction': fmt_int},
                      numeric=('threshold', 'genes', 'genes_up', 'genes_down', 'concordant', 'supported_concordant',
                               'fraction_concordant'), max_rows=20)
        if run.edge_chromatin is not None and len(run.edge_chromatin):
            doc.table(run.edge_chromatin, 'Mark change at the two ends of gained and lost interactions against shared '
                                          'interactions of the same distance (egaExplain '
                                          '<code>{}_edge_chromatin.tsv</code>).'.format(esc(run.name)),
                      formats={'n': fmt_int, 'median_change': lambda v: fmt_num(v, 3, True),
                               'shared_median_change': lambda v: fmt_num(v, 3, True), 'concordance_vs_shared': fmt_num,
                               'mannwhitney_p_expected_direction': fmt_p},
                      numeric=('n', 'median_change', 'shared_median_change', 'concordance_vs_shared',
                               'mannwhitney_p_expected_direction'))

    def select_findings(self):
        """Per cascade level the genes with the strongest support: genes first explained there; beyond the mirror first,
        then by q, P and the absolute expression change; at the gene level by the number of explaining levels and the
        absolute expression change. One promoter per gene name."""
        run = self.run
        first = run.first_level()
        levels = run.levels
        chosen = {}
        used = set()
        order = ['gene'] + ['step{}'.format(k) for k in range(1, run.max_hops + 1)]
        for level in order:
            nodes = [n for n, lv in first.items() if lv == level]
            if not nodes:
                continue
            part = levels[levels['node'].isin(nodes)]
            rows = []
            for node, gene_rows in part.groupby('node'):
                here = gene_rows[(gene_rows['cascade_level'] == level) & gene_rows['explains']]
                q = here['q'].astype(float).min() if len(here) else np.nan
                p = here['p'].astype(float).min() if len(here) else np.nan
                beyond = here['beyond_mirror'].map(truthy).any() if len(here) else False
                rows.append(dict(node=node, gene=gene_rows['gene'].iloc[0], q=q, p=p, beyond=beyond,
                                 n_levels=int(gene_rows['explains'].sum()),
                                 rna=abs(float(gene_rows['rna_log2FC'].iloc[0]))))
            frame = pd.DataFrame(rows)
            if not len(frame):
                continue
            has_q = frame['q'].notna().any()
            if level == 'gene':
                frame = frame.sort_values(['n_levels', 'rna'], ascending=[False, False])
                rule = 'by the number of explaining levels, then by the absolute expression change'
            elif has_q:
                frame = frame.assign(nq=frame['q'].fillna(2), np_=frame['p'].fillna(2))
                frame = frame.sort_values(['beyond', 'nq', 'np_', 'rna'], ascending=[False, True, True, False])
                rule = 'beyond the mirror first, then by the q value of the per-gene permutation test'
            else:
                frame = frame.sort_values(['beyond', 'n_levels', 'rna'], ascending=[False, False, False])
                rule = 'beyond the mirror first, then by the number of explaining levels and the absolute expression ' \
                       'change (no per-gene test of this step in the run outputs)'
            picked = []
            for row in frame.itertuples():
                key = row.gene or row.node
                if key in used:
                    continue
                used.add(key)
                picked.append(row)
                if len(picked) >= self.args.top:
                    break
            significant = int((frame['q'] < 0.1).sum()) if has_q else None
            chosen[level] = (picked, rule, len(frame), significant, int(frame['beyond'].sum()))
        return chosen

    def findings_section(self, doc):
        run = self.run
        doc.section('findings', 'Major findings')
        doc.p('For each level of the cascade, the {} genes with the strongest support among the genes first explained '
              'there. Each gene lists every level in rank order with what changed, by how much, the per-gene P and q, '
              'the mirror and the 3D support, followed by the interaction arcs, the neighbourhood graph with the '
              'explaining nodes and its hub.'.format(self.args.top))
        chosen = self.select_findings()
        findings = []
        for level, (picked, rule, n, significant, beyond) in chosen.items():
            title = 'Gene level (promoter)' if level == 'gene' else 'Step {}'.format(level_hops(level))
            doc.section('findings-' + level, title, level=3)
            note = '{:,} genes first explained here; ranked {}.'.format(n, rule)
            if significant is not None:
                note += ' {:,} at q &lt; 0.1{}'.format(significant, '; none reach it, the best-ranked are shown.'
                                                       if significant == 0 else '.')
            if level != 'gene':
                note += ' {:,} beyond the mirror.'.format(beyond)
            doc.p(note, cls='note')
            for row in picked:
                nodes = [row.node]
                doc.heading('{} ({} {})'.format(gene_span(row.gene or row.node), 'log2FC',
                                                fmt_num(run.classes.set_index('node').loc[row.node, 'rna_log2FC'],
                                                        signed=True)), level=4)
                self.gene_finding(doc, row.gene or row.node, nodes, hub=False)
                findings.append((row.gene or row.node, nodes))
        return findings

    def hubs_section(self, doc, findings):
        run = self.run
        doc.section('hubs', 'Hubs and communities')
        doc.p('Hubs of the genes above (the promoter, its partners and the other promoters contacting these partners) '
              'and the communities with a coordinated expression change: which member genes change in the same '
              'direction, which in the opposite direction, and which genes left or joined between {} and {}.'.format(
                  esc(run.names[0]), esc(run.names[1])))
        doc.section('hubs-genes', 'Hubs of the explained genes', level=3)
        for gene, nodes in findings:
            doc.heading('Hub of {}'.format(gene_span(gene)), level=4)
            self.hub_block(doc, gene, nodes)
            self.gene_community_text(doc, gene, nodes)
        doc.section('communities', 'Coordinated communities', level=3)
        table = run.community_expression
        if table is None or not len(table):
            doc.p('No community expression table in the run outputs.', cls='note')
            return
        coordinated = table[table['coherence_q'] < 0.1] if 'coherence_q' in table else table.head(0)
        show = coordinated if len(coordinated) else table.head(self.args.maxCommunities)
        columns = [c for c in ('community_a', 'community_b', 'transition', 'n_genes', 'edges_shared', 'edges_lost',
                               'edges_gained', 'rewired_fraction', 'mean_rna_log2FC', 'coherence_p', 'coherence_q')
                   if c in show]
        doc.table(show[columns], '{} (<code>{}_community_expression.tsv</code>): matched communities of {} and {}, '
                                 'their transition, member genes, interactions shared, lost and gained, the rewired '
                                 'fraction, the mean expression change and the coherence test against random gene sets '
                                 'of the same size.'.format(
                                     'Communities with a coordinated expression change (q &lt; 0.1)' if len(coordinated)
                                     else 'No community reaches q &lt; 0.1; the {} with the lowest P are shown'.format(
                                         len(show)), esc(run.name), esc(run.names[0]), esc(run.names[1])),
                  formats={'rewired_fraction': fmt_pct, 'mean_rna_log2FC': lambda v: fmt_num(v, signed=True),
                           'coherence_p': fmt_p, 'coherence_q': fmt_p},
                  numeric=('community_a', 'community_b', 'n_genes', 'edges_shared', 'edges_lost', 'edges_gained',
                           'rewired_fraction', 'mean_rna_log2FC', 'coherence_p', 'coherence_q'))
        summaries = {}
        if run.explain_communities is not None and len(run.explain_communities):
            summaries = dict(zip(run.explain_communities['community'], run.explain_communities['summary']))
        for k, row in enumerate(show.head(self.args.maxCommunities).itertuples()):
            ca, cb = int(row.community_a), int(row.community_b)
            doc.heading('Community {} of {} / {} of {} ({})'.format(ca, esc(run.names[0]), cb, esc(run.names[1]),
                                                                    esc(row.transition)), level=4)
            if cb in summaries:
                doc.p('egaExplain: {}'.format(esc(summaries[cb])), cls='note')
            self.community_block(doc, ca, cb, info={'transition': row.transition})

    def methods(self, doc):
        run = self.run
        doc.section('methods', 'Methods, parameters and versions')
        doc.p('The report reads the outputs of egaDifferential and egaExplain listed in Table {} and draws the figures '
              'from the union graph with the drawing functions of egaPlotDifferential (squares: baited promoters, '
              'circles: other ends; line width: interaction strength relative to the calling threshold; interaction '
              'change red-blue by log2 fold change, gained and lost dashed). Per-gene levels: {}. Thresholds of the '
              'level calls: differential genes at |log2 fold change| {:g}, explaining mark changes at |log2 fold '
              'change| {:g} ({}). The neighbourhood of a gene reaches its deepest explaining step and is reduced to at '
              'most {} nodes (the promoter, its strongest partners and the paths to the explaining nodes); explaining '
              'nodes at step k are the fragments at exactly k steps with a concordant change of at least the threshold '
              '(log2 (B + 1) / (A + 1) of the node signals in the graph), the 12 strongest per step.'.format(
                  doc.table_count + 1, esc(run.levels_source), run.rna_threshold, run.chromatin_threshold,
                  'from the run log' if run.command else 'egaDifferential defaults, no command line in the run log',
                  self.args.maxNodes))
        rows = []
        for label, path in sorted(run.used.items()):
            try:
                stamp = datetime.datetime.fromtimestamp(os.path.getmtime(path)).strftime('%Y-%m-%d %H:%M')
            except OSError:
                stamp = ''
            rows.append((label, path, stamp))
        doc.table(pd.DataFrame(rows, columns=['file', 'path', 'modified']), 'Files read for this report and their '
                                                                            'modification times.', long=('path',))
        if run.command:
            doc.p('Command line: <code>{}</code>'.format(esc(run.command)))
        info = [line for line in run.log_lines if ' INFO: ' in line or 'WARNING' in line]
        if info:
            doc.table(pd.DataFrame({'run log': info[:40]}), 'Messages of the run log (<code>{}</code>), first 40 '
                                                            'lines.'.format(esc(', '.join(run.log_paths))),
                      long=('run log',))
        else:
            doc.p('No run log found (expected at <code>{}.log</code>; see --log).'.format(esc(run.prefix)), cls='note')
        import scipy
        versions = [('epigeneticgraphanalyzer (egaReport)', __version__), ('Python', platform.python_version()),
                    ('numpy', np.__version__), ('pandas', pd.__version__), ('networkx', nx.__version__),
                    ('matplotlib', matplotlib.__version__), ('scipy', scipy.__version__),
                    ('report date', datetime.date.today().isoformat())]
        doc.table(pd.DataFrame(versions, columns=['software', 'version']), 'Software versions used to write the '
                                                                          'report.')

    # gene and community reports -------------------------------------------------------------------------------------
    def gene_report(self, gene):
        run = self.run
        nodes, differential = run.gene_nodes(gene)
        doc = Document('{}: gene report, {} vs {}'.format(gene, *run.names), self.keep_figures)
        doc.section('gene-summary', 'Summary')
        if not nodes:
            doc.p('{} is not a baited promoter of the union graph of this run.'.format(gene_span(gene)))
            return doc
        features = run.features.set_index('node') if run.features is not None else None
        rows = []
        for node in nodes:
            data = run.graph.nodes[node] if node in run.graph else {}
            counts = {'gained': 0, 'lost': 0, 'shared': 0}
            if node in run.graph:
                for w in run.graph.neighbors(node):
                    counts[run.graph.edges[node, w].get('status', 'shared')] += 1
            rna = data.get('rna_log2FC', features.loc[node, 'rna_log2FC'] if features is not None and node in
                           features.index else np.nan)
            rows.append(dict(fragment=node, genes=data.get('genes', ''), log2FC=rna, gained=counts['gained'],
                             lost=counts['lost'], shared=counts['shared'],
                             first_level=run.first_level().get(node, 'not differential')))
        doc.table(pd.DataFrame(rows), 'Promoter fragments of {g}: genes on the fragment, expression log2 fold change '
                                      '{b} / {a}, interactions gained, lost and shared, and the first explaining level.'
                  .format(g=gene_span(gene), a=esc(run.names[0]), b=esc(run.names[1])),
                  formats={'log2FC': lambda v: fmt_num(v, signed=True)},
                  numeric=('log2FC', 'gained', 'lost', 'shared'))
        if not differential:
            doc.p('{} is not differential at the thresholds of the run; the levels are not evaluated.'.format(
                gene_span(gene)), cls='note')
        doc.section('gene-levels', 'Explanation by level')
        level_nodes = [n for n in nodes if n in set(run.levels['node'])]
        if level_nodes:
            self.gene_finding(doc, gene, level_nodes, figures=False, hub=False)
            declared = '; '.join('{} {}'.format(m, r) for m, r in run.declared_roles().items())
            evidence_table(doc, run.levels[run.levels['node'].isin(level_nodes)], gene, run.names, declared)
            extra = run.levels[run.levels['node'].isin(level_nodes)]['extra'] if 'extra' in run.levels else []
            keys = sorted({k for e in extra if isinstance(e, dict) for k in e})
            if keys:
                frame = pd.DataFrame([dict(level=lv, **{k: e.get(k, '') for k in keys})
                                      for lv, e in zip(run.levels.loc[run.levels['node'].isin(level_nodes),
                                                                      'level_label'], extra)])
                doc.table(frame, 'Further columns of the all-level table for {}.'.format(gene_span(gene)))
        doc.section('gene-figures', 'Interactions and neighbourhood')
        rows = run.levels[run.levels['node'].isin(nodes)]
        self.gene_figures(doc, gene, nodes, rows)
        doc.section('gene-states', 'Chromatin states')
        self.state_block(doc, gene, nodes)
        doc.section('gene-support', '3D support per step')
        self.support_block(doc, gene, nodes)
        doc.section('gene-hub', 'Hub and community')
        self.hub_block(doc, gene, nodes)
        community = self.gene_community_text(doc, gene, nodes)
        if community is not None:
            ca, cb = int(community.get('community_a', -1)), int(community.get('community_b', -1))
            transition = community.get('community_a_transition', '')
            self.community_block(doc, ca, cb, focus=nodes,
                                 info={'transition': '' if isnan(transition) else transition})
        self.methods(doc)
        return doc

    def state_block(self, doc, gene, nodes):
        run = self.run
        classes = run.classes.set_index('node')
        rows = []
        for node in nodes:
            if node not in classes.index:
                continue
            r = classes.loc[node]
            fields = {}
            for column in classes.columns:
                if not (column.startswith('promoter_') or column.startswith('step1_')):
                    continue
                match = re.fullmatch(r'(.*?)_(mean|max1kb)', column)
                base, variant = (match.group(1), match.group(2)) if match else (column, 'value')
                fields.setdefault(base, {})[variant] = r[column]
            for base, values in fields.items():
                rows.append(dict(fragment=node, field=base, **values))
        if not rows:
            doc.p('No chromatin state calls for {} (not differential, or states not callable).'.format(
                gene_span(gene)), cls='note')
            return
        frame = pd.DataFrame(rows)
        variants = [c for c in ('value', 'mean', 'max1kb') if c in frame]

        def show(v):
            if isinstance(v, (bool, np.bool_)):
                return 'yes' if v else 'no'
            if isinstance(v, (float, np.floating)) and not isnan(v):
                return fmt_num(v)
            return '' if isnan(v) else str(v)

        fragments = list(dict.fromkeys(frame['fragment']))
        table = {}
        for row in frame.itertuples(index=False):
            data = row._asdict()
            for variant in variants:
                label = data['field'] + ('' if variant == 'value' else ' ({})'.format(
                    'mean' if variant == 'mean' else 'max 1 kb'))
                value = data.get(variant)
                if not isinstance(value, str) and isnan(value):
                    continue
                table.setdefault(label, {})[data['fragment']] = show(value)
        out = pd.DataFrame([dict(field=k, **{f: v.get(f, '') for f in fragments}) for k, v in table.items()])
        out = out[out[fragments].apply(lambda r: any(x != '' for x in r), axis=1)]
        doc.table(out, 'Chromatin state calls of {} at the promoter and at the contacting fragments (step 1) from '
                       '<code>{}_gene_classes.tsv</code>, one column per promoter fragment: mark signals per condition '
                       '(suffix _A: {}, _B: {}), states, transition and whether it explains the expression change{}.'
                  .format(gene_span(gene), esc(run.name), esc(run.names[0]), esc(run.names[1]),
                          '; calls on the mean signal over the fragment and on the maximum over 1 kb windows'
                          if len(variants) > 1 else ''), nowrap=('field',), long=tuple(fragments))

    def support_block(self, doc, gene, nodes):
        run = self.run
        rows = run.levels[run.levels['node'].isin(nodes)]
        steps = rows[rows['level'].str.fullmatch(r'step\d+', na=False)]
        if not len(steps) and not run.support:
            doc.p('The run has no 3D proximity support (egaDifferential without --proximitySupport).', cls='note')
            return
        frame = pd.DataFrame({'fragment': steps['node'], 'step': steps['level'].map(level_hops),
                              'explains': steps['explains'].map({True: 'yes', False: 'no'}),
                              'nodes': steps['nodes'], 'support': steps['support'], 'interval': steps['support_ci']})
        if frame['support'].isna().all():
            doc.p('The run has no 3D proximity support (egaDifferential without --proximitySupport); nodes per step '
                  'below.', cls='note')
            frame = frame.drop(columns=['support', 'interval'])
        doc.table(frame, '3D support of every step for {}: probability that graph neighbours at this many steps are '
                         'closer in the nucleus than distance-matched pairs (0.5 without relation; step curve of the '
                         'run), the nodes at the step and whether the step explains the gene.'.format(gene_span(gene)),
                  formats={'nodes': fmt_int, 'support': fmt_num}, numeric=('step', 'nodes', 'support'))

    def community_report(self, community, side):
        run = self.run
        doc = Document('Community {} ({}): {} vs {}'.format(community, run.names[0 if side == 'A' else 1],
                                                             *run.names), self.keep_figures)
        doc.section('community', 'Community {} of {}'.format(community, run.names[0 if side == 'A' else 1]))
        table = run.community_expression
        info, ca, cb = {}, -1, -1
        if table is not None and len(table):
            column = 'community_a' if side == 'A' else 'community_b'
            hit = table[table[column] == community]
            if len(hit):
                row = hit.iloc[0]
                ca, cb = int(row['community_a']), int(row['community_b'])
                info = {'transition': row.get('transition', '')}
                doc.table(hit, 'Expression coherence of the community (<code>{}_community_expression.tsv</code>).'
                          .format(esc(run.name)),
                          formats={'rewired_fraction': fmt_pct, 'mean_rna_log2FC': lambda v: fmt_num(v, signed=True),
                                   'coherence_p': fmt_p, 'coherence_q': fmt_p})
        if ca < 0 and cb < 0 and run.gene_communities is not None:
            rows = self.community_rows(side, community)
            if rows is not None and len(rows):
                other = 'community_b' if side == 'A' else 'community_a'
                partner = rows[other][rows[other] >= 0].mode()
                partner = int(partner.iloc[0]) if len(partner) else -1
                ca, cb = (community, partner) if side == 'A' else (partner, community)
                transition = rows['community_a_transition'].dropna()
                info = {'transition': transition.iloc[0] if len(transition) else ''}
        if ca < 0 and cb < 0:
            doc.p('Community {} not found in the run outputs.'.format(community))
            return doc
        if run.community_transitions is not None:
            ct = run.community_transitions
            hit = ct[((ct['side'] == 'A') & (ct['community_a'] == ca)) | ((ct['side'] == 'B') & (ct['community_b']
                                                                                                 == cb))]
            if len(hit):
                doc.table(hit, 'Transition records of the community (<code>{}_community_transitions.tsv</code>).'.format(
                    esc(run.name)), formats={'fraction_present_in_b': fmt_pct, 'fraction_present_in_a': fmt_pct,
                                             'jaccard': fmt_num})
        self.community_block(doc, ca, cb, info=info)
        self.methods(doc)
        return doc


# --------------------------------------------------------------------------------------------------------------------
# output

def find_chromium(explicit=None):
    for candidate in ([explicit] if explicit else []) + list(CHROMIUM_CANDIDATES):
        if candidate and (os.path.exists(candidate) or shutil.which(candidate)):
            return candidate if os.path.exists(candidate) else shutil.which(candidate)
    return None


def chromium_pdf(chromium, html_path, pdf_path, timeout=300):
    """Print an HTML file to PDF with headless Chromium; a snap Chromium cannot read or write some directories
    (e.g. /tmp), so a copy in its own directory is printed when the direct call writes nothing."""
    def run(source, target):
        command = [chromium, '--headless', '--disable-gpu', '--no-sandbox', '--no-pdf-header-footer',
                   '--run-all-compositor-stages-before-draw', '--virtual-time-budget=10000',
                   '--print-to-pdf={}'.format(target), 'file://' + os.path.abspath(source)]
        try:
            subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=timeout, check=False)
        except (OSError, subprocess.TimeoutExpired) as error:
            log.warning('chromium failed: %s', error)
        return os.path.exists(target) and os.path.getsize(target) > 0

    if os.path.exists(pdf_path):
        os.remove(pdf_path)
    if run(html_path, pdf_path):
        return True
    staging_root = os.path.join(os.path.expanduser('~'), 'snap', 'chromium', 'common')
    if not os.path.isdir(staging_root):
        return False
    staging = tempfile.mkdtemp(prefix='egaReport_', dir=staging_root)
    try:
        source = os.path.join(staging, os.path.basename(html_path))
        shutil.copy(html_path, source)
        target = os.path.join(staging, os.path.basename(pdf_path))
        if run(source, target):
            shutil.move(target, pdf_path)
            return True
        return False
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def matplotlib_pdf(doc, pdf_path, heading):
    """Plain PDF of a report with matplotlib: text pages (A4) and one page per figure."""
    from matplotlib.backends.backend_pdf import PdfPages
    width, height = 8.27, 11.69
    with PdfPages(pdf_path) as pdf:
        lines = []

        def flush():
            while lines:
                fig = new_figure(width, height)
                y = 0.95
                while lines and y > 0.05:
                    text, size, weight, family = lines.pop(0)
                    fig.text(0.07, y, text, fontsize=size, fontweight=weight, family=family, va='top')
                    y -= size / 72.0 * 1.45 / height
                pdf.savefig(fig)

        lines.append((heading, 14, 'bold', 'sans-serif'))
        lines.append(('PDF written with matplotlib (Chromium not available); the HTML report holds the same content.',
                      8, 'normal', 'sans-serif'))
        for kind, value, level in doc.text_blocks:
            if kind == 'heading':
                lines.append(('', 6, 'normal', 'sans-serif'))
                lines.append((value, 12 if level <= 2 else 10, 'bold', 'sans-serif'))
            elif kind == 'text':
                for line in textwrap.wrap(value, 110):
                    lines.append((line, 8, 'normal', 'sans-serif'))
            elif kind == 'table':
                for row in value.split('\n')[:40]:
                    for line in textwrap.wrap(row, 130) or ['']:
                        lines.append((line, 6, 'normal', 'monospace'))
            elif kind == 'figure':
                flush()
                figures, caption = doc.figures[value]
                for fig in figures:
                    pdf.savefig(fig, bbox_inches='tight')
                lines.append((textwrap.shorten(caption, 600), 7, 'normal', 'sans-serif'))
        flush()


def write_document(doc, path_stem, heading, meta, args):
    html_path = path_stem + '.html'
    with open(html_path, 'w') as handle:
        handle.write(doc.render(heading, meta))
    pdf_path = path_stem + '.pdf'
    engine = None
    if args.pdf in ('auto', 'chromium'):
        chromium = find_chromium(args.chromium)
        if chromium and chromium_pdf(chromium, html_path, pdf_path):
            engine = 'chromium'
        elif args.pdf == 'chromium':
            raise RuntimeError('Chromium could not print {}'.format(html_path))
        else:
            log.warning('Chromium not available or failed; writing the PDF with matplotlib (plain layout)')
    if engine is None and args.pdf in ('auto', 'matplotlib'):
        if doc.keep_figures:
            matplotlib_pdf(doc, pdf_path, heading)
            engine = 'matplotlib'
        else:
            log.warning('no PDF for %s: figures were not kept for the matplotlib fallback', html_path)
    log.info('wrote %s%s', html_path, ' and {} ({})'.format(pdf_path, engine) if engine else '')
    for fig_list, _ in doc.figures:
        for fig in fig_list:
            fig.clear()
    return html_path, (pdf_path if engine else None), engine


def meta_block(run):
    return ('<p>Output prefix <code>{}</code>; conditions {} (A) and {} (B).</p><p>Written {} with egaReport {}.</p>'
            .format(esc(run.prefix), esc(run.names[0]), esc(run.names[1]), datetime.date.today().isoformat(),
                    esc(__version__)))


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    os.makedirs(args.outputDir, exist_ok=True)
    run = Run(args.prefix, args.explainPrefix, args.names, args.rnaThreshold, args.chromatinThreshold, args.log,
              args.name)
    log.info('%s: conditions %s, marks %s, %d differential genes, levels from %s', run.name, run.names, run.marks,
             len(run.classes), run.levels_source)
    # keep the figures for the matplotlib PDF when Chromium is not used or may fail
    chromium = find_chromium(args.chromium) if args.pdf in ('auto', 'chromium') else None
    keep = args.pdf == 'matplotlib' or (args.pdf == 'auto' and chromium is None)
    report = Report(run, args, keep)
    outputs = []
    if args.gene:
        pages = []
        for gene in args.gene:
            doc = report.gene_report(gene)
            stem = os.path.join(args.outputDir, '{}_gene_{}'.format(slug(run.name), slug(gene)))
            heading = 'Gene report: {} ({} vs {})'.format(gene_span(gene), esc(run.names[0]), esc(run.names[1]))
            outputs.append(write_document(doc, stem, heading, meta_block(run), args))
            pages.append((gene, os.path.basename(stem)))
        if len(pages) > 1:
            index = Document('{} vs {}: gene reports'.format(*run.names), keep)
            index.section('genes', 'Gene reports')
            index.html('<ul class="index-list">{}</ul>'.format(''.join(
                '<li><a href="{s}.html">{g}</a> (<a href="{s}.pdf">PDF</a>)</li>'.format(s=esc(stem),
                                                                                         g=gene_span(gene))
                for gene, stem in pages)))
            stem = os.path.join(args.outputDir, '{}_genes_index'.format(slug(run.name)))
            path = stem + '.html'
            with open(path, 'w') as handle:
                handle.write(index.render('Gene reports: {} vs {}'.format(esc(run.names[0]), esc(run.names[1])),
                                          meta_block(run)))
            outputs.append((path, None, None))
    elif args.community is not None:
        doc = report.community_report(args.community, args.communitySide)
        stem = os.path.join(args.outputDir, '{}_community_{}{}'.format(slug(run.name), args.communitySide,
                                                                       args.community))
        outputs.append(write_document(doc, stem, 'Community report: {} ({} vs {})'.format(
            args.community, esc(run.names[0]), esc(run.names[1])), meta_block(run), args))
    else:
        doc = report.main_report()
        stem = os.path.join(args.outputDir, '{}_report'.format(slug(run.name)))
        outputs.append(write_document(doc, stem, 'Differential graph report: {} vs {}'.format(
            esc(run.names[0]), esc(run.names[1])), meta_block(run), args))
    for html_path, pdf_path, engine in outputs:
        sys.stdout.write('{}\n'.format(html_path))
        if pdf_path:
            sys.stdout.write('{} ({})\n'.format(pdf_path, engine))
    return outputs


if __name__ == '__main__':
    main()
