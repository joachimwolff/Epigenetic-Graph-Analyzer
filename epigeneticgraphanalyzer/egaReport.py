"""egaReport: a self-contained HTML and PDF report of a differential run (egaDifferential, egaExplain).

The report only reads the outputs of egaDifferential (<prefix>_*) and, when present, of egaExplain
(<explain prefix>_*). Every number shown is taken from these files or counted from them; the figures are drawn from
<prefix>_union.graphml with the drawing functions of egaPlotDifferential.

The report is written for a reader of the biology: the body names places by their genes ('the Hoxa2 promoter',
'element 1.0 Mb upstream of Hoxa3 (near Skap2)'), levels by readable labels ('Promoter: chromatin state', '2 steps
away') and numbers with at most three significant digits. File paths, parameters, versions, coordinates and the full
numeric tables are kept in the appendix 'Data and methods' at the end (a collapsible section in HTML, the last pages
of the PDF).

Main report (default): title block, key findings (sentences generated from the run), the explanation cascade, one
finding card per selected gene (the 'at a glance' paragraph, the explanation ladder and the interaction arcs, with a
link to the gene report), hubs and coordinated communities in prose with figures, appendix.

Gene report (--gene NAME [NAME ...]): title block, 'at a glance' paragraph with the explanation ladder, interactions
(arcs and neighbourhood graph), evidence by level per promoter fragment, chromatin states (with re-ChIP where
available), hub and community, appendix; an index page when several genes are given.
Community report (--community ID): one page for a community of condition B (or A with --communitySide A).

Declared mark roles are shown throughout; the result under the roles learned from the data is mentioned only where it
differs from the declared one. Upstream and downstream need the strand of the report gene (--annotation, a GTF or
BED file); without it the distance alone is given.

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
# readable names of the conditions of common data sets, and what they are (tissues, cell types)
KNOWN_CONDITIONS = {'Md': ('mandibular arch', 'tissue'), 'PA2': ('second pharyngeal arch', 'tissue'),
                    'hESC': ('human embryonic stem cells', 'cell type'),
                    'hNEC': ('hESC-derived neuroectoderm', 'cell type'),
                    'siMM': ('Kasumi-1, control siRNA', 'condition'),
                    'siRE': ('Kasumi-1, RUNX1-ETO knockdown', 'condition')}
ACCENT = '#2c5d8a'
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
    parser.add_argument('--labels', nargs=2, default=None,
                        help='Readable names of the conditions A B for the title (default: known names of common '
                             'data sets, else the condition names)')
    parser.add_argument('--annotation', default=None,
                        help='Gene annotation (GTF, GTF.gz or BED with strand) for upstream / downstream of the '
                             'report gene; without it only the distance is given')
    parser.add_argument('--rechip', default=None,
                        help='Re-ChIP signal per bait (default: <prefix dir>/rechip/rechip_baits.tsv.gz if its '
                             'columns match the condition names)')
    parser.add_argument('--top', type=int, default=2,
                        help='Genes shown per cascade level in the gene cards of the main report')
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


MINUS = '−'


def round_sig(value, digits=2):
    value = float(value)
    if value == 0 or not np.isfinite(value):
        return value
    return round(value, digits - 1 - int(np.floor(np.log10(abs(value)))))


def sig(value, digits=2):
    """A number with at most `digits` significant digits (thousands separated), no exponent above 0.001."""
    if isnan(value):
        return ''
    value = round_sig(value, digits)
    if value == 0:
        return '0'
    decimals = max(0, digits - 1 - int(np.floor(np.log10(abs(value)))))
    text = '{:,.{d}f}'.format(value, d=decimals)
    return text.replace('-', MINUS)


def fmt_log2(value, signed=True):
    """A log2 change with one decimal and a true minus sign ('+3.4', '−1.2')."""
    if isnan(value):
        return ''
    value = float(value)
    if round(value, 1) == 0:
        return '0.0'
    text = '{:+.1f}'.format(value) if signed else '{:.1f}'.format(abs(value))
    return text.replace('-', MINUS)


def fold(value):
    """Plain factor of a log2 change ('37', '1,560')."""
    factor = 2 ** abs(float(value))
    if factor < 10:
        return '{:.1f}'.format(factor)
    if factor < 100:
        return '{:.0f}'.format(factor)
    return '{:,.0f}'.format(round_sig(factor, 3))


def fmt_pq(value, markup=True):
    """P or q: two significant digits, below 0.001 as 8.5 x 10^-5 (superscript in HTML)."""
    if isnan(value):
        return ''
    value = float(value)
    if value == 0:
        return '0'
    if value < 0.001:
        mantissa, exponent = '{:.1e}'.format(value).split('e')
        exponent = int(exponent)
        power = '10<sup>{}{}</sup>'.format(MINUS if exponent < 0 else '', abs(exponent)) if markup else \
            '10^{}'.format(exponent)
        return '{} × {}'.format(mantissa, power)
    text = sig(value, 2)
    return text.rstrip('0').rstrip('.') if '.' in text else text


def join_words(items, conjunction='and'):
    items = [i for i in items if i]
    if len(items) <= 1:
        return ''.join(items)
    return '{} {} {}'.format(', '.join(items[:-1]), conjunction, items[-1])


def gene_names(genes):
    """Gene names of a fragment, informative names first (without Gm..., predicted loci and the like)."""
    names = [g for g in str(genes if not isnan(genes) else '').split(',') if g and g != 'nan']
    informative = [g for g in names if not g.startswith(IGNORED_GENE_PREFIXES)]
    return informative or names


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

SV_LEVELS = ('structural_variant', 'copy_number')


def level_order(level):
    """Rank order of a level: promoter marks, contacts, promoter state, structural variant, copy number, step 1, step 1
    state, step 2, ..."""
    fixed = {'promoter_marks': 0, 'contacts': 1, 'promoter_state': 2, 'structural_variant': 2.4, 'copy_number': 2.6,
             'step1': 3, 'step1_state': 4}
    if level in fixed:
        return fixed[level]
    match = re.fullmatch(r'step(\d+)', str(level))
    if match:
        return 3 + 2 * (int(match.group(1)) - 1)
    return 1000


def level_label(level):
    """Readable label of a level, as in the tables and the ladder."""
    labels = {'promoter_marks': 'Promoter: histone marks', 'contacts': 'Promoter: own contacts',
              'promoter_state': 'Promoter: chromatin state', 'step1': 'Direct partners (1 step): histone marks',
              'step1_state': 'Direct partners: chromatin state',
              'structural_variant': 'Promoter: structural variant nearby', 'copy_number': 'Promoter: copy number'}
    if level in labels:
        return labels[level]
    match = re.fullmatch(r'step(\d+)', str(level))
    if match:
        return '{} steps away'.format(match.group(1))
    return str(level).replace('_', ' ')


def level_phrase(level):
    """The level inside a sentence."""
    phrases = {'promoter_marks': 'the histone marks at the promoter', 'contacts': 'the own contacts of the promoter',
               'promoter_state': 'the chromatin state of the promoter',
               'step1': 'the histone marks of the direct partners',
               'step1_state': 'the chromatin state of the direct partners',
               'structural_variant': 'a structural variant near the promoter',
               'copy_number': 'the copy number of the promoter region'}
    if level in phrases:
        return phrases[level]
    match = re.fullmatch(r'step(\d+)', str(level))
    if match:
        return 'the nodes {} steps away'.format(match.group(1))
    return str(level).replace('_', ' ')


def cascade_level(level):
    """Cascade level (gene, step1, step2, ...) of an explanation level."""
    if level in ('promoter_marks', 'contacts', 'promoter_state'):
        return 'gene'
    if level in SV_LEVELS:
        return 'structural_variant'
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
               'step1marks': 'step1', 'step1state': 'step1_state', 'indirect': 'step2',
               'structuralvariant': 'structural_variant', 'sv': 'structural_variant', 'copynumber': 'copy_number',
               'cn': 'copy_number'}
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
        # roles over all nodes; the roles learned per element type of the nodes are kept apart
        self.roles_by_type = None
        if self.roles is not None and 'element_type' in self.roles:
            kind = self.roles['element_type'].fillna('all')
            self.roles_by_type = self.roles[kind != 'all'].reset_index(drop=True)
            self.roles = self.roles[kind == 'all'].reset_index(drop=True)
        self.elements = self.table('_elements.tsv')
        self.element_classes = self.table('_element_classes.tsv')
        self.element_transitions = self.table('_element_transitions.tsv')
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
        self.labels = [KNOWN_CONDITIONS.get(n, (n, None))[0] for n in self.names]
        kinds = {KNOWN_CONDITIONS.get(n, (None, None))[1] for n in self.names}
        self.kind = kinds.pop() if len(kinds) == 1 and None not in kinds else 'condition'
        self.rechip = None
        self.rechip_path = None
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

    def set_labels(self, labels):
        if labels:
            self.labels = list(labels)
            self.kind = 'condition'

    def condition_title(self):
        """'mandibular arch (Md) vs second pharyngeal arch (PA2)', or the names as given."""
        parts = ['{} ({})'.format(label, name) if label != name else name
                 for label, name in zip(self.labels, self.names)]
        return '{} vs {}'.format(*parts)

    def read_rechip(self, path=None):
        """Re-ChIP signal per bait (columns reChIP_<condition>_<mean|max1kb>); only when they match the names."""
        path = path or os.path.join(os.path.dirname(os.path.abspath(self.prefix)), 'rechip', 'rechip_baits.tsv.gz')
        if not os.path.exists(path):
            return None
        with gzip.open(path, 'rt') if path.endswith('.gz') else open(path) as handle:
            header = handle.readline().rstrip('\n').split('\t')
        if not all('reChIP_{}_mean'.format(n) in header for n in self.names) or 'node' not in header:
            return None
        self.rechip = pd.read_csv(path, sep='\t').set_index('node')
        self.rechip_path = path
        self.used[os.path.basename(path)] = path
        return self.rechip

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

    def element_lookup(self, roles='declared', presence='mean'):
        """node -> (class in A, class in B, TSS status, TSS genes) from <prefix>_elements.tsv (declared roles, fragment
        mean by default); {} without it."""
        table = self.elements
        column = 'class_{}_{}'.format(roles, presence)
        if table is None or column not in table or 'condition' not in table:
            return {}
        a = table[table['condition'].astype(str) == str(self.names[0])].set_index('node')
        b = table[table['condition'].astype(str) == str(self.names[1])].set_index('node')
        b_class = b[column].reindex(a.index)
        genes = a['tss_genes'].fillna('').astype(str) if 'tss_genes' in a else pd.Series('', index=a.index)
        return {n: (ca, cb, st, gs) for n, ca, cb, st, gs in zip(a.index, a[column], b_class, a['tss_status'], genes)}

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
                 strongest_node=col(('strongest_node',), 'strongest_node_' + roles if roles else None),
                 state_fragment=col(('state_fragment',), 'state_fragment_' + key if key else None),
                 signs=col(('signs',), 'signs_' + roles if roles else None),
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
            if level in SV_LEVELS:
                detail = value(r, 'detail')
                row['changes'] = row['all_changes'] = '' if isnan(detail) else str(detail)
                row['detail'] = '' if row['explains'] else result
            elif level == 'contacts':
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
            # structured fields for the prose of the report
            row['result'] = result
            row['mark_values'] = marks
            row['mirror_values'] = mirror
            row['concordant_marks'] = [m for m in str(value(r, 'concordant') if not isnan(value(r, 'concordant'))
                                                      else '').split(',') if m]
            row['n_concordant'] = value(r, 'n_concordant')
            row['strongest_node'] = value(r, 'strongest_node') if not isnan(value(r, 'strongest_node')) else ''
            row['strongest_text'] = value(r, 'strongest') if not isnan(value(r, 'strongest')) else ''
            row['transition'] = value(r, 'transition') if not isnan(value(r, 'transition')) else ''
            row['state_fragment'] = value(r, 'state_fragment') if not isnan(value(r, 'state_fragment')) else ''
            row['detail_raw'] = value(r, 'detail') if not isnan(value(r, 'detail')) else ''
            row['contacts'] = tuple(value(r, k) for k in ('gained', 'lost', 'shared', 'contact'))
            row['signs_declared'] = value(r, 'signs') if not isnan(value(r, 'signs')) else ''
            row['learned_explains'] = row['learned_result'] = row['learned_p'] = row['learned_q'] = np.nan
            row['signs_learned'] = ''
            if learned_key and 'explains_' + learned_key in raw:
                row['learned_explains'] = truthy(r['explains_' + learned_key])
                learned_result = r.get('result_' + learned_key, '')
                row['learned_result'] = learned_result if isinstance(learned_result, str) else ''
                row['learned_p'], row['learned_q'] = r.get('p_learned', np.nan), r.get('q_learned', np.nan)
                signs = r.get('signs_learned', '')
                row['signs_learned'] = signs if isinstance(signs, str) else ''
            row['alt_explains'] = np.nan
            row['alt_result'] = row['alt_transition'] = row['alt_detail'] = ''
            if max_key and 'explains_' + max_key in raw:
                row['alt_explains'] = truthy(r['explains_' + max_key])
                for field in ('result', 'transition', 'detail'):
                    other = r.get('{}_{}'.format(field, max_key), '')
                    row['alt_' + field] = other if isinstance(other, str) else ''
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
                    mirror='', sup=np.nan, sup_ci='', detail='', roles_level=None, contact=False, **more):
                declared, learned = self.roles_text(roles_level, contact) if roles_level else ('', '')
                fields = dict(result='concordant' if explains else 'below threshold', mark_values={},
                              mirror_values={}, concordant_marks=[], n_concordant=np.nan, strongest_node='',
                              strongest_text='', transition='', state_fragment='', detail_raw='',
                              contacts=(np.nan,) * 4, signs_declared='', signs_learned='', learned_explains=np.nan,
                              learned_result='', learned_p=np.nan, learned_q=np.nan, alt_explains=np.nan,
                              alt_result='', alt_transition='', alt_detail='')
                fields.update(more)
                rows.append(dict(base, level=level, level_label=level_label(level), cascade_level=cascade_level(level),
                                 order=level_order(level), explains=bool(explains), changes=changes,
                                 all_changes=all_changes, nodes=nodes, p=p, q=q, beyond_mirror=beyond, mirror=mirror,
                                 support=sup, support_ci=sup_ci, declared=declared, learned=learned, detail=detail,
                                 **fields))

            # promoter marks
            values = {m: r.get(m + '_promoter', np.nan) for m in self.marks}
            hit = [m for m in self.marks if self.concordant(values[m], self.signs.get(m, 1), direction)]
            add('promoter_marks', bool(hit), ', '.join('{} {}'.format(m, fmt_num(values[m], signed=True)) for m in hit),
                ', '.join('{} {}'.format(m, fmt_num(v, signed=True)) for m, v in values.items() if not isnan(v)),
                nodes=1, sup=1.0 if support else np.nan, sup_ci='NA' if support else '', roles_level='promoter',
                mark_values={m: v for m, v in values.items() if not isnan(v)}, concordant_marks=hit)
            # own interactions of the promoter
            contact = r.get('contact_direct', np.nan)
            counts = [r.get('interactions_' + s, np.nan) for s in ('gained', 'lost', 'shared')]
            text = '{} gained, {} lost, {} shared'.format(*[fmt_int(c) for c in counts]) \
                if not any(isnan(c) for c in counts) else ''
            text += '{}summed strength log2 {}'.format('; ' if text else '', fmt_num(contact, signed=True)) \
                if not isnan(contact) else ''
            add('contacts', self.concordant(contact, 1, direction), text, text,
                nodes=sum(c for c in counts if not isnan(c)) if not all(isnan(c) for c in counts) else np.nan,
                roles_level='promoter', contact=True, contacts=tuple(counts) + (contact,))
            # promoter state
            transition = r.get('promoter_transition', np.nan)
            if not isnan(transition) and transition != 'states not callable':
                explains = truthy(r.get('promoter_explains')) or truthy(r.get('promoter_unresolved'))
                kind = 'poised, trigger not measured' if truthy(r.get('promoter_unresolved')) else \
                    ('state transition' if truthy(r.get('promoter_explains')) else '')
                add('promoter_state', explains, '{}{}'.format(transition, ' ({})'.format(kind) if kind else ''),
                    str(r.get('promoter_state_summary', '') if not isnan(r.get('promoter_state_summary')) else ''),
                    nodes=1, roles_level='promoter', transition=str(transition),
                    result='poised in both' if truthy(r.get('promoter_unresolved')) else
                    ('concordant transition' if truthy(r.get('promoter_explains')) else 'unchanged'),
                    detail_raw=str(r.get('promoter_state_summary', '') if not isnan(r.get('promoter_state_summary'))
                                   else ''))
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
                    detail=detail, roles_level='step{}'.format(k),
                    mark_values={m: v for m, v in values.items() if not isnan(v)},
                    mirror_values={m: v for m, v in mirrors.items() if not isnan(v)}, concordant_marks=hit,
                    n_concordant=r.get('step1_concordant_fragments', np.nan) if k == 1 else np.nan,
                    result='concordant' if hit else 'below threshold')
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
                            roles_level='step1', transition=str(transition),
                            state_fragment=fragment if isinstance(fragment, str) else '',
                            result='poised in both' if truthy(r.get('step1_unresolved')) else
                            ('concordant transition' if truthy(r.get('step1_explains')) else 'unchanged'))
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
    """Sections, paragraphs, numbered figures and tables (captions below) and an appendix ('Data and methods') that is
    a collapsible section in HTML and the last pages of the PDF."""

    def __init__(self, title, keep_figures=False):
        self.title = title
        self.parts = []
        self.appendix_parts = []
        self.in_appendix = False
        self.toc = []
        self.figure_count = 0
        self.table_count = 0
        self.keep_figures = keep_figures
        self.figures = []      # (Figure objects, caption) for the matplotlib PDF
        self.text_blocks = []  # plain text in order, for the matplotlib PDF

    @property
    def target(self):
        return self.appendix_parts if self.in_appendix else self.parts

    def begin_appendix(self):
        self.in_appendix = True
        self.toc.append((2, 'appendix', 'Data and methods'))
        self.text_blocks.append(('heading', 'Data and methods', 2))

    def section(self, anchor, title, level=2, cls=None):
        if level == 2 and not self.in_appendix:
            self.toc.append((level, anchor, title))
        self.target.append('<h{l} id="{a}"{c}>{t}</h{l}>'.format(l=level, a=esc(anchor), t=esc(title),
                                                                  c=' class="{}"'.format(cls) if cls else ''))
        self.text_blocks.append(('heading', title, level))

    def heading(self, title, level=3, anchor=None, cls=None):
        self.target.append('<h{l}{a}{c}>{t}</h{l}>'.format(l=level, t=title, a=' id="{}"'.format(esc(anchor))
                                                           if anchor else '', c=' class="{}"'.format(cls) if cls
                                                           else ''))
        self.text_blocks.append(('heading', plain_text(title), level))

    def html(self, markup, text=None):
        self.target.append(markup)
        plain = text if text is not None else plain_text(markup)
        if plain:
            self.text_blocks.append(('text', plain, 0))

    def p(self, markup, cls=None):
        self.html('<p{}>{}</p>'.format(' class="{}"'.format(cls) if cls else '', markup))

    def figure_markup(self, figures, caption, cols=1, cls=''):
        """Markup of a numbered figure (panels as inline SVG, caption below); the caller places it."""
        self.figure_count += 1
        number = self.figure_count
        panels = ''.join('<div class="panel">{}</div>'.format(svg_markup(fig, 'f{}p{}'.format(number, k)))
                         for k, fig in enumerate(figures))
        markup = ('<figure id="figure-{n}" class="{cls}"><div class="panels cols-{c}">{p}</div><figcaption><span '
                  'class="label">Figure {n}.</span> {cap}</figcaption></figure>'.format(
                      n=number, c=cols, p=panels, cap=caption, cls=cls))
        plain = 'Figure {}. {}'.format(number, plain_text(caption))
        if self.keep_figures:
            self.figures.append((figures, plain))
            self.text_blocks.append(('figure', len(self.figures) - 1, 0))
        else:
            for fig in figures:
                fig.clear()
        return markup

    def figure(self, figures, caption, cols=1, cls=''):
        self.target.append(self.figure_markup(figures, caption, cols, cls))
        return self.figure_count

    def table(self, frame, caption, formats=None, numeric=(), long=(), escape=True, max_rows=None, nowrap=(),
              markup=(), row_classes=None, cls=''):
        """A table with its caption below; `formats` maps columns to formatting functions, `markup` names columns
        whose formatted text is HTML already."""
        self.table_count += 1
        number = self.table_count
        formats = formats or {}
        frame = frame if max_rows is None else frame.head(max_rows)
        if self.in_appendix:
            # columns without any value are left out
            keep = [c for c in frame.columns
                    if frame[c].map(lambda v: (formats[c](v) if c in formats else ('' if isnan(v) else str(v)))
                                    not in ('', 'nan')).any()]
            frame = frame[keep]
        nowrap = set(nowrap)
        head = ''.join('<th class="{}">{}</th>'.format('num' if c in numeric else '', esc(str(c).replace('_', ' ')))
                       for c in frame.columns)
        body, plain_rows = [], [' | '.join(str(c) for c in frame.columns)]
        for k, (_, row) in enumerate(frame.iterrows()):
            cells, plain = [], []
            for c in frame.columns:
                value = row[c]
                text = formats[c](value) if c in formats else ('' if isnan(value) else str(value))
                plain.append(plain_text(text))
                text = text if (not escape or c in markup) else esc(text)
                classes = ' '.join(x for x, on in (('num', c in numeric), ('long', c in long), ('key', c in nowrap))
                                   if on)
                cells.append('<td{}>{}</td>'.format(' class="{}"'.format(classes) if classes else '', text))
            row_cls = row_classes[k] if row_classes is not None else ''
            body.append('<tr{}>{}</tr>'.format(' class="{}"'.format(row_cls) if row_cls else '', ''.join(cells)))
            plain_rows.append(' | '.join(plain))
        self.target.append('<div class="table-wrap {cls}" id="table-{n}"><div class="table-scroll"><table><thead><tr>'
                           '{h}</tr></thead><tbody>{b}</tbody></table></div><p class="caption"><span class="label">'
                           'Table {n}.</span> {cap}</p></div>'.format(n=number, h=head, b=''.join(body), cap=caption,
                                                                     cls=cls))
        self.text_blocks.append(('table', '\n'.join(plain_rows), 0))
        self.text_blocks.append(('text', 'Table {}. {}'.format(number, plain_text(caption)), 0))
        return number

    def render(self, heading, nav=True, print_version=False):
        with open(os.path.join(TEMPLATE_DIR, 'report.css')) as handle:
            css = handle.read()
        with open(os.path.join(TEMPLATE_DIR, 'report.html')) as handle:
            template = string.Template(handle.read())
        links = ''.join('<a href="#{}">{}</a>'.format(esc(anchor), esc(title)) for level, anchor, title in self.toc
                        if level == 2)
        body = '\n'.join(self.parts)
        if self.appendix_parts:
            body += ('\n<section class="appendix" id="appendix"><details{}><summary><span class="appendix-title">'
                     'Data and methods</span><span class="appendix-hint">files, parameters, versions, coordinates and '
                     'the full numeric tables</span></summary>\n{}\n</details></section>'.format(
                         ' open' if print_version else '', '\n'.join(self.appendix_parts)))
        return template.substitute(version=esc(__version__), title=esc(self.title), css=css, heading=heading,
                                   nav='<nav class="contents" aria-label="Contents">{}</nav>'.format(links)
                                   if nav and len(self.toc) > 2 else '', body=body)


def plain_text(markup):
    """Text of a piece of markup (superscripts as ^)."""
    text = re.sub(r'<sup>(.*?)</sup>', r'^\1', str(markup))
    text = re.sub(r'<[^>]+>', ' ', text)
    text = html.unescape(re.sub(r'\s+', ' ', text)).strip()
    return re.sub(r'\s+([,.;:)])', r'\1', text).replace('( ', '(')


def body_outside_appendix(text):
    """The HTML body without the appendix (used by the tests and the checks of the report text)."""
    body = text.split('<body>', 1)[-1]
    return re.sub(r'<section class="appendix".*?</section>', '', body, flags=re.S)


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
    fig = new_figure(7.0, 3.3)
    grid = fig.add_gridspec(3, 1, height_ratios=[1, 1, 2.3], hspace=0.1)
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
    axes[2].legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.5), fontsize=6, frameon=False,
                   ncol=4, handlelength=2.2, columnspacing=1.2)
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
    if cascade is not None and 'structural_variant' in set(cascade['level']):
        levels.insert(1, 'structural_variant')
    labels = {'gene': 'Promoter', 'structural_variant': 'Structural variant', 'step1': 'Direct partners',
              'unexplained': 'Not explained'}
    first = cascade.set_index('level')['genes'] if cascade is not None else pd.Series(dtype=float)
    total = int(first.reindex(levels + ['unexplained']).fillna(0).sum()) if len(first) else len(run.classes)
    any_level = run.any_level_counts()
    rows = levels + ['unexplained']
    names = [labels.get(r, '{} steps away'.format(level_hops(r))) for r in rows]
    fig = new_figure(7.0, 0.42 * len(rows) + 1.0)
    axes = [fig.add_subplot(1, 2, k + 1) for k in range(2)]
    for ax, values, title in ((axes[0], first.reindex(rows).fillna(0), 'Nearest level that explains the gene'),
                              (axes[1], any_level.reindex(rows).fillna(0), 'Every level that explains the gene')):
        if title.startswith('Every'):
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

def ladder_figure(columns, levels, support=None):
    """Explanation ladder: one row per level from the promoter outwards, one column per promoter fragment. Filled
    marker: the level explains the expression change; open: it does not; ring: significant (q < 0.05) and beyond the
    distance-matched mirror loci. `support` maps levels to the 3D support written at the right."""
    n_rows, n_cols = len(levels), len(columns)
    has_support = bool(support)
    width = 2.15 + 0.62 * n_cols + (0.62 if has_support else 0)
    fig = new_figure(width, 0.255 * n_rows + 0.95)
    ax = fig.add_axes([2.0 / width, 0.62 / (0.255 * n_rows + 0.95), 1 - 2.05 / width,
                       0.255 * n_rows / (0.255 * n_rows + 0.95)])
    y_of = {lv: n_rows - 1 - i for i, lv in enumerate(levels)}
    for j, (_, rows) in enumerate(columns):
        present = [y_of[lv] for lv in levels if lv in rows]
        if present:
            ax.plot([j, j], [min(present), max(present)], color='#d9d8d3', lw=1.0, zorder=0)
        for lv in levels:
            row = rows.get(lv)
            if row is None:
                continue
            y = y_of[lv]
            if row['explains']:
                ax.scatter(j, y, s=42, color=INK, zorder=3, lw=0)
            else:
                ax.scatter(j, y, s=42, facecolor='white', edgecolor='#8c8b85', lw=1.0, zorder=3)
            if supported(row):
                ax.scatter(j, y, s=150, facecolor='none', edgecolor=ACCENT, lw=1.6, zorder=2)
    if has_support:
        for lv in levels:
            if lv in support:
                ax.text(n_cols - 0.35, y_of[lv], support[lv], va='center', ha='left', fontsize=7, color=INK)
        ax.text(n_cols - 0.35, n_rows - 0.35, '3D support', va='bottom', ha='left', fontsize=6.5, color='#5b5b57')
    ax.set_xlim(-0.5, n_cols - 0.5 + (0.95 if has_support else 0))
    ax.set_ylim(-0.6, n_rows - 0.4)
    ax.set_yticks([y_of[lv] for lv in levels])
    ax.set_yticklabels([level_label(lv) for lv in levels], fontsize=7, color=INK)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([textwrap.fill(label, 12) for label, _ in columns], fontsize=6.5, color='#5b5b57')
    ax.xaxis.tick_top()
    ax.tick_params(length=0, pad=3)
    for side in ('left', 'right', 'top', 'bottom'):
        ax.spines[side].set_visible(False)
    handles = [plt.Line2D([0], [0], marker='o', ls='', color=INK, ms=6, label='explains'),
               plt.Line2D([0], [0], marker='o', ls='', mfc='white', mec='#8c8b85', ms=6, label='does not explain'),
               plt.Line2D([0], [0], marker='o', ls='', mfc='none', mec=ACCENT, mew=1.6, ms=10,
                          label='significant, passes the distance control')]
    fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.01, 0.0), ncol=2, fontsize=6.3, frameon=False,
               handletextpad=0.3, columnspacing=1.0)
    return fig


# --------------------------------------------------------------------------------------------------------------------
# places named by genes

def read_strands(path, genes):
    """Strand of the named genes from a GTF (gene or transcript lines, gene_name, else gene_id) or a BED file with
    the strand in column 6."""
    wanted = {g for g in genes if g}
    strands = {}
    if not path or not wanted or not os.path.exists(path):
        return strands
    is_bed = re.search(r'\.bed(\.gz)?$', path) is not None
    name_pattern = re.compile(r'gene_name "([^"]+)"')
    id_pattern = re.compile(r'gene_id "([^"]+)"')
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith(('#', 'track', 'browser')):
                continue
            if is_bed:
                fields = line.rstrip('\n').split('\t')
                if len(fields) >= 6 and fields[3] in wanted:
                    strands.setdefault(fields[3], fields[5])
                continue
            if '\ttranscript\t' not in line and '\tgene\t' not in line:
                continue
            fields = line.split('\t', 8)
            match = name_pattern.search(fields[8]) or id_pattern.search(fields[8])
            if match and match.group(1) in wanted and match.group(1) not in strands:
                strands[match.group(1)] = fields[6]
                if len(strands) == len(wanted):
                    break
    return strands


def fmt_distance(bases):
    bases = abs(float(bases))
    if bases < 1000:
        return '{:.0f} bp'.format(bases)
    if bases < 10000:
        return '{:.1f} kb'.format(bases / 1000)
    if bases < 1e6:
        return '{:.0f} kb'.format(bases / 1000)
    return '{:.1f} Mb'.format(bases / 1e6)


PROMOTER_ADJECTIVE = {'active promoter': 'active', 'bivalent promoter': 'bivalent', 'poised promoter': 'poised',
                      'primed promoter': 'primed', 'repressed promoter': 'repressed',
                      'heterochromatin-like promoter': 'heterochromatin-like',
                      'architectural promoter': 'architectural', 'unmarked promoter': 'unmarked'}


ELEMENT_CLASS_ORDER = ('active promoter', 'bivalent promoter', 'poised promoter', 'primed promoter',
                       'repressed promoter', 'heterochromatin-like promoter', 'architectural promoter',
                       'promoter with marks of no role', 'unmarked promoter', 'active enhancer', 'primed enhancer',
                       'poised enhancer', 'enhancer-marked element', 'repressed element',
                       'heterochromatin-like element', 'architectural element', 'element with marks of no role',
                       'unmarked element')


def indefinite_article(phrase):
    """'an' before a vowel sound, else 'a'."""
    return 'an' if re.match(r'[aeiou]', plain_text(phrase).strip().lower()) else 'a'


def class_noun(cls):
    """Plain words of an element class for a place name: 'active enhancer', 'element' for marks without a role."""
    cls = str(cls or '')
    return 'element' if not cls or cls.endswith('marks of no role') else cls


class Places:
    """Names of graph nodes by genes: a promoter by its genes, an other end by its distance from the report gene and
    the nearest gene (baited promoters of the graph). With `elements` (node -> (class in A, class in B, TSS status,
    TSS genes), Run.element_lookup) the element class is named in plain words: 'the bivalent promoter of Hoxa2',
    'an active enhancer 120 kb upstream of Hoxa3 (near Skap2)', 'the active promoter of Skap2 (an other end)'; a
    class that differs between the conditions is given per condition in parentheses."""

    def __init__(self, graph, strands=None, elements=None, names=('A', 'B')):
        self.graph = graph
        self.strands = strands or {}
        self.elements = elements or {}
        self.names = list(names)
        index = {}
        for n, d in graph.nodes(data=True):
            if d.get('bait'):
                names = gene_names(d.get('genes', ''))
                if names:
                    index.setdefault(d.get('chrom'), []).append(((int(d['start']) + int(d['end'])) / 2, names[0]))
        self.index = {c: sorted(v) for c, v in index.items()}
        self.middles = {c: np.array([x[0] for x in v]) for c, v in self.index.items()}

    def middle(self, node):
        data = self.graph.nodes[node]
        return (int(data['start']) + int(data['end'])) / 2

    def nearest_gene(self, node, exclude=()):
        if node not in self.graph:
            return None
        chrom = self.graph.nodes[node].get('chrom')
        if chrom not in self.index:
            return None
        middle = self.middle(node)
        k = int(np.searchsorted(self.middles[chrom], middle))
        best = None
        for i in range(max(0, k - 12), min(len(self.index[chrom]), k + 12)):
            position, name = self.index[chrom][i]
            if name in exclude:
                continue
            if best is None or abs(position - middle) < best[1]:
                best = (name, abs(position - middle))
        return best

    def promoter_genes(self, node):
        data = self.graph.nodes[node] if node in self.graph else {}
        return gene_names(data.get('genes', '')) if data.get('bait') else []

    def element_class(self, node):
        """(class in A, class in B) of a node, or None without the element table."""
        info = self.elements.get(node)
        return (info[0], info[1]) if info else None

    def per_condition(self, classes, words):
        """'bivalent in Md, active in PA2' from the classes of the two conditions."""
        return ', '.join('{} in {}'.format(words(c), esc(n)) for c, n in zip(classes, self.names))

    def label(self, node, ref_node=None, ref_gene=None, markup=True, article=False, classed=True):
        """'Hoxa2 promoter' or 'element 1.0 Mb upstream of Hoxa3 (near Skap2)'; with `article` 'the' or 'an'. With
        element classes: 'bivalent promoter of Hoxa2', 'active enhancer 1.0 Mb upstream of Hoxa3 (near Skap2)';
        `classed` False leaves the class out (where the text states the state already)."""
        g = gene_span if markup else str
        a_word = article
        if node not in self.graph:
            return 'an unmapped fragment' if a_word else 'unmapped fragment'
        classes = self.element_class(node) if classed else None
        names = self.promoter_genes(node)
        if not names and classes and self.elements[node][2] == 'non-baited promoter':
            # an other end holding a TSS: named as a promoter of its genes
            names = gene_names(self.elements[node][3])
            other_end = True
        else:
            other_end = False
        if names:
            if classes is None:
                return '{}{} promoter'.format('the ' if a_word else '', g(names[0]))
            suffix = ' (an other end)' if other_end else ''
            words = [PROMOTER_ADJECTIVE.get(c) for c in classes]
            if classes[0] == classes[1] and words[0]:
                return '{}{} promoter of {}{}'.format('the ' if a_word else '', words[0], g(names[0]), suffix)
            if classes[0] == classes[1]:
                return '{}{} promoter{}'.format('the ' if a_word else '', g(names[0]), suffix)
            detail = self.per_condition(classes, lambda c: PROMOTER_ADJECTIVE.get(c) or 'marks of no role')
            return '{}{} promoter ({}{})'.format('the ' if a_word else '', g(names[0]), 'an other end; '
                                                 if other_end else '', detail)
        data = self.graph.nodes[node]
        exclude = {ref_gene} if ref_gene else set()
        near = self.nearest_gene(node, exclude)
        noun = 'element'
        extra = ''
        if classes is not None:
            if classes[0] == classes[1]:
                noun = class_noun(classes[0])
            else:
                extra = self.per_condition(classes, class_noun)
        if ref_node in self.graph and self.graph.nodes[ref_node].get('chrom') == data.get('chrom') and ref_gene:
            distance = self.middle(node) - self.middle(ref_node)
            strand = self.strands.get(ref_gene)
            if strand in ('+', '-'):
                side = 'upstream of' if (distance < 0) == (strand == '+') else 'downstream of'
            else:
                side = 'from'
            text = '{} {} {} {}'.format(noun, fmt_distance(distance), side, g(ref_gene))
            inside = []
            if near is not None and near[1] < abs(distance):
                inside.append('near {}'.format(g(near[0])))
            if extra:
                inside.append(extra)
            if inside:
                text += ' ({})'.format('; '.join(inside))
        else:
            text = '{} on {}{}'.format(noun, data.get('chrom'), ' near {}'.format(g(near[0])) if near else '')
            if extra:
                text += ' ({})'.format(extra)
        return ('{} '.format(indefinite_article(text)) if a_word else '') + text


# --------------------------------------------------------------------------------------------------------------------
# reading the levels of one gene

def significant(row):
    q, p = row.get('q', np.nan), row.get('p', np.nan)
    if not isnan(q):
        return float(q) < 0.05
    return not isnan(p) and float(p) < 0.05


def supported(row):
    """Explains, significant and beyond the distance-matched mirror loci."""
    return bool(row.get('explains')) and significant(row) and truthy(row.get('beyond_mirror'))


def is_step(level):
    return re.fullmatch(r'step\d+', str(level)) is not None


def parse_transition(text):
    """('bivalent', 'active', 'repressive mark lost') from 'bivalent -> active: repressive mark lost'."""
    text = '' if isnan(text) else str(text)
    main, _, note = text.partition(': ')
    a, arrow, b = main.partition(' -> ')
    return a.strip(), (b if arrow else a).strip(), note.strip()


def readable_note(note):
    """'gain of the repressive mark' from 'repressive mark gained'."""
    kinds = {'active': 'an activating', 'repressive': 'the repressive', 'enhancer': 'an enhancer'}
    out = []
    for part in [p.strip() for p in note.split(',') if p.strip()]:
        match = re.fullmatch(r'(active|repressive|enhancer) mark (gained|lost)', part)
        if match:
            out.append('{} of {} mark'.format('gain' if match.group(2) == 'gained' else 'loss',
                                              kinds[match.group(1)]))
        else:
            out.append(part)
    return join_words(out)


def transition_verdict(result, explains):
    """How a state transition relates to the expression change, from the result of the level."""
    result = str(result or '')
    if explains or result == 'concordant transition':
        return 'in line with the expression change'
    if 'not supported' in result:
        return 'in the direction of the expression change, but not supported by the mark changes themselves'
    if 'without a poised state' in result:
        return 'a change that involves no poised state and is not counted as an explanation'
    if 'not concordant' in result or 'discordant' in result:
        return 'against the expression change'
    return 'not counted as an explanation ({})'.format(result) if result else 'not counted as an explanation'


def state_text(transition, names):
    a, b, _ = parse_transition(transition)
    if not a:
        return ''
    if a == b:
        return '{} in {} and {}'.format(a, names[0], names[1])
    return '{} in {}, {} in {}'.format(a, names[0], b, names[1])


def strongest_change(row):
    """(mark, log2 change) of the strongest node from the all-level text ('... H3K27ac log2 +5.2')."""
    match = re.search(r'(\S+) log2 ([+-]?[0-9.]+)\s*$', str(row.get('strongest_text') or ''))
    return (match.group(1), float(match.group(2))) if match else (None, np.nan)


def fold_values(detail):
    """'H3K27ac 0.6 / 0.6, H3K4me2 4.9 / 4.9' from '... (H3K27ac 0.6/0.6, H3K4me2 4.9/4.9; mean) ...'."""
    match = re.search(r'\(([^()]*?); (?:mean|max1kb)\)', str(detail or ''))
    if not match:
        return ''
    return ', '.join(item.strip().replace('/', ' / ') for item in match.group(1).split(','))


def transition_count(detail):
    match = re.search(r'\((\d+) of (\d+) with this transition\)', str(detail or ''))
    return (int(match.group(1)), int(match.group(2))) if match else None


def sign_map(text):
    out = {}
    for item in str(text or '').split(';'):
        if ':' in item:
            mark, value = item.rsplit(':', 1)
            try:
                out[mark] = int(float(value))
            except ValueError:
                continue
    return out


def marks_phrase(values, marks, first_prefix='log2 '):
    """'H3K27ac (log2 +3.4) and H3K4me2 (+1.6) rise and H3K27me3 (−1.4) falls'."""
    marks = [m for m in marks if m in values and not isnan(values[m])]
    up = [m for m in marks if values[m] > 0]
    down = [m for m in marks if values[m] < 0]
    state = {'first': True}

    def item(m):
        text = '{} ({}{})'.format(m, first_prefix if state['first'] else '', fmt_log2(values[m]))
        state['first'] = False
        return text
    parts = []
    if up:
        parts.append('{} {}'.format(join_words([item(m) for m in up]), 'rises' if len(up) == 1 else 'rise'))
    if down:
        parts.append('{} {}'.format(join_words([item(m) for m in down]), 'falls' if len(down) == 1 else 'fall'))
    return join_words(parts)


def number_word(n):
    words = {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six', 7: 'seven', 8: 'eight', 9: 'nine'}
    return words.get(n, '{:,}'.format(n))


def capitalise(text):
    return text[:1].upper() + text[1:] if text else text


class GeneStory:
    """Levels of one gene (all its promoter fragments) and their reading in sentences, tables and the ladder."""

    def __init__(self, run, places, gene, nodes):
        self.run = run
        self.places = places
        self.gene = gene
        self.nodes = list(nodes)
        self.fragments = []  # (node, {level: row})
        for node in self.nodes:
            part = run.levels[run.levels['node'] == node]
            if len(part):
                self.fragments.append((node, {r['level']: r for _, r in part.sort_values('order').iterrows()}))
        self.fragments.sort(key=lambda f: self.fragment_key(f[1]))
        rna = [rows[next(iter(rows))]['rna_log2FC'] for _, rows in self.fragments]
        if not rna:
            rna = [run.graph.nodes[n].get('rna_log2FC', np.nan) for n in self.nodes if n in run.graph]
        self.rna = next((float(v) for v in rna if not isnan(v)), np.nan)
        self.differential = bool(self.fragments)

    def fragment_key(self, rows):
        best = self.best(rows)
        return (best is None, float(best['q']) if best is not None and not isnan(best['q']) else 2.0)

    @staticmethod
    def best(rows, require_beyond=True):
        candidates = [r for lv, r in rows.items() if is_step(lv) and r['explains'] and significant(r) and
                      (truthy(r.get('beyond_mirror')) or not require_beyond)]
        if not candidates:
            return None
        return min(candidates, key=lambda r: (float(r['q']) if not isnan(r['q']) else 2.0,
                                              float(r['p']) if not isnan(r['p']) else 2.0, r['order']))

    @staticmethod
    def nearest(rows):
        explaining = [r for r in rows.values() if r['explains']]
        return min(explaining, key=lambda r: r['order']) if explaining else None

    # names ----------------------------------------------------------------------------------------------------------
    def g(self):
        return gene_span(self.gene)

    def fragment_others(self, node):
        data = self.run.graph.nodes[node] if node in self.run.graph else {}
        genes = data.get('genes', '') or self.run.classes.set_index('node')['genes'].get(node, '')
        return [x for x in gene_names(genes) if x != self.gene]

    def fragment_name(self, node, markup=True):
        """'the fragment shared with Hoxa7 and Mira' (several fragments) or 'its promoter'."""
        others = self.fragment_others(node)
        g = gene_span if markup else str
        if others:
            return 'the fragment shared with {}'.format(join_words([g(x) for x in others[:3]]) +
                                                        (' and others' if len(others) > 3 else ''))
        return 'the fragment with {} alone'.format(g(self.gene))

    def fragment_heading(self, node):
        others = self.fragment_others(node)
        if len(self.fragments) <= 1 and len(self.nodes) <= 1:
            return 'Promoter of {}'.format(self.g())
        if others:
            return 'Promoter fragment shared with {}'.format(join_words([gene_span(x) for x in others[:3]]))
        return 'Promoter fragment of {} alone'.format(self.g())

    def column_label(self, node):
        if len(self.fragments) <= 1:
            return ''
        others = self.fragment_others(node)
        return 'with {}'.format(', '.join(others[:2])) if others else '{} alone'.format(self.gene)

    def place(self, node, ref_node, article=True, markup=True, classed=True):
        return self.places.label(node, ref_node, self.gene, markup=markup, article=article, classed=classed)

    # sentences ------------------------------------------------------------------------------------------------------
    def expression_sentence(self):
        a, b = self.run.names
        if isnan(self.rna):
            return 'No expression change is recorded for {}.'.format(self.g())
        if not self.differential:
            return ('{} is not differential between {} and {} at the thresholds of this comparison (log2 {}); the '
                    'levels are not evaluated.'.format(self.g(), esc(a), esc(b), fmt_log2(self.rna)))
        return 'Expression of {} is about {}-fold {} in {} than in {} (log2 {}).'.format(
            self.g(), fold(self.rna), 'higher' if self.rna > 0 else 'lower', esc(b), esc(a), fmt_log2(self.rna))

    def promoter_sentence(self, node, rows, nearest=None):
        threshold = self.run.chromatin_threshold
        nearest_level = nearest['level'] if nearest is not None else None
        nearest_tag = ', the nearest level that explains it'
        where = 'its promoter' if len(self.fragments) == 1 else 'its promoter fragment {}'.format(
            self.fragment_name(node)[len('the fragment '):])
        where = where.replace('promoter fragment with', 'promoter fragment that carries')
        clauses = []
        marks = rows.get('promoter_marks')
        if marks is not None:
            values = marks['mark_values'] or {}
            if marks['explains'] and marks['concordant_marks']:
                clauses.append('{}, in line with the expression change{}'.format(
                    marks_phrase(values, marks['concordant_marks']),
                    nearest_tag if nearest_level == 'promoter_marks' else ''))
            else:
                large = [m for m, v in values.items() if not isnan(v) and abs(v) >= threshold]
                if large:
                    clauses.append('{}, against the expression change'.format(marks_phrase(values, large)))
                elif values:
                    mark = max(values, key=lambda m: abs(values[m]) if not isnan(values[m]) else -1)
                    clauses.append('the histone marks change by at most log2 {} ({})'.format(
                        fmt_log2(abs(values[mark]), signed=False), mark))
                else:
                    clauses.append('no histone mark changes by log2 {:g} or more'.format(threshold))
        contacts = rows.get('contacts')
        if contacts is not None:
            gained, lost, shared, strength = contacts['contacts']
            counts = '{} gained, {} lost{}; summed strength log2 {}'.format(
                fmt_int(gained), fmt_int(lost), ', {} shared'.format(fmt_int(shared)) if not isnan(shared) else '',
                fmt_log2(strength)) if not isnan(gained) else 'summed strength log2 {}'.format(fmt_log2(strength))
            if contacts['explains']:
                verb = 'strengthen' if not isnan(strength) and strength > 0 else 'weaken'
                clauses.append('its own contacts {} ({}){}'.format(verb, counts, nearest_tag if nearest_level ==
                                                                  'contacts' else ''))
            elif not isnan(strength) and abs(strength) >= threshold:
                clauses.append('its own contacts change against the expression change ({})'.format(counts))
            else:
                clauses.append('its own contacts change little ({})'.format(counts))
        if not clauses:
            return ''
        return 'At {}, {}.'.format(where, '; '.join(clauses))

    def state_sentence(self, rows, nearest):
        state = rows.get('promoter_state')
        if state is None:
            return ''
        a_state, b_state, note = parse_transition(state['transition'])
        if not a_state:
            return ''
        a, b = self.run.names
        kind = '{}s'.format(self.run.kind) if self.run.kind != 'condition' else 'conditions'
        subject = 'The promoter' if len(self.fragments) == 1 else 'This fragment'
        nearest_here = nearest is not None and nearest['level'] == 'promoter_state'
        if state['result'] == 'poised in both' or (state['explains'] and a_state == b_state):
            trigger = 'the trigger of its {} is not among the measured layers'.format(
                'activation' if self.rna > 0 else 'silencing')
            if nearest_here:
                trigger = 'this poised state is the nearest level that accounts for the change, and ' + trigger
            if a_state == 'bivalent':
                return '{} carries both an activating and a repressive mark in both {} (bivalent); {}.'.format(
                    subject, kind, trigger)
            return '{} is {} in both {}; {}.'.format(subject, a_state, kind, trigger)
        if a_state != b_state:
            return '{} changes from {} in {} to {} in {}{}, {}{}.'.format(
                subject, a_state, esc(a), b_state, esc(b), ' ({})'.format(readable_note(note)) if note else '',
                transition_verdict(state['result'], state['explains']),
                '; this is the nearest level that explains it' if nearest_here else '')
        return '{} is {} in both {}.'.format(subject, a_state, kind)

    def structural_sentences(self, rows, nearest):
        """A structural variant near the promoter and the copy number of the promoter region, in plain words:
        'The promoter lies 0.4 Mb from the t(8;21) breakpoint; the rearrangement is present in t821 only and the
        promoter's contacts to chromosome 21 near RUNX1 are lost there; this is the nearest level that explains the
        change.'"""
        out = []
        subject = 'The promoter' if len(self.fragments) == 1 else 'This promoter fragment'
        a, b = self.run.names
        sv = rows.get('structural_variant')
        if sv is not None and str(sv.get('result') or '') not in ('', 'no breakpoint nearby'):
            extra = sv.get('extra') or {}
            name = str(extra.get('sv_name', '') or 'a')
            distance = extra.get('sv_distance', np.nan)
            status = str(extra.get('sv_status', '') or '')
            partner = str(extra.get('sv_partner', '') or '')
            genes = [g for g in str(extra.get('sv_partner_genes', '') or '').split(',') if g and g != 'nan'][:2]
            chrom = partner.split(':')[0]
            where = 'chromosome {}'.format(chrom[3:] if chrom.lower().startswith('chr') else chrom) \
                if chrom and not name.startswith(chrom + ' ') else 'the far side of the join'
            near = ' near {}'.format(join_words([gene_span(g) for g in genes])) if genes else ''
            text = '{} lies {} from the {} breakpoint'.format(
                subject, fmt_distance(distance) if not isnan(distance) else 'close to', esc(name))
            if status in ('present in both', 'similar support in both'):
                text += '; the rearrangement is present in both {} and {}'.format(esc(a), esc(b))
            elif status == 'gained':
                text += '; the rearrangement is present in {} only'.format(esc(b))
            elif status == 'lost':
                text += '; the rearrangement is present in {} only'.format(esc(a))
            elif status:
                text += '; the rearrangement is {}'.format(esc(status))
            junction = extra.get('sv_junction_log2', np.nan)
            if not isnan(junction) and abs(junction) >= 0.5:
                text += ', and its contacts to {}{} are {} in {} (log2 {})'.format(
                    where, near, 'gained' if junction > 0 else 'lost', esc(b), fmt_log2(junction))
            elif not isnan(junction):
                text += ', and its contacts to {}{} are similar in both'.format(where, near)
            else:
                text += ' (partner side: {}{})'.format(where, near)
            if sv['explains']:
                text += '; this accounts for the expression change{}'.format(
                    ', the nearest level that explains it'
                    if nearest is not None and nearest['level'] == 'structural_variant' else '')
            else:
                text += {'present in both': '; it does not explain the change',
                         'discordant': '; the change of the rearrangement runs against the expression change',
                         'no active partner element': '; the partner side carries no active element',
                         'junction contacts discordant': '; the promoter\'s own contacts across the junction change '
                                                         'against the expression change'}.get(
                    str(sv.get('result') or ''), '')
            out.append(text + '.')
        cn = rows.get('copy_number')
        if cn is not None and str(cn.get('result') or '') not in ('', 'no value'):
            extra = cn.get('extra') or {}
            va, vb = extra.get('copy_number_A', np.nan), extra.get('copy_number_B', np.nan)
            value = extra.get('copy_number_log2', np.nan)
            if not isnan(va) and not isnan(vb):
                verdict = 'in line with the expression change' if cn['explains'] else \
                    'against the expression change' if cn.get('result') == 'discordant' else \
                    ('practically unchanged' if not isnan(value) and abs(value) < 0.1 else
                     'a change below the threshold')
                out.append('The relative copy number of the promoter region is {} in {} and {} in {} (log2 {}), '
                           '{}.'.format(sig(va, 2), esc(a), sig(vb, 2), esc(b), fmt_log2(value), verdict))
        return out

    def step_sentence(self, node, rows, nearest):
        """The nearest explaining level (when not at the promoter) and the strongest supported explanation."""
        best = self.best(rows)
        out = []
        nearest_said = nearest is not None and (nearest['level'] in ('promoter_marks', 'contacts', 'promoter_state',
                                                                     'structural_variant', 'copy_number') or
                                                cascade_level(nearest['level']) == getattr(self, 'known_level', None))
        if nearest is None:
            return ['None of the measured levels explains the change.']
        if best is not None:
            lead = 'The strongest supported explanation'
            if best['level'] == nearest['level']:
                lead = 'The nearest level that explains the change, and the strongest supported one,'
            elif not nearest_said:
                out.append('The nearest level that explains the change is {}.'.format(level_phrase(nearest['level'])))
            where = 'at the direct partners' if best['level'] == 'step1' else \
                '{} steps away'.format(level_hops(best['level']))
            values = best['mark_values'] or {}
            marks = best['concordant_marks'] or [m for m in values if abs(values[m]) >= self.run.chromatin_threshold]
            across = 'across the {} nodes there'.format(fmt_int(best['nodes'])) if not isnan(best['nodes']) else \
                'across these nodes'
            sentence = '{} lies {}: {} {}'.format(lead, where, marks_phrase(values, marks), across)
            mark, value = strongest_change(best)
            if best.get('strongest_node') and mark and not isnan(value) and \
                    abs(value) >= self.run.chromatin_threshold:
                sentence += ', most strongly at {}, where {} {} by log2 {} (about {}-fold)'.format(
                    self.place(best['strongest_node'], node), mark, 'rises' if value > 0 else 'falls',
                    fmt_log2(value), fold(value))
            sentence += ', beyond the change at distance-matched mirror loci (P = {}, q = {}).'.format(
                fmt_pq(best['p']), fmt_pq(best['q']))
            out.append(sentence)
            if not isnan(best.get('support')):
                ci = str(best.get('support_ci') or '')
                out.append('In the nucleus, nodes {} apart in the graph are closer than distance-matched pairs '
                           '(3D support {}{}; 0.5 means no relation).'.format(
                               'one step' if best['level'] == 'step1' else '{} steps'.format(level_hops(best['level'])),
                               sig(best['support'], 2), ', interval {}'.format(ci.replace('-', ' to '))
                               if ci and ci != 'NA' else ''))
            learned = self.learned_note(best)
            if learned:
                out.append('Under the mark roles learned from these data, this step {}.'.format(learned))
            return out
        if not nearest_said:
            out.append('The nearest level that explains the change is {}.'.format(level_phrase(nearest['level'])))
        weak = self.best(rows, require_beyond=False)
        if weak is not None:
            out.append('{} explain the change significantly (P = {}, q = {}), but distance-matched mirror loci change '
                       'the same way, so the graph adds no support beyond genomic distance.'.format(
                           capitalise(level_phrase(weak['level'])), fmt_pq(weak['p']), fmt_pq(weak['q'])))
        else:
            steps = [r for lv, r in rows.items() if is_step(lv) and r['explains']]
            tested = [r for r in steps if not isnan(r['p']) or not isnan(r['q'])]
            if tested:
                top = min(tested, key=lambda r: (float(r['q']) if not isnan(r['q']) else 2.0, r['order']))
                out.append('The strongest of the explaining steps away from the promoter, {}, does not reach '
                           'significance (P = {}, q = {}{}).'.format(
                               level_phrase(top['level']), fmt_pq(top['p']), fmt_pq(top['q']),
                               '; the change also exceeds that at distance-matched mirror loci'
                               if truthy(top.get('beyond_mirror')) else ''))
            elif steps:
                out.append('{} explain the change without a per-gene test.'.format(
                    capitalise(join_words([level_phrase(r['level']) for r in steps]))))
            else:
                out.append('No step away from the promoter explains the change.')
        return out

    def learned_note(self, row, plural=False):
        """How the result of a level differs under the learned roles, as a verb phrase ('' when it does not)."""
        learned = row.get('learned_explains')
        if isnan(learned) or learned is None:
            return ''
        if bool(learned) != bool(row['explains']):
            if learned:
                return 'also explain{} the gene'.format('' if plural else 's')
            return '{} not explain the gene because {}'.format('do' if plural else 'does', self.learned_reason(row))
        if row['explains'] and significant(row) and not isnan(row.get('learned_q')) and \
                float(row['learned_q']) >= 0.05:
            return '{} not significant (q = {})'.format('are' if plural else 'is', fmt_pq(row['learned_q']))
        return ''

    def learned_reason(self, row):
        result = str(row.get('learned_result') or '')
        if row['level'] == 'contacts':
            return 'the own contacts show no consistent association with expression in these data'
        if result == 'not callable':
            return 'the chromatin states cannot be called with these roles'
        declared, learned = sign_map(row.get('signs_declared')), sign_map(row.get('signs_learned'))
        reasons = []
        for mark in declared:
            if mark in learned and learned[mark] != declared[mark]:
                if learned[mark] == 0:
                    reasons.append('{} shows no consistent association with expression at this level'.format(mark))
                else:
                    reasons.append('{} behaves as {} mark at this level'.format(
                        mark, 'an activating' if learned[mark] > 0 else 'a repressive'))
        if reasons:
            return join_words(reasons)
        return 'the marks change {}'.format('against the expression change' if result == 'discordant' else
                                           'below the threshold' if result == 'below threshold' else result)

    def fragments_sentence(self):
        if len(self.fragments) < 2:
            return ''
        items = []
        for node, rows in self.fragments:
            nearest = self.nearest(rows)
            best = self.best(rows)
            first = 'is first explained by {}'.format(level_phrase(nearest['level'])) if nearest is not None else \
                'is not explained'
            control = 'passes the distance control {}'.format(
                'at the direct partners' if best['level'] == 'step1' else '{} steps away'.format(
                    level_hops(best['level']))) if best is not None else 'passes no distance control'
            items.append((self.fragment_name(node), first, control))
        if len({(f, c) for _, f, c in items}) == 1:
            return 'The {} promoter fragments of {} give the same result.'.format(number_word(len(items)), self.g())
        parts = ['{} {} and {}'.format(name, first, control) for name, first, control in items]
        return 'The {} promoter fragments of {} differ: {}.'.format(number_word(len(items)), self.g(),
                                                                   '; '.join(parts))

    def glance(self, known_level=None):
        """The 'at a glance' paragraph; `known_level` (cascade level of a card group) drops the sentence that would
        repeat the group heading."""
        self.known_level = known_level
        sentences = [self.expression_sentence()]
        if not self.fragments:
            return ' '.join(sentences)
        node, rows = self.fragments[0]
        nearest = self.nearest(rows)
        sentences.append(self.promoter_sentence(node, rows, nearest))
        sentences.append(self.state_sentence(rows, nearest))
        sentences += self.structural_sentences(rows, nearest)
        sentences += self.step_sentence(node, rows, nearest)
        sentences.append(self.fragments_sentence())
        return ' '.join(s for s in sentences if s)

    # ladder ---------------------------------------------------------------------------------------------------------
    def ladder(self):
        if not self.fragments:
            return None
        levels = sorted({lv for _, rows in self.fragments for lv in rows}, key=level_order)
        columns = [(self.column_label(node), rows) for node, rows in self.fragments]
        support = {}
        node, rows = self.fragments[0]
        for lv, row in rows.items():
            if is_step(lv) and not isnan(row.get('support')):
                support[lv] = sig(row['support'], 2)
        return ladder_figure(columns, levels, support)

    # evidence -------------------------------------------------------------------------------------------------------
    def what_changed(self, row, node):
        level = row['level']
        threshold = self.run.chromatin_threshold
        a, b = self.run.names
        if level in SV_LEVELS:
            return esc(row.get('detail_raw') or row.get('all_changes') or '')
        if level == 'contacts':
            gained, lost, shared, strength = row['contacts']
            if isnan(gained):
                return esc(row.get('all_changes') or '')
            return '{} gained, {} lost, {} shared; summed strength log2 {}'.format(
                fmt_int(gained), fmt_int(lost), fmt_int(shared), fmt_log2(strength))
        if level in ('promoter_state', 'step1_state'):
            transition = row['transition']
            if not transition or transition in ('no poised fragment', 'states not callable'):
                return 'no poised partner' if 'poised' in str(transition or row['result']) else 'not callable'
            text = state_text(transition, (esc(a), esc(b)))
            _, _, note = parse_transition(transition)
            if note:
                text += ' ({})'.format(readable_note(note))
            if level == 'step1_state' and row.get('state_fragment'):
                text += ' at {}'.format(self.place(row['state_fragment'], node, classed=False))
                count = transition_count(row.get('detail_raw'))
                if count:
                    text += '; {} of {} partners'.format(*count)
            if row['result'] == 'poised in both':
                text += '; poised, trigger not measured'
            elif not row['explains'] and parse_transition(transition)[0] != parse_transition(transition)[1]:
                result = str(row['result'] or '')
                text += '; {}'.format('not supported by the mark changes' if 'not supported' in result else
                                      'no poised state' if 'without a poised state' in result else
                                      'against the expression change' if 'not concordant' in result else result)
            return text
        values = row['mark_values'] or {}
        if row['explains'] and row['concordant_marks']:
            text = ', '.join('{} {}'.format(m, fmt_log2(values[m])) for m in row['concordant_marks'] if m in values)
            if level != 'promoter_marks':
                n = row.get('n_concordant')
                if row.get('strongest_node') and not isnan(n) and int(n) >= 1:
                    others = int(n) - 1
                    text += ' at {}{}'.format(self.place(row['strongest_node'], node),
                                              ' and {} other node{}'.format(others, '' if others == 1 else 's')
                                              if others else '')
                elif not isnan(row.get('nodes')):
                    text += ' over {} nodes'.format(fmt_int(row['nodes']))
            return text
        listing = ', '.join('{} {}'.format(m, fmt_log2(v)) for m, v in values.items() if not isnan(v))
        large = [m for m, v in values.items() if not isnan(v) and abs(v) >= threshold]
        if large:
            return 'against the expression change: {}'.format(', '.join('{} {}'.format(m, fmt_log2(values[m]))
                                                                       for m in large))
        return 'all below log2 {:g}{}'.format(threshold, ' ({})'.format(listing) if listing else '')

    def evidence_frame(self, node, rows, with_support):
        out = []
        classes = []
        for lv in sorted(rows, key=level_order):
            row = rows[lv]
            control = 'n.a.'
            if is_step(lv) and not isnan(row.get('beyond_mirror')):
                control = 'passed' if truthy(row['beyond_mirror']) else 'not passed'
            entry = {'Level': level_label(lv), 'Explains': 'yes' if row['explains'] else 'no',
                     'What changed': self.what_changed(row, node), 'P': fmt_pq(row['p']), 'q': fmt_pq(row['q']),
                     'Distance control': control}
            if with_support:
                support = ''
                if level_hops(lv) > 0 and not isnan(row.get('support')):
                    ci = str(row.get('support_ci') or '')
                    support = '{}{}'.format(sig(row['support'], 2), ' ({})'.format(ci.replace('-', ' to '))
                                            if ci and ci != 'NA' else '')
                entry['3D support'] = support
            out.append(entry)
            classes.append(('explains' if row['explains'] else 'muted') + (' supported' if supported(row) else ''))
        return pd.DataFrame(out), classes

    def learned_items(self, rows):
        """Sentences on the levels whose result differs under the learned roles, levels with the same difference
        grouped ('the nodes 2 and 3 steps away do not explain the gene because ...')."""
        groups = {}
        for lv in sorted(rows, key=level_order):
            key = self.learned_note(rows[lv], plural=True)
            if key:
                groups.setdefault(re.sub(r'\(q = [^)]*\)', '', key), []).append(lv)
        items = []
        for key, levels in groups.items():
            steps = [level_hops(lv) for lv in levels if is_step(lv)]
            subjects = [level_phrase(lv) for lv in levels if not is_step(lv)]
            if steps:
                subjects.append('the nodes {} {} away'.format(join_words([str(k) for k in steps]),
                                                             'step' if steps == [1] else 'steps'))
            plural = len(levels) > 1 or not levels[0].endswith('state')
            note = self.learned_note(rows[levels[0]], plural=plural)
            if len(levels) > 1:
                note = re.sub(r' \(q = [^)]*\)', ' (q = {})'.format(', '.join(fmt_pq(rows[lv]['learned_q'])
                                                                         for lv in levels)), note)
                note = note.replace('at this level', 'at these levels')
            items.append('{} {}'.format(join_words(subjects), note))
        return items


# --------------------------------------------------------------------------------------------------------------------
# report builder

class Report:
    def __init__(self, run, args, keep_figures, places):
        self.run = run
        self.args = args
        self.keep_figures = keep_figures
        self.places = places

    def story(self, gene, nodes=None):
        if nodes is None:
            nodes, _ = self.run.gene_nodes(gene)
        return GeneStory(self.run, self.places, gene, nodes)

    def title_block(self, title, subtitle, kicker):
        return ('<p class="kicker">{}</p><h1>{}</h1><p class="subtitle">{}</p><p class="date">{}</p>'.format(
            kicker, title, subtitle, datetime.date.today().isoformat()))

    def glance_block(self, story, doc, caption_extra='', arcs=None, known_level=None):
        """The 'at a glance' paragraph with the ladder next to it; in a card of the main report (`arcs`: a function
        returning the markup of the arcs figure) the paragraph above, the ladder and the arcs below side by side."""
        fig = story.ladder()
        text = '<p class="glance">{}</p>'.format(story.glance(known_level))
        if fig is None:
            doc.html('<div class="glance-grid single">{}</div>'.format(text))
            if arcs is not None:
                doc.html(arcs(), text='')
            return
        cols = len(story.fragments)
        first = getattr(doc, 'first_ladder', None)
        if not first:
            doc.first_ladder = doc.figure_count + 1
        caption = 'Explanation ladder of {}, read as Figure {}.'.format(story.g(), first) if first else ('Explanation ladder of {g}{f}: each row is a level, from the promoter (top) outwards. Filled '
                   'marker: the level explains the expression change; open: it does not; blue ring: significant '
                   '(q &lt; 0.05) and beyond distance-matched mirror loci.{s}{x}'.format(
                       g=story.g(), f=', one column per promoter fragment' if cols > 1 else '',
                       s=' 3D support at the right: how much closer graph neighbours at this step are in the nucleus '
                         'than distance-matched pairs (0.5: no relation).' if any(
                           is_step(lv) and not isnan(r.get('support')) for _, rows in story.fragments
                           for lv, r in rows.items()) else '', x=caption_extra))
        figure = doc.figure_markup([fig], caption, cls='ladder')
        if arcs is not None:
            doc.html('{}<div class="card-figs">{}{}</div>'.format(text, figure, arcs()), text=plain_text(text))
            return
        doc.html('<div class="glance-grid">{}{}</div>'.format(text, figure), text=plain_text(text))

    def arcs_markup(self, doc, gene, nodes, story):
        hood = self.gene_neighbourhood(gene, nodes, story)
        if hood is None or hood['sub'].number_of_edges() == 0:
            return '<p class="note">No interactions of {} to draw.</p>'.format(gene_span(gene))
        return doc.figure_markup([arcs_figure(self.run, hood['sub'], hood['seeds'], hood['explaining'])],
                                 self.arcs_caption(gene, hood['sub'], doc), cls='arcs')

    # gene figures -----------------------------------------------------------------------------------------------
    def gene_neighbourhood(self, gene, nodes, story):
        run = self.run
        seeds = [n for n in nodes if n in run.graph]
        if not seeds:
            return None
        rows = run.levels[run.levels['node'].isin(nodes)]
        direction = np.sign(story.rna) if not isnan(story.rna) else 1.0
        levels = rows.loc[rows['explains'], 'level'].tolist() if len(rows) else []
        found, which = explaining_nodes(run, seeds, levels, direction)
        explaining = [n for k in sorted(found) for n in found[k]]
        hops = max([level_hops(lv) for lv in levels] + [1])
        sub = neighbourhood_subgraph(run, seeds, hops, explaining, self.args.maxNodes)
        explaining = [n for n in explaining if n in sub]
        return dict(seeds=seeds, sub=sub, explaining=explaining, found=found, which=which, hops=hops, levels=levels,
                    rows=rows)

    def arcs_caption(self, gene, sub, doc=None):
        a, b = (esc(n) for n in self.run.names)
        first = getattr(doc, 'first_arcs', None) if doc is not None else None
        if first:
            return 'Interactions around {g}, drawn as in Figure {f}; {n} nodes and {e} interactions.'.format(
                g=gene_span(gene), f=first, n=sub.number_of_nodes(), e=sub.number_of_edges())
        if doc is not None:
            doc.first_arcs = doc.figure_count + 1
        return ('Interactions around {g} in {a} (top), in {b} (middle) and their change (bottom: gained in {b} above '
                'the axis in red, lost below in blue, shared in grey). Black squares: the promoter of {g}; open purple '
                'circles: nodes that carry the explanation. Line width scales with interaction strength; {n} nodes and '
                '{e} interactions{r}.'.format(g=gene_span(gene), a=a, b=b, n=sub.number_of_nodes(),
                                              e=sub.number_of_edges(), r=', reduced to the strongest partners and the '
                                              'paths to the explaining nodes' if sub.number_of_nodes() >=
                                              self.args.maxNodes else ''))

    def gene_figures(self, doc, gene, nodes, story, arcs_only=False):
        hood = self.gene_neighbourhood(gene, nodes, story)
        if hood is None:
            doc.p('The promoter of {} is not in the graph.'.format(gene_span(gene)), cls='note')
            return
        sub, seeds, explaining = hood['sub'], hood['seeds'], hood['explaining']
        if sub.number_of_edges() == 0:
            doc.p('{} has no interactions in either {}.'.format(gene_span(gene), self.run.kind), cls='note')
            return
        doc.figure([arcs_figure(self.run, sub, seeds, explaining)], self.arcs_caption(gene, sub))
        if arcs_only:
            return
        mark = None
        for lv in sorted(hood['levels'], key=level_order, reverse=True):
            row = hood['rows'][hood['rows']['level'] == lv].iloc[0]
            if row.get('concordant_marks'):
                mark = row['concordant_marks'][0]
                break
        mark = mark or (self.run.marks[0] if self.run.marks else None)
        named = []
        for k in sorted(hood['found']):
            items = ['{} ({} {})'.format(self.places.label(n, seeds[0], gene), hood['which'][n][1],
                                         fmt_log2(hood['which'][n][2])) for n in hood['found'][k] if n in sub][:3]
            if items:
                named.append('{}: {}'.format('direct partners' if k == 1 else '{} steps away'.format(k),
                                             join_words(items)))
        a, b = (esc(n) for n in self.run.names)
        hops = hood['hops']
        doc.figure(neighbourhood_figures(self.run, sub, seeds, explaining, mark, self.args.maxLog2),
                   'Neighbourhood of {g} up to {h} step{s}, drawn on one layout: interactions in {a} and in {b}, their '
                   'change, and the expression change of promoters (squares) with the change of {m} at other ends '
                   '(circles). Purple outline: nodes that carry the explanation{w}.'.format(
                       g=gene_span(gene), h=hops, s='' if hops == 1 else 's', a=a, b=b, m=esc(mark or 'the mark'),
                       w=' ({})'.format('; '.join(named)) if named else ''), cols=2, cls='compact')

    # hub and community ------------------------------------------------------------------------------------------
    def gene_list(self, nodes, values=True, limit=8):
        graph = self.run.graph
        seen, items = set(), []
        for n in nodes:
            names = gene_names(graph.nodes[n].get('genes', ''))
            if not names or names[0] in seen:
                continue
            seen.add(names[0])
            value = graph.nodes[n].get('rna_log2FC')
            if values and value is not None and not isnan(value):
                items.append('{} {}'.format(gene_span(names[0]), fmt_log2(value)))
            else:
                items.append(gene_span(names[0]))
        more = len(items) - limit
        items = items[:limit]
        if more > 0:
            items.append('{} more'.format(more))
        return join_words(items), len(seen)

    def hub_members(self, seeds):
        graph = self.run.graph
        members = {}
        for which, excluded in (('A', 'gained'), ('B', 'lost')):
            partners = {w for s in seeds for w in graph.neighbors(s)
                        if graph.edges[s, w].get('status', 'shared') != excluded}
            promoters = {p for w in partners for p in graph.neighbors(w)
                         if graph.nodes[p].get('bait') and graph.edges[w, p].get('status', 'shared') != excluded}
            members[which] = (partners, (promoters | {w for w in partners if graph.nodes[w].get('bait')}) -
                              set(seeds))
        return members

    def hub_block(self, doc, gene, nodes, draw=True, prose=True):
        """Hub of a gene in prose: its partners, the other promoters contacting them, which change with it, which
        left or joined; then the side-by-side figure."""
        run = self.run
        graph = run.graph
        seeds = [n for n in nodes if n in graph]
        if not seeds:
            return
        members = self.hub_members(seeds)
        (partners_a, promoters_a), (partners_b, promoters_b) = members['A'], members['B']
        all_promoters = promoters_a | promoters_b
        differential = set(run.classes['node'])
        value = graph.nodes[seeds[0]].get('rna_log2FC')
        direction = np.sign(float(value)) if value is not None and not isnan(value) else 1.0

        def rna(n):
            v = graph.nodes[n].get('rna_log2FC')
            return float(v) if v is not None and not isnan(v) else np.nan
        same = sorted((n for n in all_promoters if n in differential and not isnan(rna(n)) and
                       np.sign(rna(n)) == direction), key=lambda n: -abs(rna(n)))
        opposite = sorted((n for n in all_promoters if n in differential and not isnan(rna(n)) and
                           np.sign(rna(n)) == -direction), key=lambda n: -abs(rna(n)))
        left = sorted(promoters_a - promoters_b, key=lambda n: -abs(rna(n)) if not isnan(rna(n)) else 0)
        joined = sorted(promoters_b - promoters_a, key=lambda n: -abs(rna(n)) if not isnan(rna(n)) else 0)
        a, b = (esc(n) for n in run.names)
        g = gene_span(gene)
        text = ['In {a}, the promoter of {g} contacts {pa:,} fragments, which in turn contact {xa:,} other promoters; '
                'in {b}, {pb:,} fragments and {xb:,} other promoters.'.format(
                    a=a, b=b, g=g, pa=len(partners_a), xa=len(promoters_a), pb=len(partners_b), xb=len(promoters_b))]
        same_text, n_same = self.gene_list(same)
        opposite_text, n_opposite = self.gene_list(opposite)
        if n_same or n_opposite:
            parts = []
            if n_same:
                parts.append('{} {} in the same direction as {} ({})'.format(
                    number_word(n_same), 'changes' if n_same == 1 else 'change', g, same_text))
            if n_opposite:
                parts.append('{} in the opposite direction ({})'.format(number_word(n_opposite), opposite_text))
            text.append('Among the differential promoters of this hub, {}.'.format('; '.join(parts)))
        else:
            text.append('No other promoter of the hub is differential.')
        left_text, n_left = self.gene_list(left, values=False)
        joined_text, n_joined = self.gene_list(joined, values=False)
        moves = []
        if n_left:
            moves.append('{} promoter{} {} the hub in {} ({})'.format(number_word(n_left), '' if n_left == 1 else 's',
                                                                   'leaves' if n_left == 1 else 'leave', b, left_text))
        if n_joined:
            moves.append('{} {} {}'.format(number_word(n_joined), 'joins' if n_joined == 1 else 'join',
                                           'it ({})'.format(joined_text) if n_left else
                                           'the hub in {} ({})'.format(b, joined_text)))
        if moves:
            text.append(capitalise(join_words(moves)) + '.')
        if prose:
            doc.p(' '.join(text))
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
                            'outline), its partners and the other promoters contacting them, coloured by the '
                            'expression change. Orange dotted outline: left the hub in {b}; green outline: joined; '
                            'faded: present in the other {k} only.{r}'.format(
                                g=g, a=a, b=b, k=run.kind, r=' Reduced to the named promoters and their partners.'
                                if sub.number_of_nodes() >= self.args.maxNodes else ''), cols=2, cls='compact')

    def community_name(self, members):
        """A community named by its differential gene with the largest expression change."""
        table = self.run.gene_communities
        rows = table[table['node'].isin(members)]
        if 'differential' in rows:
            diff = rows[rows['differential'].map(truthy)]
            rows = diff if len(diff) else rows
        rows = rows.assign(abs_fc=rows['rna_log2FC'].abs()).sort_values('abs_fc', ascending=False)
        for genes in rows['genes']:
            names = gene_names(genes)
            if names:
                return names[0]
        return ''

    def coherence(self, community_a, community_b):
        table = self.run.community_expression
        if table is None:
            return None
        hit = table[(table['community_a'] == community_a) & (table['community_b'] == community_b)]
        return hit.iloc[0] if len(hit) else None

    def transition_phrase(self, transition):
        b = esc(self.run.names[1])
        return {'stable': 'stays stable in {}'.format(b), 'split': 'splits in {}'.format(b),
                'merged': 'merges with other communities in {}'.format(b),
                'remodelled': 'is remodelled in {}'.format(b), 'dissolved': 'dissolves in {}'.format(b),
                'formed': 'forms in {}'.format(b)}.get(str(transition), '')

    def gene_community_text(self, doc, gene, nodes):
        run = self.run
        if run.gene_communities is None:
            return None
        rows = run.gene_communities[run.gene_communities['node'].isin(nodes)]
        if not len(rows):
            doc.p('{} is not assigned to a community.'.format(gene_span(gene)), cls='note')
            return None
        r = rows.iloc[0]
        a, b = (esc(n) for n in run.names)
        g = gene_span(gene)
        status = {'kept': 'stays in its community between {} and {}'.format(a, b),
                  'switch': 'moves to another community between {} and {}'.format(a, b),
                  'lost': 'belongs to a community in {} but to none in {}'.format(a, b),
                  'gained': 'belongs to no community in {} and joins one in {}'.format(a, b)}.get(
                      str(r.get('status')), 'has community status {}'.format(esc(r.get('status'))))
        sizes = []
        if not isnan(r.get('size_a')) and int(r['size_a']) > 0:
            sizes.append('{:,} promoters in {}'.format(int(r['size_a']), a))
        if not isnan(r.get('size_b')) and int(r['size_b']) > 0:
            sizes.append('{:,} in {}'.format(int(r['size_b']), b))
        text = '{} {}{}'.format(g, status, ' ({}'.format(', '.join(sizes)) if sizes else '')
        if not isnan(r.get('member_jaccard')) and sizes:
            text += '; {:.0f}% of the members shared)'.format(100 * float(r['member_jaccard']))
        elif sizes:
            text += ')'
        text += '.'
        transition = r.get('community_a_transition')
        phrase = self.transition_phrase(transition) if not isnan(transition) else ''
        if phrase:
            text += ' The {} community {}.'.format(a, phrase)
        coherent = self.coherence(r.get('community_a'), r.get('community_b'))
        if coherent is not None and not isnan(coherent.get('coherence_q')):
            q = float(coherent['coherence_q'])
            text += ' The expression of its {:,} measured member genes changes {} (mean log2 {}, q = {}).'.format(
                int(coherent['n_genes']), 'coherently' if q < 0.1 else 'no more coherently than in random gene sets',
                fmt_log2(coherent['mean_rna_log2FC']), fmt_pq(q))

        def listing(value, limit=12):
            names = [x for x in str(value).split(',') if x and x != 'nan'] if not isnan(value) else []
            names = [x for x in names if not x.startswith(IGNORED_GENE_PREFIXES)] or names
            if not names:
                return ''
            shown = [gene_span(x) for x in names[:limit]]
            if len(names) > limit:
                shown.append('{} more'.format(len(names) - limit))
            return join_words(shown)
        joined, left = listing(r.get('joined_genes')), listing(r.get('left_genes'))
        if joined:
            text += ' Joined in {}: {}.'.format(b, joined)
        if left:
            text += ' Left in {}: {}.'.format(b, left)
        doc.p(text)
        return r

    def community_members(self, community_a, community_b):
        table = self.run.gene_communities
        members_a = set(table.loc[table['community_a'] == community_a, 'node']) if community_a >= 0 else set()
        members_b = set(table.loc[table['community_b'] == community_b, 'node']) if community_b >= 0 else set()
        rows = table[table['node'].isin(members_a | members_b)].copy()
        rows['in'] = np.where(rows['node'].isin(members_a) & rows['node'].isin(members_b), 'both',
                              np.where(rows['node'].isin(members_a), 'left', 'joined'))
        rows['gene'] = rows['genes'].map(first_gene)
        return members_a, members_b, rows

    def community_block(self, doc, community_a, community_b, focus=(), info=None, draw=True):
        """A community in prose (members, coherent change, same and opposite direction, left and joined) and the
        side-by-side figure; the member table goes to the appendix."""
        run = self.run
        if run.gene_communities is None:
            doc.p('No community table in the run outputs.', cls='note')
            return
        members_a, members_b, rows = self.community_members(community_a, community_b)
        a, b = (esc(n) for n in run.names)
        mean = rows['rna_log2FC'].mean()
        direction = np.sign(mean) if not isnan(mean) else 1.0
        diff = rows['differential'].map(truthy) if 'differential' in rows else \
            rows['rna_log2FC'].abs() >= run.rna_threshold
        same = rows[diff & (np.sign(rows['rna_log2FC']) == direction)].sort_values('rna_log2FC', key=np.abs,
                                                                                   ascending=False)
        opposite = rows[diff & (np.sign(rows['rna_log2FC']) == -direction)].sort_values('rna_log2FC', key=np.abs,
                                                                                        ascending=False)

        def listing(frame, values=True, limit=10):
            names, seen = [], set()
            for r in frame.itertuples():
                label = gene_names(r.genes)
                if not label or label[0] in seen:
                    continue
                seen.add(label[0])
                names.append('{} {}'.format(gene_span(label[0]), fmt_log2(r.rna_log2FC)) if values and
                             not isnan(r.rna_log2FC) else gene_span(label[0]))
            more = len(names) - limit
            return join_words(names[:limit] + (['{} more'.format(more)] if more > 0 else [])), len(names)
        transition = (info or {}).get('transition', '')
        phrase = self.transition_phrase(transition)
        text = ['The community has {:,} promoters in {} and {:,} in {} ({:,} in both){}.'.format(
            len(members_a), a, len(members_b), b, len(members_a & members_b),
            '; it {}'.format(phrase) if phrase else '')]
        coherent = self.coherence(community_a, community_b)
        if coherent is not None and not isnan(coherent.get('coherence_q')):
            q = float(coherent['coherence_q'])
            text.append('The expression of its {:,} measured member genes changes {} (mean log2 {}, P = {}, q = {}); '
                        '{:.0f}% of its interactions are rewired.'.format(
                            int(coherent['n_genes']), 'coherently' if q < 0.1 else
                            'no more coherently than random gene sets of the same size',
                            fmt_log2(coherent['mean_rna_log2FC']), fmt_pq(coherent['coherence_p']), fmt_pq(q),
                            100 * float(coherent['rewired_fraction'])))
        same_text, n_same = listing(same)
        opposite_text, n_opposite = listing(opposite)
        if n_same or n_opposite:
            text.append('Differential members change {} ({}){}.'.format(
                'up' if direction > 0 else 'down', same_text or 'none',
                '; {:,} change the other way ({})'.format(n_opposite, opposite_text) if n_opposite else
                '; none changes the other way'))
        left_text, n_left = listing(rows[rows['in'] == 'left'], values=False)
        joined_text, n_joined = listing(rows[rows['in'] == 'joined'], values=False)
        moves = []
        if n_left:
            moves.append('{} left in {}'.format(left_text, b))
        if n_joined:
            moves.append('{} joined'.format(joined_text))
        if moves:
            text.append(capitalise(join_words(moves, 'and')) + '.')
        doc.p(' '.join(text))
        if draw:
            self.community_figure(doc, members_a, members_b, same, opposite, rows, focus, community_a, community_b)

    def community_figure(self, doc, members_a, members_b, same, opposite, rows, focus, community_a, community_b):
        run = self.run
        graph = run.graph
        baits = [n for n in members_a | members_b if n in graph]
        if not baits:
            return
        counts = {}
        for n in baits:
            for w in graph.neighbors(n):
                if not graph.nodes[w].get('bait'):
                    counts[w] = counts.get(w, 0) + 1
        connectors = sorted((w for w, c in counts.items() if c >= 2), key=lambda w: (-counts[w], w))
        nodes = set(baits) | set(connectors[:max(0, self.args.maxNodes - len(baits))])
        if len(baits) > self.args.maxNodes:
            keep = set(focus) | set(same['node'].head(12)) | set(opposite['node'].head(12)) | \
                set(rows.loc[rows['in'] != 'both', 'node'].head(24))
            nodes = keep | {w for w in connectors if sum(graph.has_edge(w, p) for p in keep) >= 2}
        sub = ordered_subgraph(graph, nodes)
        sub = ordered_subgraph(sub, [n for n in sub.nodes if sub.degree[n] > 0])
        if sub.number_of_edges() == 0:
            doc.p('The promoters of the community share no interacting fragment.', cls='note')
            return
        positions = separate(spread_layout(sub, aspect=1.2), sub, min_distance=0.06)
        others = set(n for n in sub if not sub.nodes[n].get('bait'))
        name = self.community_name(members_a | members_b)
        a, b = run.names
        figures = hub_figures(run, sub, positions, (members_a | others) & set(sub.nodes),
                              (members_b | others) & set(sub.nodes), set(focus),
                              '{}: community of {}'.format(a, name), '{}: community of {}'.format(b, name),
                              self.args.maxLog2, members_label='community')
        first = getattr(doc, 'first_community', None)
        if first:
            doc.figure(figures, 'Community of {} in {} (left) and {} (right), drawn as in Figure {}.'.format(
                gene_span(name), esc(a), esc(b), first), cols=2, cls='compact')
            return
        doc.first_community = doc.figure_count + 1
        doc.figure(figures, 'Community of {g} in {a} (left) and {b} (right) on one layout: member promoters (squares, '
                            'coloured by the expression change) and the other ends shared by at least two of them. '
                            'Orange dotted outline: left the community in {b}; green outline: joined; faded: member '
                            'in the other {k} only.{r}'.format(g=gene_span(name), a=esc(a), b=esc(b), k=run.kind,
                                                               r=' Reduced to the named members and their shared '
                                                                 'other ends.' if len(baits) > self.args.maxNodes
                                                               else ''), cols=2, cls='compact')

    # gene report ------------------------------------------------------------------------------------------------
    def gene_report(self, gene):
        run = self.run
        nodes, differential = run.gene_nodes(gene)
        doc = Document('{}: {} vs {}'.format(gene, *run.names), self.keep_figures)
        heading = self.title_block(gene_span(gene), esc(run.condition_title()), 'Gene report')
        story = self.story(gene, nodes)
        doc.section('glance', 'At a glance')
        if not nodes:
            doc.p('{} is not a baited promoter of the graph of this comparison.'.format(gene_span(gene)))
            self.appendix(doc, genes=[])
            return doc, heading
        self.glance_block(story, doc)
        doc.section('interactions', 'Interactions')
        self.gene_figures(doc, gene, nodes, story)
        if story.fragments:
            doc.section('evidence', 'Evidence by level')
            with_support = any(not isnan(r.get('support')) and level_hops(lv) > 0 for _, rows in story.fragments
                               for lv, r in rows.items())
            first_table = None
            learned = {}
            for node, rows in story.fragments:
                if len(story.fragments) > 1:
                    doc.heading(story.fragment_heading(node), level=3)
                frame, classes = story.evidence_frame(node, rows, with_support)
                a, b = (esc(n) for n in run.names)
                if first_table is None:
                    caption = ('Evidence for {g}{f}, from the promoter outwards. Mark changes are log2 fold changes '
                               '{b} / {a}, at the promoter itself or taken over the nodes of a step; explains: the '
                               'change agrees with the expression change by at least log2 {t:g}. P and q: per-gene '
                               'permutation test of the step; distance control: the change exceeds that at '
                               'distance-matched mirror loci.{s}'.format(
                                   g=gene_span(gene), f=' at the {}'.format(story.fragment_name(node)[4:])
                                   if len(story.fragments) > 1 else '', a=a, b=b, t=run.chromatin_threshold,
                                   s=' 3D support: probability that graph neighbours at this step are closer in the '
                                     'nucleus than distance-matched pairs.' if with_support else ''))
                else:
                    caption = 'As Table {}, for the {}.'.format(first_table, story.fragment_name(node)[4:])
                number = doc.table(frame, caption, markup=tuple(frame.columns), escape=False,
                                   nowrap=('Level', 'P', 'q', 'Distance control', 'Explains'),
                                   long=('What changed',), row_classes=classes, cls='evidence')
                first_table = first_table or number
                for item in story.learned_items(rows):
                    learned.setdefault(item, []).append(node)
            if learned:
                notes = []
                for item, where in learned.items():
                    if len(story.fragments) > 1:
                        item = item[:-1] + '; {})' if item.endswith(')') else item + ' ({})'
                        item = item.format('both fragments' if len(where) == len(story.fragments) == 2 else
                                               'all fragments' if len(where) == len(story.fragments) else
                                               join_words([story.fragment_name(n)[4:] for n in where]))
                    notes.append(item)
                doc.p('Under the mark roles learned from these data, {}.'.format('; '.join(notes)), cls='note')
        doc.section('states', 'Chromatin states')
        self.state_block(doc, story)
        self.partner_elements(doc, gene, nodes)
        doc.section('hub', 'Hub and community')
        self.hub_block(doc, gene, nodes, draw=False)
        community = self.gene_community_text(doc, gene, nodes)
        self.hub_block(doc, gene, nodes, prose=False)
        self.appendix(doc, genes=[(gene, nodes, story)], community=community)
        return doc, heading

    def state_block(self, doc, story):
        run = self.run
        if not story.fragments:
            doc.p('No chromatin states are called for {}: the gene is not differential.'.format(story.g()),
                  cls='note')
            return
        a, b = (esc(n) for n in run.names)
        items = []
        described = set()
        for node, rows in story.fragments:
            sentences = []
            for lv, where in (('promoter_state', 'Promoter'), ('step1_state', 'Direct partners')):
                row = rows.get(lv)
                if row is None:
                    continue
                transition = row['transition']
                if not transition or transition in ('no poised fragment', 'states not callable'):
                    sentences.append('{}: {}.'.format(where, 'no poised partner' if 'poised' in str(
                        transition or row['result']) else 'states not callable'))
                    continue
                count = transition_count(row.get('detail_raw'))
                fragment = row.get('state_fragment') if lv == 'step1_state' else ''
                if fragment and (fragment, transition) in described:
                    sentences.append('{}: the same {} element{}.'.format(
                        where, parse_transition(transition)[0], ' ({} of {} partners)'.format(*count) if count else ''))
                    continue
                text = state_text(transition, (a, b))
                _, _, note = parse_transition(transition)
                if note:
                    text += ' ({})'.format(readable_note(note))
                if fragment:
                    described.add((fragment, transition))
                    text += ' at {}'.format(story.place(fragment, node, classed=False))
                    if count:
                        text += ' ({} of {} partners)'.format(*count)
                values = fold_values(row.get('detail_raw'))
                if values:
                    text += '; fold enrichment {} / {}: {}'.format(a, b, values)
                if lv == 'promoter_state':
                    text += self.rechip_text(node, transition)
                sentence = '{}: {}.'.format(where, text)
                alt = row.get('alt_transition')
                if alt and alt != transition:
                    sentence += ' On the strongest 1 kb window of the fragment the call is {}{}.'.format(
                        state_text(alt, (a, b)), ', which explains the change' if truthy(row.get('alt_explains'))
                        and not row['explains'] else '')
                sentences.append(sentence)
            if sentences:
                title = story.fragment_heading(node) if len(story.fragments) > 1 else ''
                items.append('<p>{}{}</p>'.format('<strong>{}.</strong> '.format(title) if title else '',
                                                   ' '.join(sentences)))
        doc.html('<div class="states">{}</div>'.format(''.join(items)))
        doc.p('States are called from the presence of each mark on the fragment (fold enrichment over input of at '
              'least {:g}, fragment mean); bivalent: an activating and a repressive mark together; poised enhancer: '
              'an enhancer mark with the repressive mark. A poised state present in both {}s explains a change '
              'whose trigger is not among the measured layers.{}'.format(
                  run.command_value('--stateThreshold', 2.0), run.kind,
                  ' Re-ChIP: sequential ChIP signal of both marks on the same fragment; a signal of at least 2 '
                  'confirms the bivalent call.' if run.rechip is not None else ''), cls='note')

    def class_counts_text(self, classes):
        """'five active enhancers, three unmarked elements and one active promoter' from a list of classes."""
        counts = pd.Series(list(classes), dtype=object).value_counts()
        parts = []
        for cls, n in counts.items():
            noun = class_noun(cls)
            if noun == 'element':
                noun = 'element with marks of no declared role' if n == 1 else \
                    'elements with marks of no declared role'
            elif n != 1:
                noun = re.sub(r'(enhancer|element|promoter)$', r'\1s', noun)
            parts.append('{} {}'.format(number_word(int(n)), noun))
        return join_words(parts)

    def partner_elements(self, doc, gene, nodes):
        """Element classes of the promoter fragments of the gene and of their direct partners in both conditions,
        and the partners whose class changes (declared roles, fragment mean)."""
        run = self.run
        places = self.places
        if not places.elements:
            return
        a, b = (esc(n) for n in run.names)
        seeds = [n for n in nodes if n in run.graph]
        partners = sorted({m for n in seeds for m in run.graph.neighbors(n)} - set(seeds))
        info = [(m, places.elements.get(m)) for m in partners]
        info = [(m, i) for m, i in info if i]
        sentences = []
        def adjective(c):
            return PROMOTER_ADJECTIVE.get(c, 'with marks of no declared role')

        def state(c):
            return '{} in both {}s'.format(adjective(c[0]), run.kind) if c[0] == c[1] else \
                '{} in {} and {} in {}'.format(adjective(c[0]), a, adjective(c[1]), b)
        own = [(n, places.element_class(n)) for n in seeds if places.element_class(n)]
        story = self.story(gene, nodes)
        if len(own) == 1:
            sentences.append('The promoter of {} is {}.'.format(gene_span(gene), state(own[0][1])))
        elif own:
            items = ['{} is {}'.format(story.fragment_name(n), state(c)) for n, c in own]
            sentences.append('Of the promoter fragments of {}, {}.'.format(gene_span(gene), join_words(items)))
        if info:
            promoters = [i for _, i in info if i[2] in ('baited promoter', 'non-baited promoter')]
            uncaptured = sum(1 for i in promoters if i[2] == 'non-baited promoter')
            sentences.append('Of its {} direct partner{}, {} {} {}{}.'.format(
                number_word(len(info)), '' if len(info) == 1 else 's',
                number_word(len(promoters)) if promoters else 'none', 'is' if len(promoters) == 1 else 'are',
                'the promoter of another gene' if len(promoters) == 1 else 'promoters of other genes',
                '' if not uncaptured else ' ({} not captured)'.format(
                    'it is' if len(promoters) == 1 else
                    'all of them' if uncaptured == len(promoters) else '{} of them'.format(number_word(uncaptured)))))
            texts = [self.class_counts_text([i[k] for _, i in info]) for k in (0, 1)]
            if texts[0] == texts[1]:
                sentences.append('In both {}s, the partners are {}.'.format(run.kind, texts[0]))
            else:
                for name, text in zip((a, b), texts):
                    sentences.append('In {}, the partners are {}.'.format(name, text))
            changed = [i for _, i in info if i[0] != i[1]]
            if changed:
                moves = pd.Series(['{} to {}'.format(class_noun(i[0]), class_noun(i[1])) for i in changed]) \
                    .value_counts()
                items = ['{} from {}'.format(number_word(int(n)), move) for move, n in moves.head(4).items()]
                more = len(moves) - 4
                sentences.append('{} partner{} change{} class between {} and {}: {}{}.'.format(
                    capitalise(number_word(len(changed))), '' if len(changed) == 1 else 's',
                    's' if len(changed) == 1 else '', a, b, join_words(items),
                    ' and {} other change{}'.format(number_word(more), '' if more == 1 else 's') if more > 0 else ''))
            else:
                sentences.append('No partner changes class between {} and {}.'.format(a, b))
        if not sentences:
            return
        doc.heading('Element classes of the promoter and its partners', level=3)
        doc.p(' '.join(sentences))
        doc.p('Element classes follow from the marks present on a fragment (fold enrichment of at least {:g}, fragment '
              'mean) and the declared roles of the marks: active promoter or enhancer (an activating mark), primed '
              'enhancer (the enhancer mark without an activating one), bivalent promoter or poised enhancer (with the '
              'repressive mark), repressed (the repressive mark only), unmarked (none). A fragment is a promoter when '
              'it holds an annotated transcription start site, whether it is captured or not.'.format(
                  run.command_value('--stateThreshold', 2.0)), cls='note')

    def rechip_text(self, node, transition):
        run = self.run
        if run.rechip is None or node not in run.rechip.index:
            return ''
        row = run.rechip.loc[node]
        values = [row.get('reChIP_{}_mean'.format(n), np.nan) for n in run.names]
        if any(isnan(v) for v in values):
            return ''
        a_state, b_state, _ = parse_transition(transition)
        confirmed = [n for n, v, s in zip(run.names, values, (a_state, b_state)) if s == 'bivalent' and v >= 2]
        bivalent = [n for s, n in zip((a_state, b_state), run.names) if s == 'bivalent']
        tag = ''
        if bivalent:
            if len(confirmed) == len(bivalent):
                tag = ' (bivalent confirmed)'
            elif confirmed:
                tag = ' (bivalent confirmed in {})'.format(esc(confirmed[0]))
            else:
                tag = ' (bivalent not confirmed)'
        def value(v):
            return '{:.0f}'.format(v) if v >= 10 else sig(v, 2)
        return '; re-ChIP signal {} in {} and {} in {}{}'.format(value(values[0]), esc(run.names[0]), value(values[1]),
                                                                 esc(run.names[1]), tag)

    # main report ------------------------------------------------------------------------------------------------
    def main_report(self):
        run = self.run
        doc = Document('{} vs {}: differential graph report'.format(*run.names), self.keep_figures)
        heading = self.title_block(esc(run.condition_title()), 'Differential chromatin interaction graph',
                                   'Comparison report')
        self.key_findings(doc)
        self.cascade_section(doc)
        self.elements_section(doc)
        cards = self.gene_cards(doc)
        self.communities_section(doc)
        self.appendix(doc, genes=cards, main=True)
        return doc, heading

    def level_table(self):
        if run_has(self.run.level_counts):
            return self.run.variant_rows(self.run.level_counts).set_index('level')
        return None

    def key_findings(self, doc):
        run = self.run
        a, b = (esc(n) for n in run.names)
        doc.section('findings', 'Key findings')
        bullets = []
        n = len(run.classes)
        up = int((run.classes['rna_log2FC'] > 0).sum())
        bullets.append('{:,} genes differ in expression between {} and {} by at least {:g}-fold (|log2| ≥ {:g}): '
                       '{:,} are higher and {:,} lower in {}.'.format(n, a, b, 2 ** run.rna_threshold,
                                                                     run.rna_threshold, up, n - up, b))
        levels = self.level_table()
        if levels is not None and 'genes_first' in levels:
            first = levels['genes_first']

            def count(lv):
                return int(first.get(lv, 0)) if not isnan(first.get(lv, np.nan)) else 0
            promoter = sum(count(lv) for lv in ('promoter_marks', 'promoter_contacts', 'promoter_state'))
            direct = count('step1_marks') + count('step1_state')
            poised = levels['first_poised_in_both'].get('promoter_state', np.nan) if 'first_poised_in_both' in \
                levels else np.nan
            steps = [lv for lv in levels.index if re.fullmatch(r'step\d+', str(lv)) and lv != 'step1']
            far = join_words(['{} {:,}'.format('the nodes {} steps away for'.format(level_hops(lv)) if k == 0 else
                                               '{} steps away for'.format(level_hops(lv)), count(lv))
                              for k, lv in enumerate(steps)])
            bullets.append('Taking the nearest level that explains each gene, the promoter accounts for {p:,} genes '
                           '({pp:.0f}%): {m:,} through its histone marks, {c:,} through its own contacts and {s:,} '
                           'through its chromatin state{po}. The direct partners account for {d:,} ({dp:.0f}%){f}; '
                           '{u:,} genes remain unexplained.'.format(
                               p=promoter, pp=100 * promoter / n if n else 0, m=count('promoter_marks'),
                               c=count('promoter_contacts'), s=count('promoter_state'),
                               po=' ({:,} of them bivalent or poised in both {}s)'.format(int(poised), run.kind)
                               if not isnan(poised) else '', d=direct, dp=100 * direct / n if n else 0,
                               f=', ' + far if far else '', u=count('none')))
            if 'structural_variant' in levels.index or 'copy_number' in levels.index:
                sv_any = levels['genes_any'].get('structural_variant', np.nan)
                cn_any = levels['genes_any'].get('copy_number', np.nan)
                items = []
                if not isnan(sv_any):
                    items.append('a structural variant near the promoter (a breakpoint within the window, its change '
                                 'between the {}s in line with the expression change and active elements on the '
                                 'partner side) explains {:,} genes, {:,} of them first'.format(
                                     run.kind, int(sv_any), count('structural_variant')))
                if not isnan(cn_any):
                    items.append('a copy number change of the promoter region explains {:,}, {:,} of them first'.format(
                        int(cn_any), count('copy_number')))
                bullets.append(capitalise(join_words(items)) + '.')
            parts = []
            for lv in ['step1_marks'] + steps:
                if lv not in levels.index:
                    continue
                r = levels.loc[lv]
                parts.append('{} {:,} ({:,} passing the distance control, {} at q &lt; 0.05)'.format(
                    'the histone marks of the direct partners explain' if lv == 'step1_marks' else
                    'nodes {} steps away'.format(level_hops(lv)), int(r['genes_any']),
                    int(r.get('any_beyond_mirror', 0) or 0),
                    '{:,}'.format(int(r['any_q_lt_0.05'])) if not isnan(r.get('any_q_lt_0.05')) and
                    int(r['any_q_lt_0.05']) else 'none'))
            if parts:
                bullets.append('Counting every level that explains a gene, {}.'.format(join_words(parts)))
        strongest = self.strongest_cases()
        if strongest:
            bullets.append('The strongest cases supported by the graph are {}.'.format(join_words(
                ['{} ({}, q = {})'.format(gene_span(g), 'direct partners' if lv == 'step1' else
                                          '{} steps away'.format(level_hops(lv)), fmt_pq(q))
                 for g, lv, q in strongest])))
        bivalent = self.bivalent_sentence()
        if bivalent:
            bullets.append(bivalent)
        communities = self.community_sentence()
        if communities:
            bullets.append(communities)
        perm = run.permutation or {}
        if perm.get('observed_indirect_delta_r2') is not None:
            bullets.append('Across {:,} genes, adding the nodes 2 steps away to a regression of the expression change '
                           'on the mark changes raises the explained variance by {} percentage points, against {} in '
                           'permuted graphs (P = {}).'.format(
                               int(perm.get('n_genes', 0)), sig(100 * perm['observed_indirect_delta_r2'], 2),
                               sig(100 * perm.get('null_indirect_delta_r2_mean', np.nan), 2),
                               fmt_pq(perm.get('empirical_p_delta_r2'))))
        doc.html('<ul class="findings">{}</ul>'.format(''.join('<li>{}</li>'.format(x) for x in bullets)))

    def strongest_cases(self, limit=3):
        levels = self.run.levels
        if not len(levels):
            return []
        part = levels[levels['level'].map(is_step) & levels['explains']]
        part = part[part.apply(supported, axis=1)] if len(part) else part
        if not len(part):
            return []
        part = part.assign(qf=part['q'].astype(float)).sort_values(['qf', 'order'])
        out, seen = [], set()
        for r in part.itertuples():
            if r.gene in seen:
                continue
            seen.add(r.gene)
            out.append((r.gene, r.level, r.qf))
            if len(out) >= limit:
                break
        return out

    def bivalent_sentence(self):
        table = self.run.bivalent
        if table is None or not len(table):
            return ''
        threshold = self.run.command_value('--stateThreshold', 2.0)
        if 'threshold' in table and (table['threshold'] == threshold).any():
            table = table[table['threshold'] == threshold]
        table = table[(table['level'] == 'promoter') & table['transition'].str.contains('bivalent')]
        b = esc(self.run.names[1])
        stay = table[table['direction'] == 0]
        changed = table[table['direction'] != 0].sort_values('genes', ascending=False).head(3)
        parts = []
        for r in changed.itertuples():
            a_state, b_state, note = parse_transition(r.transition)
            parts.append('{:,} change from {} to {}{} ({} with the expected expression change)'.format(
                int(r.genes), a_state, b_state, ' by {}'.format(readable_note(note)) if note else '',
                fmt_int(r.concordant) if not isnan(r.concordant) else 'none'))
        stay_n = int(stay.loc[stay['transition'] == 'bivalent -> bivalent', 'genes'].sum()) if len(stay) else 0
        if not parts and not stay_n:
            return ''
        text = 'At the promoters of differential genes, '
        if stay_n:
            text += '{:,} are bivalent in both {}s'.format(stay_n, self.run.kind)
        if parts:
            text += '{}in {} {}'.format('; ' if stay_n else '', b, join_words(parts))
        return text + '.'

    def community_sentence(self):
        run = self.run
        table = run.community_expression
        summary = run.community_summary or {}
        if table is None or not len(table) or 'coherence_q' not in table:
            return ''
        coordinated = table[table['coherence_q'] < 0.1]
        tested = summary.get('communities_tested', len(table))
        if not len(coordinated):
            return 'None of {:,} matched community pairs changes expression coherently (q &lt; 0.1).'.format(tested)
        items = []
        for r in coordinated.head(4).itertuples():
            members_a, members_b, _ = self.community_members(int(r.community_a), int(r.community_b))
            name = self.community_name(members_a | members_b)
            items.append('the community of {} ({:,} genes, mean log2 {}, {})'.format(
                gene_span(name), int(r.n_genes), fmt_log2(r.mean_rna_log2FC),
                self.transition_phrase(r.transition) or esc(r.transition)))
        return 'Of {:,} matched community pairs, {} change{} expression coherently (q &lt; 0.1): {}.'.format(
            tested, number_word(len(coordinated)), 's' if len(coordinated) == 1 else '', join_words(items))

    def cascade_section(self, doc):
        run = self.run
        doc.section('cascade', 'Explanation cascade')
        if run.cascade is None:
            doc.p('The run has no cascade table.', cls='note')
            return
        doc.p('Each differential gene is assigned the nearest level that explains its expression change, from the '
              'promoter outwards: the histone marks, own contacts and chromatin state of the promoter, then its direct '
              'partners, then the nodes 2 and 3 steps away in the graph (left). The right panel counts every level '
              'that explains a gene, so one gene can count at several levels. A level explains a gene when its marks, '
              'contacts or chromatin state change in the direction expected from the expression change by at least '
              'log2 {:g}; away from the promoter the change is taken over the nodes at that step.'.format(
                  run.chromatin_threshold))
        fig, total = cascade_figure(run)
        doc.figure([fig], 'Differential genes ({:,}) by the nearest level that explains their expression change '
                          '(left) and by every level that explains it (right), with the percentage of all '
                          'differential genes.'.format(total))

    def elements_section(self, doc):
        """Element classes of the fragments (declared roles, fragment mean), promoters among the other ends, class
        transitions with the expression change of the genes in contact at step 1, and the roles learned per element
        type."""
        run = self.run
        counts = run.element_classes
        if counts is None or not len(counts):
            return
        a, b = run.names
        doc.section('elements', 'Elements')
        part = counts[(counts['roles'] == 'declared') & (counts['presence'] == 'mean')]
        sentences = []
        if 'other_end_fragments' in part:
            ends = part[part['condition'].astype(str) == str(a)]
            total = int(ends['other_end_fragments'].sum())
            captured = int(ends.loc[ends['tss_status'] == 'baited promoter', 'other_end_fragments'].sum())
            uncaptured = int(ends.loc[ends['tss_status'] == 'non-baited promoter', 'other_end_fragments'].sum())
            if total:
                sentences.append('Of the {:,} fragments at the other end of an interaction, {:,} ({:.0f}%) hold an '
                                 'annotated transcription start site: {:,} are captured promoters themselves and '
                                 '{:,} are promoters that were not captured.'.format(
                                     total, captured + uncaptured, 100 * (captured + uncaptured) / total, captured,
                                     uncaptured))
        sentences.append('Every fragment is classed in each {} from the marks present on it and the declared roles '
                         'of the marks (Table {}).'.format(run.kind, doc.table_count + 1))
        doc.p(' '.join(sentences), cls='lead')
        pivot = part.pivot_table(index='element_class', columns='condition', values='fragments', aggfunc='sum',
                                 fill_value=0)
        order = {c: k for k, c in enumerate(ELEMENT_CLASS_ORDER)}
        rows = []
        for cls in sorted(pivot.index, key=lambda c: (order.get(c, 99), c)):
            rows.append({'Element class': cls, a: int(pivot.loc[cls].get(a, 0)), b: int(pivot.loc[cls].get(b, 0))})
        doc.table(pd.DataFrame(rows), 'Fragments of the union graph by element class in {} and {} (declared mark '
                                      'roles, presence at a fold enrichment of at least {:g}, fragment mean). '
                                      'Promoters hold an annotated transcription start site, captured or '
                                      'not.'.format(esc(a), esc(b), run.command_value('--stateThreshold', 2.0)),
                  formats={a: fmt_int, b: fmt_int}, numeric=(a, b))
        transitions = run.element_transitions
        if transitions is not None and len(transitions):
            t = transitions[(transitions['roles'] == 'declared') & (transitions['presence'] == 'mean') &
                            (transitions['step'] == 1) & transitions['changed'].map(truthy)]
            t = t[t['n_genes'] >= 20].sort_values('n_genes', ascending=False).head(8)
            if len(t):
                baseline = float(t['baseline_mean'].iloc[0])

                def l2(value):
                    return fmt_num(value, 2, signed=True).replace('-', '−')
                frame = pd.DataFrame({
                    'Class in {}, in {}'.format(a, b): ['{} to {}'.format(x, y) for x, y in zip(t['class_A'],
                                                                                              t['class_B'])],
                    'Fragments': t['n_fragments'].map(fmt_int).to_numpy(),
                    'Genes in contact': t['n_genes'].map(fmt_int).to_numpy(),
                    'Mean log2 expression change': ['{} ({} to {})'.format(l2(m), l2(lo), l2(hi))
                                                    if not isnan(lo) else l2(m) for m, lo, hi in
                                                    zip(t['mean_rna_log2FC'], t['ci_low'], t['ci_high'])]})
                doc.table(frame, 'Element classes that change between {a} and {b} at the direct partners of '
                                 'promoters, the most frequent first (at least 20 genes): the genes whose promoter '
                                 'contacts at least one fragment with that change, and their mean log2 expression '
                                 'change {b} / {a} with a 95% chromosome block-bootstrap interval. All genes with a '
                                 'direct partner: mean log2 {m}. Steps 0 to {k} and all role sets are in the '
                                 'element transition table of the run.'.format(
                                     a=esc(a), b=esc(b), m=l2(baseline), k=run.max_hops),
                          numeric=('Fragments', 'Genes in contact'))
        roles = run.roles_by_type
        if roles is not None and len(roles):
            r = roles[roles['element_type'].isin(['promoter', 'non-promoter element']) &
                      roles['level'].astype(str).str.fullmatch(r'step\d+')]
            if len(r):
                rows = []
                for (mark, level), group in r.groupby(['mark', 'level'], sort=False):
                    row = {'Mark': mark, 'Step': level_hops(level)}
                    for kind, label in (('promoter', 'Promoter nodes'), ('non-promoter element', 'Other elements')):
                        g = group[group['element_type'] == kind]
                        rho = fmt_num(g['rho'].iloc[0], signed=True).replace('-', '\u2212') if len(g) else ''
                        row[label] = '{} (rho {})'.format(g['learned_role'].iloc[0], rho) \
                            if len(g) else 'no nodes'
                    rows.append(row)
                frame = pd.DataFrame(rows).sort_values(['Mark', 'Step'], kind='mergesort')
                doc.table(frame, 'Roles of the marks learned from the data at each step, separately for the nodes '
                                 'that are promoters (captured or not) and for the other elements: the association '
                                 '(Spearman rho) of the mark change summed over those nodes with the expression '
                                 'change; active or repressive when the 95% chromosome block-bootstrap interval '
                                 'excludes 0, otherwise no consistent effect. Roles per chromatin class are in the '
                                 'mark roles table of the run.')

    def select_findings(self):
        """Per cascade level the genes with the strongest support among those first explained there: beyond the
        mirror first, then by q, P and the absolute expression change; at the promoter by the number of explaining
        levels and the absolute expression change. One promoter per gene name."""
        run = self.run
        first = run.first_level()
        levels = run.levels
        chosen = {}
        used = set()
        order = ['gene', 'structural_variant'] + ['step{}'.format(k) for k in range(1, run.max_hops + 1)]
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
            if level in ('gene', 'structural_variant'):
                frame = frame.sort_values(['n_levels', 'rna'], ascending=[False, False])
            elif frame['q'].notna().any():
                frame = frame.assign(nq=frame['q'].fillna(2), np_=frame['p'].fillna(2))
                frame = frame.sort_values(['beyond', 'nq', 'np_', 'rna'], ascending=[False, True, True, False])
            else:
                frame = frame.sort_values(['beyond', 'n_levels', 'rna'], ascending=[False, False, False])
            picked = []
            for row in frame.itertuples():
                key = row.gene or row.node
                if key in used:
                    continue
                used.add(key)
                picked.append(row)
                if len(picked) >= self.args.top:
                    break
            chosen[level] = picked
        return chosen

    def gene_cards(self, doc):
        run = self.run
        doc.section('genes', 'Genes')
        chosen = self.select_findings()
        doc.p('For each level of the cascade, the {} gene{} with the strongest support among the genes first explained '
              'there: beyond distance-matched mirror loci first, then by q value; at the promoter by the number of '
              'levels that explain the gene and the size of the expression change. Each card links to the full gene '
              'report.'.format(number_word(self.args.top), '' if self.args.top == 1 else 's'), cls='lead')
        cards = []
        for level, picked in chosen.items():
            if not picked:
                continue
            title = 'Explained at the promoter' if level == 'gene' else \
                'Explained first by a structural variant or copy number change near the promoter' \
                if level == 'structural_variant' else ('Explained first at the direct partners' if level == 'step1'
                                                       else 'Explained first {} steps away'.format(level_hops(level)))
            doc.heading(title, level=3, cls='group')
            for row in picked:
                gene = row.gene or row.node
                nodes, _ = run.gene_nodes(gene)
                nodes = [n for n in nodes if n in set(run.levels['node'])] or [row.node]
                story = self.story(gene, nodes)
                link = '{}_gene_{}.html'.format(slug(run.name), slug(gene))
                doc.html('<article class="card"><h4>{} <a class="report-link" href="{}">gene report</a></h4>'.format(
                    gene_span(gene), esc(link)), text=gene)
                self.glance_block(story, doc, arcs=lambda: self.arcs_markup(doc, gene, nodes, story),
                                  known_level=level if level != 'gene' else None)
                doc.html('</article>', text='')
                cards.append((gene, nodes, story))
        return cards

    def communities_section(self, doc):
        run = self.run
        doc.section('communities', 'Hubs and communities')
        table = run.community_expression
        if table is None or not len(table) or run.gene_communities is None:
            doc.p('The run has no community tables.', cls='note')
            return
        coordinated = table[table['coherence_q'] < 0.1] if 'coherence_q' in table else table.head(0)
        show = coordinated if len(coordinated) else table.head(self.args.maxCommunities)
        doc.p('Communities are groups of promoters that are densely connected through shared interacting fragments, '
              'matched between {} and {}. {}'.format(
                  esc(run.names[0]), esc(run.names[1]),
                  'Shown are the communities whose genes change expression more coherently than random gene sets of '
                  'the same size (q &lt; 0.1).' if len(coordinated) else 'No community reaches q &lt; 0.1 in the '
                  'coherence test; the {} with the lowest P are shown.'.format(number_word(len(show)))), cls='lead')
        for row in show.head(self.args.maxCommunities).itertuples():
            ca, cb = int(row.community_a), int(row.community_b)
            members_a, members_b, _ = self.community_members(ca, cb)
            name = self.community_name(members_a | members_b)
            doc.html('<article class="card">', text='')
            doc.heading('Community of {}'.format(gene_span(name)), level=3)
            self.community_block(doc, ca, cb, info={'transition': row.transition})
            doc.html('</article>', text='')

    # community report -------------------------------------------------------------------------------------------
    def community_report(self, community, side):
        run = self.run
        doc = Document('Community {} ({}): {} vs {}'.format(community, run.names[0 if side == 'A' else 1],
                                                             *run.names), self.keep_figures)
        table = run.community_expression
        info, ca, cb = {}, -1, -1
        if table is not None and len(table):
            column = 'community_a' if side == 'A' else 'community_b'
            hit = table[table[column] == community]
            if len(hit):
                row = hit.iloc[0]
                ca, cb = int(row['community_a']), int(row['community_b'])
                info = {'transition': row.get('transition', '')}
        if ca < 0 and cb < 0 and run.gene_communities is not None:
            rows = self.community_rows(side, community)
            if rows is not None and len(rows):
                other = 'community_b' if side == 'A' else 'community_a'
                partner = rows[other][rows[other] >= 0].mode()
                partner = int(partner.iloc[0]) if len(partner) else -1
                ca, cb = (community, partner) if side == 'A' else (partner, community)
                transition = rows['community_a_transition'].dropna()
                info = {'transition': transition.iloc[0] if len(transition) else ''}
        name = ''
        if (ca >= 0 or cb >= 0) and run.gene_communities is not None:
            members_a, members_b, _ = self.community_members(ca, cb)
            name = self.community_name(members_a | members_b)
        heading = self.title_block('Community {}'.format(community) if not name else
                                   'Community of {}'.format(gene_span(name)), esc(run.condition_title()),
                                   'Community report')
        doc.section('community', 'Community')
        if ca < 0 and cb < 0:
            doc.p('Community {} is not found in the run outputs.'.format(community))
            self.appendix(doc, genes=[])
            return doc, heading
        self.community_block(doc, ca, cb, info=info)
        self.appendix(doc, genes=[], communities=[(ca, cb)])
        return doc, heading

    def community_rows(self, side, community):
        table = self.run.gene_communities
        if table is None:
            return None
        column = 'community_a' if side == 'A' else 'community_b'
        return table[table[column] == community]

    # appendix ---------------------------------------------------------------------------------------------------
    def appendix(self, doc, genes=(), community=None, communities=(), main=False):
        run = self.run
        doc.begin_appendix()
        a, b = run.names
        doc.section('methods', 'Methods', level=3)
        doc.p('Comparison of {a} (condition A) with {b} (condition B), outputs of egaDifferential under the prefix '
              '<code>{p}</code>{e}. Differential genes: |log2 fold change| of at least {r:g}; a level explains a gene '
              'when its mark change (log2 (B + 1) / (A + 1)), contact change or chromatin state transition agrees '
              'with the expression change by at least log2 {c:g} ({src}). Mark roles declared: {roles}. The result '
              'under the roles learned from the data is mentioned in the body only where it differs. P and q: '
              'per-gene permutation test of a step against random nodes; distance control: the explaining marks '
              'change more than at mirrored loci of the same genomic distance. Places are named by the gene of the '
              'nearest baited promoter{s}. The neighbourhood of a gene reaches its deepest explaining step and is '
              'reduced to at most {m} nodes; the nodes that carry the explanation at step k are the fragments at '
              'exactly k steps with a concordant change of at least the threshold, the 12 strongest per step.{rc}'
              .format(a=esc(a), b=esc(b), p=esc(run.prefix), e=', explanations of egaExplain under <code>{}</code>'
                      .format(esc(run.explain_prefix)) if run.explain_prefix else '', r=run.rna_threshold,
                      c=run.chromatin_threshold, src='thresholds from the run log' if run.command else
                      'egaDifferential defaults; no command line in the run log',
                      roles=esc(', '.join('{} {}'.format(k, v) for k, v in run.declared_roles().items()
                                          if k != 'contact')),
                      s='; upstream and downstream from the strand in <code>{}</code>'.format(
                          esc(self.args.annotation)) if self.places.strands else
                      '; without a gene annotation (--annotation) only the distance is given',
                      m=self.args.maxNodes,
                      rc=' Re-ChIP signal from <code>{}</code> (fragment mean; at least 2 confirms a bivalent call).'
                      .format(esc(run.rechip_path)) if run.rechip is not None else ''))
        for gene, nodes, story in genes:
            self.appendix_gene(doc, gene, nodes, story)
        if community is not None:
            self.appendix_community(doc, int(community.get('community_a', -1)), int(community.get('community_b', -1)),
                                    nodes=genes[0][1] if genes else ())
        for ca, cb in communities:
            self.appendix_community(doc, ca, cb)
        if main:
            self.appendix_run_tables(doc)
        self.appendix_files(doc)

    def appendix_gene(self, doc, gene, nodes, story):
        run = self.run
        doc.section('appendix-{}'.format(slug(gene)), 'Numbers for {}'.format(gene), level=3)
        rows = []
        first = run.first_level()
        for node in nodes:
            data = run.graph.nodes[node] if node in run.graph else {}
            counts = {'gained': 0, 'lost': 0, 'shared': 0}
            if node in run.graph:
                for w in run.graph.neighbors(node):
                    counts[run.graph.edges[node, w].get('status', 'shared')] += 1
            rows.append({'fragment': node, 'genes': ', '.join(gene_names(data.get('genes', ''))),
                         'expression log2': fmt_num(data.get('rna_log2FC', np.nan), signed=True),
                         'gained': counts['gained'], 'lost': counts['lost'], 'shared': counts['shared'],
                         'first level': first.get(node, 'not differential')})
        doc.table(pd.DataFrame(rows), 'Promoter fragments of {}: coordinates, genes on the fragment, expression log2 '
                                      'fold change {} / {}, interactions gained, lost and shared, and the first '
                                      'explaining level (cascade name).'.format(gene_span(gene), esc(b_name(run)),
                                                                                esc(run.names[0])),
                  numeric=('expression log2', 'gained', 'lost', 'shared'), nowrap=('fragment',))
        for node, levels in story.fragments:
            frame = []
            for lv in sorted(levels, key=level_order):
                r = levels[lv]
                frame.append({'level': lv, 'explains': 'yes' if r['explains'] else 'no',
                              'changes (mirror)': r.get('all_changes') or '', 'nodes': fmt_int(r.get('nodes')),
                              'concordant': fmt_int(r.get('n_concordant')),
                              'strongest': '{}{}'.format(r.get('strongest_node') or '',
                                                         ' ({})'.format(r.get('strongest_text'))
                                                         if r.get('strongest_text') else ''),
                              'P': fmt_p(r.get('p')), 'q': fmt_p(r.get('q')),
                              'beyond mirror': '' if isnan(r.get('beyond_mirror')) else
                              ('yes' if truthy(r['beyond_mirror']) else 'no'),
                              'learned roles': '' if isnan(r.get('learned_explains')) else '{} ({}{})'.format(
                                  'yes' if r['learned_explains'] else 'no',
                                  r.get('learned_result') or '', ', q {}'.format(fmt_p(r['learned_q']))
                                  if not isnan(r.get('learned_q')) else ''),
                              'max 1 kb': '' if isnan(r.get('alt_explains')) or lv not in
                              ('promoter_state', 'step1_state') else '{} ({})'.format(
                                  'yes' if r['alt_explains'] else 'no',
                                  r.get('alt_transition') or r.get('alt_result') or ''),
                              '3D support': '{}{}'.format(fmt_num(r.get('support')), ' ({})'.format(r['support_ci'])
                                                          if r.get('support_ci') and r['support_ci'] != 'NA' else '')
                              if not isnan(r.get('support')) else ''})
            frame = pd.DataFrame(frame)
            frame = frame[[c for c in frame.columns if (frame[c] != '').any()]]
            doc.table(frame, 'All levels of {} at fragment {} (declared roles, fragment mean), with the mark changes '
                             'at the mirrored loci in brackets (states: calls with the fold enrichment {} / {} and '
                             'the state fragment), the coordinates of the strongest node, the result under the '
                             'learned roles and on the strongest 1 kb window, and the 3D support with its '
                             'interval.'.format(gene_span(gene), esc(node), esc(run.names[0]), esc(run.names[1])),
                      long=('changes (mirror)', 'strongest'), numeric=('nodes', 'concordant', 'P', 'q'),
                      cls='full')

    def appendix_community(self, doc, community_a, community_b, nodes=()):
        run = self.run
        if run.gene_communities is None or (community_a < 0 and community_b < 0):
            return
        members_a, members_b, rows = self.community_members(community_a, community_b)
        doc.section('appendix-community-{}-{}'.format(community_a, community_b),
                    'Community {} of {} and {} of {}'.format(community_a, run.names[0], community_b, run.names[1]),
                    level=3)
        show = rows.sort_values('rna_log2FC', key=np.abs, ascending=False)
        frame = pd.DataFrame({'gene': show['genes'].map(lambda g: ', '.join(gene_names(g))),
                              'fragment': show['node'], 'membership': show['in'],
                              'expression log2': show['rna_log2FC'].map(lambda v: fmt_num(v, signed=True)),
                              'differential': show['differential'].map(lambda v: 'yes' if truthy(v) else 'no')
                              if 'differential' in show else ''})
        features = run.features.set_index('node') if run.features is not None else None
        if features is not None:
            for mark in run.marks:
                if mark + '_promoter' in features:
                    frame['{} promoter'.format(mark)] = features[mark + '_promoter'].reindex(show['node']).map(
                        lambda v: fmt_num(v, signed=True)).to_numpy()
        doc.table(frame, 'Promoters of the community by absolute expression change (up to 40): membership (both, left '
                         'in {b}, joined in {b}), expression log2 fold change {b} / {a} and the log2 change of the '
                         'marks at the promoter.'.format(a=esc(run.names[0]), b=esc(run.names[1])), max_rows=40,
                  numeric=tuple(c for c in frame.columns if c.endswith(('log2', 'promoter'))), nowrap=('fragment',))

    def appendix_run_tables(self, doc):
        run = self.run
        doc.section('appendix-run', 'Run tables', level=3)
        counts = run.edge_counts()
        rows = [('differential genes', '{:,}'.format(len(run.classes))),
                ('promoters with an expression change', fmt_int(len(run.features)) if run.features is not None else ''),
                ('interactions gained in {}'.format(run.names[1]), '{:,}'.format(counts.get('gained', 0))),
                ('interactions lost in {}'.format(run.names[1]), '{:,}'.format(counts.get('lost', 0))),
                ('interactions shared', '{:,}'.format(counts.get('shared', 0)))]
        perm = run.permutation or {}
        for key, label in (('n_genes', 'genes in the models'), ('observed_indirect_delta_r2', 'R2 gain, 2 steps'),
                           ('null_indirect_delta_r2_mean', 'R2 gain, 2 steps, permutations (mean)'),
                           ('empirical_p_delta_r2', 'empirical P of the R2 gain'),
                           ('n_permutations', 'permutations')):
            if key in perm:
                rows.append((label, fmt_num(perm[key], 5) if isinstance(perm[key], float) else str(perm[key])))
        summary = run.community_summary or {}
        if summary.get('transitions'):
            rows.append(('community transitions', ', '.join('{} {:,}'.format(k, v)
                                                             for k, v in summary['transitions'].items())))
        if summary.get('gene_status'):
            rows.append(('promoters by community fate', ', '.join('{} {:,}'.format(k, v)
                                                                 for k, v in summary['gene_status'].items())))
        doc.table(pd.DataFrame(rows, columns=['quantity', 'value']), 'Counts and graph-level tests of the run.',
                  long=('value',))
        if run.cascade is not None:
            cascade = run.cascade
            columns = [c for c in ('level', 'genes', 'genes_reaching', 'fraction_of_differential',
                                   'fraction_of_remaining', 'beyond_mirror', 'p_lt_0.05', 'q_lt_0.1', 'q_lt_0.05',
                                   'proximity_support', 'proximity_ci', 'median_nodes_tested') if c in cascade]
            doc.table(cascade[columns], 'Explanation cascade (first explaining level).',
                      formats={'fraction_of_differential': fmt_pct, 'fraction_of_remaining': fmt_pct,
                               'beyond_mirror': fmt_int, 'p_lt_0.05': fmt_int, 'q_lt_0.1': fmt_int,
                               'q_lt_0.05': fmt_int,
                               'proximity_support': fmt_num, 'median_nodes_tested': fmt_num},
                      numeric=tuple(c for c in columns if c != 'level'))
        if run_has(run.level_counts):
            counts = run.variant_rows(run.level_counts)
            columns = [c for c in ('level', 'genes_any', 'genes_first', 'fraction_any', 'fraction_first',
                                   'any_transition', 'any_poised_in_both', 'any_beyond_mirror', 'any_p_lt_0.05',
                                   'any_q_lt_0.05', 'proximity_support') if c in counts]
            ints = [c for c in columns if c.startswith(('genes_', 'any_', 'first_'))]
            doc.table(counts[columns], 'Genes explained at every level (any) and first explained there (first), '
                                       'declared roles, fragment mean.',
                      formats=dict({c: fmt_int for c in ints}, fraction_any=fmt_pct, fraction_first=fmt_pct,
                                   proximity_support=fmt_num, median_nodes_tested=fmt_num),
                      numeric=tuple(c for c in columns if c != 'level'))
        combos, from_file = run.level_combinations()
        doc.table(combos, 'Most frequent combinations of explaining levels per gene (first 15).', max_rows=15,
                  formats={'genes': fmt_int, 'fraction_of_differential': fmt_pct, 'genes_up': fmt_int,
                           'genes_down': fmt_int},
                  numeric=('rank', 'n_levels', 'genes', 'genes_up', 'genes_down', 'fraction_of_differential'))
        if run_has(run.types):
            doc.table(run.types, 'Explanation types at the promoter and at the direct partners.',
                      formats={c: fmt_int for c in ('genes', 'genes_up', 'genes_down', 'mark_gained', 'mark_lost')},
                      numeric=('genes', 'genes_up', 'genes_down', 'mark_gained', 'mark_lost'))
        if run_has(run.roles):
            doc.figure([roles_figure(run)], 'Association of every mark change with the expression change at the '
                                            'promoter and at each step (Spearman rho with block-bootstrap interval). '
                                            'Filled markers: learned role equals the declared role; open: it '
                                            'differs.')
            roles = run.roles
            frame = pd.DataFrame({'mark': roles['mark'], 'level': roles['level'], 'versus': roles.get('versus', ''),
                                  'declared': roles.get('declared_role', ''), 'learned': roles.get('learned_role', ''),
                                  'rho': roles['rho'].map(fmt_num) if 'rho' in roles else '',
                                  'genes': roles.get('n_genes', np.nan)})
            doc.table(frame, 'Declared and learned role of every mark per level.', formats={'genes': fmt_int},
                      numeric=('genes', 'rho'))
        if run_has(run.state_transitions):
            table = run.state_transitions
            if 'changed' in table:
                table = table[table['changed'].map(truthy)]
            if 'n_genes' in table:
                table = table.sort_values('n_genes', ascending=False)
            columns = [c for c in ('level', 'transition', 'n_fragments', 'n_genes', 'mean_rna_log2FC', 'n_up',
                                   'n_down', 'concordance') if c in table]
            doc.table(table[columns], 'Changed mark combinations at the promoter (step 0) and the direct partners '
                                      '(step 1), first 20 rows.', max_rows=20,
                      formats={'n_fragments': fmt_int, 'n_genes': fmt_int, 'n_up': fmt_int, 'n_down': fmt_int,
                               'mean_rna_log2FC': lambda v: fmt_num(v, signed=True), 'concordance': fmt_num},
                      numeric=tuple(c for c in columns if c not in ('level', 'transition')))
        if run_has(run.bivalent):
            threshold = run.command_value('--stateThreshold', 2.0)
            table = run.bivalent
            if 'threshold' in table and (table['threshold'] == threshold).any():
                table = table[table['threshold'] == threshold]
            table = table.sort_values('genes', ascending=False) if 'genes' in table else table
            doc.table(table, 'Declared states (bivalent promoter, poised enhancer) at the presence threshold {:g}, '
                             'first 20 rows.'.format(threshold), max_rows=20,
                      formats={'threshold': fmt_num, 'genes': fmt_int, 'genes_up': fmt_int, 'genes_down': fmt_int,
                               'concordant': fmt_int, 'supported_concordant': fmt_int,
                               'fraction_concordant': fmt_pct, 'direction': fmt_int},
                      numeric=('threshold', 'genes', 'genes_up', 'genes_down', 'concordant', 'supported_concordant',
                               'fraction_concordant'))
        if run_has(run.community_expression):
            table = run.community_expression.head(15)
            columns = [c for c in ('community_a', 'community_b', 'transition', 'n_genes', 'edges_shared',
                                   'edges_lost', 'edges_gained', 'rewired_fraction', 'mean_rna_log2FC', 'coherence_p',
                                   'coherence_q') if c in table]
            doc.table(table[columns], 'Matched communities with the lowest coherence P (first 15).',
                      formats={'rewired_fraction': fmt_pct, 'mean_rna_log2FC': lambda v: fmt_num(v, signed=True),
                               'coherence_p': fmt_p, 'coherence_q': fmt_p},
                      numeric=tuple(c for c in columns if c != 'transition'))
        if run_has(run.edge_chromatin):
            doc.table(run.edge_chromatin, 'Mark change at the two ends of gained and lost interactions against shared '
                                          'interactions of the same distance (egaExplain).',
                      formats={'n': fmt_int, 'median_change': lambda v: fmt_num(v, 3, True),
                               'shared_median_change': lambda v: fmt_num(v, 3, True),
                               'concordance_vs_shared': fmt_num, 'mannwhitney_p_expected_direction': fmt_p},
                      numeric=('n', 'median_change', 'shared_median_change', 'concordance_vs_shared',
                               'mannwhitney_p_expected_direction'))

    def appendix_files(self, doc):
        run = self.run
        doc.section('appendix-files', 'Files, parameters and versions', level=3)
        rows = []
        for label, path in sorted(run.used.items()):
            try:
                stamp = datetime.datetime.fromtimestamp(os.path.getmtime(path)).strftime('%Y-%m-%d %H:%M')
            except OSError:
                stamp = ''
            rows.append((label, path, stamp))
        doc.table(pd.DataFrame(rows, columns=['file', 'path', 'modified']), 'Files read for this report.',
                  long=('path',))
        parameters = [('levels read from', run.levels_source), ('expression threshold (|log2|)', run.rna_threshold),
                      ('mark change threshold (|log2|)', run.chromatin_threshold),
                      ('state presence threshold', run.command_value('--stateThreshold', 2.0)),
                      ('largest neighbourhood drawn (nodes)', self.args.maxNodes),
                      ('colour scale limit (log2)', self.args.maxLog2), ('genes per cascade level', self.args.top),
                      ('gene annotation', self.args.annotation or 'none')]
        doc.table(pd.DataFrame([(k, str(v)) for k, v in parameters], columns=['parameter', 'value']),
                  'Parameters of the report.', long=('value',))
        if run.command:
            doc.p('Command line: <code>{}</code>'.format(esc(run.command)))
        info = [line for line in run.log_lines if ' INFO: ' in line or 'WARNING' in line]
        if info:
            doc.table(pd.DataFrame({'run log': info[:40]}), 'Messages of the run log, first 40 lines.',
                      long=('run log',))
        import scipy
        versions = [('epigeneticgraphanalyzer (egaReport)', __version__), ('Python', platform.python_version()),
                    ('numpy', np.__version__), ('pandas', pd.__version__), ('networkx', nx.__version__),
                    ('matplotlib', matplotlib.__version__), ('scipy', scipy.__version__),
                    ('report date', datetime.date.today().isoformat())]
        doc.table(pd.DataFrame(versions, columns=['software', 'version']), 'Software versions.')


def run_has(table):
    return table is not None and len(table) > 0


def b_name(run):
    return run.names[1]


def level_short(level):
    short = {'promoter_marks': 'promoter marks', 'contacts': 'contacts', 'promoter_state': 'promoter state',
             'structural_variant': 'structural variant', 'copy_number': 'copy number', 'step1': 'step 1',
             'step1_state': 'step 1 state'}
    return short.get(level, 'step {}'.format(level_hops(level)) if re.fullmatch(r'step\d+', str(level)) else level)


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


def write_document(doc, path_stem, heading, args, nav=True):
    """HTML (appendix collapsed) and PDF (printed from a copy with the appendix open)."""
    html_path = path_stem + '.html'
    with open(html_path, 'w') as handle:
        handle.write(doc.render(heading, nav=nav))
    pdf_path = path_stem + '.pdf'
    engine = None
    if args.pdf in ('auto', 'chromium'):
        chromium = find_chromium(args.chromium)
        printed = False
        if chromium:
            print_path = path_stem + '.print.html'
            with open(print_path, 'w') as handle:
                handle.write(doc.render(heading, nav=nav, print_version=True))
            try:
                printed = chromium_pdf(chromium, print_path, pdf_path)
            finally:
                os.remove(print_path)
        if printed:
            engine = 'chromium'
        elif args.pdf == 'chromium':
            raise RuntimeError('Chromium could not print {}'.format(html_path))
        else:
            log.warning('Chromium not available or failed; writing the PDF with matplotlib (plain layout)')
    if engine is None and args.pdf in ('auto', 'matplotlib'):
        if doc.keep_figures:
            matplotlib_pdf(doc, pdf_path, plain_text(heading))
            engine = 'matplotlib'
        else:
            log.warning('no PDF for %s: figures were not kept for the matplotlib fallback', html_path)
    log.info('wrote %s%s', html_path, ' and {} ({})'.format(pdf_path, engine) if engine else '')
    for fig_list, _ in doc.figures:
        for fig in fig_list:
            fig.clear()
    return html_path, (pdf_path if engine else None), engine


def report_genes(run, args):
    """Genes whose strand is needed: the genes of the gene reports (the main report reads the strands of the
    selected genes later)."""
    return list(args.gene or [])


def main(args=None):
    args = parse_arguments(args)
    setup_logging(args.verbose)
    os.makedirs(args.outputDir, exist_ok=True)
    run = Run(args.prefix, args.explainPrefix, args.names, args.rnaThreshold, args.chromatinThreshold, args.log,
              args.name)
    run.set_labels(args.labels)
    run.read_rechip(args.rechip)
    log.info('%s: conditions %s, marks %s, %d differential genes, levels from %s', run.name, run.names, run.marks,
             len(run.classes), run.levels_source)
    # keep the figures for the matplotlib PDF when Chromium is not used or may fail
    chromium = find_chromium(args.chromium) if args.pdf in ('auto', 'chromium') else None
    keep = args.pdf == 'matplotlib' or (args.pdf == 'auto' and chromium is None)
    places = Places(run.graph, elements=run.element_lookup(), names=run.names)
    report = Report(run, args, keep, places)
    if args.annotation:
        genes = report_genes(run, args)
        if not args.gene and args.community is None:
            genes = [row.gene or row.node for picked in report.select_findings().values() for row in picked]
        places.strands = read_strands(args.annotation, genes)
        log.info('strand of %d of %d genes from %s', len(places.strands), len(genes), args.annotation)
    outputs = []
    if args.gene:
        pages = []
        for gene in args.gene:
            doc, heading = report.gene_report(gene)
            stem = os.path.join(args.outputDir, '{}_gene_{}'.format(slug(run.name), slug(gene)))
            outputs.append(write_document(doc, stem, heading, args, nav=False))
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
                handle.write(index.render(report.title_block('Gene reports', esc(run.condition_title()),
                                                             'Index'), nav=False))
            outputs.append((path, None, None))
    elif args.community is not None:
        doc, heading = report.community_report(args.community, args.communitySide)
        stem = os.path.join(args.outputDir, '{}_community_{}{}'.format(slug(run.name), args.communitySide,
                                                                       args.community))
        outputs.append(write_document(doc, stem, heading, args, nav=False))
    else:
        doc, heading = report.main_report()
        stem = os.path.join(args.outputDir, '{}_report'.format(slug(run.name)))
        outputs.append(write_document(doc, stem, heading, args))
    for html_path, pdf_path, engine in outputs:
        sys.stdout.write('{}\n'.format(html_path))
        if pdf_path:
            sys.stdout.write('{} ({})\n'.format(pdf_path, engine))
    return outputs


if __name__ == '__main__':
    main()
