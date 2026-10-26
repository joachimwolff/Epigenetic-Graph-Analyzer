# Configuration file for the Sphinx documentation builder.

project = 'Epigenetic Graph Analyzer'
author = 'Joachim Wolff'

extensions = [
    'myst_parser',
    'sphinxcontrib.mermaid',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

master_doc = 'index'

# keep '--option' as written (no en dash)
smartquotes = False

# $...$ in the Markdown pages
myst_enable_extensions = ['dollarmath']
# anchors for the section headings, so that links of the form page.md#heading resolve
myst_heading_anchors = 3

html_static_path = ['_static']
html_css_files = ['wrap_tables.css']
