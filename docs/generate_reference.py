"""Command reference of the documentation: one Markdown page per console script of pyproject.toml, built from the
argument parser of the command (description, usage and every option with its help and default).

The parser of each command is captured by intercepting ArgumentParser.parse_args when the module's
parse_arguments() is called, so the pages always follow the code. Run after changing a command line:

    python docs/generate_reference.py
"""

import argparse
import importlib
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

GROUPS = [
    ('Graphs', ['egaCreateGraph', 'egaIndirectContacts', 'egaAnnotateGraph']),
    ('Multi-way data', ['egaBuildMultiway', 'egaValidateMultiway']),
    ('Two conditions and reports', ['egaDifferential', 'egaExplain', 'egaReport', 'egaExplorer']),
    ('Plots', ['egaPlotNeighbourhood', 'egaPlotLocus', 'egaPlotRegion', 'egaPlotDifferential']),
    ('Further modules', ['egaCompartments', 'egaBoundaries', 'egaTadCrossing', 'egaEqtl']),
]


class Captured(Exception):
    pass


def scripts():
    text = open(os.path.join(ROOT, 'pyproject.toml')).read()
    block = text.split('[project.scripts]', 1)[1].split('\n[', 1)[0]
    return dict(re.findall(r'^(\w+)\s*=\s*"([\w.]+):\w+"', block, re.M))


def capture(module_name):
    module = importlib.import_module(module_name)
    holder = {}
    original = argparse.ArgumentParser.parse_args

    def intercept(self, args=None, namespace=None):
        holder['parser'] = self
        raise Captured()
    argparse.ArgumentParser.parse_args = intercept
    try:
        module.parse_arguments([])
    except Captured:
        pass
    finally:
        argparse.ArgumentParser.parse_args = original
    return holder['parser']


def option_rows(parser):
    rows = []
    for action in parser._actions:
        if isinstance(action, argparse._HelpAction):
            continue
        names = ', '.join('`{}`'.format(o) for o in action.option_strings) or '`{}`'.format(action.dest)
        if action.metavar:
            value = action.metavar if isinstance(action.metavar, str) else ' '.join(action.metavar)
        elif action.choices is not None:
            value = '{' + ', '.join(str(c) for c in action.choices) + '}'
        elif action.nargs == 0:
            value = ''
        else:
            value = action.dest.upper()
        if action.nargs in ('+', '*'):
            value += ' ...'
        default = ''
        # machine-dependent defaults (the number of cores) are described in the help text instead
        if action.default not in (None, argparse.SUPPRESS, False) and action.nargs != 0 and action.dest != 'threads':
            default = '`{}`'.format(action.default)
        help_text = (action.help or '').replace('%(default)s', str(action.default)).replace('|', '\\|')
        help_text = ' '.join(help_text.split())
        option = names + (' `{}`'.format(value.replace('|', '\\|')) if value else '')
        if action.required:
            option += ' (required)'
        rows.append('| {} | {} | {} |'.format(option, default, help_text))
    return rows


def page(command, parser):
    lines = ['# {}'.format(command), '']
    if parser.description:
        lines += [' '.join(parser.description.split()), '']
    usage = parser.format_usage().replace('usage: ', '').strip()
    usage = re.sub(r'^\S+', command, usage)
    lines += ['```text', usage, '```', '', '| Option | Default | Description |', '|---|---|---|']
    lines += option_rows(parser)
    return '\n'.join(lines) + '\n'


def main():
    out = os.path.join(HERE, 'reference')
    os.makedirs(out, exist_ok=True)
    available = scripts()
    listed = [c for _, commands in GROUPS for c in commands]
    missing = sorted(set(available) - set(listed))
    if missing:
        raise SystemExit('commands without a group: {}'.format(', '.join(missing)))
    index = ['# Command reference', '',
             'Every command prints the same information with `--help`. These pages are generated from the argument '
             'parsers by `docs/generate_reference.py`.', '']
    index += ['```{toctree}', ':hidden:', ''] + [c for _, commands in GROUPS for c in commands] + ['```', '']
    for group, commands in GROUPS:
        index += ['## {}'.format(group), '']
        for command in commands:
            parser = capture(available[command])
            with open(os.path.join(out, command + '.md'), 'w') as handle:
                handle.write(page(command, parser))
            description = ' '.join((parser.description or '').split())
            experimental = description.startswith('[Experimental]')
            if experimental:
                description = description.split('exploratory. ', 1)[-1]
            first = description.split('. ')[0].rstrip('.') + (' (experimental)' if experimental else '')
            index.append('* [{}]({}.md): {}'.format(command, command, first))
        index.append('')
    with open(os.path.join(out, 'index.md'), 'w') as handle:
        handle.write('\n'.join(index))
    print('{} command pages written to {}'.format(len(listed), out))


if __name__ == '__main__':
    main()
