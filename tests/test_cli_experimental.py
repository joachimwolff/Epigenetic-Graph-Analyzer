import pytest

from epigeneticgraphanalyzer import (egaBoundaries, egaCompartments, egaTimecourse, egaTriggers,
                                     switching)
from epigeneticgraphanalyzer._cli import EXPERIMENTAL, base_parser


def test_base_parser_marks_experimental_modules():
    assert base_parser('x', experimental=True).description.startswith(EXPERIMENTAL)
    assert base_parser('x').description == 'x'


@pytest.mark.parametrize('module', [switching, egaTriggers, egaCompartments, egaTimecourse, egaBoundaries])
def test_experimental_help(module, capsys):
    with pytest.raises(SystemExit):
        module.parse_args(['--help']) if hasattr(module, 'parse_args') else module.main(['--help'])
    assert '[Experimental]' in capsys.readouterr().out
