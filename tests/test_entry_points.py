import importlib
import re
import os


def test_console_scripts_resolve_and_wrappers_return_none(monkeypatch):
    """Every console script of pyproject.toml points to a callable; the cli() wrappers discard the return value of
    main(), which setuptools would otherwise print and turn into exit status 1."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    block = open(os.path.join(root, 'pyproject.toml')).read().split('[project.scripts]', 1)[1].split('\n[', 1)[0]
    targets = re.findall(r'^\w+\s*=\s*"([\w.]+):(\w+)"', block, re.M)
    assert len(targets) >= 17
    for module_name, function in targets:
        module = importlib.import_module(module_name)
        assert callable(getattr(module, function))
        if function == 'cli':
            monkeypatch.setattr(module, 'main', lambda args=None: ['an output'])
            assert module.cli([]) is None
