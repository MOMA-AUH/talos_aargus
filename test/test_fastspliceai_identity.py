"""FastSpliceAI predictions are reused only for identical inputs and scorer."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from talos import fastspliceai_identity as module


@pytest.fixture
def identity_inputs(tmp_path):
    paths = {}
    for name in ('gff', 'annotation', 'reference', 'reference_fai', 'selection'):
        paths[name] = tmp_path / name
        paths[name].write_text(name)
    return paths


def identity(paths, **changes: object):
    options = {'distance': 500, 'mask': 0, 'scorer': {'version': 'stub', 'weights': 'abc'}}
    options.update(changes)
    return module.create_identity(**paths, **options)


def test_same_content_and_settings_have_stable_key(identity_inputs, tmp_path):
    first = identity(identity_inputs)
    assert first == identity(identity_inputs)
    assert first['settings'] == {'distance': 500, 'mask': 0}
    moved = {}
    for name, path in identity_inputs.items():
        moved[name] = tmp_path / f'moved-{name}'
        moved[name].write_bytes(path.read_bytes())
    assert identity(moved)['key'] == first['key']


@pytest.mark.parametrize('name', ['gff', 'annotation', 'reference', 'reference_fai', 'selection'])
def test_resource_content_changes_invalidate_cache(identity_inputs, name):
    first = identity(identity_inputs)
    identity_inputs[name].write_text('changed bytes')
    assert identity(identity_inputs)['key'] != first['key']


@pytest.mark.parametrize(
    'settings',
    [
        {'distance': 100},
        {'mask': 1},
        {'scorer': {'version': 'new', 'weights': 'abc'}},
    ],
)
def test_scoring_changes_invalidate_cache(identity_inputs, settings):
    assert identity(identity_inputs, **settings)['key'] != identity(identity_inputs)['key']


def test_resource_mutation_during_capture_is_rejected(identity_inputs, monkeypatch):
    original = module._hash_file  # noqa: SLF001 - exercise mutation during the hash pass

    def mutate(path) -> str:
        digest = original(path)
        if path == identity_inputs['annotation']:
            identity_inputs['gff'].write_text('modified after hashing')
        return digest

    monkeypatch.setattr(module, '_hash_file', mutate)
    with pytest.raises(ValueError, match='changed'):
        identity(identity_inputs)


def test_scorer_source_and_weights_are_content_bound(tmp_path, monkeypatch):
    source = tmp_path / 'spliceai/run.py'
    weights = tmp_path / 'spliceai/models/default.h5'
    weights.parent.mkdir(parents=True)
    source.write_text('source')
    weights.write_text('weights')
    distribution = SimpleNamespace(
        version='2.0.0',
        files=[Path('spliceai/run.py'), Path('spliceai/models/default.h5')],
        locate_file=lambda path: tmp_path / path,
        read_text=lambda _name: json.dumps({'vcs_info': {'commit_id': '1234'}}),
    )
    monkeypatch.setattr(module, 'distribution', lambda _name: distribution)
    first = module.scorer_identity()
    assert first['version'] == '2.0.0'
    assert first['source_revision'] == '1234'
    with monkeypatch.context() as runtime:
        runtime.setattr(module, 'version', lambda _name: 'changed-runtime')
        assert module.scorer_identity() != first
    weights.write_text('different weights')
    assert module.scorer_identity() != first
    weights.write_text('weights')
    source.write_text('different source')
    assert module.scorer_identity() != first


def test_unidentified_model_weights_are_rejected(monkeypatch):
    monkeypatch.setattr(module, 'distribution', lambda _name: SimpleNamespace(version='2', files=[]))
    with pytest.raises(ValueError, match='weights'):
        module.scorer_identity()


def test_cli_identity_includes_selection_provenance(identity_inputs, tmp_path, monkeypatch):
    scorer = {'version': 'stub', 'weights': 'abc'}
    monkeypatch.setattr(module, 'scorer_identity', lambda: scorer)
    output = tmp_path / 'identity.json'
    arguments = ['fastspliceai_identity']
    for name, path in identity_inputs.items():
        arguments.extend(['--' + name.replace('_', '-'), str(path)])
    arguments.extend(['--distance', '500', '--mask', '0', '--output', str(output)])
    monkeypatch.setattr('sys.argv', arguments)
    module.main()
    assert json.loads(output.read_text()) == identity(identity_inputs)
