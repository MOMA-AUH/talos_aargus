"""AARgus captures contain only the supported contract and verified artifacts."""

import hashlib
import json
import os
from typing import Any

import pytest

from talos import workflow_metadata as module
from talos.version import __version__
from talos.workflow_metadata import _artifact, hash_artifact, main, prepare_capture, verify_capture


def inputs(tmp_path):
    config = tmp_path / 'config.toml'
    config.write_text('[RunHailFiltering]\nspliceai = 0.5\n')
    resource = tmp_path / 'genes.gff3'
    resource.write_text('##gff-version 3\n')
    spec = tmp_path / 'spec.json'
    spec.write_text(
        json.dumps(
            {
                'assembly': 'GRCh38',
                'resources': [
                    {
                        'name': 'Ensembl',
                        'version': '115',
                        'path': str(resource),
                    }
                ],
            }
        )
    )
    return config, resource, spec


def test_rejects_duplicate_resource_names_and_missing_versions(tmp_path):
    config, _, spec = inputs(tmp_path)
    data = json.loads(spec.read_text())
    data['resources'].append(dict(data['resources'][0]))
    spec.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='Duplicate'):
        prepare_capture(spec, config, tmp_path / 'result.json', tmp_path / 'capture')
    data['resources'] = [data['resources'][0]]
    data['resources'][0]['version'] = ''
    spec.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='version'):
        prepare_capture(spec, config, tmp_path / 'result.json', tmp_path / 'capture')


@pytest.mark.parametrize('installed_version', ['2.0.0', '2.1.0'])
def test_fastspliceai_records_software_version_and_annotation_hash(tmp_path, monkeypatch, installed_version):
    config, annotation, spec = inputs(tmp_path)
    data = json.loads(spec.read_text())
    data['resources'][0].update(name='FastSpliceAI', version='1.5', parameter='fastspliceai_annotation')
    spec.write_text(json.dumps(data))

    def installed(distribution: str) -> str:
        assert distribution == 'spliceai'
        return installed_version

    monkeypatch.setattr(module, 'distribution_version', installed)
    capture = prepare_capture(
        spec,
        config,
        tmp_path / 'result.json',
        tmp_path / 'capture',
        parameters={'resource_paths': {'fastspliceai_annotation': str(annotation)}},
    )
    assert capture['resources'] == [
        {
            'name': 'FastSpliceAI',
            'version': installed_version,
            'sha256': hashlib.sha256(annotation.read_bytes()).hexdigest(),
        }
    ]


def test_result_schema_and_software_are_checked(tmp_path):
    config, _, spec = inputs(tmp_path)
    result = tmp_path / 'results.json'
    capture = prepare_capture(spec, config, result, tmp_path / 'capture')
    result.write_text(json.dumps({'version': '2.3.0', 'metadata': {'version': capture['software_version']}}))
    verified = verify_capture(tmp_path / 'capture/capture.json', result)
    assert set(verified['result']) == {'path', 'sha256'}
    result.write_text(json.dumps({'version': '2.3.0', 'metadata': {'version': 'other'}}))
    with pytest.raises(ValueError, match='software'):
        verify_capture(tmp_path / 'capture/capture.json', result)


def test_bundle_is_immutable_and_uses_published_paths(tmp_path):
    config, _, spec = inputs(tmp_path)
    destination = tmp_path / 'published/cohort.metadata'
    bundle = tmp_path / 'task/cohort.metadata'
    capture = prepare_capture(spec, config, tmp_path / 'result.json', bundle, published_directory=destination)
    assert capture['configuration']['path'] == str(destination / 'configuration.toml')
    assert {path.name for path in bundle.iterdir()} == {'capture.json', 'configuration.toml'}
    assert (bundle / 'configuration.toml').read_bytes() == config.read_bytes()
    with pytest.raises(FileExistsError):
        prepare_capture(spec, config, tmp_path / 'result.json', bundle)


def test_cli_prepare_and_verify(tmp_path, monkeypatch):
    config, _, spec = inputs(tmp_path)
    result = tmp_path / 'results.json'

    result.write_text(json.dumps({'version': '2.3.0', 'metadata': {'version': __version__}}))
    monkeypatch.setattr(
        'sys.argv',
        [
            'metadata',
            'prepare',
            '--spec',
            str(spec),
            '--configuration',
            str(config),
            '--result-path',
            str(result),
            '--output-directory',
            str(tmp_path / 'capture'),
            '--verify-result',
            str(result),
            '--parameters-json',
            '{"delta": 0.2}',
        ],
    )
    main()
    capture = json.loads((tmp_path / 'capture/capture.json').read_text())
    assert capture['result_path'] == str(result)
    assert set(capture) == {
        'software_version',
        'assembly',
        'resources',
        'configuration',
        'result_path',
    }
    assert all(set(resource) == {'name', 'version', 'sha256'} for resource in capture['resources'])
    assert set(capture['configuration']) == {'path', 'sha256'}
    assert all(len(resource['sha256']) == 64 for resource in capture['resources'])
    monkeypatch.setattr(
        'sys.argv',
        [
            'metadata',
            'verify',
            '--capture',
            str(tmp_path / 'capture/capture.json'),
            '--result',
            str(result),
            '--output',
            str(tmp_path / 'verified.json'),
        ],
    )
    main()
    assert json.loads((tmp_path / 'verified.json').read_text())['result']['path'] == str(result)


def test_unsupported_result_is_rejected_without_a_capture(tmp_path):
    config, _, spec = inputs(tmp_path)
    result = tmp_path / 'results.json'
    result.write_text('{"version":"2.2.0","metadata":{}}')
    with pytest.raises(ValueError, match=r'schema 2\.3\.0'):
        prepare_capture(spec, config, result, tmp_path / 'capture', observed_result=result)
    assert not (tmp_path / 'capture').exists()


def test_declared_resource_must_match_bound_workflow_input(tmp_path):
    config, resource, spec = inputs(tmp_path)
    data = json.loads(spec.read_text())
    data['resources'][0]['parameter'] = 'ensembl_gff'
    spec.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='actual workflow resource path'):
        prepare_capture(spec, config, tmp_path / 'results.json', tmp_path / 'capture')
    other = tmp_path / 'other.gff3'
    other.write_text('different annotation bytes')
    with pytest.raises(ValueError, match='checksum differs'):
        prepare_capture(
            spec,
            config,
            tmp_path / 'results.json',
            tmp_path / 'capture',
            parameters={'resource_paths': {'ensembl_gff': str(other)}},
        )
    prepare_capture(
        spec,
        config,
        tmp_path / 'results.json',
        tmp_path / 'capture',
        parameters={'resource_paths': {'ensembl_gff': str(resource)}},
    )


def test_receipts_allow_resource_assembly_without_reading_unchanged_resources(tmp_path, monkeypatch):
    config, resource, spec = inputs(tmp_path)
    receipt = tmp_path / 'resource-hash.json'
    receipt.write_text(json.dumps(module.hash_artifact(resource)))
    original_open = type(resource).open

    def guarded_open(path, *args: Any, **kwargs: Any) -> Any:
        if path.resolve() == resource.resolve():
            raise AssertionError('metadata assembly reopened resource content')
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(type(resource), 'open', guarded_open)
    capture = prepare_capture(spec, config, tmp_path / 'result.json', tmp_path / 'capture', resource_hashes=(receipt,))
    assert capture['resources'][0]['sha256'] == json.loads(receipt.read_text())['sha256']


def test_same_size_resource_rewrite_with_restored_mtime_invalidates_receipt(tmp_path):
    config, resource, spec = inputs(tmp_path)
    before = resource.stat()
    receipt = tmp_path / 'resource-hash.json'
    receipt.write_text(json.dumps(hash_artifact(resource)))
    resource.write_bytes(b'x' * before.st_size)
    os.utime(resource, ns=(before.st_atime_ns, before.st_mtime_ns))
    with pytest.raises(ValueError, match='stale'):
        prepare_capture(spec, config, tmp_path / 'result.json', tmp_path / 'capture', resource_hashes=(receipt,))
    assert not (tmp_path / 'capture').exists()


def test_mutation_during_resource_hash_is_rejected(tmp_path, monkeypatch):
    resource = tmp_path / 'changing.txt'
    resource.write_bytes(b'abc')
    original_artifact = _artifact

    def mutate_after_read(path, published_path=None) -> dict:
        result = original_artifact(path, published_path)
        path.write_bytes(b'changed while hashing')
        return result

    monkeypatch.setattr(module, '_artifact', mutate_after_read)
    with pytest.raises(ValueError, match='changed'):
        module.hash_artifact(resource)


@pytest.mark.parametrize(
    'update',
    [
        {'sha256': 'not-a-digest'},
        {'path': 'relative.txt'},
        {'profile': 'unknown/receipt'},
        {'fingerprint': {'size': True, 'mtime_ns': 1, 'ctime_ns': 1}},
    ],
)
def test_invalid_hash_receipt_is_rejected(tmp_path, update):
    config, resource, spec = inputs(tmp_path)
    receipt = tmp_path / 'resource-hash.json'
    receipt.write_text(json.dumps(hash_artifact(resource) | update))
    with pytest.raises(ValueError):
        prepare_capture(spec, config, tmp_path / 'result.json', tmp_path / 'capture', resource_hashes=(receipt,))


def test_missing_receipt_and_conflicting_path_records_are_rejected(tmp_path):
    config, resource, spec = inputs(tmp_path)
    other = tmp_path / 'other.txt'
    other.write_bytes(b'other')
    receipt = tmp_path / 'other-hash.json'
    receipt.write_text(json.dumps(hash_artifact(other)))
    with pytest.raises(ValueError, match='Missing resource hash'):
        prepare_capture(spec, config, tmp_path / 'result.json', tmp_path / 'capture', resource_hashes=(receipt,))
    first = tmp_path / 'first.json'
    second = tmp_path / 'second.json'
    record = hash_artifact(resource)
    first.write_text(json.dumps(record))
    second.write_text(json.dumps(record | {'sha256': 'a' * 64}))
    with pytest.raises(ValueError, match='Conflicting'):
        prepare_capture(spec, config, tmp_path / 'result.json', tmp_path / 'capture', resource_hashes=(first, second))


def test_verify_does_not_require_original_resource_files(tmp_path):
    config, resource, spec = inputs(tmp_path)
    receipt = tmp_path / 'hash.json'
    receipt.write_text(json.dumps(module.hash_artifact(resource)))
    prepare_capture(spec, config, tmp_path / 'result.json', tmp_path / 'capture', resource_hashes=(receipt,))
    resource.unlink()
    verify_capture(tmp_path / 'capture/capture.json')


def test_standalone_preparation_hashes_shared_resource_only_once(tmp_path, monkeypatch):
    config, resource, spec = inputs(tmp_path)
    data = json.loads(spec.read_text())
    data['resources'][0]['parameter'] = 'ensembl_gff'
    spec.write_text(json.dumps(data))
    calls = []

    def record_hash(path, published_path=None) -> dict:
        if path.resolve() == resource.resolve():
            calls.append(path)
        return _artifact(path, published_path)

    monkeypatch.setattr(module, '_artifact', record_hash)
    prepare_capture(
        spec,
        config,
        tmp_path / 'result.json',
        tmp_path / 'capture',
        parameters={'resource_paths': {'ensembl_gff': str(resource)}},
    )
    assert len(calls) == 1


def test_hash_receipt_cli_checks_expected_source_fingerprint(tmp_path, monkeypatch):
    source = tmp_path / 'source.txt'
    source.write_bytes(b'original bytes')
    staged = tmp_path / 'staged.txt'
    staged.symlink_to(source)
    receipt = tmp_path / 'receipt.json'
    monkeypatch.setattr(
        'sys.argv',
        [
            'metadata',
            'hash-artifact',
            '--input',
            str(staged),
            '--source-path',
            str(source),
            '--fingerprint-json',
            json.dumps(module.fingerprint(source)),
            '--output',
            str(receipt),
        ],
    )
    main()
    record = json.loads(receipt.read_text())
    assert record['path'] == str(source)
    assert record['fingerprint']['size'] == len(b'original bytes')
    source.write_bytes(b'changed')
    with pytest.raises(ValueError, match='changed before hashing'):
        hash_artifact(staged, source, record['fingerprint'])


def test_one_changed_source_uses_new_receipt_and_reuses_other_receipt_without_reads(tmp_path, monkeypatch):
    config, first, spec = inputs(tmp_path)
    second = tmp_path / 'second.txt'
    second.write_bytes(b'before')
    data = json.loads(spec.read_text())
    data['resources'].append({'name': 'second', 'version': 'test', 'path': str(second)})
    spec.write_text(json.dumps(data))
    first_receipt = tmp_path / 'first-hash.json'
    first_receipt.write_text(json.dumps(hash_artifact(first)))
    second.write_bytes(b'after')
    second_receipt = tmp_path / 'second-hash.json'
    second_receipt.write_text(json.dumps(hash_artifact(second)))
    original_open = type(first).open

    def guarded_open(path, *args: Any, **kwargs: Any) -> Any:
        if path.resolve() in {first.resolve(), second.resolve()}:
            raise AssertionError('assembly must consume cached hash records')
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(type(first), 'open', guarded_open)
    prepare_capture(
        spec,
        config,
        tmp_path / 'result.json',
        tmp_path / 'capture',
        resource_hashes=(first_receipt, second_receipt, first_receipt),
    )


@pytest.mark.parametrize('completed', [False, True])
def test_pipeline_release_uses_software_version(tmp_path, completed):
    config, _, spec = inputs(tmp_path)
    result = tmp_path / 'result.json'
    if completed:
        result.write_text(json.dumps({'version': '2.3.0', 'metadata': {'version': __version__}}))
    capture = prepare_capture(spec, config, result, tmp_path / 'capture', observed_result=result if completed else None)
    assert 'release_key' not in capture
    assert capture['software_version'] == __version__
