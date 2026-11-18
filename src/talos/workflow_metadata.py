"""Generate and verify the preparation capture consumed by AARgus."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from importlib.metadata import version as distribution_version
from pathlib import Path
from typing import Any

from talos.version import __version__

RESOURCE_HASH_PROFILE = 'talos/resource-hash/1'


def _text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{label} must be a nonempty string')
    return value


def _artifact(path: Path, published_path: Path | None = None) -> dict:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'path': str((published_path or path).resolve()), 'sha256': digest.hexdigest()}


def fingerprint(path: Path) -> dict[str, int]:
    status = path.stat()
    return {'size': status.st_size, 'mtime_ns': status.st_mtime_ns, 'ctime_ns': status.st_ctime_ns}


def hash_artifact(input_path: Path, source_path: Path | None = None, expected_fingerprint: dict | None = None) -> dict:
    """Hash one stable source; record its original identity through staged symlinks."""
    source = (source_path or input_path).resolve()
    before = fingerprint(source)
    if expected_fingerprint is not None and before != expected_fingerprint:
        raise ValueError(f'Resource changed before hashing: {source}')
    if fingerprint(input_path) != before:
        raise ValueError(f'Staged resource differs from source fingerprint: {source}')
    artifact = _artifact(input_path, source)
    if fingerprint(source) != before or fingerprint(input_path) != before:
        raise ValueError(f'Resource changed during hashing: {source}')
    return artifact | {'profile': RESOURCE_HASH_PROFILE, 'fingerprint': before}


class ResourceHasher:
    """Use validated task receipts, or hash each source once within this preparation."""

    def __init__(self, receipts: tuple[Path, ...]) -> None:
        self.required_receipts = bool(receipts)
        self.records: dict[Path, dict] = {}
        for path in receipts:
            record = _read_json(path)
            source = self._validate(record)
            if source in self.records and self.records[source] != record:
                raise ValueError(f'Conflicting resource hash receipts: {source}')
            self.records[source] = record

    @staticmethod
    def _validate(record: dict) -> Path:
        if set(record) != {'profile', 'path', 'sha256', 'fingerprint'} or record['profile'] != RESOURCE_HASH_PROFILE:
            raise ValueError('Invalid resource hash receipt schema')
        source = Path(_text(record['path'], 'resource hash path'))
        if not source.is_absolute() or source != source.resolve():
            raise ValueError('Resource hash receipt path must be canonical and absolute')
        if not isinstance(record['sha256'], str) or not re.fullmatch(r'[0-9a-f]{64}', record['sha256']):
            raise ValueError('Invalid resource hash SHA256')
        stamp = record['fingerprint']
        if (
            not isinstance(stamp, dict)
            or set(stamp) != {'size', 'mtime_ns', 'ctime_ns'}
            or any(type(value) is not int or value < 0 for value in stamp.values())
        ):
            raise ValueError('Invalid resource hash fingerprint')
        if fingerprint(source) != stamp:
            raise ValueError(f'Resource hash receipt is stale: {source}')
        return source

    def artifact(self, path: Path) -> dict:
        source = path.resolve()
        record = self.records.get(source)
        if record is None:
            if self.required_receipts:
                raise ValueError(f'Missing resource hash receipt: {source}')
            record = hash_artifact(source)
            self.records[source] = record
        self._validate(record)
        return {'path': str(source), 'sha256': record['sha256']}

    def validate_all(self) -> None:
        for record in self.records.values():
            self._validate(record)


def _write_json(path: Path, data: dict) -> None:
    with path.open('x', encoding='utf-8') as handle:
        json.dump(data, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def _read_json(path: Path) -> dict:
    data = json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(data, dict):
        raise ValueError(f'{path} must contain a JSON object')
    return data


def _sha256(value: Any) -> str:
    if not isinstance(value, str) or not re.fullmatch(r'[0-9a-f]{64}', value):
        raise ValueError('Resource SHA256 must contain 64 lowercase hexadecimal characters')
    return value


def _resources(spec: dict, base: Path, parameters: dict, hasher: ResourceHasher) -> list[dict]:
    captures = []
    names = set()
    runtime_paths = parameters.get('resource_paths', {})
    if not isinstance(runtime_paths, dict):
        raise ValueError('resource_paths must map actual workflow input parameter names to paths')
    bound_parameters = set()
    for resource in spec.get('resources', []):
        name = _text(resource.get('name'), 'resource name')
        if name in names:
            raise ValueError(f'Duplicate resource name: {name}')
        names.add(name)
        version = _text(resource.get('version'), 'resource version')
        prepared = 'sha256' in resource
        digest = (
            _sha256(resource['sha256'])
            if prepared
            else hasher.artifact(base / _text(resource.get('path'), 'resource path'))['sha256']
        )
        capture = {
            'name': name,
            'version': version,
            'sha256': digest,
        }
        captures.append(capture)
        parameter = resource.get('parameter')
        if parameter == 'fastspliceai_annotation':
            capture['version'] = distribution_version('spliceai')
        if parameter is not None:
            _text(parameter, 'resource parameter')
            if not prepared:
                actual_path = runtime_paths.get(parameter)
                if not isinstance(actual_path, str) or not actual_path:
                    raise ValueError(f'Missing actual workflow resource path for {parameter}')
                observed = hasher.artifact(Path(actual_path))
                if observed['sha256'] != digest:
                    raise ValueError(f'Resource checksum differs from actual workflow input: {parameter}')
            bound_parameters.add(parameter)
    missing = set(runtime_paths) - bound_parameters
    if missing:
        raise ValueError(f'Actual workflow resources need versioned declarations and bindings: {sorted(missing)}')
    return captures


def _check_result(path: Path, version: str) -> str:
    result = _read_json(path)
    if result.get('version') != '2.3.0':
        raise ValueError('AARgus requires Talos result schema 2.3.0')
    metadata = result.get('metadata')
    if not isinstance(metadata, dict) or metadata.get('version') != version:
        raise ValueError('Result software version differs from captured software version')
    return metadata['version']


def prepare_capture(
    specification: Path,
    configuration: Path,
    result_path: Path,
    output_directory: Path,
    *,
    published_directory: Path | None = None,
    parameters: dict | None = None,
    observed_result: Path | None = None,
    resource_hashes: tuple[Path, ...] = (),
) -> dict:
    """Write one immutable capture bundle, optionally for later Nextflow publishing."""
    spec = _read_json(specification)
    software_version = _check_result(observed_result, __version__) if observed_result is not None else __version__
    assembly = _text(spec.get('assembly'), 'assembly')
    actual = {} if parameters is None else parameters
    if not isinstance(actual, dict):
        raise ValueError('parameters must be a JSON object')
    hasher = ResourceHasher(resource_hashes)
    captures = _resources(spec, specification.resolve().parent, actual, hasher)
    published = (published_directory or output_directory).resolve()
    config_bytes = configuration.read_bytes()
    hasher.validate_all()
    output_directory.mkdir(parents=True, exist_ok=False)
    copied_config = output_directory / 'configuration.toml'
    copied_config.write_bytes(config_bytes)
    config = _artifact(copied_config, published / copied_config.name)
    capture = {
        'software_version': software_version,
        'assembly': assembly,
        'resources': captures,
        'configuration': config,
        'result_path': str(result_path.resolve()),
    }
    _write_json(output_directory / 'capture.json', capture)
    return capture


def _verify_artifact(artifact: dict) -> None:
    path = Path(artifact['path'])
    if _artifact(path)['sha256'] != artifact['sha256']:
        raise ValueError(f'SHA256 mismatch: {path}')


def verify_capture(capture_path: Path, result: Path | None = None) -> dict:
    """Verify configuration and recorded digests, then optionally hash a result."""
    capture = _read_json(capture_path)
    if result is not None and result.resolve() != Path(capture['result_path']).resolve():
        raise ValueError('Result path differs from captured planned result path')
    _verify_artifact(capture['configuration'])
    for resource in capture['resources']:
        _sha256(resource['sha256'])
    if result is not None:
        _check_result(result, capture['software_version'])
        return {
            'preparation': _artifact(capture_path),
            'result': _artifact(result),
        }
    return {'preparation': _artifact(capture_path)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    prepare = commands.add_parser('prepare')
    prepare.add_argument('--spec', type=Path, required=True)
    prepare.add_argument('--configuration', type=Path, required=True)
    prepare.add_argument('--result-path', type=Path, required=True)
    prepare.add_argument('--output-directory', type=Path, required=True)
    prepare.add_argument('--published-directory', type=Path)
    prepare.add_argument('--parameters-json', type=json.loads, default={})
    prepare.add_argument('--verify-result', type=Path)
    prepare.add_argument('--resource-hash', type=Path, action='append', default=[])
    hashing = commands.add_parser('hash-artifact')
    hashing.add_argument('--input', type=Path, required=True)
    hashing.add_argument('--source-path', type=Path)
    hashing.add_argument('--fingerprint-json', type=json.loads)
    hashing.add_argument('--output', type=Path, required=True)
    verify = commands.add_parser('verify')
    verify.add_argument('--capture', type=Path, required=True)
    verify.add_argument('--result', type=Path)
    verify.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'prepare':
        prepare_capture(
            args.spec,
            args.configuration,
            args.result_path,
            args.output_directory,
            published_directory=args.published_directory,
            parameters=args.parameters_json,
            observed_result=args.verify_result,
            resource_hashes=tuple(args.resource_hash),
        )
    elif args.command == 'hash-artifact':
        _write_json(args.output, hash_artifact(args.input, args.source_path, args.fingerprint_json))
    else:
        _write_json(args.output, verify_capture(args.capture, args.result))


if __name__ == '__main__':
    main()
