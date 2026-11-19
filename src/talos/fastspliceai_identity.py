"""Bind a FastSpliceAI cache namespace to resources, settings, and scorer bytes."""

from __future__ import annotations

import argparse
import hashlib
import json
from importlib.metadata import PackageNotFoundError, distribution, version
from pathlib import Path

PROFILE = 'talos/fastspliceai-cache/2'
WEIGHT_SUFFIXES = {'.h5', '.hdf5', '.pt', '.pth', '.ckpt', '.safetensors', '.onnx'}
SCORING_DEPENDENCIES = ('torch', 'numpy', 'pysam', 'pyfaidx', 'h5py')


def _fingerprint(path: Path) -> tuple[int, int, int, int]:
    status = path.stat()
    return status.st_size, status.st_mtime_ns, status.st_ctime_ns, status.st_ino


def _hash_file(path: Path) -> str:
    before = _fingerprint(path)
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    if before != _fingerprint(path):
        raise ValueError(f'FastSpliceAI resource changed during hashing: {path}')
    return digest.hexdigest()


def scorer_identity() -> dict:
    """Identify installed source and packaged weights without importing the model."""
    package = distribution('spliceai')
    files = package.files or []
    sources = sorted(path for path in files if str(path).startswith('spliceai/') and path.suffix == '.py')
    weights = sorted(path for path in files if str(path).startswith('spliceai/') and path.suffix in WEIGHT_SUFFIXES)
    if not weights:
        raise ValueError('Cannot identify packaged FastSpliceAI model weights')
    if not sources:
        raise ValueError('Cannot identify packaged FastSpliceAI source files')
    paths = {str(path): Path(package.locate_file(path)) for path in [*sources, *weights]}
    stamps = {name: _fingerprint(path) for name, path in paths.items()}
    hashes = {name: _hash_file(path) for name, path in paths.items()}
    if any(_fingerprint(paths[name]) != stamp for name, stamp in stamps.items()):
        raise ValueError('FastSpliceAI source or model weights changed during identity capture')
    dependencies = {}
    for name in SCORING_DEPENDENCIES:
        try:
            dependencies[name] = version(name)
        except PackageNotFoundError:
            dependencies[name] = None
    identity = {'version': package.version, 'files': hashes, 'dependencies': dependencies}
    direct_url = package.read_text('direct_url.json')
    if direct_url:
        provenance = json.loads(direct_url)
        revision = provenance.get('vcs_info', {}).get('commit_id')
        archive_hashes = provenance.get('archive_info', {}).get('hashes')
        if revision:
            identity['source_revision'] = revision
        if archive_hashes:
            identity['source_archive'] = archive_hashes
    return identity


def create_identity(
    gff: Path,
    annotation: Path,
    reference: Path,
    reference_fai: Path,
    selection: Path,
    *,
    distance: int,
    mask: int,
    scorer: dict | None = None,
) -> dict:
    """Return a reproducible namespace key after checking stable resource reads."""
    if type(distance) is not int or not 0 <= distance <= 4999:
        raise ValueError('FastSpliceAI distance must be an integer between 0 and 4999')
    if type(mask) is not int or mask not in {0, 1}:
        raise ValueError('FastSpliceAI mask must be 0 or 1')
    paths = {
        'gff': gff,
        'annotation': annotation,
        'reference': reference,
        'reference_fai': reference_fai,
        'selection': selection,
    }
    stamps = {name: _fingerprint(path) for name, path in paths.items()}
    resources = {name: _hash_file(path) for name, path in paths.items()}
    payload = {
        'profile': PROFILE,
        'resources': resources,
        'settings': {'distance': distance, 'mask': mask},
        'scorer': scorer_identity() if scorer is None else scorer,
    }
    if any(_fingerprint(paths[name]) != stamp for name, stamp in stamps.items()):
        raise ValueError('FastSpliceAI resource changed during identity capture')
    encoded = json.dumps(payload, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    return payload | {'key': hashlib.sha256(encoded).hexdigest()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('gff', 'annotation', 'reference', 'reference-fai', 'selection', 'output'):
        parser.add_argument(f'--{name}', type=Path, required=True)
    parser.add_argument('--distance', type=int, required=True)
    parser.add_argument('--mask', type=int, required=True)
    args = parser.parse_args()
    identity = create_identity(
        args.gff,
        args.annotation,
        args.reference,
        args.reference_fai,
        args.selection,
        distance=args.distance,
        mask=args.mask,
    )
    args.output.write_text(json.dumps(identity, sort_keys=True, indent=2) + '\n')


if __name__ == '__main__':
    main()
