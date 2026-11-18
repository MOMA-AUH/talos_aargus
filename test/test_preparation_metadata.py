"""Preparation publishes hashes and normal outputs without extra copying tasks."""

import csv
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import NoReturn

import pytest

from talos import workflow_metadata
from talos.version import __version__

ROOT = Path(__file__).resolve().parents[1]


def _prepare_project(project: Path) -> tuple[Path, Path]:
    shutil.copytree(ROOT / 'nextflow', project / 'nextflow')
    (project / 'preparation.nf').write_text(
        (ROOT / 'preparation.nf')
        .read_text()
        .replace('nextflow.enable.dsl=2', 'nextflow.enable.dsl=2\nnextflow.preview.output=true', 1)
    )
    outputs = {
        'EncodeAlphaMissense': ['alphamissense.vcf.gz', 'alphamissense.zip'],
        'EncodeMitimpact': ['mitimpact.vcf.gz', 'mitimpact.zip'],
        'EncodeMitotip': ['mitotip.vcf.gz', 'mitotip.zip'],
        'EncodeNapogee': ['napogee.vcf.gz', 'napogee.zip'],
        'ResummariseRawSubmissions': [
            'clinvarbitration_${timestamp}.vcf.bgz',
            'clinvarbitration_${timestamp}.vcf.bgz.tbi',
        ],
        'AnnotateClinvarWithBcftools': ['clinvarbitration.annotated.tsv'],
        'MakeClinvarbitrationPm5': ['clinvarbitration_${timestamp}.pm5.tsv'],
        'CreateRoiFromGff3': ['GRCh38.bed', 'GRCh38_merged.bed', 'GRCh38_symbol_to_ensg.json'],
        'ParseManeIntoJson': ['mane.json'],
        'DownloadPanelApp': ['panelapp_${timestamp}.json'],
    }
    for name, files in outputs.items():
        module = project / 'nextflow/modules/prep' / name / 'main.nf'
        header = module.read_text().split('script:', 1)[0]
        commands = [f'printf \'%s\\n\' "prepared {filename}" > "{filename}"' for filename in files]
        module.write_text(header + 'script:\n    """\n    ' + '\n    '.join(commands) + '\n    """\n}\n')
    raw, prepared = project / 'raw', project / 'prepared'
    raw.mkdir()
    month = datetime.now().strftime('%Y-%m')  # noqa: DTZ005
    for name in [
        'genes.gff',
        'gnomad.zip',
        'am.tsv',
        'mitimpact.tsv',
        'mitotip.tsv',
        'napogee.tsv',
        'mane.tsv',
        'hp.obo',
        'gen2phen.txt',
        'phenio.db',
        f'submissions_{month}.txt.gz',
        f'variants_{month}.txt.gz',
    ]:
        (raw / name).write_text('resource-a\n')
    (raw / 'ref.fa').write_text('>chr1\nACGT\n')
    (raw / 'ref.dict').write_text('@SQ\tSN:chr1\tLN:4\tM5:f1f8f4bf413b16ad135722aa4591043e\n')
    values = {
        'large_files': raw,
        'processed_annotations': prepared,
        'outdir': prepared,
        'ref_genome': raw / 'ref.fa',
        'ensembl_gff': raw / 'genes.gff',
        'gnomad_zip': raw / 'gnomad.zip',
        'alphamissense_tsv': raw / 'am.tsv',
        'alphamissense_zip': prepared / 'alphamissense.zip',
        'alphamissense_version': 'alpha-test',
        'mitimpact_tsv': raw / 'mitimpact.tsv',
        'mitimpact_zip': prepared / 'mitimpact.zip',
        'mitotip_tsv': raw / 'mitotip.tsv',
        'mitotip_zip': prepared / 'mitotip.zip',
        'mitotip_version': 'mitotip-test',
        'napogee_tsv': raw / 'napogee.tsv',
        'napogee_zip': prepared / 'napogee.zip',
        'mane': raw / 'mane.tsv',
        'mane_json': prepared / 'mane.json',
        'ensembl_bed': prepared / 'GRCh38.bed',
        'ensembl_merged_bed': prepared / 'GRCh38_merged.bed',
        'ensembl_symbol_lookup': prepared / 'GRCh38_symbol_to_ensg.json',
        'hpo': raw / 'hp.obo',
        'gen2phen': raw / 'gen2phen.txt',
        'phenio_db': raw / 'phenio.db',
        'phenio_version': 'phenio-test',
    }
    parameters = '\n'.join(f"    {key} = '{value}'" for key, value in values.items())
    (project / 'nextflow.config').write_text(f"""includeConfig '{ROOT}/nextflow.config'
params {{
{parameters}
}}
docker.enabled = false
singularity.enabled = false
apptainer.enabled = false
outputDir = params.outdir
workflow.output.mode = 'symlink'
process {{
    executor = 'local'
    cpus = 1
    memory = 1.GB
    withName: '.*' {{
        cpus = 1
        memory = 1.GB
    }}
    beforeScript = 'export PYTHONPATH="{ROOT}/src"; export PATH="{Path(sys.executable).parent}:$PATH"'
}}
""")
    return raw, prepared


def _run_preparation(project: Path, trace: str) -> list[dict[str, str]]:
    executable = os.environ.get('NEXTFLOW_BIN') or shutil.which('nextflow')
    if executable is None:
        pytest.skip('Set NEXTFLOW_BIN to run the preparation integration test')
    result = subprocess.run(  # noqa: S603 - explicitly selected local Nextflow executable
        [executable, 'run', 'preparation.nf', '-ansi-log', 'false', '-with-trace', trace],
        cwd=project,
        env={**os.environ, 'NXF_OFFLINE': 'true', 'NXF_SYNTAX_PARSER': 'v2'},
        capture_output=True,
        text=True,
        timeout=90,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    with (project / trace).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def test_preparation_publishes_hashes_without_copying_jobs_or_resource_locations(tmp_path, monkeypatch):
    raw, prepared = _prepare_project(tmp_path)
    trace = _run_preparation(tmp_path, 'cold.tsv')
    assert any('HashResource' in row['name'] for row in trace)
    assert not any('StorePreparedResource' in row['name'] for row in trace)
    catalog_path = prepared / 'resources.metadata/resources.json'
    catalog = json.loads(catalog_path.read_text())
    assert catalog['profile'] == 'talos/prepared-resources/2'
    assert {resource['name'] for resource in catalog['resources']} == {
        'MANE',
        'Ensembl',
        'gnomAD',
        'AlphaMissense',
        'ClinVar',
        'PM5',
        'HPO',
        'HPO-G2P',
        'Phenio',
        'MitImpact',
        'MitoTIP',
        'nAPOGEE',
    }
    assert {path.name for path in catalog_path.parent.iterdir()} == {'resources.json'}
    assert all(set(resource) == {'name', 'version', 'parameter', 'sha256'} for resource in catalog['resources'])
    ensembl = next(resource for resource in catalog['resources'] if resource['name'] == 'Ensembl')
    assert ensembl['sha256'] == workflow_metadata.hash_artifact(raw / 'genes.gff')['sha256']
    assert all(len(resource['sha256']) == 64 for resource in catalog['resources'])
    shutil.rmtree(tmp_path / 'work')
    assert (prepared / 'mane.json').is_file() and not (prepared / 'mane.json').is_symlink()
    config = tmp_path / 'analysis.toml'
    config.write_text('[RunHailFiltering]\nspliceai = 0.5\n')
    specification = tmp_path / 'analysis-spec.json'
    specification.write_text(json.dumps({'assembly': catalog['assembly'], 'resources': catalog['resources']}))
    report = tmp_path / 'result.json'
    report.write_text(json.dumps({'version': '2.3.0', 'metadata': {'version': __version__}}))

    def reject_hash(_self: object, _path: Path) -> NoReturn:
        raise AssertionError('Analysis must use prepared hashes without reopening resources')

    monkeypatch.setattr(workflow_metadata.ResourceHasher, 'artifact', reject_hash)
    capture = workflow_metadata.prepare_capture(
        specification,
        config,
        report,
        tmp_path / 'analysis.metadata',
        observed_result=report,
    )
    assert all(set(resource) == {'name', 'version', 'sha256'} for resource in capture['resources'])
    workflow_metadata.verify_capture(tmp_path / 'analysis.metadata/capture.json')


def test_preparation_accepts_an_empty_resource_override(tmp_path):
    _, prepared = _prepare_project(tmp_path)
    override = tmp_path / 'override.json'
    override.write_text(json.dumps({'assembly': 'GRCh38', 'resources': []}))
    with (tmp_path / 'nextflow.config').open('a') as handle:
        handle.write(f"\nparams.aargus_metadata_spec = '{override}'\n")
    rows = _run_preparation(tmp_path, 'empty.tsv')
    assert not any('HashResource' in row['name'] for row in rows)
    catalog = json.loads((prepared / 'resources.metadata/resources.json').read_text())
    assert catalog == {'profile': 'talos/prepared-resources/2', 'assembly': 'GRCh38', 'resources': []}


def test_preparation_override_hashes_only_selected_resources(tmp_path):
    raw, prepared = _prepare_project(tmp_path)
    extra = raw / 'source.txt'
    extra.write_text('historical input\n')
    specification = {
        'assembly': 'GRCh38',
        'resources': [
            {'name': 'Ensembl', 'version': '116', 'path': str(raw / 'genes.gff'), 'parameter': 'ensembl_gff'},
            {'name': 'additional resource', 'version': 'test', 'path': str(extra)},
        ],
    }
    override = tmp_path / 'override.json'
    override.write_text(json.dumps(specification))
    with (tmp_path / 'nextflow.config').open('a') as handle:
        handle.write(f"\nparams.aargus_metadata_spec = '{override}'\n")
    rows = _run_preparation(tmp_path, 'override.tsv')
    catalog = json.loads((prepared / 'resources.metadata/resources.json').read_text())
    assert {entry['name'] for entry in catalog['resources']} == {'Ensembl', 'additional resource'}
    assert len([row for row in rows if 'HashResource' in row['name']]) == 2
    captured = next(entry for entry in catalog['resources'] if entry['name'] == 'additional resource')
    assert captured['sha256'] == workflow_metadata.hash_artifact(extra)['sha256']
    assert 'path' not in captured


def test_resource_hashes_resume_independently_and_detect_restored_mtime(tmp_path):
    """Unchanged resources reuse receipts; same-size edits still invalidate one task."""
    nextflow = os.environ.get('NEXTFLOW_BIN') or shutil.which('nextflow')
    if nextflow is None:
        pytest.skip('Set NEXTFLOW_BIN to run the Nextflow integration test')
    first, second = tmp_path / 'first.txt', tmp_path / 'second.txt'
    first.write_text('AAAA')
    second.write_text('BBBB')
    (tmp_path / 'main.nf').write_text(f"""nextflow.enable.dsl=2
include {{ HashResource }} from '{ROOT}/nextflow/modules/prep/HashResource/main'
workflow {{
    def inputs = ['first.txt', 'second.txt'].collect {{ name ->
        def source = file(name).toRealPath()
        def attributes = java.nio.file.Files.readAttributes(source, 'unix:size,lastModifiedTime,ctime')
        def fingerprint = [size: attributes.size,
            mtime_ns: attributes.lastModifiedTime.to(java.util.concurrent.TimeUnit.NANOSECONDS),
            ctime_ns: attributes.ctime.to(java.util.concurrent.TimeUnit.NANOSECONDS)]
        tuple(source.toString(), source, fingerprint)
    }}
    HashResource(channel.fromList(inputs))
}}
""")
    (tmp_path / 'nextflow.config').write_text(f"""process {{
    cpus = 1
    memory = 1.GB
    beforeScript = 'export PYTHONPATH="{ROOT}/src"; export PATH="{Path(sys.executable).parent}:$PATH"'
}}
""")

    def run_hashes(label, *, resume) -> dict[str, str]:
        trace = tmp_path / f'{label}.tsv'
        command = [nextflow, 'run', 'main.nf', '-ansi-log', 'false', '-with-trace', str(trace)]
        if resume:
            command.append('-resume')
        result = subprocess.run(  # noqa: S603 - explicitly selected local Nextflow executable
            command,
            cwd=tmp_path,
            env={**os.environ, 'NXF_OFFLINE': 'true', 'NXF_SYNTAX_PARSER': 'v2'},
            capture_output=True,
            text=True,
            timeout=90,
            check=False,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        with trace.open() as handle:
            records = list(csv.DictReader(handle, delimiter='\t'))
        return {record['name']: record['status'] for record in records}

    expected_first = f'HashResource ({first})'
    expected_second = f'HashResource ({second})'
    assert run_hashes('cold', resume=False) == {expected_first: 'COMPLETED', expected_second: 'COMPLETED'}
    assert run_hashes('warm', resume=True) == {expected_first: 'CACHED', expected_second: 'CACHED'}
    before = first.stat()
    first.write_text('CCCC')
    os.utime(first, ns=(before.st_atime_ns, before.st_mtime_ns))
    assert first.stat().st_size == before.st_size
    assert first.stat().st_mtime_ns == before.st_mtime_ns
    assert first.stat().st_ctime_ns != before.st_ctime_ns
    assert run_hashes('changed', resume=True) == {expected_first: 'COMPLETED', expected_second: 'CACHED'}
