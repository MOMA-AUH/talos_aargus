"""Public option and metadata behavior through the real TALOS Nextflow graph."""

import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pytest

from talos.version import __version__

ROOT = Path(__file__).resolve().parents[1]


def test_metadata_spec_override_is_optional_and_named_for_aargus():
    schema = json.loads((ROOT / 'nextflow_schema.json').read_text())
    parameters = schema['definitions']['default']
    assert parameters['properties']['aargus_metadata_spec']['type'] == 'string'
    assert 'aargus_metadata_spec' not in parameters['required']


def _prepare_workflow(tmp_path, annotation, mito) -> tuple[Path, dict[str, str]]:
    shutil.copytree(ROOT / 'nextflow', tmp_path / 'nextflow')
    # Retain real process contracts while replacing expensive biological work.
    outputs = {
        'StartupChecks': '${cohort}_checked',
        'UnifiedPanelAppParser': '${cohort}_panelapp.json',
        'RunHailFiltering': '${cohort}_small_variants_labelled.vcf.bgz ${cohort}_small_variants_labelled.vcf.bgz.tbi',
        'AnnotateMitoVcf': '${cohort}_mito_labelled.vcf.bgz',
        'ValidateMOI': '${cohort}_results_${timestamp}.json',
        'HPOFlagging': '${cohort}_full_results_${timestamp}.json',
        'CreateTalosHTML': '${cohort}_report_${timestamp}.html',
    }
    for name, output in outputs.items():
        module = tmp_path / 'nextflow/modules/talos' / name / 'main.nf'
        header = module.read_text().split('script:', 1)[0]
        command = f'touch {output}'
        if name == 'HPOFlagging':
            report = json.dumps({'version': '2.3.0', 'metadata': {'version': __version__}})
            command = f"printf '%s\\n' '{report}' > {output}"
        module.write_text(header + 'script:\n    """\n    ' + command + '\n    """\n}\n')
    (tmp_path / 'main.nf').write_text("""nextflow.enable.dsl=2
nextflow.preview.output=true
include { TALOS } from './nextflow/talos'
workflow {
    main:
    def base = file(params.resources)
    def inputs = channel.of(tuple(
        'cohort', [base.resolve('input.vcf')], base.resolve('pedigree.ped'),
        base.resolve('configuration.toml'), base.resolve('NO_HISTORY'), base.resolve('NO_FILE'),
        base.resolve('NO_SEQR_FILE'), base.resolve(params.mito ? 'mito.vcf' : 'NO_MITO')))
    TALOS(channel.value(base.resolve('mane.json')), channel.value(base.resolve('genes.gff')),
        channel.value(base.resolve('ref.fa')), inputs, params.annotation,
        channel.value(base.resolve('spliceai_transcripts.tsv')))
    publish:
        metadata = TALOS.out.metadata
}
output {
    metadata { path { cohort, bundle -> "${cohort}/outputs" } }
}
""")
    resources = tmp_path / 'resources'
    resources.mkdir()
    month = datetime.now().strftime('%Y-%m')  # noqa: DTZ005
    bindings = {
        'mane_json': 'mane.json',
        'clinvar': f'clinvarbitration_{month}.vcf.bgz',
        'pm5': f'clinvarbitration_{month}.pm5.tsv',
        'hpo': 'hp.obo',
        'gen2phen': 'gen2phen.txt',
        'phenio_db': 'phenio.db',
    }
    if annotation:
        bindings.update(
            ensembl_gff='genes.gff',
            gnomad_zip='gnomad.zip',
            alphamissense_zip='alphamissense.zip',
            fastspliceai_annotation='spliceai_transcripts.tsv',
        )
    if mito:
        bindings.update(
            ensembl_gff='genes.gff',
            mitimpact='mitimpact.zip',
            mitotip='mitotip.zip',
            napogee='napogee.zip',
        )
    for name in {
        *bindings.values(),
        'ref.fa',
        'genes.gff',
        'GRCh38.bed',
        'GRCh38_merged.bed',
        'input.vcf',
        'pedigree.ped',
        'NO_HISTORY',
        'NO_FILE',
        'NO_SEQR_FILE',
        'NO_MITO',
        'mito.vcf',
        f'clinvarbitration_{month}.vcf.bgz.tbi',
        f'panelapp_{month}.json',
    }:
        (resources / name).write_text(name + '\n')
    (resources / 'ref.dict').write_text('@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\tM5:f1f8f4bf413b16ad135722aa4591043e\n')
    (resources / 'configuration.toml').write_text('[RunHailFiltering]\nspliceai = 0.5\n')
    return resources, bindings


@pytest.mark.parametrize(
    ('annotation', 'mito', 'metadata_source'),
    [
        (False, False, 'override'),
        (True, True, 'override'),
        (False, False, 'config'),
        (True, True, 'config'),
        (False, False, 'missing-catalog'),
        (False, False, 'invalid-hash'),
        (False, False, 'missing-hash'),
    ],
)
def test_talos_publishes_verifiable_capture_with_actual_resources(  # noqa: PLR0915
    tmp_path, annotation, mito, metadata_source
):
    """Stub analysis tools; execute the real metadata/hash processes and channels."""
    nextflow = os.environ.get('NEXTFLOW_BIN') or shutil.which('nextflow')
    if nextflow is None:
        pytest.skip('Set NEXTFLOW_BIN to run the Nextflow integration test')
    resources, bindings = _prepare_workflow(tmp_path, annotation, mito)
    specification = {
        'assembly': 'GRCh38',
        'resources': [
            {
                'name': key,
                'version': 'test',
                'path': str(resources / name),
                'parameter': key,
            }
            for key, name in bindings.items()
        ],
    }
    (tmp_path / 'specification.json').write_text(json.dumps(specification))
    versions = {
        'mane_version': 'test-mane',
        'ensembl_version': 'test-ensembl',
        'gnomad_version': 'test-gnomad',
        'alphamissense_version': 'test-alphamissense',
        'hpo_version': 'test-hpo',
        'gen2phen_version': 'test-gen2phen',
        'phenio_version': 'test-phenio',
        'mitimpact_version': '3.1.3',
        'mitotip_version': 'test-mitotip',
        'napogee_version': '1.0.0',
    }
    metadata_root = resources / 'resources.metadata'
    metadata_root.mkdir()
    month = datetime.now().strftime('%Y-%m')  # noqa: DTZ005
    version_by_binding = {
        'mane_json': versions['mane_version'],
        'ensembl_gff': versions['ensembl_version'],
        'gnomad_zip': versions['gnomad_version'],
        'alphamissense_zip': versions['alphamissense_version'],
        'clinvar': month,
        'pm5': month,
        'hpo': versions['hpo_version'],
        'gen2phen': versions['gen2phen_version'],
        'phenio_db': versions['phenio_version'],
        'mitimpact': versions['mitimpact_version'],
        'mitotip': versions['mitotip_version'],
        'napogee': versions['napogee_version'],
        'fastspliceai_annotation': versions['mane_version'],
    }
    captured_resources = [
        {key: value for key, value in resource.items() if key != 'path'}
        | {
            'version': version_by_binding[resource['parameter']],
            'sha256': hashlib.sha256(Path(resource['path']).read_bytes()).hexdigest(),
        }
        for resource in specification['resources']
    ]
    catalog = {
        'profile': 'talos/prepared-resources/2',
        'assembly': 'GRCh38',
        'resources': [
            resource for resource in captured_resources if resource['parameter'] != 'fastspliceai_annotation'
        ],
    }
    if metadata_source == 'override':
        # A stale declared annotation hash must not replace the table generated
        # for this analysis; ordinary resources still use prepared identities.
        specification['resources'] = [
            resource | ({'sha256': '0' * 64} if resource['parameter'] == 'fastspliceai_annotation' else {})
            for resource in captured_resources
        ]
        (tmp_path / 'specification.json').write_text(json.dumps(specification))
    catalog_path = metadata_root / 'resources.json'
    catalog_path.write_text(json.dumps(catalog))
    if metadata_source == 'missing-catalog':
        catalog_path.unlink()
    if metadata_source == 'missing-hash':
        catalog['resources'][0].pop('sha256')
        catalog_path.write_text(json.dumps(catalog))
    if metadata_source == 'invalid-hash':
        catalog['resources'][0]['sha256'] = 'invalid'
        catalog_path.write_text(json.dumps(catalog))
    version_options = '\n'.join(f"        {key} = 'changed-after-preparation'" for key in versions)
    metadata_option = (
        f"aargus_metadata_spec = '{tmp_path}/specification.json'\n"
        if metadata_source == 'override'
        else (f"aargus_metadata_spec = null\n        ref_genome_assembly = 'GRCh38'\n{version_options}")
    )
    executables = tmp_path / 'bin'
    executables.mkdir()
    distribution = executables / 'spliceai-2.0.0.dist-info'
    distribution.mkdir()
    (distribution / 'METADATA').write_text('Metadata-Version: 2.1\nName: spliceai\nVersion: 2.0.0\n')
    for name in ['bgzip', 'tabix', 'bcftools', 'echtvar', 'fastspliceai']:
        executable = executables / name
        executable.write_text('#!/bin/sh\nprintf test-version\n')
        executable.chmod(0o755)
    task_path = f'{executables}:{Path(sys.executable).parent}'
    config = f"""params {{
        annotation = {str(annotation).lower()}
        mito = {str(mito).lower()}
        resources = '{resources}'
        processed_annotations = '{resources}'
        hpo = '{resources}/hp.obo'
        gen2phen = '{resources}/gen2phen.txt'
        phenio_db = '{resources}/phenio.db'
        ref_genome = '{resources}/ref.fa'
        ensembl_gff = '{resources}/genes.gff'
        gnomad_zip = '{resources}/gnomad.zip'
        vcf_split_n = 1000
        input_vcf_extension = 'vcf.gz'
        {metadata_option}
        outdir = '{tmp_path}/published'
    }}
    docker.enabled = false
    singularity.enabled = false
    apptainer.enabled = false
    outputDir = params.outdir
    workflow.output.mode = 'copy'
    process {{
        executor = 'local'
        cpus = 1
        memory = 1.GB
        beforeScript = 'export PYTHONPATH="{ROOT}/src:{executables}"; export PATH="{task_path}:$PATH"'
    }}
    """
    (tmp_path / 'nextflow.config').write_text(f"includeConfig '{ROOT}/nextflow.config'\n" + config)
    result = subprocess.run(  # noqa: S603 - explicitly selected local Nextflow executable
        [nextflow, 'run', 'main.nf', '-ansi-log', 'false', '-with-trace', str(tmp_path / 'trace.tsv')],
        cwd=tmp_path,
        env={**os.environ, 'NXF_OFFLINE': 'true', 'NXF_SYNTAX_PARSER': 'v2'},
        capture_output=True,
        text=True,
        timeout=90,
        check=False,
    )
    if metadata_source in {'missing-catalog', 'invalid-hash', 'missing-hash'}:
        assert result.returncode != 0
        assert 'run preparation' in result.stdout + result.stderr
        assert 'Submitted process > TALOS:HashResource' not in result.stdout
        return
    assert result.returncode == 0, result.stdout + result.stderr
    bundles = list((tmp_path / 'published/cohort/outputs').glob('cohort.metadata_*'))
    assert len(bundles) == 1
    bundle = bundles[0]
    assert not bundle.is_symlink()
    assert all(not path.is_symlink() for path in bundle.rglob('*'))
    assert {path.name for path in bundle.iterdir()} == {'capture.json', 'configuration.toml'}
    capture = json.loads((bundle / 'capture.json').read_text())
    assert set(capture) == {
        'software_version',
        'assembly',
        'resources',
        'configuration',
        'result_path',
    }
    assert capture['software_version'] == __version__
    assert capture['assembly'] == 'GRCh38'
    assert capture['result_path'].startswith(str(tmp_path / 'published/cohort/outputs/cohort_full_results_'))
    result_date = Path(capture['result_path']).stem.rsplit('_', 1)[1]
    assert bundle.name == f'cohort.metadata_{result_date}'
    assert capture['configuration']['path'] == str(bundle / 'configuration.toml')
    expected_names = set(bindings)
    if annotation and metadata_source == 'config':
        expected_names.remove('fastspliceai_annotation')
        expected_names.add('FastSpliceAI')
    assert {resource['name'] for resource in capture['resources']} == expected_names
    for resource in capture['resources']:
        assert set(resource) == {'name', 'version', 'sha256'}
        expected_version = (
            '2.0.0'
            if resource['name'] in {'fastspliceai_annotation', 'FastSpliceAI'}
            else version_by_binding[resource['name']]
        )
        assert resource['version'] == expected_version
        binding = 'fastspliceai_annotation' if resource['name'] == 'FastSpliceAI' else resource['name']
        source = resources / bindings[binding]
        assert resource['sha256'] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert set(capture['configuration']) == {'path', 'sha256'}
    assert (bundle / 'configuration.toml').read_bytes() == (resources / 'configuration.toml').read_bytes()
    with (tmp_path / 'trace.tsv').open() as handle:
        tasks = list(csv.DictReader(handle, delimiter='\t'))
    assert not any('HashResource' in row['name'] for row in tasks)
    if annotation:
        task = next(row for row in tasks if row['name'].startswith('TALOS:WorkflowMetadata '))
        prefix, suffix = task['hash'].split('/')
        task_directory = next((tmp_path / 'work' / prefix).glob(f'{suffix}*'))
        assert (task_directory / 'annotation/spliceai_transcripts.tsv').read_bytes() == (
            resources / 'spliceai_transcripts.tsv'
        ).read_bytes()
    from talos.workflow_metadata import verify_capture  # noqa: PLC0415

    verify_capture(bundle / 'capture.json')
