"""MANE identifies transcripts; Ensembl supplies validated exon models and symbols."""

import csv
import json

import pytest

from talos.annotation_scripts.generate_fastspliceai_annotation import generate_annotation
from talos.annotation_scripts.generate_fastspliceai_annotation import main as annotation_main

MANE_COLUMNS = (
    'Ensembl_Gene',
    'Ensembl_nuc',
    'MANE_status',
    'GRCh38_chr',
    'chr_start',
    'chr_end',
    'chr_strand',
    'symbol',
)


def transcript(
    identifier, *, gene='ENSG1', symbol='GENE', chrom='1', strand='+', exons=((10, 100),), cds=(), tags='', version=1
):
    return {
        'identifier': identifier,
        'gene': gene,
        'symbol': symbol,
        'chrom': chrom,
        'strand': strand,
        'exons': exons,
        'cds': cds,
        'tags': tags,
        'version': version,
    }


def mane_row(model, status='MANE Select', **changes: object):
    return {
        'Ensembl_Gene': model['gene'] + '.7',
        'Ensembl_nuc': f'{model["identifier"]}.{model["version"]}',
        'MANE_status': status,
        'GRCh38_chr': model['chrom'],
        'chr_start': min(first for first, _last in model['exons']) + 1,
        'chr_end': max(last for _first, last in model['exons']),
        'chr_strand': model['strand'],
        'symbol': 'OLD_SYMBOL',
    } | changes


def run_annotation(tmp_path, models, mane=(), *, canonical=None, reverse=False):
    records = []
    for model in models:
        chrom, strand = model['chrom'], model['strand']
        first = min(first for first, _last in model['exons'])
        last = max(last for _first, last in model['exons'])
        gene_attributes = f'ID=gene:{model["gene"]};Name={model["symbol"]}'
        canonical_id = canonical.get((model['gene'], chrom, strand)) if isinstance(canonical, dict) else canonical
        if canonical_id:
            gene_attributes += f';canonical_transcript={canonical_id}'
        records.append(f'{chrom}\tEnsembl\tgene\t{first + 1}\t{last}\t.\t{strand}\t.\t{gene_attributes}')
        attributes = f'ID=transcript:{model["identifier"]};Parent=gene:{model["gene"]};{model["tags"]}'
        if model['version'] is not None:
            attributes += f';version={model["version"]}'
        records.append(f'{chrom}\tEnsembl\ttranscript\t{first + 1}\t{last}\t.\t{strand}\t.\t{attributes}')
        for feature, intervals in [('exon', model['exons']), ('CDS', model['cds'])]:
            records.extend(
                f'{chrom}\tEnsembl\t{feature}\t{start + 1}\t{end}\t.\t{strand}\t.'
                f'\tParent=transcript:{model["identifier"]}'
                for start, end in intervals
            )
    gff = tmp_path / 'genes.gff'
    gff.write_text('\n'.join(reversed(records) if reverse else records) + '\n')
    fai = tmp_path / 'ref.fa.fai'
    fai.write_text('chr1\t1000\t0\t1000\t1001\nchr2\t1000\t0\t1000\t1001\n')
    summary = tmp_path / 'mane.tsv'
    with summary.open('w') as stream:
        writer = csv.DictWriter(stream, MANE_COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(mane)
    output = tmp_path / 'annotation.tsv'
    report = tmp_path / 'selection.json'
    generate_annotation(gff, fai, output, mane_summary=summary, selection_report=report)
    with output.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    return rows, json.loads(report.read_text())


def test_mane_select_and_plus_clinical_keep_ensembl_names_and_exon_pairs(tmp_path):
    selected = transcript('ENST1', exons=((10, 30), (80, 100)))
    clinical = transcript('ENST2', exons=((10, 40), (80, 100)))
    other = transcript('ENST3', exons=((0, 150),), tags='tag=Ensembl_canonical')
    rows, report = run_annotation(
        tmp_path,
        [selected, clinical, other],
        [mane_row(selected, GRCh38_chr='NC_000001.11'), mane_row(clinical, 'MANE Plus Clinical')],
    )
    assert [row['#NAME'] for row in rows] == ['GENE', 'GENE']
    assert [row['EXON_END'] for row in rows] == ['30,100,', '40,100,']
    assert [(row['transcript_id'], row['source']) for row in report['selected_transcripts']] == [
        ('ENST1.1', 'MANE Select'),
        ('ENST2.1', 'MANE Plus Clinical'),
    ]
    assert report['unusable_mane'] == []


def test_identical_selected_models_are_scored_once(tmp_path):
    selected, clinical = transcript('ENST1'), transcript('ENST2')
    rows, report = run_annotation(
        tmp_path, [selected, clinical], [mane_row(selected), mane_row(clinical, 'MANE Plus Clinical')]
    )
    assert len(rows) == 1
    assert len(report['selected_transcripts']) == 2
    assert {row['annotation_row'] for row in report['selected_transcripts']} == {1}


@pytest.mark.parametrize(
    ('first', 'second', 'canonical', 'expected_source'),
    [
        ({'tags': 'tag=Ensembl_canonical'}, {}, None, 'Ensembl canonical'),
        ({}, {'tags': 'tag=Ensembl_canonical'}, 'ENST1.1', 'Ensembl canonical'),
        (
            {'tags': 'tag=basic;tag=appris_principal_1', 'cds': ((10, 40),)},
            {'tags': 'tag=appris_principal_2', 'cds': ((10, 90),)},
            None,
            'APPRIS',
        ),
        ({'cds': ((10, 90),)}, {'cds': ((10, 40),)}, None, 'longest CDS'),
        ({'exons': ((0, 100),)}, {}, None, 'longest transcript'),
    ],
)
def test_fallback_policy(tmp_path, first, second, canonical, expected_source):
    rows, report = run_annotation(
        tmp_path, [transcript('ENST1', **first), transcript('ENST2', **second)], canonical=canonical
    )
    assert len(rows) == 1
    assert report['selected_transcripts'][0]['transcript_id'] == 'ENST1.1'
    assert report['selected_transcripts'][0]['source'] == expected_source


@pytest.mark.parametrize(
    ('changes', 'reason'),
    [
        ({'Ensembl_nuc': 'ENST1.2'}, 'version_mismatch'),
        ({'Ensembl_Gene': 'ENSG99.1'}, 'gene_mismatch'),
        ({'GRCh38_chr': '2'}, 'chromosome_mismatch'),
        ({'GRCh38_chr': 'NW_015148969.2'}, 'reference_contig_absent'),
        ({'chr_strand': '-'}, 'strand_mismatch'),
        ({'chr_start': 12}, 'bounds_mismatch'),
        ({'chr_end': 99}, 'bounds_mismatch'),
        ({'Ensembl_nuc': 'ENST99.1'}, 'transcript_absent'),
    ],
)
def test_unusable_mane_is_reported_and_uses_fallback(tmp_path, changes, reason):
    model = transcript('ENST1')
    fallback = transcript('ENST2', exons=((0, 200),), tags='tag=Ensembl_canonical')
    rows, report = run_annotation(tmp_path, [model, fallback], [mane_row(model, **changes)])
    assert len(rows) == 1
    assert report['selected_transcripts'][0]['transcript_id'] == 'ENST2.1'
    assert report['selected_transcripts'][0]['source'] == 'Ensembl canonical'
    assert report['unusable_mane'][0]['reason'] == reason


def test_plus_clinical_does_not_replace_a_missing_select_fallback(tmp_path):
    clinical = transcript('ENST1')
    fallback = transcript('ENST2', exons=((0, 200),), tags='tag=Ensembl_canonical')
    rows, report = run_annotation(tmp_path, [clinical, fallback], [mane_row(clinical, 'MANE Plus Clinical')])
    assert len(rows) == 2
    assert {row['source'] for row in report['selected_transcripts']} == {'MANE Plus Clinical', 'Ensembl canonical'}


def test_shared_symbols_do_not_merge_different_ensembl_genes(tmp_path):
    rows, report = run_annotation(tmp_path, [transcript('ENST1'), transcript('ENST2', gene='ENSG2')])
    assert len(rows) == 2
    assert {row['gene_id'] for row in report['selected_transcripts']} == {'ENSG1', 'ENSG2'}


def test_selection_and_report_are_independent_of_gff_record_order(tmp_path):
    models = [transcript('ENST1'), transcript('ENST2')]
    original = run_annotation(tmp_path, models)
    reordered = run_annotation(tmp_path, models, reverse=True)
    assert reordered == original


def test_cds_length_sums_coding_exons_without_counting_introns(tmp_path):
    short_cds = transcript('ENST1', exons=((0, 10), (900, 910)), cds=((0, 10), (900, 910)))
    long_cds = transcript('ENST2', exons=((20, 100),), cds=((20, 100),))
    _, report = run_annotation(tmp_path, [short_cds, long_cds])
    assert report['selected_transcripts'][0]['transcript_id'] == 'ENST2.1'


def test_gene_symbol_is_taken_from_gene_record_even_if_transcript_has_an_alias(tmp_path):
    model = transcript('ENST1', symbol='CURRENT', tags='gene_name=OLD')
    rows, report = run_annotation(tmp_path, [model], [mane_row(model, symbol='OTHER_ALIAS')])
    assert rows[0]['#NAME'] == report['selected_transcripts'][0]['symbol'] == 'CURRENT'


@pytest.mark.parametrize('canonical', ['ENST1', 'transcript:ENST1.1'])
def test_gene_canonical_identifier_forms(tmp_path, canonical):
    _, report = run_annotation(tmp_path, [transcript('ENST1'), transcript('ENST2')], canonical=canonical)
    assert report['selected_transcripts'][0]['transcript_id'] == 'ENST1.1'
    assert report['selected_transcripts'][0]['source'] == 'Ensembl canonical'


def test_stale_canonical_version_uses_ranked_fallback(tmp_path):
    first = transcript('ENST1')
    second = transcript('ENST2', cds=((10, 90),))
    _, report = run_annotation(tmp_path, [first, second], canonical='ENST1.2')
    assert report['selected_transcripts'][0]['transcript_id'] == 'ENST2.1'
    assert report['selected_transcripts'][0]['source'] == 'longest CDS'


def test_mane_requires_a_known_ensembl_transcript_version(tmp_path):
    model = transcript('ENST1', version=None)
    _, report = run_annotation(tmp_path, [model], [mane_row(model, Ensembl_nuc='ENST1.1')])
    assert report['unusable_mane'][0]['reason'] == 'version_mismatch'
    assert report['selected_transcripts'][0]['transcript_id'] == 'ENST1'


@pytest.mark.parametrize(
    ('changes', 'message'),
    [
        ({'Ensembl_nuc': 'ENST1'}, 'versioned transcript'),
        ({'Ensembl_nuc': 'NM_000001.1'}, 'versioned transcript'),
        ({'Ensembl_Gene': 'SYMBOL'}, 'Ensembl gene ID'),
        ({'MANE_status': 'Unknown'}, 'Unsupported MANE status'),
        ({'chr_start': 0}, 'coordinates'),
        ({'chr_end': 1}, 'coordinates'),
        ({'chr_strand': '.'}, 'coordinates'),
    ],
)
def test_invalid_mane_records_fail_without_writing_annotation(tmp_path, changes, message):
    model = transcript('ENST1')
    with pytest.raises(ValueError, match=message):
        run_annotation(tmp_path, [model], [mane_row(model, **changes)])
    assert not (tmp_path / 'annotation.tsv').exists()


def test_conflicting_mane_records_fail(tmp_path):
    model = transcript('ENST1')
    with pytest.raises(ValueError, match='Conflicting MANE records'):
        run_annotation(tmp_path, [model], [mane_row(model), mane_row(model, 'MANE Plus Clinical')])


def test_identical_mane_records_do_not_add_repeated_inference(tmp_path):
    model = transcript('ENST1')
    rows, report = run_annotation(tmp_path, [model], [mane_row(model)] * 2)
    assert len(rows) == len(report['selected_transcripts']) == 1


def test_multiple_select_transcripts_for_same_gene_locus_fail(tmp_path):
    first, second = transcript('ENST1'), transcript('ENST2')
    with pytest.raises(ValueError, match='Multiple MANE Select'):
        run_annotation(tmp_path, [first, second], [mane_row(first), mane_row(second)])


def test_clinical_model_also_chosen_as_fallback_is_not_duplicated(tmp_path):
    model = transcript('ENST1')
    rows, report = run_annotation(tmp_path, [model], [mane_row(model, 'MANE Plus Clinical')])
    assert len(rows) == len(report['selected_transcripts']) == 1
    assert report['fallback_genes'] == 1


@pytest.mark.parametrize(
    ('exons', 'cds', 'message'),
    [
        (((10, 40), (60, 100)), ((30, 70),), 'not contained in an exon'),
        (((10, 100),), ((10, 50), (40, 90)), 'Overlapping CDS'),
    ],
)
def test_invalid_coding_geometry_is_rejected(tmp_path, exons, cds, message):
    with pytest.raises(ValueError, match=message):
        run_annotation(tmp_path, [transcript('ENST1', exons=exons, cds=cds)])


def test_same_gene_id_on_distinct_reference_loci_is_not_collapsed(tmp_path):
    rows, report = run_annotation(tmp_path, [transcript('ENST1'), transcript('ENST2', chrom='2')])
    assert len(rows) == 2
    assert {row['chromosome'] for row in report['selected_transcripts']} == {'chr1', 'chr2'}


def test_gene_canonical_declarations_remain_specific_to_each_reference_locus(tmp_path):
    models = [transcript('ENST1'), transcript('ENST2', chrom='2')]
    _, report = run_annotation(
        tmp_path,
        models,
        canonical={
            ('ENSG1', '1', '+'): 'ENST1.1',
            ('ENSG1', '2', '+'): 'ENST2.1',
        },
    )
    assert {row['transcript_id'] for row in report['selected_transcripts']} == {'ENST1.1', 'ENST2.1'}
    assert {row['source'] for row in report['selected_transcripts']} == {'Ensembl canonical'}


def test_mane_provenance_changes_even_if_selected_exons_are_identical(tmp_path):
    model = transcript('ENST1')
    rows, original = run_annotation(tmp_path, [model], [mane_row(model)])
    same_rows, renamed = run_annotation(tmp_path, [model], [mane_row(model, symbol='RENAMED')])
    assert rows == same_rows
    assert original['mane_summary_sha256'] != renamed['mane_summary_sha256']


def test_missing_mane_columns_are_rejected(tmp_path):
    run_annotation(tmp_path, [transcript('ENST1')])
    summary = tmp_path / 'mane.tsv'
    summary.write_text('Ensembl_nuc\nENST1.1\n')
    with pytest.raises(ValueError, match='missing required columns'):
        generate_annotation(
            tmp_path / 'genes.gff', tmp_path / 'ref.fa.fai', tmp_path / 'invalid.tsv', mane_summary=summary
        )
    assert not (tmp_path / 'invalid.tsv').exists()


def test_generator_cli_writes_annotation_and_selection_report(tmp_path, monkeypatch):
    model = transcript('ENST1')
    expected_rows, expected_report = run_annotation(tmp_path, [model], [mane_row(model)])
    output, report = tmp_path / 'cli.tsv', tmp_path / 'cli.json'
    monkeypatch.setattr(
        'sys.argv',
        [
            'generate_fastspliceai_annotation',
            '--gff',
            str(tmp_path / 'genes.gff'),
            '--reference-fai',
            str(tmp_path / 'ref.fa.fai'),
            '--mane-summary',
            str(tmp_path / 'mane.tsv'),
            '--output',
            str(output),
            '--selection-report',
            str(report),
        ],
    )
    annotation_main()
    assert json.loads(report.read_text()) == expected_report
    with output.open() as stream:
        assert list(csv.DictReader(stream, delimiter='\t')) == expected_rows
