"""FastSpliceAI models use every Ensembl transcript present in the reference."""

import gzip

import pytest

from talos.annotation_scripts.generate_fastspliceai_annotation import generate_annotation


def annotation_inputs(tmp_path):
    gff = tmp_path / 'genes.gff.gz'
    fai = tmp_path / 'ref.fa.fai'
    fai.write_text('chr2\t1000\t0\t1000\t1001\nchr10\t1000\t0\t1000\t1001\nchrM\t1000\t0\t1000\t1001\n')
    records = [
        '10\tEnsembl\tgene\t11\t30\t.\t-\t.\tID=gene:ENSG10;Name=TEN',
        '10\tEnsembl\ttranscript\t11\t30\t.\t-\t.\tID=transcript:ENST10;Parent=gene:ENSG10;version=2',
        '10\tEnsembl\texon\t21\t30\t.\t-\t.\tParent=transcript:ENST10',
        '10\tEnsembl\texon\t11\t15\t.\t-\t.\tParent=transcript:ENST10',
        '2\tEnsembl\tgene\t11\t200\t.\t+\t.\tID=gene:ENSG2;Name=TWO',
        '2\tEnsembl\tmRNA\t101\t200\t.\t+\t.\tID=transcript:ENST2.1;Parent=gene:ENSG2',
        '2\tEnsembl\texon\t101\t120\t.\t+\t.\tParent=transcript:ENST2.1',
        '2\tEnsembl\texon\t181\t200\t.\t+\t.\tParent=transcript:ENST2.1',
        '2\tEnsembl\tlnc_RNA\t11\t20\t.\t+\t.\tID=transcript:ENST1;Parent=gene:ENSG2',
        '2\tEnsembl\texon\t11\t20\t.\t+\t.\tParent=transcript:ENST1',
        'MT\tEnsembl\tgene\t5\t10\t.\t+\t.\tID=gene:ENSG25',
        'MT\tEnsembl\ttRNA\t5\t10\t.\t+\t.\tID=transcript:ENST25;Parent=gene:ENSG25',
        'MT\tEnsembl\texon\t5\t10\t.\t+\t.\tParent=transcript:ENST25',
        'PATCH\tEnsembl\ttranscript\t1\t10\t.\t+\t.\tID=transcript:ENST99;gene_name=PATCH',
        'PATCH\tEnsembl\texon\t1\t10\t.\t+\t.\tParent=transcript:ENST99',
    ]
    with gzip.open(gff, 'wt') as stream:
        stream.write('##gff-version 3\n' + '\n'.join(records) + '\n')
    return gff, fai, tmp_path / 'annotation.tsv'


def test_all_isoforms_noncoding_and_mt_use_reference_order_and_paired_exons(tmp_path):
    gff, fai, output = annotation_inputs(tmp_path)
    generate_annotation(gff, fai, output)
    assert output.read_text().splitlines() == [
        '#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END',
        'TWO\tchr2\t+\t10\t20\t10,\t20,',
        'TWO\tchr2\t+\t100\t200\t100,180,\t120,200,',
        'TEN\tchr10\t-\t10\t30\t10,20,\t15,30,',
        'ENSG25\tchrM\t+\t4\t10\t4,\t10,',
    ]


@pytest.mark.parametrize(
    ('old', 'new', 'message'),
    [
        ('ENST10;Parent=gene:ENSG10;version=2', 'ENST10.3;Parent=gene:ENSG10;version=2', 'version'),
        ('exon\t21\t30', 'exon\t21\t29', 'span'),
        ('exon\t21\t30', 'exon\t0\t30', 'coordinates'),
        ('10\tEnsembl\texon', '2\tEnsembl\texon', 'contig'),
        ('21\t30\t.\t-\t.', '21\t30\t.\t+\t.', 'strand'),
        ('Parent=gene:ENSG10', 'Parent=gene:MISSING', 'gene'),
        ('Parent=transcript:ENST10', 'Parent=transcript:ENST404', 'transcript'),
    ],
)
def test_invalid_models_fail_without_output(tmp_path, old, new, message):
    gff, fai, output = annotation_inputs(tmp_path)
    with gzip.open(gff, 'rt') as stream:
        content = stream.read()
    with gzip.open(gff, 'wt') as stream:
        stream.write(content.replace(old, new))
    with pytest.raises(ValueError, match=message):
        generate_annotation(gff, fai, output)
    assert not output.exists()


def test_transcript_without_exons_is_not_invented(tmp_path):
    gff, fai, output = annotation_inputs(tmp_path)
    with gzip.open(gff, 'rt') as stream:
        content = stream.read()
    with gzip.open(gff, 'wt') as stream:
        stream.write(content + '2\tEnsembl\ttranscript\t201\t250\t.\t+\t.\tID=transcript:ENST5;gene_name=EMPTY\n')
    generate_annotation(gff, fai, output)
    assert 'EMPTY' not in output.read_text()


def test_plain_gff_with_gene_name_and_multi_parent_exons(tmp_path):
    gff = tmp_path / 'genes.gff'
    gff.write_text(
        '1\tEnsembl\tsnRNA\t1\t20\t.\t+\t.\tID=transcript:ENST1;gene_name=SN%2BRNA\n'
        '1\tEnsembl\tsnRNA\t1\t20\t.\t+\t.\tID=transcript:ENST2;gene_name=SN%2BRNA\n'
        '1\tEnsembl\texon\t1\t20\t.\t+\t.\tParent=transcript:ENST1,transcript:ENST2\n'
    )
    fai = tmp_path / 'ref.fai'
    fai.write_text('1\t100\t0\t100\t101\n')
    output = tmp_path / 'annotation.tsv'
    generate_annotation(gff, fai, output)
    assert output.read_text().splitlines()[1:] == ['SN+RNA\t1\t+\t0\t20\t0,\t20,'] * 2


def test_exon_version_attribute_does_not_version_its_parent(tmp_path):
    gff, fai, output = annotation_inputs(tmp_path)
    with gzip.open(gff, 'rt') as stream:
        content = stream.read()
    with gzip.open(gff, 'wt') as stream:
        stream.write(content.replace('Parent=transcript:ENST10', 'Parent=transcript:ENST10;version=99'))
    generate_annotation(gff, fai, output)
    assert 'TEN\tchr10' in output.read_text()


def test_unversioned_transcript_cannot_merge_different_parent_versions(tmp_path):
    gff = tmp_path / 'genes.gff'
    gff.write_text(
        '1\tEnsembl\ttranscript\t1\t20\t.\t+\t.\tID=transcript:ENST1;gene_name=GENE\n'
        '1\tEnsembl\texon\t1\t5\t.\t+\t.\tParent=transcript:ENST1.1\n'
        '1\tEnsembl\texon\t11\t20\t.\t+\t.\tParent=transcript:ENST1.2\n'
    )
    fai = tmp_path / 'ref.fai'
    fai.write_text('1\t100\t0\t100\t101\n')
    output = tmp_path / 'annotation.tsv'
    with pytest.raises(ValueError, match='version'):
        generate_annotation(gff, fai, output)
    assert not output.exists()


def test_single_parent_version_is_accepted_for_unversioned_transcript(tmp_path):
    gff = tmp_path / 'genes.gff'
    gff.write_text(
        '1\tEnsembl\ttranscript\t1\t20\t.\t+\t.\tID=transcript:ENST1;gene_name=GENE\n'
        '1\tEnsembl\texon\t1\t5\t.\t+\t.\tParent=transcript:ENST1.1\n'
        '1\tEnsembl\texon\t11\t20\t.\t+\t.\tParent=transcript:ENST1.1\n'
    )
    fai = tmp_path / 'ref.fai'
    fai.write_text('1\t100\t0\t100\t101\n')
    output = tmp_path / 'annotation.tsv'
    generate_annotation(gff, fai, output)
    assert output.read_text().splitlines()[1] == 'GENE\t1\t+\t0\t20\t0,10,\t5,20,'
