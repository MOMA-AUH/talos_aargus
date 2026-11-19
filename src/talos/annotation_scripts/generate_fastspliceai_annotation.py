"""Generate FastSpliceAI exon models from the workflow Ensembl GFF and reference.

MANE Select and Plus Clinical identify preferred Ensembl exon models. Genes
without a usable Select model use an Ensembl canonical/APPRIS/length fallback.
GFF coordinates are one-based closed; output starts are zero-based and ends exclusive.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import TextIO
from urllib.parse import unquote

HEADER = '#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n'
POLICY = 'talos/fastspliceai-transcript-selection/1'
MANE_STATUSES = {'MANE Select', 'MANE Plus Clinical'}


def _open_text(path: Path) -> TextIO:
    return gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()


def _contig_key(value: str) -> str:
    bare = value.removeprefix('chr')
    if bare in {'M', 'MT'} or value.split('.', 1)[0] == 'NC_012920':
        return 'MT'
    accession = re.fullmatch(r'NC_0*(\d+)(?:\.\d+)?', value)
    if accession and 1 <= int(accession[1]) <= 24:
        number = int(accession[1])
        return str(number) if number <= 22 else {23: 'X', 24: 'Y'}[number]
    return bare


def _reference(reference_fai: Path) -> dict[str, int]:
    lengths = {}
    with reference_fai.open() as stream:
        for line in stream:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2 or not fields[0] or int(fields[1]) < 1 or fields[0] in lengths:
                raise ValueError('Reference FAI must contain unique contigs and positive lengths')
            lengths[fields[0]] = int(fields[1])
    if not lengths:
        raise ValueError('Reference FAI is empty')
    return lengths


def _reference_contig(value: str, lengths: dict[str, int]) -> str | None:
    if value in lengths:
        return value
    candidates = [name for name in lengths if _contig_key(name) == _contig_key(value)]
    if len(candidates) > 1:
        raise ValueError(f'Ambiguous reference contig alias: {value}')
    return candidates[0] if candidates else None


def _identifier(value: str) -> tuple[str, str | None] | None:
    matched = re.fullmatch(r'(?:transcript:|rna-)?(ENST\d+)(?:\.(\d+))?', unquote(value))
    return (matched[1], matched[2]) if matched else None


def _gene_key(value: str) -> str:
    return unquote(value).removeprefix('gene:').split('.')[0]


def _coordinates(start: str, end: str, strand: str, length: int) -> tuple[int, int]:
    first, last = int(start) - 1, int(end)
    if first < 0 or last <= first or last > length or strand not in {'+', '-'}:
        raise ValueError('Invalid transcript/exon coordinates or strand')
    return first, last


@dataclass(frozen=True)
class Gene:
    symbol: str
    canonical: tuple[str, str | None] | None


@dataclass(frozen=True)
class Transcript:
    chrom: str
    strand: str
    start: int
    end: int
    version: str | None
    parent: str
    symbol: str
    gene_id: str
    tags: frozenset[str]


@dataclass(frozen=True)
class AnnotationModel:
    identifier: str
    transcript: Transcript
    gene_id: str
    symbol: str
    version: str | None
    exons: tuple[tuple[int, int], ...]
    cds_length: int

    @property
    def group(self) -> tuple[str, str, str]:
        return self.gene_id, self.transcript.chrom, self.transcript.strand

    @property
    def fields(self) -> tuple[str, ...]:
        model = self.transcript
        return (
            self.symbol,
            model.chrom,
            model.strand,
            str(model.start),
            str(model.end),
            ''.join(f'{first},' for first, _last in self.exons),
            ''.join(f'{last},' for _first, last in self.exons),
        )


@dataclass(frozen=True)
class ManeTranscript:
    gene_id: str
    identifier: str
    version: str
    status: str
    contig: str
    chrom: str | None
    start: int
    end: int
    strand: str


def _version(identifier_version: str | None, explicit_version: str | None) -> str | None:
    if identifier_version is not None and explicit_version is not None and identifier_version != explicit_version:
        raise ValueError('Conflicting transcript version in ID and version attribute')
    return identifier_version or explicit_version


def _models(gff: Path, lengths: dict[str, int]) -> tuple[dict, dict, dict, dict]:
    genes = {}
    transcripts = {}
    exons = defaultdict(set)
    coding = defaultdict(set)
    resolved_contigs = {}
    with _open_text(gff) as stream:
        for line in stream:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError('GFF records must contain nine fields')
            raw_chrom, _source, feature, start, end, _score, strand, _phase, raw_attributes = fields
            if raw_chrom not in resolved_contigs:
                resolved_contigs[raw_chrom] = _reference_contig(raw_chrom, lengths)
            chrom = resolved_contigs[raw_chrom]
            if chrom is None:
                continue
            attributes = dict(item.split('=', 1) for item in raw_attributes.split(';') if '=' in item)
            identifier = attributes.get('ID', '')
            if feature in {'gene', 'pseudogene'} or identifier.startswith('gene:'):
                gene_id = _gene_key(identifier or attributes.get('gene_id', ''))
                key = (gene_id, chrom, strand)
                symbol = unquote(attributes.get('Name', attributes.get('gene_name', gene_id)))
                gene = Gene(symbol, _identifier(attributes.get('canonical_transcript', '')))
                if key in genes and genes[key] != gene:
                    raise ValueError(f'Conflicting gene records: {key}')
                genes[key] = gene
            elif feature in {'exon', 'CDS'}:
                first, last = _coordinates(start, end, strand, lengths[chrom])
                for parent in attributes.get('Parent', '').split(','):
                    identity = _identifier(parent)
                    if identity:
                        # An exon version attribute versions the exon, not its Parent transcript.
                        target = exons if feature == 'exon' else coding
                        target[identity[0]].add((chrom, strand, first, last, identity[1]))
            elif (identity := _identifier(identifier)) is not None:
                first, last = _coordinates(start, end, strand, lengths[chrom])
                model = Transcript(
                    chrom,
                    strand,
                    first,
                    last,
                    _version(identity[1], attributes.get('version')),
                    _gene_key(attributes.get('Parent', '')),
                    unquote(attributes.get('gene_name', '')),
                    unquote(attributes.get('gene_id', '')),
                    frozenset(
                        unquote(tag)
                        for item in raw_attributes.split(';')
                        if item.startswith('tag=')
                        for tag in item.split('=', 1)[1].split(',')
                    ),
                )
                if identity[0] in transcripts:
                    raise ValueError(f'Duplicate transcript ID: {identity[0]}')
                transcripts[identity[0]] = model
    return genes, transcripts, exons, coding


def _validated_models(genes: dict, transcripts: dict, exons: dict, coding: dict) -> dict[str, AnnotationModel]:
    models = {}
    for identifier, records in exons.items():
        if identifier not in transcripts:
            raise ValueError(f'Exon references missing transcript: {identifier}')
        model = transcripts[identifier]
        coding_records = coding.get(identifier, set())
        observed_versions = {record[4] for record in records | coding_records if record[4] is not None}
        if model.version is not None:
            observed_versions.add(model.version)
        if len(observed_versions) > 1:
            raise ValueError(f'Conflicting exon/transcript versions: {identifier}')
        coordinates = {(record[2], record[3]) for record in records}
        for chrom, strand, first, last, _record_version in records | coding_records:
            if chrom != model.chrom:
                raise ValueError(f'Exon contig differs from transcript: {identifier}')
            if strand != model.strand:
                raise ValueError(f'Exon strand differs from transcript: {identifier}')
            if first < model.start or last > model.end:
                raise ValueError(f'Exon coordinates exceed transcript: {identifier}')
        ordered = sorted(coordinates)
        if (ordered[0][0], ordered[-1][1]) != (model.start, model.end):
            raise ValueError(f'Exon span differs from transcript: {identifier}')
        if any(left[1] > right[0] for left, right in pairwise(ordered)):
            raise ValueError(f'Overlapping exons for transcript: {identifier}')
        gene = genes.get((model.parent or _gene_key(model.gene_id), model.chrom, model.strand))
        symbol = (gene.symbol if gene else '') or model.symbol or model.gene_id
        if not symbol:
            raise ValueError(f'Transcript references missing gene: {identifier}')
        if any(character in symbol for character in '\t\n\r|,'):
            raise ValueError(f'Invalid gene symbol: {identifier}')
        coding_intervals = sorted({(record[2], record[3]) for record in coding_records})
        if any(not any(start <= first < last <= end for start, end in ordered) for first, last in coding_intervals):
            raise ValueError(f'CDS coordinates are not contained in an exon: {identifier}')
        if any(left[1] > right[0] for left, right in pairwise(coding_intervals)):
            raise ValueError(f'Overlapping CDS intervals: {identifier}')
        models[identifier] = AnnotationModel(
            identifier,
            model,
            model.parent or _gene_key(model.gene_id) or symbol,
            symbol,
            next(iter(observed_versions), None),
            tuple(ordered),
            sum(last - first for first, last in coding_intervals),
        )
    if not models:
        raise ValueError('GFF contains no exon-bearing Ensembl transcripts on the reference')
    return models


def _mane_records(summary: Path, lengths: dict[str, int]) -> list[ManeTranscript]:
    required = {'Ensembl_Gene', 'Ensembl_nuc', 'MANE_status', 'GRCh38_chr', 'chr_start', 'chr_end', 'chr_strand'}
    records = {}
    with _open_text(summary) as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('MANE summary is missing required columns')
        for row in reader:
            identity = _identifier(row['Ensembl_nuc'])
            gene_id = _gene_key(row['Ensembl_Gene'])
            if identity is None or identity[1] is None or not re.fullmatch(r'ENSG\d+', gene_id):
                raise ValueError('MANE summary requires an Ensembl gene ID and a versioned transcript ID')
            if row['MANE_status'] not in MANE_STATUSES:
                raise ValueError(f'Unsupported MANE status: {row["MANE_status"]}')
            start, end = _coordinates(row['chr_start'], row['chr_end'], row['chr_strand'], int(row['chr_end']))
            record = ManeTranscript(
                gene_id,
                identity[0],
                identity[1],
                row['MANE_status'],
                row['GRCh38_chr'],
                _reference_contig(row['GRCh38_chr'], lengths),
                start,
                end,
                row['chr_strand'],
            )
            if identity[0] in records and records[identity[0]] != record:
                raise ValueError(f'Conflicting MANE records: {identity[0]}')
            records[identity[0]] = record
    return sorted(records.values(), key=lambda record: (record.gene_id, record.status, record.identifier))


def _mane_mismatch(record: ManeTranscript, candidate: AnnotationModel | None) -> str | None:
    if record.chrom is None:
        return 'reference_contig_absent'
    if candidate is None:
        return 'transcript_absent'
    checks = (
        (record.gene_id, candidate.gene_id, 'gene_mismatch'),
        (record.chrom, candidate.transcript.chrom, 'chromosome_mismatch'),
        (record.strand, candidate.transcript.strand, 'strand_mismatch'),
        (record.version, candidate.version, 'version_mismatch'),
        ((record.start, record.end), (candidate.transcript.start, candidate.transcript.end), 'bounds_mismatch'),
    )
    return next((reason for expected, observed, reason in checks if expected != observed), None)


def _appris_rank(model: AnnotationModel) -> int:
    ranks = []
    for tag in model.transcript.tags:
        match = re.fullmatch(r'appris_(principal|alternative|alternate)(?:_(\d+))?', tag)
        if match:
            ranks.append((0 if match[1] == 'principal' else 10) + int(match[2] or 1))
    return min(ranks, default=100)


def _canonical_rank(model: AnnotationModel, genes: dict[tuple[str, str, str], Gene]) -> int:
    gene = genes.get(model.group)
    if (
        gene
        and gene.canonical
        and gene.canonical[0] == model.identifier
        and (gene.canonical[1] is None or gene.canonical[1] == model.version)
    ):
        return 0
    return 1 if 'Ensembl_canonical' in model.transcript.tags else 2


def _fallback(
    candidates: list[AnnotationModel], genes: dict[tuple[str, str, str], Gene]
) -> tuple[AnnotationModel, str]:
    chosen = min(
        candidates,
        key=lambda model: (
            _canonical_rank(model, genes),
            not bool(model.cds_length),
            _appris_rank(model) if model.cds_length else 100,
            -model.cds_length,
            -sum(last - first for first, last in model.exons),
            model.identifier,
        ),
    )
    if _canonical_rank(chosen, genes) < 2:
        source = 'Ensembl canonical'
    elif chosen.cds_length and _appris_rank(chosen) < 100:
        source = 'APPRIS'
    else:
        source = 'longest CDS' if chosen.cds_length else 'longest transcript'
    return chosen, source


def _select_models(
    models: dict[str, AnnotationModel],
    genes: dict[tuple[str, str, str], Gene],
    mane: list[ManeTranscript],
) -> tuple[list, list, int]:
    groups = defaultdict(list)
    preferred = defaultdict(list)
    unusable = []
    for model in models.values():
        groups[model.group].append(model)
    for record in mane:
        candidate = models.get(record.identifier)
        reason = _mane_mismatch(record, candidate)
        if reason:
            unusable.append(
                {
                    'gene_id': record.gene_id,
                    'transcript_id': f'{record.identifier}.{record.version}',
                    'status': record.status,
                    'contig': record.contig,
                    'reason': reason,
                }
            )
        elif candidate is not None:
            preferred[candidate.group].append((candidate, record.status))
    selected = []
    fallbacks = 0
    for group, candidates in sorted(groups.items()):
        choices = preferred[group]
        select_count = sum(source == 'MANE Select' for _model, source in choices)
        if select_count > 1:
            raise ValueError(f'Multiple MANE Select transcripts for gene/locus: {group}')
        if not select_count:
            fallbacks += 1
            fallback, source = _fallback(candidates, genes)
            if all(model.identifier != fallback.identifier for model, _source in choices):
                choices.append((fallback, source))
        selected.extend(choices)
    return selected, unusable, fallbacks


def generate_annotation(
    gff: Path,
    reference_fai: Path,
    output: Path,
    *,
    mane_summary: Path,
    selection_report: Path | None = None,
) -> None:
    """Write MANE-led Ensembl models and deterministic selection provenance."""
    lengths = _reference(reference_fai)
    genes, transcripts, exons, coding = _models(gff, lengths)
    models = _validated_models(genes, transcripts, exons, coding)
    selected, unusable, fallbacks = _select_models(models, genes, _mane_records(mane_summary, lengths))
    rank = {chrom: index for index, chrom in enumerate(lengths)}
    selected.sort(
        key=lambda item: (
            rank[item[0].transcript.chrom],
            item[0].transcript.start,
            item[0].transcript.end,
            item[0].symbol,
            item[0].identifier,
            item[1],
        )
    )
    rows = {}
    choices = []
    for model, source in selected:
        key = (model.group, model.fields)
        if key not in rows:
            rows[key] = len(rows) + 1
        choices.append(
            {
                'gene_id': model.gene_id,
                'symbol': model.symbol,
                'transcript_id': model.identifier + (f'.{model.version}' if model.version is not None else ''),
                'chromosome': model.transcript.chrom,
                'strand': model.transcript.strand,
                'source': source,
                'annotation_row': rows[key],
            }
        )
    with output.open('w') as stream:
        stream.write(HEADER)
        for _group, fields in rows:
            stream.write('\t'.join(fields) + '\n')
    if selection_report is not None:
        report = {
            'policy': POLICY,
            'mane_summary_sha256': hashlib.sha256(mane_summary.read_bytes()).hexdigest(),
            'annotation_rows': len(rows),
            'fallback_genes': fallbacks,
            'selected_transcripts': choices,
            'unusable_mane': unusable,
        }
        selection_report.write_text(json.dumps(report, sort_keys=True, indent=2) + '\n')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff', type=Path, required=True)
    parser.add_argument('--reference-fai', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mane-summary', type=Path, required=True)
    parser.add_argument('--selection-report', type=Path, required=True)
    args = parser.parse_args()
    generate_annotation(
        args.gff,
        args.reference_fai,
        args.output,
        mane_summary=args.mane_summary,
        selection_report=args.selection_report,
    )


if __name__ == '__main__':
    main()
