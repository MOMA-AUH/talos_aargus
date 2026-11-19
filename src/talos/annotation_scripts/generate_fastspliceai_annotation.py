"""Generate FastSpliceAI exon models from the workflow Ensembl GFF and reference.

Each exon-bearing Ensembl transcript is retained separately. GFF coordinates
are one-based closed; the generated starts are zero-based and ends exclusive.
"""

from __future__ import annotations

import argparse
import gzip
import re
from collections import defaultdict
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import TextIO
from urllib.parse import unquote

HEADER = '#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n'


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
class Transcript:
    chrom: str
    strand: str
    start: int
    end: int
    version: str | None
    parent: str
    symbol: str
    gene_id: str


def _version(identifier_version: str | None, explicit_version: str | None) -> str | None:
    if identifier_version is not None and explicit_version is not None and identifier_version != explicit_version:
        raise ValueError('Conflicting transcript version in ID and version attribute')
    return identifier_version or explicit_version


def _models(gff: Path, lengths: dict[str, int]) -> tuple[dict, dict, dict]:
    genes = {}
    transcripts = {}
    exons = defaultdict(set)
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
                key = _gene_key(identifier or attributes.get('gene_id', ''))
                symbol = unquote(attributes.get('Name', attributes.get('gene_name', key)))
                if key in genes and genes[key] != symbol:
                    raise ValueError(f'Conflicting gene records: {key}')
                genes[key] = symbol
            elif feature == 'exon':
                first, last = _coordinates(start, end, strand, lengths[chrom])
                for parent in attributes.get('Parent', '').split(','):
                    identity = _identifier(parent)
                    if identity:
                        # An exon version attribute versions the exon, not its Parent transcript.
                        exons[identity[0]].add((chrom, strand, first, last, identity[1]))
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
                )
                if identity[0] in transcripts:
                    raise ValueError(f'Duplicate transcript ID: {identity[0]}')
                transcripts[identity[0]] = model
    return genes, transcripts, exons


def _validated_rows(genes: dict, transcripts: dict, exons: dict, lengths: dict) -> list:
    rows = []
    rank = {chrom: index for index, chrom in enumerate(lengths)}
    for identifier, records in exons.items():
        if identifier not in transcripts:
            raise ValueError(f'Exon references missing transcript: {identifier}')
        model = transcripts[identifier]
        observed_versions = {record[4] for record in records if record[4] is not None}
        if model.version is not None:
            observed_versions.add(model.version)
        if len(observed_versions) > 1:
            raise ValueError(f'Conflicting exon/transcript versions: {identifier}')
        coordinates = set()
        for chrom, strand, first, last, version in records:
            if chrom != model.chrom:
                raise ValueError(f'Exon contig differs from transcript: {identifier}')
            if strand != model.strand:
                raise ValueError(f'Exon strand differs from transcript: {identifier}')
            if version is not None and model.version is not None and version != model.version:
                raise ValueError(f'Exon version differs from transcript: {identifier}')
            if first < model.start or last > model.end:
                raise ValueError(f'Exon coordinates exceed transcript: {identifier}')
            coordinates.add((first, last))
        ordered = sorted(coordinates)
        if (ordered[0][0], ordered[-1][1]) != (model.start, model.end):
            raise ValueError(f'Exon span differs from transcript: {identifier}')
        if any(left[1] > right[0] for left, right in pairwise(ordered)):
            raise ValueError(f'Overlapping exons for transcript: {identifier}')
        symbol = model.symbol or genes.get(model.parent) or model.gene_id
        if not symbol:
            raise ValueError(f'Transcript references missing gene: {identifier}')
        if any(character in symbol for character in '\t\n\r|,'):
            raise ValueError(f'Invalid gene symbol: {identifier}')
        fields = [
            symbol,
            model.chrom,
            model.strand,
            str(model.start),
            str(model.end),
            ''.join(f'{first},' for first, _last in ordered),
            ''.join(f'{last},' for _first, last in ordered),
        ]
        rows.append(((rank[model.chrom], model.start, model.end, symbol, identifier), fields))
    if not rows:
        raise ValueError('GFF contains no exon-bearing Ensembl transcripts on the reference')
    return sorted(rows)


def generate_annotation(gff: Path, reference_fai: Path, output: Path) -> None:
    """Write every reference-compatible Ensembl transcript without merging isoforms."""
    lengths = _reference(reference_fai)
    genes, transcripts, exons = _models(gff, lengths)
    rows = _validated_rows(genes, transcripts, exons, lengths)
    with output.open('w') as stream:
        stream.write(HEADER)
        for _order, fields in rows:
            stream.write('\t'.join(fields) + '\n')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff', type=Path, required=True)
    parser.add_argument('--reference-fai', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    generate_annotation(args.gff, args.reference_fai, args.output)


if __name__ == '__main__':
    main()
