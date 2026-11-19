# Talos Aargus

[![Docs](https://img.shields.io/badge/docs-populationgenomics.github.io%2Ftalos-blue)](https://populationgenomics.github.io/talos/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
![Test](https://github.com/populationgenomics/automated-interpretation-pipeline/actions/workflows/test.yaml/badge.svg)

This repository is an Aarhus University Hospital extension of the Talos workflow.

It keeps the Talos v11 TSV-driven analysis model, but adds an AUH-oriented runtime
and preflight layer:

- Docker execution with FastSpliceAI included in the Talos image
- Apptainer support for cluster deployments
- standalone `preflight/` workflow for singleton, duo, and trio VCF preparation
- In-house resource wiring for common variants, variant artefacts, and splicing variants

## Part 1. AUH Extended Pipeline

### Overview

The AUH extension is intended to support:

- local or cluster execution with explicit runtime paths
- preflight preparation of family VCFs before Talos
- Slurm-friendly bundle generation for full end-to-end runs
- optional Apptainer execution with the tested Talos runtime baked into the image

The main components are:

- `preparation.nf`
  Creates processed Talos annotations such as ClinVar, PM5, PanelApp, and encoded annotation resources
- `main.nf`
  Runs annotation plus Talos prioritisation from a TSV input file
- `talos_only.nf`
  Runs Talos only on previously annotated inputs
- `preflight/main.nf`
  Standalone AUH preflight workflow for manifest-driven family preparation

### Runtime Layout

Expected local structure:

- `large_files/`
  Downloaded static reference resources
- `processed_annotations/`
  Prepared Talos annotation outputs from `preparation.nf`
- `resources/`
  Small versioned resources shipped with this checkout, including the
  FastSpliceAI-compatible MANE v1.5 transcript table
- `/faststorage/project/reanalyses_auh/env/talos2_env/bin`
  Shared AUH runtime path used by the preflight layer and local example scripts,
  including `nextflow`
- `../../../resources/`
  Shared AUH sibling resource folder for common-dbSNP and legacy precomputed
  SpliceAI resources

### Install

Requirements:

- [Nextflow](https://www.nextflow.io/docs/latest/install.html)
- Java
- Docker or Apptainer for main-workflow tasks
- Python 3.10 or 3.11
- `bcftools`, `tabix`, `bgzip`, `samtools`, `wget`

For a shared conda-based setup:

```bash
bash scripts/install_talos_env.sh
```

### Prepare Static Resources

Download the larger static resource set:

```bash
bash large_files/gather_files.sh
```

On Slurm systems, a batch wrapper is also provided:

```bash
sbatch large_files/download_resources_slurm.sh
```

Then prepare Talos processed annotations:

```bash
nextflow -c nextflow.config run preparation.nf \
  --processed_annotations processed_annotations \
  --large_files large_files
```

### Workflow Smoke Test

The repository includes a smoke-test TSV for the Talos v11 input model:

```bash
nextflow -c nextflow.config run main.nf \
  --input_tsv nextflow/inputs/test.tsv \
  --processed_annotations processed_annotations \
  -output-dir nextflow_output
```

If you want to test the Apptainer path instead:

```bash
bash build_apptainer_frontend.sh
bash run_nextflow_example_apptainer.sh
```

### Cached FastSpliceAI

The annotation workflow always runs FastSpliceAI through the bcftools VCF
cache. There are no precomputed, uncached, or disabled SpliceAI modes.
The default cache root is `<processed_annotations>/spliceai_cache`; no mode
or cache path needs to be supplied. `fastspliceai_cache` optionally overrides
this **directory**, not a VCF filename.

All input splits from all cohorts are normalized first, including existing
`shards` inputs and whole VCFs split by the workflow. FastSpliceAI first gathers
the splits within each cohort into one indexed, site-only VCF, then combines
those cohort VCFs into the distinct global allele set. Both gathers stream
through bcftools, avoiding temporary site-only files for every split.
The global sites are intersected with the existing cache. Only cache misses
reach one FastSpliceAI scoring task. The predictions
are then applied back to each cohort's splits without combining samples or
changing genotypes. Fully cached runs do not submit a scoring job.

The existing per-split `SplitVcf`, `NormaliseAndRegionFilterVcf`,
`AnnotateWithEchtvar`, and `AnnotateCsqWithBcftools` processes are unchanged.
FastSpliceAI gathers a side branch from the normalized splits while the ordinary
annotation path continues. Its completed cache is applied to each split after
the existing Echtvar and CSQ steps.

The `FASTSPLICEAI` subworkflow in `nextflow/fastspliceai.nf` gathers splits per
cohort, gathers cohorts, finds cache misses, scores missing sites, extends the
cache, and applies it to the annotated splits. It receives the normalized VCFs,
the per-split annotated VCFs, and the reference channel from `ANNOTATION`.
Its processes live in `nextflow/modules/fastspliceai/<Process>/main.nf`, following
the same layout as the Talos process modules. The published cache is
`<processed_annotations>/spliceai_cache/<identity>/cache.vcf.gz` with a `.tbi` index.
New annotations are combined with all existing cached annotations, including
alleles absent from the current cohorts, and replace the published cache.
One Talos pipeline runs at a time; no locking is used.

Only sites with complete SpliceAI predictions enter the cache; numeric zero
scores remain valid predictions. Each namespace is keyed by the GFF, generated
annotation, reference FASTA/index, scorer source and model weights, runtime
versions, prediction distance, and masking. Device and mixed precision do not
change the identity. Changed resources select a new namespace;
legacy cache files directly under the cache root are not reused. Each namespace
also contains `cache_identity.json` explaining its identity. Resources are hashed
once per run, including on resume, to detect changes inside a mutable image.
The transcript selection policy, MANE summary hash, and selected models are
recorded in `transcript_selection.json` beside the cache. This report is part of
the cache identity, so changes to MANE or the selection policy select a new namespace.

### Talos image with FastSpliceAI

FastSpliceAI is declared as the `spliceai` dependency in `pyproject.toml`,
using the fork at revision `ca5444366c53a3695a0711ebac3b5e526aa8235b`.
The same `uv.lock` and `uv sync` install Talos, FastSpliceAI, PyTorch, and CUDA
libraries into `/talos/.venv`. No separate scoring image or installation step
is required. The bundled default SpliceAI weights are used directly.

After changing dependencies, regenerate the lockfile before building:

```bash
uv lock
docker build -t talos:11.1.0 .
docker run --rm talos:11.1.0 spliceai --help
```

The Dockerfile uses `uv sync --locked` to require a matching lockfile, and
checks the installed FastSpliceAI CLI. uv HTTP requests have a 300-second
timeout; increase it with `--build-arg UV_HTTP_TIMEOUT=600` if needed.

Scoring uses the main `container` image. For Apptainer, convert
the Talos image to a SIF and select it in your deployment configuration:

```bash
apptainer build talos.sif docker-daemon://talos:11.1.0
```

```groovy
// cached-spliceai.config
params.container = '/absolute/path/talos.sif'
```

```bash
nextflow -c nextflow.config \
  -c cached-spliceai.config run main.nf -with-apptainer \
  --input_tsv /path/to/talos_input.tsv \
  --processed_annotations /path/to/processed_annotations \
  -output-dir /path/to/talos_outputs
```

The scoring process enables GPU passthrough (`--gpus all` for Docker, `--nv`
for Apptainer). The Slurm configuration allocates GPU resources only for
scoring. Reference FASTA/index, transcript tables, and the cache must be
accessible on compute nodes and inside containers. Scoring runs `spliceai` from
the Talos image PATH and always uses `ref_genome`. `GenerateFastSpliceAiAnnotation`
builds exon models from the same `params.ensembl_gff` channel used by consequence
annotation. It selects MANE Select and MANE Plus Clinical transcripts using the
existing `params.mane` summary. MANE is matched by stable Ensembl gene/transcript
IDs, with transcript version, chromosome, strand, and outer-bound checks against
the GFF and reference. Gene symbols come from Ensembl gene records; aliases and
renames in the MANE summary do not affect matching. The original summary is
required during analysis because the prepared `mane.json` omits versions and loci.

Each gene/reference-chromosome/strand group without a usable MANE Select model
receives one fallback model, while usable MANE Plus Clinical models are retained.
Fallback order is the gene's `canonical_transcript` declaration, a transcript
tagged `Ensembl_canonical`, then APPRIS rank and longest CDS. APPRIS principal
ranks 1–5 precede alternative ranks 1–2. CDS length sums coding intervals rather
than including introns. Without a canonical or coding model, the longest exonic
transcript is used. Remaining ties use longest exonic length then transcript ID.
Stable gene IDs distinguish genes sharing a symbol; a symbol is used only to
group parentless GFF models that provide no gene ID. Identical selected exon
models are emitted once per gene/locus. Different MANE clinical models may
intentionally produce multiple rows with the same gene symbol.

`transcript_selection.json` lists the selected versioned transcript IDs, their
gene IDs and symbols, their annotation row numbers, selection sources, and MANE
models that could not be used, including patch-contig or version mismatches.
This policy uses intact Ensembl transcript models; it does not merge splice sites
across isoforms or reproduce the original SpliceAI GENCODE V24 selection policy.

Rows retain zero-based starts, exclusive ends, and paired exon boundaries, sorted
in reference contig order and genomic coordinates. Chromosome aliases, including
M/MT/chrM, are resolved to reference names. Generated table bytes are hashed during
analysis and recorded in cohort metadata. Rebuild the Talos container after
updating the generator or cache identity code.

Full parameter details are in
[`docs/NextflowConfiguration.md`](docs/NextflowConfiguration.md).

### Web Runner

The optional Reflex runner is a web interface only: it submits the same
standalone preflight, preparation, and TALOS scripts that can be run from a
terminal. Its project defaults, shared resources, runtime paths, and GPU
settings are configured in [`runner/runner.config`](runner/runner.config).

```bash
REFLEX_CONDA_ENV=/faststorage/project/reanalyses_auh/env/talos2_env \
  bash runner/run_reflex_excel_runner.sh
```

See [`runner/README.md`](runner/README.md) for runner usage and configuration.

### Standalone Preflight

The AUH preflight workflow is kept separate from Talos proper under [`preflight/`](preflight/).

It covers:

- direct manifest input
- optional Excel-to-manifest mapping
- single-sample VCF preparation
- optional common-dbSNP subtraction
- family merge
- cohort merge
- optional SpliceAI annotation

Generic manifest-mode example:

```bash
nextflow run preflight/main.nf \
  -c preflight/nextflow.config \
  -c preflight/local.config \
  --sample_manifest preflight/examples/first_trio/sample_manifest.tsv \
  --pedigree preflight/examples/first_trio/pedigree.ped \
  -output-dir preflight/results
```

Synthetic preflight smoke test:

```bash
bash preflight/example_preflight.sh
```

Synthetic preflight plus Talos Apptainer smoke test:

```bash
bash preflight/example_with_talos_apptainer.sh
```

### Two-Family Direct Example

The simplest direct layout for many samples is:

- one `sample_manifest.tsv` with one row per single-sample VCF
- one `pedigree.ped` with all families
- one preflight run producing a merged cohort VCF
- one Talos TSV row pointing at that merged cohort VCF

Example `sample_manifest.tsv` for two trios:

```tsv
sample_id	vcf_path
PROBAND_A	/data/vcfs/PROBAND_A.vcf.gz
FATHER_A	/data/vcfs/FATHER_A.vcf.gz
MOTHER_A	/data/vcfs/MOTHER_A.vcf.gz
PROBAND_B	/data/vcfs/PROBAND_B.vcf.gz
FATHER_B	/data/vcfs/FATHER_B.vcf.gz
MOTHER_B	/data/vcfs/MOTHER_B.vcf.gz
```

Example `pedigree.ped`:

```tsv
FAM_A	PROBAND_A	FATHER_A	MOTHER_A	1	2
FAM_A	FATHER_A	0	0	1	1
FAM_A	MOTHER_A	0	0	2	1
FAM_B	PROBAND_B	FATHER_B	MOTHER_B	2	2
FAM_B	FATHER_B	0	0	1	1
FAM_B	MOTHER_B	0	0	2	1
```

Native preflight:

```bash
export PATH=/path/to/talos2_env/bin:$PATH
nextflow run preflight/main.nf \
  -c preflight/nextflow.config \
  -c preflight/local.config \
  --sample_manifest /path/to/sample_manifest.tsv \
  --pedigree /path/to/pedigree.ped \
  -output-dir /path/to/preflight_results
```

Native Talos:

```bash
nextflow -c nextflow.config run main.nf \
  --input_tsv /path/to/talos_input.tsv \
  --processed_annotations /path/to/processed_annotations \
  -output-dir /path/to/talos_outputs
```

Apptainer preflight:

```bash
apptainer exec \
  --bind /path/to/repo:/path/to/repo \
  --bind /path/to/data:/path/to/data \
  --pwd /path/to/repo \
  /path/to/talos-nextflow.sif \
  nextflow run preflight/main.nf \
    -c preflight/nextflow.config \
    -c preflight/local.config \
    --sample_manifest /path/to/sample_manifest.tsv \
    --pedigree /path/to/pedigree.ped \
    -output-dir /path/to/preflight_results
```

Apptainer Talos:

```bash
apptainer exec \
  --bind /path/to/repo:/path/to/repo \
  --bind /path/to/data:/path/to/data \
  --pwd /path/to/repo \
  /path/to/talos-nextflow.sif \
  nextflow -c nextflow.config run main.nf \
    --input_tsv /path/to/talos_input.tsv \
    --processed_annotations /path/to/processed_annotations \
    -output-dir /path/to/talos_outputs
```

For this AUH checkout, a local ignored runner script can also be created and run
in-place for two real trios. It generates its own `sample_manifest.tsv`,
`pedigree.ped`, `preflight.local.config`, and `talos_input.tsv` in the folder
where the script lives.

Local ignored example runners live under `local_examples/`, for example:

```bash
bash local_examples/run_single_trio_example.sh apptainer
bash local_examples/run_two_trios_example.sh native
```

### Family Bundle Generator

For real inputs, generate a timestamped ignored bundle with generic family roles:

```bash
bash preflight/create_family_bundle.sh \
  --proband /path/to/PROBAND.vcf.gz \
  --mother /path/to/MOTHER.vcf.gz \
  --father /path/to/FATHER.vcf.gz
```

Fast smoke mode keeps only the first variant from each provided VCF:

```bash
bash preflight/create_family_bundle.sh \
  --proband /path/to/PROBAND.vcf.gz \
  --mother /path/to/MOTHER.vcf.gz \
  --father /path/to/FATHER.vcf.gz \
  --fast true
```

This creates a local ignored bundle under `trio_runs/` with:

- `sample_manifest.tsv`
- `pedigree.ped`
- `preflight.local.config`
- `run_full_slurm.sh`

The generated script can be run either directly or through Slurm:

```bash
bash trio_runs/<bundle_name>/run_full_slurm.sh
```

```bash
sbatch trio_runs/<bundle_name>/run_full_slurm.sh
```

### Apptainer Build

The Apptainer tooling lives in [`apptainer/`](apptainer/).

Typical flow:

```bash
bash apptainer/fetch_env_to_folder.sh
bash apptainer/build_sif.sh
```

See [apptainer/README.md](apptainer/README.md) for details.

### Output Summary

Talos produces:

- JSON results for each cohort
- optional HTML reports
- labelled VCF outputs
- PanelApp cohort JSON
- reanalysis metadata such as `first_seen` and `evidence_last_updated`
- immutable `.metadata` source identity bundles with
  resource releases specified in `nextflow.config` and hashed during preparation;
  see [AARgus workflow metadata](docs/AargusMetadata.md)

The preflight layer produces:

- `prepared_manifest.tsv`
- family merge plan
- cohort merged VCF
- optional `cohort_merged_spliceai.vcf.gz`

## Part 2. Original Talos Basis

This repository is based on **Talos**, developed by the Population Genomics team.
The AUH extension keeps the core Talos method and overall workflow structure while
modifying runtime packaging, preflight handling, and local deployment patterns.

Upstream Talos overview:

- scalable variant prioritisation in known disease genes
- ACMG/AMP-aligned rule-based logic modules
- ClinVar and PanelApp-driven reanalysis support
- cohort-scale reruns with prior-history tracking

Primary upstream workflows:

- `preparation.nf`
- `main.nf`
- `talos_only.nf`

Primary upstream documentation:

- <https://populationgenomics.github.io/talos/>
- [docs/getting-started.md](docs/getting-started.md)
- [docs/Configuration.md](docs/Configuration.md)
- [docs/NextflowConfiguration.md](docs/NextflowConfiguration.md)

### Acknowledgement

This AUH fork builds on the original Talos project and its workflow design,
documentation base, and variant prioritisation framework.

If you use Talos in research or clinical work, cite the original publication:

Welland MJ, Ahlquist KD, De Fazio P, et al. Scalable automated reanalysis of
genomic data in research and clinical rare disease cohorts. *Nature Medicine*
(2026). <https://doi.org/10.1038/s41591-026-04477-5>
