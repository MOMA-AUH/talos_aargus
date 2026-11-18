# AARgus metadata

Both `main.nf` and `talos_only.nf` publish a cohort metadata bundle at
`<output-directory>/<cohort>/outputs/<cohort>.metadata_<date>/`, using the same
`yyyy-MM-dd` date as the results. The bundle contains:

- `capture.json`: software version, assembly, resource hashes, configuration,
  and result path.
- `configuration.toml`: the configuration bytes referenced by the capture.

Each resource contains only `name`, `version`, and `sha256`:

```json
{
  "name": "FastSpliceAI",
  "version": "2.0.0",
  "sha256": "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"
}
```

The `FastSpliceAI` version comes from the installed `spliceai` package inside
Talos's image/SIF. Its hash identifies the transcript annotation table generated from the configured Ensembl GFF and reference during analysis. The generated artifact is passed directly from
annotation to metadata capture; it is not taken from the preparation catalog. Reference FASTAs, generated BED files, model
weights, and container images are not hashed for this metadata. PanelApp versions
are already in the Talos result.

Resource paths are not exported. The configuration reference still contains its
absolute path and SHA256; `result_path` identifies the cohort result. The exporter
checks completed results use schema 2.3.0 and the installed Talos software version.

The AARgus importer requires a corresponding update to accept hash-only resources;
this format is currently implemented only in `talos_aargus`.

## Prepare resources

Set release identifiers beside the resource paths in `nextflow.config`. Replace
`null` versions with the release or snapshot being prepared. ClinVar and PM5 use
`clinvar_version`, which defaults to the preparation month (`yyyy-MM`).

Preparation hashes task outputs or existing inputs directly, then publishes the
normal resources once. There is no separate resource-copying process. The catalog
at `<processed_annotations>/resources.metadata/resources.json` contains resource
names, versions, parameter bindings, and SHA256 hashes; it contains no resource
locations or references to work directories.

Checksum tasks use Nextflow's normal `-resume` caching. A source fingerprint in
the task inputs invalidates the checksum cache when that source changes, including
same-size edits whose modification time was restored. Retain the Nextflow work
and cache directories to reuse those tasks. Generated resources that are later
reused from their published locations can require a new checksum task because
the task input identity has changed.

Analysis consumes the catalog without reopening resource files for metadata
verification. Rerun preparation whenever a resource changes. To regenerate a
prepared output from changed raw inputs, remove that output or use a fresh
preparation destination. Configuration changes made only during analysis do not
replace prepared resource versions or hashes.

Rerun preparation to migrate catalogs from `talos/prepared-resources/1` to
`talos/prepared-resources/2`. Old per-resource receipt directories are no longer
used. Metadata remains usable after deleting preparation work directories.

Default names are `MANE`, `Ensembl`, `gnomAD`, `AlphaMissense`, `ClinVar`, `PM5`,
`HPO`, `HPO-G2P`, `Phenio`, `MitImpact`, `MitoTIP`, `nAPOGEE`, and `FastSpliceAI`.

## Override resource declarations

`--aargus_metadata_spec /absolute/path/specification.json` overrides preparation
resource declarations. Input paths are needed to calculate hashes but are not
included in the published catalog or capture. For example:

```json
{
  "assembly": "GRCh38",
  "resources": [
    {
      "name": "Ensembl",
      "version": "116",
      "path": "/resources/genes.gff3.gz",
      "parameter": "ensembl_gff"
    }
  ]
}
```

Use the same resource selection during preparation and analysis. `parameter`
binds a declaration to its workflow input; it is not exported in the cohort
capture. Names must be unique and versions nonempty. Paths may be absolute or
relative to the input specification. Additional preparation resources may omit
`parameter` and are selected by name during analysis.

## Standalone capture and verification

```bash
python -m talos.workflow_metadata prepare \
  --spec /resources/specification.json \
  --configuration /run/config.toml \
  --result-path /run/cohort_full_results.json \
  --output-directory /run/cohort.metadata

python -m talos.workflow_metadata verify \
  --capture /run/cohort.metadata/capture.json \
  --result /run/cohort_full_results.json \
  --output /run/verified.json
```

Standalone specifications can supply `sha256` directly instead of `path`.
Verification validates resource digest syntax and checks configuration and result
bytes. It does not require the original resource files.
