# Nextflow Configuration

This README documents the parameters used in the nextflow configuration file. Some default values for configuration variables are initially set up to facilitate the demonstration test-run.

| **Parameter**             | **Description**                                                                                                                                                                                                                        |
| ------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `large_files`             | Directory containing large external resources (e.g. gnomAD, MANE, AlphaMissense). See [large_files README](https://github.com/populationgenomics/talos/blob/main/large_files/README.md) for instructions on how to source this content |
| `processed_annotations`   | Output directory for all processed/reformatted data. Cohort/callset independent, to be used across multiple runs.                                                                                                                      |
| `outputDir`               | Defaults to `processed_annotations`, required for the 'workflow outputs' NextFlow functionality. Override for Annotation/Talos workflow using `-output-dir <path>`                                                                     |
| `input_vcf_extension`     | Linked to `shards`. If a sharded VCF is provided, this argument is the file extension to glob for. By default this is `vcf.bgz`, but could be altered to `vcf.gz`.                                                                     |
| `alphamissense_tsv`       | Path in `large_files` to the AlphaMissense raw data, reformatted into an echtvar archive during preparation                                                                                                                            |
| `alphamissense_zip`       | Path in `processed_annotations` to the reformatted AlphaMissense data, ready to be used by Echtvar for annotation                                                                                                                      |
| `ensembl_bed`             | Path in `processed_annotations` to the Ensembl BED file, one line per gene                                                                                                                                                             |
| `ensembl_merged_bed`      | Path in `processed_annotations` to the Ensembl BED file, overlapping regions merged                                                                                                                                                    |
| `mane_json`               | Path in `processed_annotations` to pre-processed MANE data                                                                                                                                                                             |
| `mitimpact_zip`           | Path in `processed_annotations` to the reformatted MitImpact reference data                                                                                                                                                            |
| `mitotip_zip`             | Path in `processed_annotations` to the reformatted MitoTip reference data                                                                                                                                                              |
| `napogee_zip`             | Path in `processed_annotations` to the reformatted nAPOGEE reference data                                                                                                                                                              |
| `ensembl_gff`             | Path in `large_files` to the Ensembl GFF3 file                                                                                                                                                                                         |
| `gnomad_zip`              | Path in `large_files` to the echtvar annotation file from [Zenodo](https://zenodo.org/records/15222100)                                                                                                                                |
| `fastspliceai_cache` | Cache root, default `<processed_annotations>/spliceai_cache`. Resource/scorer/selection identities automatically select subdirectories containing `cache.vcf.gz`, its index, `cache_identity.json` and `transcript_selection.json`. |
| `fastspliceai_distance` | Prediction distance, default 500, allowed 0–4999. Changes automatically select a new cache namespace. |
| `fastspliceai_mask` | Prediction mask, 0 or 1, default 0. Changes automatically select a new cache namespace. |
| `fastspliceai_mixed_precision` | Enable `--mixed-precision`, default `false`. Does not change cache identity. |
| `fastspliceai_device` | `cuda` (default) or `cpu`. Does not change cache identity. |
| `fastspliceai_threads` | Scoring threads, default 6.  |
| `mane`                    | Path in `large_files` to the MANE transcript data                                                                                                                                                                                      |
| `ref_genome`              | Path in `large_files` to the reference genome FASTA file                                                                                                                                                                               |
| `hpo`                     | Path in `large_files` to the HPO ontology file in `.obo` form                                                                                                                                                                          |
| `gen2phen`                | Path in `large_files` to the Geno~Phenotype databset in `.txt` form                                                                                                                                                                    |
| `phenio_db`               | Path in `large_files` to the Monarch Phenotype DB in `.db.gz` form                                                                                                                                                                     |
| `mitimpact_tsv`           | Path in `large_files` to the raw Mitimpact file                                                                                                                                                                                        |
| `mitotip_tsv`             | Path in `large_files` to the raw MitoTip file                                                                                                                                                                                          |
| `napogee_tsv`             | Path in `large_files` to the raw Mitochondrial nAPOGEE file                                                                                                                                                                            |
| `submission_summary`      | FTP URL for weekly ClinVar submission dump                                                                                                                                                                                             |
| `variant_summary`         | FTP URL for weekly ClinVar variant dump                                                                                                                                                                                                |
| `clinvar_blacklist`       | A string, containing quoted, space-delimited entries for all ClinVar submitter sites to blacklist (ignore).                                                                                                                            |
| `container`               | Talos image/SIF for tasks. Run with `-with-apptainer`; Docker is disabled by default. |
| `aargus_metadata_spec` | Optional JSON source specification used during preparation and analysis; metadata exports are mandatory. |
| `docker.enabled`          | Whether the Docker runtime is enabled in the active Nextflow profile                                                                                                                                                                   |

FastSpliceAI analysis also requires the original `mane` summary TSV (plain or
gzipped). It selects compatible MANE Select and Plus Clinical Ensembl models;
genes without a usable Select model receive an Ensembl fallback. The prepared
`mane_json` remains in use for downstream annotation. Selection provenance and
unusable MANE models are recorded in `transcript_selection.json` in the cache.

All processes use Slurm. Export deployment settings before starting Nextflow or
the runner:

```bash
export SLURM_ACCOUNT=your_account
export SLURM_GPU_PARTITIONS=gpu-l40s,gpu-h200
```

`SLURM_ACCOUNT` applies to CPU and GPU tasks and runner `sbatch` submissions.
`SLURM_GPU_PARTITIONS` selects the partitions for `RunFastSpliceAi`, defaulting to
`gpu-l40s,gpu-h200` when unset or empty. An unset `SLURM_ACCOUNT` uses Slurm's site
default account. The shared account and `--gpus=1` are included automatically. Other
resource and scoring settings are configured in Nextflow.

Run the workflow with `-with-apptainer`. FastSpliceAI scoring defaults to 6 CPU
cores, 16 GB RAM, and 8 hours. Apptainer scoring uses `--nv`. These GPU settings
apply only to `RunFastSpliceAi`; fully cached runs skip that process.

For project-specific scoring settings, use the project's Nextflow configuration:

```groovy
params.fastspliceai_threads = 8
params.fastspliceai_mixed_precision = true
process {
    withName: RunFastSpliceAi {
        memory = 24.GB
        time = '12:00:00'
    }
}
```
