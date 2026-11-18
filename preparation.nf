#!/usr/bin/env nextflow

/*
This workflow takes ownership over all the input data preparation steps across both Talos workflows (annotation + talos)

The data prepared in this workflow is:

- AlphaMissense data, reformatted into an echtvar archive
- PanelApp data, freshly downloaded
- ClinVar data, freshly downloaded and reformatted into a VCF plus PM5 TSV
- Mane data, downloaded and reformatted
- BED file representing ROI for all Ensembl genes
*/

nextflow.enable.dsl=2

include { preparedResourceFingerprint } from './nextflow/modules/talos/WorkflowMetadata/main'

include { AnnotateClinvarWithBcftools } from './nextflow/modules/prep/AnnotateClinvarWithBcftools/main'
include { CreateRoiFromGff3 } from './nextflow/modules/prep/CreateRoiFromGff3/main'
include { DownloadClinVarFiles } from './nextflow/modules/prep/DownloadClinVarFiles/main'
include { DownloadPanelApp } from './nextflow/modules/prep/DownloadPanelApp/main'
include { EncodeAlphaMissense } from './nextflow/modules/prep/EncodeAlphaMissense/main'
include { EncodeMitimpact } from './nextflow/modules/prep/EncodeMitimpact/main'
include { EncodeMitotip } from './nextflow/modules/prep/EncodeMitotip/main'
include { EncodeNapogee } from './nextflow/modules/prep/EncodeNapogee/main'
include { HashResource } from './nextflow/modules/prep/HashResource/main'
include { MakeClinvarbitrationPm5 } from './nextflow/modules/prep/MakeClinvarbitrationPm5/main'
include { ResummariseRawSubmissions } from './nextflow/modules/prep/ResummariseRawSubmissions/main'
include { ParseManeIntoJson } from './nextflow/modules/prep/ParseManeIntoJson/main'
include { WritePreparedResourceCatalog } from './nextflow/modules/prep/WritePreparedResourceCatalog/main'


def preparedResourceDeclarations(paths) {
    def requiredText = { value, label ->
        if (!(value instanceof CharSequence) || !value.toString().trim()) {
            throw new IllegalArgumentException("Set params.${label} to a nonempty string for AARgus metadata")
        }
        value.toString()
    }
    // binding: [version parameter, display name]
    def resources = [
        mane_json: ['mane_version', 'MANE'],
        ensembl_gff: ['ensembl_version', 'Ensembl'],
        gnomad_zip: ['gnomad_version', 'gnomAD'],
        alphamissense_zip: ['alphamissense_version', 'AlphaMissense'],
        clinvar: ['clinvar_version', 'ClinVar'],
        pm5: ['clinvar_version', 'PM5'],
        hpo: ['hpo_version', 'HPO'],
        gen2phen: ['gen2phen_version', 'HPO-G2P'],
        phenio_db: ['phenio_version', 'Phenio'],
        mitimpact: ['mitimpact_version', 'MitImpact'],
        mitotip: ['mitotip_version', 'MitoTIP'],
        napogee: ['napogee_version', 'nAPOGEE']
    ]
    paths.collect { binding, observed ->
        def resource = resources[binding]
        if (resource == null) throw new IllegalArgumentException("No AARgus metadata resource mapping for ${binding}")
        def declaration = [parameter: binding, path: file(observed, checkIfExists: true).toRealPath().toString(),
            name: resource[1], version: requiredText.call(params[resource[0]], resource[0])]
        declaration
    }
}

def preparationMetadataPlan(stored, destination, timestamp) {
    def paths = [
        mane_json: stored['mane.json'],
        ensembl_gff: file(params.ensembl_gff, checkIfExists: true),
        gnomad_zip: file(params.gnomad_zip, checkIfExists: true),
        alphamissense_zip: stored['alphamissense.zip'],
        clinvar: stored["clinvarbitration_${timestamp}.vcf.bgz".toString()],
        pm5: stored["clinvarbitration_${timestamp}.pm5.tsv".toString()],
        hpo: file(params.hpo, checkIfExists: true),
        gen2phen: file(params.gen2phen, checkIfExists: true),
        phenio_db: file(params.phenio_db, checkIfExists: true),
        mitimpact: stored['mitimpact.zip'],
        mitotip: stored['mitotip.zip'],
        napogee: stored['napogee.zip']
    ].collectEntries { binding, source -> [binding, source.toRealPath().toString()] }
    def declarations
    if (params.aargus_metadata_spec) {
        def specification = file(params.aargus_metadata_spec, checkIfExists: true).toRealPath()
        def explicit = new groovy.json.JsonSlurper().parse(specification.toFile())
        if (explicit.assembly != params.ref_genome_assembly) throw new IllegalArgumentException('Metadata specification assembly differs from ref_genome_assembly')
        declarations = (explicit.resources ?: []).collect { resource ->
            def captured = resource.subMap(['name', 'version', 'path', 'parameter'])
            if (!(captured.version instanceof CharSequence) || !captured.version.toString().trim()) {
                throw new IllegalArgumentException("Set resource ${captured.name} version before preparation")
            }
            def source = java.nio.file.Paths.get(captured.path)
            captured.path = file(source.isAbsolute() ? source : specification.parent.resolve(source), checkIfExists: true).toRealPath().toString()
            captured
        }
    } else {
        declarations = preparedResourceDeclarations(paths)
    }
    def catalog = [profile: 'talos/prepared-resources/2', assembly: params.ref_genome_assembly, resources: declarations]
    def inputs = declarations.collect { file(it.path, checkIfExists: true) }
    inputs = inputs.collect { it.toRealPath() }.unique().sort { it.toString() }
    def to_hash = inputs.collect { source ->
        tuple(source.toString(), source, preparedResourceFingerprint(source))
    }
    [catalog: catalog, hash_inputs: to_hash]
}

workflow {
    main:
    def timestamp = new java.util.Date().format('yyyy-MM')

    ch_gff = channel.fromPath(params.ensembl_gff, checkIfExists: true)

    // generate the AlphaMissense HT - long running, stored in a separate folder
    if (!file(params.alphamissense_zip).exists()) {
	ch_alphamissense_tsv = channel.fromPath(params.alphamissense_tsv, checkIfExists: true)
        EncodeAlphaMissense(ch_alphamissense_tsv)
        ch_alphamissense_zip = EncodeAlphaMissense.out
    } else {
        ch_alphamissense_zip = channel.fromPath(params.alphamissense_zip)
    }

    // generate/format the various Mitochondrial annotation sources
    if (!file(params.mitimpact_zip).exists()) {
	ch_mitimpact_tsv = channel.fromPath(params.mitimpact_tsv, checkIfExists: true)
        EncodeMitimpact(ch_mitimpact_tsv)
        ch_mitimpact_zip = EncodeMitimpact.out
    } else {
        ch_mitimpact_zip = channel.fromPath(params.mitimpact_zip)
    }

    // generate/format the various Mitochondrial annotation sources
    if (!file(params.mitotip_zip).exists()) {
	ch_mitotip_tsv = channel.fromPath(params.mitotip_tsv, checkIfExists: true)
        EncodeMitotip(ch_mitotip_tsv)
        ch_mitotip_zip = EncodeMitotip.out
    } else {
        ch_mitotip_zip = channel.fromPath(params.mitotip_zip)
    }

    // generate/format the various Mitochondrial annotation sources
    if (!file(params.napogee_zip).exists()) {
	ch_napogee_tsv = channel.fromPath(params.napogee_tsv, checkIfExists: true)
        EncodeNapogee(ch_napogee_tsv)
        ch_napogee_zip = EncodeNapogee.out
    } else {
        ch_napogee_zip = channel.fromPath(params.napogee_zip)
    }

    // does this month's clinvarbitration data exist?
    String current_clinvarbitration_all = "${params.processed_annotations}/clinvarbitration_${timestamp}.vcf.bgz"
    String current_clinvarbitration_pm5 = "${params.processed_annotations}/clinvarbitration_${timestamp}.pm5.tsv"

    if (file(current_clinvarbitration_pm5).exists()) {
        ch_clinvar_all = channel.of(
            tuple(
                file(current_clinvarbitration_all, checkIfExists: true),
                file("${current_clinvarbitration_all}.tbi", checkIfExists: true),
            ),
        )
        ch_clinvar_pm5 = channel.fromPath(current_clinvarbitration_pm5)
    } else {
        // new workflow elements to go and create it from raw data
        String subfile = "${params.large_files}/submissions_${timestamp}.txt.gz"
        String varfile = "${params.large_files}/variants_${timestamp}.txt.gz"

        if (file(subfile).exists() && file(varfile).exists()) {
            ch_clinvar_sub = channel.fromPath(subfile)
            ch_clinvar_var = channel.fromPath(varfile)
        } else {
            println "Attempting to download ClinVar raw data, requires internet connection."
            println "If this step fails, try re-running gather_files.sh in the `large_files` directory."
            // this step requires an internet connection, which may be problematic at some sites
            DownloadClinVarFiles(timestamp)
            ch_clinvar_sub = DownloadClinVarFiles.out.submissions
            ch_clinvar_var = DownloadClinVarFiles.out.variants
        }

        // reinterpret the results using altered heuristics
        ResummariseRawSubmissions(
            ch_clinvar_var,
            ch_clinvar_sub,
            timestamp,
        )

        ch_ref_fa = channel.fromPath(params.ref_genome, checkIfExists: true)

        // annotate the SNV VCF using BCFtools
        AnnotateClinvarWithBcftools(
            ResummariseRawSubmissions.out.vcf,
            ch_ref_fa,
            ch_gff,
        )

        MakeClinvarbitrationPm5(
            AnnotateClinvarWithBcftools.out,
            timestamp,
        )

        ch_clinvar_all = ResummariseRawSubmissions.out.vcf
        ch_clinvar_pm5 = MakeClinvarbitrationPm5.out
    }

    // generate the Region-of-interest BED file from Ensembl GFF3
    // generates a per-gene BED file with ID annotations
    // and a overlap-merged version of the same for more efficient region filtering
    if (!file(params.ensembl_bed).exists() || !file(params.ensembl_merged_bed).exists() || !file(params.ensembl_symbol_lookup).exists()) {
        CreateRoiFromGff3(ch_gff)
        ch_bed = CreateRoiFromGff3.out.bed
        ch_merged_bed = CreateRoiFromGff3.out.merged_bed
        ch_symbol_lookup = CreateRoiFromGff3.out.json_lookup
    } else {
        ch_merged_bed = channel.fromPath(params.ensembl_merged_bed, checkIfExists: true)
        ch_bed = channel.fromPath(params.ensembl_bed, checkIfExists: true)
        ch_symbol_lookup = channel.fromPath(params.ensembl_symbol_lookup, checkIfExists: true)
    }

    // pull and parse the MANE data into JSON for downstream annotation
    if (!file(params.mane_json).exists()) {
        ch_mane_summary = channel.fromPath(params.mane, checkIfExists: true)
        ParseManeIntoJson(ch_mane_summary)
        ch_mane_json = ParseManeIntoJson.out
    } else {
        ch_mane_json = channel.fromPath(params.mane_json, checkIfExists: true)
    }

    String current_panelapp = "${params.processed_annotations}/panelapp_${timestamp}.json"

    if (!file(current_panelapp).exists()) {
        DownloadPanelApp(
            ch_mane_json,
            timestamp,
        )
        panelapp_out = DownloadPanelApp.out
    } else {
        panelapp_out = channel.fromPath(current_panelapp, checkIfExists: true)
    }

    // Hash the existing task outputs; normal publishing copies resources once.
    def destination = file(workflow.outputDir).toAbsolutePath().toString()
    ch_prepared_files = ch_alphamissense_zip.mix(ch_bed, ch_merged_bed, ch_symbol_lookup,
        ch_clinvar_all, ch_clinvar_pm5, ch_mitimpact_zip, ch_mitotip_zip, ch_mane_json, ch_napogee_zip, panelapp_out)
        .flatMap { value -> (value instanceof Collection) ? value : [value] }
        .map { source -> tuple(java.nio.file.Paths.get(source.toString()).getFileName().toString(), source) }
    ch_metadata_plan = ch_prepared_files.toList().map { entries ->
        def stored = entries.collectEntries { name, source -> [name, source] }
        preparationMetadataPlan(stored, destination, timestamp)
    }
    HashResource(ch_metadata_plan.flatMap { it.hash_inputs })
    ch_receipts = HashResource.out
        .map { source, receipt -> [source: source, receipt: receipt] }.toList()
    ch_catalog_inputs = ch_metadata_plan.combine(ch_receipts.map { [receipts: it] }).map { plan, received ->
        def catalog = new LinkedHashMap(plan.catalog)
        def hashes = received.receipts.collectEntries { entry ->
            def receipt = new groovy.json.JsonSlurper().parse(entry.receipt.toFile())
            [entry.source, receipt.sha256]
        }
        catalog.resources = catalog.resources.collect { resource ->
            def captured = resource.subMap(['name', 'version', 'parameter'])
            captured.sha256 = hashes[resource.path]
            captured
        }
        catalog
    }
    WritePreparedResourceCatalog(ch_catalog_inputs)

    publish:
        resources = ch_prepared_files.map { name, source -> source }
        metadata = WritePreparedResourceCatalog.out
}

output {
    resources {
        mode 'copy'
    }
    metadata {
        mode 'copy'
    }
}
