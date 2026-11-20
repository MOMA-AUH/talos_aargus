#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { AnnotateMitoVcf } from './modules/talos/AnnotateMitoVcf/main'
include { UnifiedPanelAppParser } from './modules/talos/UnifiedPanelAppParser/main'
include { RunHailFiltering } from './modules/talos/RunHailFiltering/main'
include { ValidateMOI } from './modules/talos/ValidateMOI/main'
include { HPOFlagging } from './modules/talos/HPOFlagging/main'
include { CreateTalosHTML } from './modules/talos/CreateTalosHTML/main'
include { StartupChecks } from './modules/talos/StartupChecks/main'
include { WorkflowMetadata; WriteWorkflowMetadataSpecification } from './modules/talos/WorkflowMetadata/main'

def processedAnnotationPath(String name) {
    return file("${params.processed_annotations}/${name}")
}

def preparedMetadataCatalog() {
    try {
        def path = processedAnnotationPath('resources.metadata/resources.json')
        def catalog = new groovy.json.JsonSlurper().parse(path.toFile())
        if (catalog.profile != 'talos/prepared-resources/2' || !(catalog.resources instanceof List)) {
            throw new IllegalArgumentException('Invalid prepared resource catalog')
        }
        if (catalog.assembly != params.ref_genome_assembly) {
            throw new IllegalArgumentException('Prepared assembly differs from ref_genome_assembly')
        }
        if (catalog.resources.any { !(it.sha256 instanceof String) || !(it.sha256 ==~ /[0-9a-f]{64}/) }) {
            throw new IllegalArgumentException('Prepared resources require SHA256 hashes')
        }
        catalog
    } catch (Exception error) {
        throw new IllegalArgumentException("Prepared resource catalog is missing or incompatible (${error.message}); run preparation.nf", error)
    }
}

def selectedPreparedMetadata(actual, catalog, explicit) {
    def requested = explicit != null ? (explicit.resources ?: []) : actual.resource_paths.collect { binding, path ->
        [parameter: binding]
    }
    if (explicit != null && explicit.assembly != catalog.assembly) {
        throw new IllegalArgumentException('Metadata override assembly differs from prepared resources; run preparation.nf')
    }
    def resources = requested.collect { resource ->
        if (resource.parameter == 'fastspliceai_annotation' && actual.resource_paths.fastspliceai_annotation) {
            // Generated during this analysis: hash the table used by the scorer,
            // rather than reusing a prepared identity for a bundled transcript table.
            return [name: resource.name ?: 'FastSpliceAI', version: params.ensembl_version,
                path: actual.resource_paths.fastspliceai_annotation, parameter: 'fastspliceai_annotation']
        }
        def matches = catalog.resources.findAll {
            resource.parameter != null ? it.parameter == resource.parameter : it.name == resource.name
        }
        if (matches.size() != 1 || (resource.version != null && matches[0].version != resource.version)) {
            throw new IllegalArgumentException("No prepared resource identity for ${resource.parameter ?: resource.name}; run preparation.nf")
        }
        def captured = new LinkedHashMap(matches[0])
        if (resource.name != null) captured.name = resource.name
        captured
    }
    actual.resource_paths.each { binding, path ->
        if (!resources.any { it.parameter == binding }) {
            throw new IllegalArgumentException("Metadata must select the prepared ${binding} input; run preparation.nf")
        }
    }
    def specification = [assembly: catalog.assembly, resources: resources]
    specification
}

workflow TALOS {
	take:
		ch_mane
		ch_gff
		ch_ref_genome
		ch_vcfs
        annotation_performed
        ch_fastspliceai_annotation

    main:
    // existence of these files is necessary for starting the workflow
    // we open them as a channel, and pass the channel through to the method
    ch_hpo_file = channel.fromPath(params.hpo, checkIfExists: true).first()
    ch_gen2phen = channel.fromPath(params.gen2phen, checkIfExists: true).first()
    ch_phenio = channel.fromPath(params.phenio_db, checkIfExists: true).first()

    // current year-month as a String, used to prompt for up to date resource updates
    def current_month = new java.util.Date().format('yyyy-MM')
    def timestamp = new java.util.Date().format('yyyy-MM-dd')

    // check if clinvar and panelapp data exist using the timestamp
    String current_clinvarbitration_all = "${params.processed_annotations}/clinvarbitration_${current_month}.vcf.bgz"
    String current_clinvarbitration_pm5 = "${params.processed_annotations}/clinvarbitration_${current_month}.pm5.tsv"

    if (!file(current_clinvarbitration_pm5).exists()) {
        println "ClinvArbitration data for this month (${current_clinvarbitration_pm5}) doesn't exist, run the Talos Prep workflow"
        exit 1
    }

    // read in each Clinvar input source as channel
    ch_clinvar_all = channel.fromPath(current_clinvarbitration_all, checkIfExists: true).first()
    ch_clinvar_all_tbi = channel.fromPath("${current_clinvarbitration_all}.tbi", checkIfExists: true).first()
    ch_clinvar_pm5 = channel.fromPath(current_clinvarbitration_pm5, checkIfExists: true).first()
    ch_bed = channel.fromPath(processedAnnotationPath('GRCh38.bed'), checkIfExists: true).first()

    String panelapp_path = "${params.processed_annotations}/panelapp_${current_month}.json"

    if (!file(panelapp_path).exists()) {
        println "PanelApp data for this month (${panelapp_path}) doesn't exist, run the Talos Prep workflow"
        exit 1
    }
    ch_panelapp = channel.fromPath(panelapp_path, checkIfExists: true).first()

    // run pre-Talos startup checks
    StartupChecks(
        ch_vcfs,
        ch_clinvar_all,
    )

    // UnifiedPanelAppParser
    ch_panel_app_inputs = StartupChecks.out
        .join(ch_vcfs)
        .map { cohort, check_file, _vcfs, pedigree, config, _history, _ext, _seqr, _mito ->
            tuple(cohort, check_file, config, pedigree)
        }

    UnifiedPanelAppParser(
        ch_panel_app_inputs,
	ch_panelapp,
	ch_hpo_file,
    )

    ch_run_hail_inputs = ch_vcfs
        .join(UnifiedPanelAppParser.out)
        .join(StartupChecks.out)
        .map { cohort, vcfs, pedigree, config, _history, _ext, _seqr, _mito, panelapp_data, check_file ->
            tuple(cohort, vcfs, panelapp_data, check_file, pedigree, config)
        }

    RunHailFiltering(
        ch_run_hail_inputs,
        ch_clinvar_all,
        ch_clinvar_all_tbi,
        ch_clinvar_pm5,
        ch_bed,
        ch_mane,
    )

    // surprise! It's Mito data!
    ch_mito_joined = ch_vcfs
        .join(UnifiedPanelAppParser.out)
        .join(StartupChecks.out)
        .map { cohort, _vcfs, pedigree, config, _history, _ext, _seqr, mito, panelapp_data, check_file ->
          tuple(cohort, mito, panelapp_data, pedigree, config)
    }

    ch_mito_branched = ch_mito_joined.branch {
        real:     it[1].name != 'NO_MITO'
        sentinel: it[1].name == 'NO_MITO'
    }

    ch_mito_for_annotation = ch_mito_branched.real
        .map { cohort, mito, panelapp, ped, config ->
            tuple(cohort, mito, panelapp, ped, config,
                  file(processedAnnotationPath('mitimpact.zip'), checkIfExists: true),
                  file(processedAnnotationPath('mitotip.zip'), checkIfExists: true),
                  file(processedAnnotationPath('napogee.zip'), checkIfExists: true))
        }

    AnnotateMitoVcf(
        ch_mito_for_annotation,
        ch_ref_genome,
        ch_gff,
        ch_clinvar_all,
        ch_clinvar_all_tbi,
    )

    ch_mito_resolved = AnnotateMitoVcf.out
        .mix(ch_mito_branched.sentinel.map { cohort, mito, _pa, _ped, _cfg -> tuple(cohort, mito) })

    // Validate MOI of all variants
    ch_validate_moi_inputs = RunHailFiltering.out
        .join(UnifiedPanelAppParser.out)
        .join(ch_vcfs)
        .join(ch_mito_resolved)
        .map { cohort, labelled_vcf, labelled_vcf_index, panelapp_out, _vcfs, pedigree, config, history, _ext, _seqr, _mito, anno_mito ->
            tuple(cohort, labelled_vcf, labelled_vcf_index, anno_mito, panelapp_out, pedigree, config, history)
        }

    ValidateMOI(
        ch_validate_moi_inputs,
        timestamp,
    )

    // Flag any relevant HPO terms
    ch_hpo_inputs = ValidateMOI.out
        .join(ch_vcfs)
        .map { cohort, talos_result_json, _vcfs, _pedigree, config, _history, _ext, _seqr, _mito ->
            tuple(cohort, talos_result_json, config)
        }

    HPOFlagging(
        ch_hpo_inputs,
        ch_mane,
        ch_gen2phen,
        ch_phenio,
        timestamp,
    )

    // Generate HTML report
    ch_create_html_inputs = HPOFlagging.out
        .join(UnifiedPanelAppParser.out)
        .join(ch_vcfs)
        .map { cohort, result_json, panelapp_data, _vcfs, _pedigree, config, _history, ext, seqr, _mito ->
            tuple(cohort, result_json, panelapp_data, config, ext, seqr)
        }

    CreateTalosHTML(
        ch_create_html_inputs,
        timestamp,
    )

    def prepared_catalog = preparedMetadataCatalog()
    ch_metadata_inputs = HPOFlagging.out.join(ch_vcfs)
        .map { cohort, result, _vcfs, _pedigree, config, _history, _ext, _seqr, _mito ->
            tuple(cohort, result, config)
        }
    def effective = [resource_paths: [
        mane_json: processedAnnotationPath('mane.json').toAbsolutePath().toString(),
        clinvar: file(current_clinvarbitration_all).toAbsolutePath().toString(),
        pm5: file(current_clinvarbitration_pm5).toAbsolutePath().toString(),
        hpo: file(params.hpo).toAbsolutePath().toString(),
        gen2phen: file(params.gen2phen).toAbsolutePath().toString(),
        phenio_db: file(params.phenio_db).toAbsolutePath().toString()
    ]]
    if (annotation_performed) {
        effective.resource_paths += [
            ensembl_gff: file(params.ensembl_gff).toAbsolutePath().toString(),
            gnomad_zip: file(params.gnomad_zip).toAbsolutePath().toString(),
            alphamissense_zip: processedAnnotationPath('alphamissense.zip').toAbsolutePath().toString()
        ]
    }
    // Select the resources used by this execution.
    ch_effective_base = ch_mito_branched.real.map { true }.toList().map { mitochondrial_cohorts ->
        if (mitochondrial_cohorts) {
            effective.resource_paths += [
                ensembl_gff: file(params.ensembl_gff).toAbsolutePath().toString(),
                mitimpact: processedAnnotationPath('mitimpact.zip').toAbsolutePath().toString(),
                mitotip: processedAnnotationPath('mitotip.zip').toAbsolutePath().toString(),
                napogee: processedAnnotationPath('napogee.zip').toAbsolutePath().toString()
            ]
        }
        effective
    }
    ch_effective = annotation_performed
        ? ch_effective_base.combine(ch_fastspliceai_annotation).map { actual, annotation ->
            actual.resource_paths.fastspliceai_annotation = annotation.toAbsolutePath().toString()
            actual
        }.first()
        : ch_effective_base
    def declarations = null
    if (params.aargus_metadata_spec) {
        def specification = file(params.aargus_metadata_spec, checkIfExists: true)
        declarations = new groovy.json.JsonSlurper().parse(specification.toFile())
    }
    ch_prepared_specification = ch_effective.map { actual ->
        selectedPreparedMetadata(actual, prepared_catalog, declarations)
    }
    WriteWorkflowMetadataSpecification(ch_prepared_specification)
    WorkflowMetadata(ch_metadata_inputs, WriteWorkflowMetadataSpecification.out,
        file(workflow.outputDir).toAbsolutePath().toString(), ch_effective, timestamp,
        ch_effective.map { actual ->
            actual.resource_paths.fastspliceai_annotation ? [file(actual.resource_paths.fastspliceai_annotation)] : []
        })
    ch_metadata = WorkflowMetadata.out.metadata

    emit:
        metadata = ch_metadata
	json = HPOFlagging.out
	html = CreateTalosHTML.out
	labelled = RunHailFiltering.out
	panelapp = UnifiedPanelAppParser.out
}
