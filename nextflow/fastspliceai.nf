#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CollectCohortFastSpliceAiSites } from './modules/fastspliceai/CollectCohortFastSpliceAiSites/main'
include { CollectFastSpliceAiSites } from './modules/fastspliceai/CollectFastSpliceAiSites/main'
include { CheckFastSpliceAiCache } from './modules/fastspliceai/CheckFastSpliceAiCache/main'
include { RunFastSpliceAi } from './modules/fastspliceai/RunFastSpliceAi/main'
include { CompleteFastSpliceAi } from './modules/fastspliceai/CompleteFastSpliceAi/main'
include { AnnotateWithFastSpliceAi } from './modules/fastspliceai/AnnotateWithFastSpliceAi/main'
include { FastSpliceAiCacheIdentity } from './modules/fastspliceai/FastSpliceAiCacheIdentity/main'
include { GenerateFastSpliceAiAnnotation } from './modules/fastspliceai/GenerateFastSpliceAiAnnotation/main'

workflow FASTSPLICEAI {
    take:
        ch_normalized
        ch_vcfs
        ch_ref_genome
        ch_gff

    main:
    def cache = file(params.fastspliceai_cache).toFile().canonicalPath
    log.info "FastSpliceAI cache: ${cache}"
    ch_scoring_reference = ch_ref_genome.first()
    ch_fai = ch_scoring_reference.map { reference -> file("${reference}.fai", checkIfExists: true) }
    ch_mane_summary = channel.fromPath(params.mane, checkIfExists: true).first()
    GenerateFastSpliceAiAnnotation(ch_gff, ch_fai, ch_mane_summary)
    ch_annotation = GenerateFastSpliceAiAnnotation.out.annotation.first()
    ch_selection = GenerateFastSpliceAiAnnotation.out.selection.first()
    FastSpliceAiCacheIdentity(ch_gff, ch_annotation, ch_scoring_reference, ch_fai, ch_selection,
        params.fastspliceai_distance, params.fastspliceai_mask)
    ch_cache = FastSpliceAiCacheIdentity.out.map { identity ->
        def record = new groovy.json.JsonSlurper().parse(identity.toFile())
        def directory = "${cache}/${record.key}"
        log.info "FastSpliceAI cache identity: ${directory}"
        directory
    }.first()

    // Gather normalized sites across splits and cohorts without combining samples.
    ch_cohort_inputs = ch_normalized.groupTuple(by: 0).map { cohort, vcfs, indices ->
        tuple(cohort, vcfs.sort { it.name }, indices.sort { it.name })
    }
    CollectCohortFastSpliceAiSites(ch_cohort_inputs)
    ch_union_inputs = CollectCohortFastSpliceAiSites.out.collect(flat: false).filter { !it.isEmpty() }.map { cohorts ->
        cohorts.sort { a, b -> a[0] <=> b[0] }
        tuple(cohorts.collect { it[1] }, cohorts.collect { it[2] })
    }
    CollectFastSpliceAiSites(ch_union_inputs)
    CheckFastSpliceAiCache(CollectFastSpliceAiSites.out, ch_cache)
    ch_checked = CheckFastSpliceAiCache.out.branch {
        misses: it[4].text.trim().toInteger() > 0
        hits: true
    }
    RunFastSpliceAi(ch_checked.misses, ch_scoring_reference, ch_fai, ch_annotation)
    ch_cached = ch_checked.hits.map { old_cache, old_index, misses, misses_index, count ->
        tuple(old_cache, old_index, misses, misses_index)
    }
    CompleteFastSpliceAi(RunFastSpliceAi.out.mix(ch_cached), ch_cache, FastSpliceAiCacheIdentity.out.first(), ch_selection)

    // Apply the completed cache to each original per-split annotated VCF.
    ch_completed_cache = CompleteFastSpliceAi.out.cache.first()
    AnnotateWithFastSpliceAi(
        ch_vcfs,
        ch_completed_cache.map { cache_vcf, index -> cache_vcf },
        ch_completed_cache.map { cache_vcf, index -> index },
    )

    emit:
        vcfs = AnnotateWithFastSpliceAi.out
        annotation = ch_annotation
        selection = ch_selection
}
