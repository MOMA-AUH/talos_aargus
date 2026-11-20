include { shellQuote; quoteVcfInputs; completeFastSpliceAiPredictions } from '../common'

process CompleteFastSpliceAi {
    publishDir { cache_directory }, mode: 'copy', overwrite: true
    input:
        tuple path(old_cache), path(old_index), path(scored), path(scored_index)
        val cache_directory
        path identity
        path selection, stageAs: 'selection_source.json'
    output:
        tuple path('cache.vcf.gz'), path('cache.vcf.gz.tbi'), emit: cache
        path 'cache_identity.json', emit: receipt
        path 'transcript_selection.json', emit: selection
    script:
        """
        set -euo pipefail
        bcftools concat --no-version -a -Ov ${shellQuote(old_cache)} ${quoteVcfInputs(scored)} \\
            | ${completeFastSpliceAiPredictions()} \\
            | bcftools sort -Ou | bcftools norm --no-version -d exact -Oz -o cache.vcf.gz
        bcftools index -t cache.vcf.gz
        cp ${shellQuote(identity)} cache_identity.json
        cp ${shellQuote(selection)} transcript_selection.json
        """
}
