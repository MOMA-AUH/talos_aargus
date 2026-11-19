include { shellQuote; completeFastSpliceAiPredictions } from '../common'

process CheckFastSpliceAiCache {
    cache false
    input:
        tuple path(vcf), path(vcf_index)
        val cache_directory
    output:
        tuple path('old_cache.vcf.gz'), path('old_cache.vcf.gz.tbi'), path('cache_misses.vcf.gz'), path('cache_misses.vcf.gz.tbi'), path('missing_count.txt')
    script:
        """
        set -euo pipefail
        if [[ -f ${shellQuote(cache_directory + '/cache.vcf.gz')} ]]; then
            bcftools view --no-version -G -Ov ${shellQuote(cache_directory + '/cache.vcf.gz')} \\
                | ${completeFastSpliceAiPredictions()} \\
                | bcftools view --no-version -Oz -o old_cache.vcf.gz
        else
            bcftools view --no-version -h ${shellQuote(vcf)} | bcftools view --no-version -Oz -o old_cache.vcf.gz
        fi
        bcftools index -t old_cache.vcf.gz
        bcftools isec --no-version -c none -C -w1 -Oz -o cache_misses.vcf.gz ${shellQuote(vcf)} old_cache.vcf.gz
        bcftools index -t cache_misses.vcf.gz
        bcftools index -n cache_misses.vcf.gz > missing_count.txt
        """
}
