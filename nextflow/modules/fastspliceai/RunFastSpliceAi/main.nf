include { shellQuote } from '../common'

process RunFastSpliceAi {
    input:
        tuple path(old_cache), path(old_index), path(misses), path(misses_index), path(missing_count)
        path reference
        path reference_fai
        path annotation
    output:
        tuple path(old_cache), path(old_index), path('scored.vcf.gz'), path('scored.vcf.gz.tbi')
    script:
        def precision = params.fastspliceai_mixed_precision ? '--mixed-precision' : ''
        """
        set -euo pipefail
        spliceai \\
            -I ${shellQuote(misses)} \\
            -O scored.vcf.gz \\
            -R ${shellQuote(reference)} \\
            -A ${shellQuote(annotation)} \\
            -D ${shellQuote(params.fastspliceai_distance)} \\
            -M ${shellQuote(params.fastspliceai_mask)} \\
            --device ${shellQuote(params.fastspliceai_device)} \\
            --batch-size 16 \\
            --threads ${task.cpus} \\
            --output-type z \\
            --overwrite-existing \\
            ${precision}
        bcftools index -t scored.vcf.gz
        """
}
