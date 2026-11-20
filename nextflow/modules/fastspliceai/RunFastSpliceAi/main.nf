include { shellQuote } from '../common'

process RunFastSpliceAi {
    tag "${shard_id}"

    input:
        tuple val(shard_id), path(misses), path(misses_index)
        path reference
        path reference_fai
        path annotation
    output:
        tuple val(shard_id), path("scored_${shard_id}.vcf.gz"), path("scored_${shard_id}.vcf.gz.tbi")
    script:
        def precision = params.fastspliceai_mixed_precision ? '--mixed-precision' : ''
        def scored = "scored_${shard_id}.vcf.gz"
        """
        set -euo pipefail
        spliceai \\
            -I ${shellQuote(misses)} \\
            -O ${shellQuote(scored)} \\
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
        bcftools index -t ${shellQuote(scored)}
        """
}
