include { shellQuote } from '../common'

process FastSpliceAiCacheIdentity {
    cache false

    input:
        path gff
        path annotation
        path reference
        path reference_fai
        val distance
        val mask

    output:
        path 'identity.json'

    script:
        """
        set -euo pipefail
        python -m talos.fastspliceai_identity \\
            --gff ${shellQuote(gff)} \\
            --annotation ${shellQuote(annotation)} \\
            --reference ${shellQuote(reference)} \\
            --reference-fai ${shellQuote(reference_fai)} \\
            --distance ${shellQuote(distance)} \\
            --mask ${shellQuote(mask)} \\
            --output identity.json
        """
}
