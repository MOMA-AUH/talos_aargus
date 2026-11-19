include { shellQuote } from '../common'

process GenerateFastSpliceAiAnnotation {
    input:
        path gff
        path reference_fai

    output:
        path 'fastspliceai_annotation.tsv'

    script:
        """
        set -euo pipefail
        python -m talos.annotation_scripts.generate_fastspliceai_annotation \\
            --reference-fai ${shellQuote(reference_fai)} \\
            --gff ${shellQuote(gff)} \\
            --output fastspliceai_annotation.tsv
        """
}
