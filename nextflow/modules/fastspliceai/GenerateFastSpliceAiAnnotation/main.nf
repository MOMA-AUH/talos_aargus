include { shellQuote } from '../common'

process GenerateFastSpliceAiAnnotation {
    input:
        path gff
        path reference_fai
        path mane_summary

    output:
        path 'fastspliceai_annotation.tsv', emit: annotation
        path 'transcript_selection.json', emit: selection

    script:
        """
        set -euo pipefail
        python -m talos.annotation_scripts.generate_fastspliceai_annotation \\
            --reference-fai ${shellQuote(reference_fai)} \\
            --gff ${shellQuote(gff)} \\
            --mane-summary ${shellQuote(mane_summary)} \\
            --output fastspliceai_annotation.tsv \\
            --selection-report transcript_selection.json
        """
}
