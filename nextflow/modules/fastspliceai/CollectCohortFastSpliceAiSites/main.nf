include { quoteVcfInputs } from '../common'

process CollectCohortFastSpliceAiSites {
    tag "$cohort"
    input:
        tuple val(cohort), path(vcfs, stageAs: 'parts/??/*'), path(vcf_indices, stageAs: 'parts/??/*')
    output:
        tuple val(cohort), path('cohort_sites.vcf.gz'), path('cohort_sites.vcf.gz.tbi')
    script:
        def inputs = quoteVcfInputs(vcfs)
        """
        set -euo pipefail
        bcftools concat --no-version -a -G -Ov ${inputs} \\
            | awk -v info='##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI predictions">' '
                BEGIN { FS=OFS="\\t" }
                /^##fileformat=|^##contig=/ { print; next }
                /^#CHROM/ { print info; print; next }
                /^#/ { next }
                { print \$1,\$2,".",\$4,\$5,".",".","." }
            ' | bcftools sort -Ou \\
            | bcftools norm --no-version -d exact -Oz -o cohort_sites.vcf.gz
        bcftools index -t cohort_sites.vcf.gz
        """
}
