include { quoteVcfInputs } from '../common'

process CollectFastSpliceAiSites {
    input:
        tuple path(vcfs, stageAs: 'cohorts/??/*'), path(vcf_indices, stageAs: 'cohorts/??/*')
    output:
        tuple path('union_sites.vcf.gz'), path('union_sites.vcf.gz.tbi')
    script:
        def inputs = quoteVcfInputs(vcfs)
        """
        set -euo pipefail
        bcftools concat --no-version -a -Ou ${inputs} | bcftools sort -Ou \\
            | bcftools norm --no-version -d exact -Oz -o union_sites.vcf.gz
        bcftools index -t union_sites.vcf.gz
        """
}
