process EncodeAlphaMissense {
    memory 16.GB

    input:
        path tsv

    output:
        tuple path("alphamissense.vcf.gz"), path("alphamissense.zip")

    script:
        """
        set -euo pipefail

        python -m talos.annotation_scripts.parse_alphamissense \
            --input ${tsv} \
            --output alphamissense.vcf.gz

        echtvar encode alphamissense.zip ${projectDir}/echtvar/am_config.json alphamissense.vcf.gz
        """
}
