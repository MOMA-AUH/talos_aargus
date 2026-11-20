include { shellQuote } from '../common'

process SplitFastSpliceAiMisses {
    memory 2.GB
    cpus 1

    input:
        tuple path(misses), path(misses_index), path(missing_count)
    output:
        tuple path('shard_*.vcf.gz', arity: '1..*'), path('shard_*.vcf.gz.tbi', arity: '1..*')
    script:
        if (!(params.fastspliceai_split_n.toString() ==~ /[1-9][0-9]*/)) {
            throw new IllegalArgumentException('fastspliceai_split_n must be a positive integer')
        }
        def splitSize = params.fastspliceai_split_n as long
        """
        set -euo pipefail
        expected_records=\$(cat ${shellQuote(missing_count)})
        bcftools view --no-version -h ${shellQuote(misses)} > header.txt
        bcftools view --no-version -H ${shellQuote(misses)} \\
            | split -l ${splitSize} -d -a 8 - shard_ \\
                --filter='cat header.txt - | bcftools view --no-version -Oz -o "\${FILE}.vcf.gz"'

        shopt -s nullglob
        shards=(shard_*.vcf.gz)
        total_records=0
        for shard in "\${shards[@]}"; do
            bcftools index -t "\$shard"
            records=\$(bcftools index -n "\$shard")
            total_records=\$((total_records + records))
        done
        if [[ \$total_records -ne \$expected_records ]]; then
            echo 'FastSpliceAI shards do not contain every missing record' >&2
            exit 1
        fi
        """
}
