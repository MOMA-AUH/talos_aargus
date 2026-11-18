// One independently cached task per resource, irrespective of cohort or resource version labels.
// Keep ordinary Nextflow file caching: deep caching would hash every input on each launch.
process HashResource {
    tag "${source_path}"
    stageInMode 'symlink'

    input:
        tuple val(source_path), path(artifact), val(fingerprint)

    output:
        tuple val(source_path), path('hash-*.json')

    script:
        def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
        def stamp = groovy.json.JsonOutput.toJson(fingerprint)
        def key = java.security.MessageDigest.getInstance('SHA-256').digest(source_path.bytes).encodeHex().toString()
        """
        set -euo pipefail
        python -m talos.workflow_metadata hash-artifact \
            --input ${quote.call(artifact)} --source-path ${quote.call(source_path)} \
            --fingerprint-json ${quote.call(stamp)} --output 'hash-${key}.json'
        """
}
