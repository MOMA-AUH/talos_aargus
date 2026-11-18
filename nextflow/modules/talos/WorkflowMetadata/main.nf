// Export the preparation capture consumed by AARgus.
process WorkflowMetadata {
    stageInMode 'symlink'
    input:
        tuple val(cohort), path(result, stageAs: 'result/*'), path(configuration, stageAs: 'configuration/*')
        path specification, stageAs: 'specification/*'
        val published_directory
        val effective_parameters
        val timestamp
        path generated_annotations, stageAs: 'annotation/*'

    output:
        tuple val(cohort), path("${cohort}.metadata_${timestamp}"), emit: metadata

    script:
        if (!(cohort ==~ /[A-Za-z0-9_.-]+/)) {
            throw new IllegalArgumentException('Metadata cohort must be a safe filename component')
        }
        def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
        def parameters = groovy.json.JsonOutput.toJson(effective_parameters)
        def bundle_name = "${cohort}.metadata_${timestamp}"
        def destination = "${published_directory}/${cohort}/outputs/${bundle_name}"
        def result_name = java.nio.file.Paths.get(result.toString()).getFileName().toString()
        def result_destination = "${published_directory}/${cohort}/outputs/${result_name}"
        """
        set -euo pipefail
        python -m talos.workflow_metadata prepare \
            --spec ${quote.call(specification)} \
            --configuration ${quote.call(configuration)} \
            --result-path ${quote.call(result_destination)} \
            --output-directory ${quote.call(bundle_name)} \
            --published-directory ${quote.call(destination)} \
            --parameters-json ${quote.call(parameters)} \
            --verify-result ${quote.call(result)}
        """
}

// Materialize the selected configuration as an immutable, reviewable input.
process WriteWorkflowMetadataSpecification {
    input:
        val specification

    output:
        path 'aargus-specification.json'

    script:
        def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
        def content = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(specification))
        """
        set -euo pipefail
        printf '%s\\n' ${quote.call(content)} > aargus-specification.json
        """
}

// Small receipt checks are shared by preparation reuse and analysis validation.
def preparedResourceFingerprint(source) {
    def attributes = java.nio.file.Files.readAttributes(source, 'unix:size,lastModifiedTime,ctime')
    [size: attributes.size,
        mtime_ns: attributes.lastModifiedTime.to(java.util.concurrent.TimeUnit.NANOSECONDS),
        ctime_ns: attributes.ctime.to(java.util.concurrent.TimeUnit.NANOSECONDS)]
}
