process WritePreparedResourceCatalog {
    input:
        val catalog

    output:
        path 'resources.metadata'

    script:
        def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
        def content = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(catalog))
        """
        set -euo pipefail
        mkdir -p resources.metadata
        printf '%s\\n' ${quote.call(content)} > resources.metadata/resources.json
        """
}
