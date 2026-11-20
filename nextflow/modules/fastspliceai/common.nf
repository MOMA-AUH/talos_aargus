def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

def quoteVcfInputs(vcfs) {
    // Nextflow unwraps a singleton path input; Path.collect iterates its components.
    def files = vcfs instanceof Collection ? vcfs : [vcfs]
    return files.collect { shellQuote(it) }.join(' ')
}
