def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

def quoteVcfInputs(vcfs) {
    // Nextflow unwraps a singleton path input; Path.collect iterates its components.
    def files = vcfs instanceof Collection ? vcfs : [vcfs]
    return files.collect { shellQuote(it) }.join(' ')
}

def completeFastSpliceAiPredictions() {
    // A populated INFO string can still contain missing scores or positions.
    // Numeric zero predictions are complete and must remain cache hits.
    return 'awk ' + shellQuote('''
        BEGIN { FS=OFS="\\t" }
        /^#/ { print; next }
        {
            count = split($8, info, ";")
            for (i=1; i<=count; i++) {
                if (info[i] !~ /^SpliceAI=/) continue
                entries = split(substr(info[i], 10), predictions, ",")
                complete = (entries > 0)
                for (j=1; j<=entries; j++) {
                    fields = split(predictions[j], values, "[|]")
                    if (fields != 10 || values[1] != $5) complete = 0
                    for (k=2; k<=fields; k++) {
                        if (values[k] == "." || values[k] == "") complete = 0
                    }
                }
                if (complete) print
                break
            }
        }
    ''')
}
