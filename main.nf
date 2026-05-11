#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { HONEYPI } from './workflows/honeypi'

// Returns a list of [meta, r1, r2] maps from a FASTQ directory (Groovy only, no Channel).
def discover_read_pairs(dir_path) {
    def dir = file(dir_path, checkIfExists: true)
    if (!dir.isDirectory()) error "--input_dir '${dir_path}' is not a directory."

    def r1_re = /(?i)(.+?)(?:_R1(?:[._\-]\S*)?|_1)\.(fastq|fq)(\.gz)?$/
    def found  = []

    dir.listFiles().sort { it.name }.each { f ->
        def m = (f.name =~ r1_re)
        if (!m) return
        def stem   = m[0][1]
        def r2_name = f.name.replaceFirst(/(?i)(?<=_)R1(?=[._\-])/, 'R2')
                             .replaceFirst(/(?i)_R1\./, '_R2.')
                             .replaceFirst(/_1\.(fastq|fq)/, '_2.$1')
        def r2 = new File(dir.toString(), r2_name)
        if (!r2.exists()) return
        def sid = stem.replaceAll('_', '-').replaceAll(/[-]+$/, '')
        found << [[id: sid], file(f.toString()), file(r2.toString())]
    }
    if (found.isEmpty()) error "No paired FASTQ files found in --input_dir '${dir_path}'."
    log.info "Auto-discovered ${found.size()} sample(s) from ${dir_path}"
    return found
}

workflow {
    if (params.input_dir) {
        def pairs = discover_read_pairs(params.input_dir)
        pairs.each { meta, r1, r2 ->
            if ((meta.id as String).contains('_'))
                error "Sample ID '${meta.id}' (from --input_dir) contains underscore '_'. Replace with hyphen '-'."
        }
        ch_reads = Channel.fromList(pairs)
    } else if (params.input) {
        ch_reads = Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true, strip: true)
            .map { row ->
                if (!row.containsKey('sample') || !row.containsKey('fastq_1') || !row.containsKey('fastq_2'))
                    error "Samplesheet must contain columns: sample, fastq_1, fastq_2"
                if (row.sample.contains('_'))
                    error "Sample ID '${row.sample}' contains underscore '_'. Replace with hyphen '-' to avoid RDP output parsing issues."
                def meta = [id: row.sample]
                def r1   = file(row.fastq_1, checkIfExists: true)
                def r2   = file(row.fastq_2, checkIfExists: true)
                [meta, r1, r2]
            }
    } else {
        error "Provide either --input <samplesheet.csv> or --input_dir <directory_of_fastqs>"
    }

    HONEYPI(ch_reads)
}
