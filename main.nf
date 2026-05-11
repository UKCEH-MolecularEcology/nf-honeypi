#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { HONEYPI } from './workflows/honeypi'

workflow {
    Channel
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
        .set { ch_reads }

    HONEYPI(ch_reads)
}
