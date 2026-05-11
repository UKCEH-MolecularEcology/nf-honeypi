process TRIM_GALORE {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'

    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.fq.gz') ? "trimmed/${fn}" : "logs/${fn}" }

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1_val_1.fq.gz"), path("${meta.id}_R2_val_2.fq.gz"), emit: reads
    path "${meta.id}_trim_stats.txt",                                                         emit: stats
    path 'versions.yml',                                                                       emit: versions

    script:
    """
    ln -sf "${r1}" "${meta.id}_R1.fastq.gz"
    ln -sf "${r2}" "${meta.id}_R2.fastq.gz"

    trim_galore \\
        --paired \\
        --cores ${task.cpus} \\
        --quality 20 \\
        --length 50 \\
        "${meta.id}_R1.fastq.gz" "${meta.id}_R2.fastq.gz" \\
        2>&1 | tee "${meta.id}_trim_stats.txt"

    # Trim Galore appends _val_1 / _val_2 to the basename without extension
    # Input: {id}_R1.fastq.gz → output: {id}_R1_val_1.fq.gz
    # Rename if Trim Galore used a different suffix
    for suf in "_val_1.fq.gz" "_val_1.fastq.gz"; do
        [ -f "${meta.id}_R1\${suf}" ] && [ "\${suf}" != "_val_1.fq.gz" ] && \\
            mv "${meta.id}_R1\${suf}" "${meta.id}_R1_val_1.fq.gz"
    done
    for suf in "_val_2.fq.gz" "_val_2.fastq.gz"; do
        [ -f "${meta.id}_R2\${suf}" ] && [ "\${suf}" != "_val_2.fq.gz" ] && \\
            mv "${meta.id}_R2\${suf}" "${meta.id}_R2_val_2.fq.gz"
    done

    trim_galore --version | grep 'Trim Galore' | sed 's/Trim Galore version //' | \\
        awk '{print "\\"TRIM_GALORE\\":\\n    trim_galore: " \$1}' > versions.yml
    """

    stub:
    """
    touch "${meta.id}_R1_val_1.fq.gz" "${meta.id}_R2_val_2.fq.gz"
    touch "${meta.id}_trim_stats.txt"
    printf '"TRIM_GALORE":\n    trim_galore: 0.6.10\n' > versions.yml
    """
}
