process ITSX {
    tag "itsx_${its_region}"
    label 'process_medium'

    container 'quay.io/biocontainers/itsx:1.1.3--hdfd78af_1'

    publishDir "${params.outdir}/itsx", mode: 'copy'

    input:
    path asvs_fasta
    val  its_region
    val  its_min_len
    val  its_max_len

    output:
    path 'ASVs_ITSxed.fasta',  emit: itsxed
    path 'itsx_out.summary.txt', emit: summary
    path 'versions.yml',       emit: versions

    script:
    """
    # Count input sequences
    n_seqs=\$(grep -c '^>' "${asvs_fasta}" || echo 0)
    echo "ITSx input: \${n_seqs} sequences"

    ITSx \\
        -i "${asvs_fasta}" \\
        -o itsx_out \\
        --complement T \\
        --cpu ${task.cpus} \\
        --graphical F \\
        --save_regions ${its_region} \\
        --partial 50

    # Collect the extracted region
    ITS_FILE="itsx_out.${its_region}.fasta"
    if [ ! -f "\${ITS_FILE}" ] || [ ! -s "\${ITS_FILE}" ]; then
        echo "WARNING: ITSx produced no ${its_region} sequences. Creating empty file."
        touch ASVs_ITSxed.fasta
    else
        # Length filter using awk (replaces vsearch; the ITSx container has no python)
        awk -v min=${its_min_len} -v max=${its_max_len} '
            function flush() {
                if (hdr != "") {
                    if (length(seq) >= min && length(seq) <= max) { print hdr; print seq; kept++ }
                    else removed++
                }
            }
            /^>/ { flush(); hdr = \$0; seq = ""; next }
            { seq = seq \$0 }
            END {
                flush()
                printf "Length filter (%d-%d bp): kept %d, removed %d\\n", min, max, kept+0, removed+0 > "/dev/stderr"
            }' "itsx_out.${its_region}.fasta" > ASVs_ITSxed.fasta
    fi

    n_out=\$(grep -c '^>' ASVs_ITSxed.fasta 2>/dev/null || echo 0)
    echo "ITSx output: \${n_out} sequences after length filter"

    [ -f itsx_out.summary.txt ] || touch itsx_out.summary.txt

    ITSx --version 2>&1 | head -1 | sed 's/ITSx//' | tr -d 'v' | \\
        awk '{print "\\"ITSX\\":\\n    ITSx: " \$1}' > versions.yml || \\
        printf '"ITSX":\\n    ITSx: 1.1.3\\n' > versions.yml
    """

    stub:
    """
    printf '>ASV_0000000001\nACGTACGT\n' > ASVs_ITSxed.fasta
    touch itsx_out.summary.txt
    printf '"ITSX":\n    ITSx: 1.1.3\n' > versions.yml
    """
}
