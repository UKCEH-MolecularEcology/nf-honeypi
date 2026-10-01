process RDP_CLASSIFIER {
    tag "rdp_classifier"
    label 'process_medium'

    container 'quay.io/biocontainers/rdp_classifier:2.13--hdfd78af_1'

    publishDir "${params.outdir}/rdp_classifier", mode: 'copy'

    input:
    path  fasta
    path  db_dir
    val   confidence_threshold

    output:
    path 'assigned_taxonomy_prelim.txt', emit: prelim
    path 'ASVs_taxonomy.txt',            emit: taxonomy
    path 'versions.yml',                 emit: versions

    script:
    def conf     = confidence_threshold ?: 0.5
    def n_chunks = task.cpus
    def mem_gb   = Math.max(1, task.memory.toGiga().intdiv(task.cpus))
    """
    db_properties=\$(find -L "${db_dir}" -name '*.properties' | head -1)
    if [ -z "\${db_properties}" ]; then
        echo "ERROR: No .properties file found in ${db_dir}" >&2; exit 1
    fi

    n_seqs=\$(grep -c '^>' "${fasta}" || echo 0)
    if [ "\${n_seqs}" -eq 0 ]; then echo "ERROR: no sequences in ${fasta}" >&2; exit 1; fi
    echo "Classifying \${n_seqs} sequences against \${db_properties} in ${n_chunks} parallel chunks (${mem_gb} GB each)"

    # rdp_classifier is single-threaded: split into contiguous chunks and run them in parallel
    awk -v n=\${n_seqs} -v k=${n_chunks} '
        BEGIN { size = int((n + k - 1) / k); if (size < 1) size = 1 }
        /^>/  { i++; f = sprintf("chunk_%03d.fa", int((i - 1) / size) + 1) }
        { print > f }' "${fasta}"

    pids=()
    for chunk in chunk_*.fa; do
        rdp_classifier classify \\
            -Xmx${mem_gb}g \\
            -t "\${db_properties}" \\
            -o "\${chunk}.rdp" \\
            "\${chunk}" &
        pids+=(\$!)
    done
    rc=0
    for pid in "\${pids[@]}"; do wait "\${pid}" || rc=1; done
    [ \${rc} -eq 0 ] || { echo "ERROR: a classifier chunk failed" >&2; exit 1; }

    # chunk names are zero-padded, so the glob keeps the original ASV order
    cat chunk_*.fa.rdp > assigned_taxonomy_prelim.txt

    # Reformat RDP output to TSV taxonomy table (awk only — no python3 in this container)
    awk -v threshold=${conf} -f "${projectDir}/bin/parse_rdp.awk" assigned_taxonomy_prelim.txt > ASVs_taxonomy.txt

    rdp_classifier 2>&1 | head -1 | \\
        awk '{print "\\"RDP_CLASSIFIER\\":\\n    rdp_classifier: 2.13"}' > versions.yml || \\
        printf '"RDP_CLASSIFIER":\\n    rdp_classifier: 2.13\\n' > versions.yml
    """

    stub:
    """
    printf 'ASV_0000000001\t+\troot\trootrank\t1.0\tViridiplantae\tdomain\t0.99\tMagnoliophyta\tphylum\t0.95\tLiliopsida\tclass\t0.80\tPoales\torder\t0.75\tPoaceae\tfamily\t0.70\tHordeum\tgenus\t0.65\n' > assigned_taxonomy_prelim.txt
    printf 'ASV_0000000001\tk__Viridiplantae; p__Magnoliophyta; c__Liliopsida; o__Poales; f__Poaceae; g__Hordeum\t0.65\nASV_0000000002\tUnassignable\t1.0\n' > ASVs_taxonomy.txt
    printf '"RDP_CLASSIFIER":\n    rdptools: 2.0.2\n' > versions.yml
    """
}
