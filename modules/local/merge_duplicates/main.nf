process MERGE_DUPLICATES {
    tag "merge_duplicates"
    label 'process_low'

    container 'python:3.10-slim'

    publishDir "${params.outdir}/counts", mode: 'copy'

    // Inputs are staged under distinct names: the outputs below are called ASVs.fasta /
    // ASVs_taxonomy.txt, and writing them must not overwrite an upstream file through a symlink.
    input:
    path counts_txt,    stageAs: 'in_counts_filtered.txt'    // ASVs_counts_filtered.txt
    path asvs_fasta,    stageAs: 'in_consolidated.fasta'     // consolidated ASVs.fasta
    path taxonomy_txt,  stageAs: 'in_rdp_taxonomy.txt'       // native-format taxonomy (no header)
    path sample_ids,    stageAs: 'in_sample_ids.txt'         // original sample IDs, one per line

    output:
    path 'ASVs_counts.txt',        emit: counts          // native honeypi final table
    path 'ASVs_taxonomy.txt',      emit: taxonomy        // native honeypi final taxonomy
    path 'ASVs.fasta',             emit: fasta           // native honeypi final FASTA
    path 'ASVs_counts_merged.txt', emit: counts_by_taxon // extra: counts per identical taxonomy
    path 'versions.yml',           emit: versions

    script:
    """
    # Native honeypi (honeypi_mergeDuplicateASV): merge ASVs with identical sequences
    python3 "${projectDir}/bin/merge_duplicate_asvs.py" \\
        --fasta    in_consolidated.fasta \\
        --counts   in_counts_filtered.txt \\
        --taxonomy in_rdp_taxonomy.txt \\
        --sample-ids in_sample_ids.txt \\
        --strip-sample-regex='${params.sample_strip_regex}' \\
        --outdir .

    python3 -c "import platform; print('\\"MERGE_DUPLICATES\\":\\n    python: ' + platform.python_version())" > versions.yml
    """

    stub:
    """
    printf 'sample-A\\tsample-B\nASV_0000000001\\t100\\t80\nASV_0000000002\\t50\\t120\n' > ASVs_counts.txt
    printf 'ASV_0000000001\\tk__Viridiplantae; p__Streptophyta; c__Liliopsida; o__Poales; f__Poaceae; g__Hordeum\\t0.65\nASV_0000000002\\tUnassignable\\t1.0\n' > ASVs_taxonomy.txt
    printf '>ASV_0000000001\nACGTACGT\n>ASV_0000000002\nTGCATGCA\n' > ASVs.fasta
    printf 'taxonomy\\tsample-A\\tsample-B\nUnassignable\\t150\\t200\n' > ASVs_counts_merged.txt
    printf '"MERGE_DUPLICATES":\\n    python: 3.10\\n' > versions.yml
    """
}
