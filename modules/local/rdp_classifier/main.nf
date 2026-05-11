process RDP_CLASSIFIER {
    tag "rdp_classifier"
    label 'process_medium'

    container 'quay.io/biocontainers/rdptools:2.0.2--hdfd78af_1'

    publishDir "${params.outdir}/rdp_classifier", mode: 'copy'

    input:
    path  fasta
    path  db_properties
    val   confidence_threshold

    output:
    path 'assigned_taxonomy_prelim.txt', emit: prelim
    path 'ASVs_taxonomy.txt',            emit: taxonomy
    path 'versions.yml',                 emit: versions

    script:
    def conf = confidence_threshold ?: 0.5
    """
    n_seqs=\$(grep -c '^>' "${fasta}" || echo 0)
    echo "Classifying \${n_seqs} sequences against ${db_properties}"

    classifier classify \\
        -Xmx${task.memory.toGiga()}g \\
        -t "${db_properties}" \\
        -o assigned_taxonomy_prelim.txt \\
        "${fasta}"

    # Reformat RDP output to TSV taxonomy table
    python3 << 'PYEOF'
import sys
import os

conf_threshold = ${conf}

rank_prefix = {
    'domain':       'k__',
    'superkingdom': 'k__',
    'phylum':       'p__',
    'class':        'c__',
    'order':        'o__',
    'family':       'f__',
    'genus':        'g__',
    'species':      's__',
}

out_lines = ['ASV_ID\\ttaxonomy\\tconfidence']

with open('assigned_taxonomy_prelim.txt') as fin:
    for line in fin:
        line = line.rstrip('\\n')
        if not line:
            continue
        parts = line.split('\\t')
        if len(parts) < 3:
            continue
        seq_id = parts[0].strip()
        if not seq_id or seq_id.startswith('#'):
            continue

        # RDP format: seq_id, orientation, (name, rank, conf)*
        # Orientation field may or may not be present; detect by checking
        # if parts[1] is '+' or '-' or 'no rank'
        start = 1
        if len(parts) > 1 and parts[1] in ('+', '-'):
            start = 2

        taxa_parts   = []
        last_conf    = 1.0
        i = start
        while i + 2 < len(parts):
            name  = parts[i].strip()
            rank  = parts[i + 1].strip().lower()
            try:
                conf = float(parts[i + 2])
            except (ValueError, IndexError):
                conf = 0.0
            i += 3

            if rank == 'rootrank' or not name or name.lower() == 'root':
                continue
            if rank not in rank_prefix:
                continue
            if conf < conf_threshold:
                break
            taxa_parts.append(rank_prefix[rank] + name)
            last_conf = conf

        taxonomy = ';'.join(taxa_parts) if taxa_parts else 'k__unclassified'
        out_lines.append(f'{seq_id}\\t{taxonomy}\\t{last_conf:.4f}')

with open('ASVs_taxonomy.txt', 'w') as fout:
    fout.write('\\n'.join(out_lines) + '\\n')

n = len(out_lines) - 1
classified = sum(1 for l in out_lines[1:] if 'unclassified' not in l)
print(f"Taxonomy: {classified}/{n} sequences classified above threshold {conf_threshold}",
      flush=True)
PYEOF

    classifier --version 2>&1 | head -1 | \\
        awk '{print "\\"RDP_CLASSIFIER\\":\\n    rdptools: " \$NF}' > versions.yml || \\
        printf '"RDP_CLASSIFIER":\\n    rdptools: 2.0.2\\n' > versions.yml
    """
}
