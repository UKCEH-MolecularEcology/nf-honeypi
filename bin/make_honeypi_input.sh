#!/usr/bin/env bash
# Prepare already-trimmed (Trim Galore) reads for the ORIGINAL honeypi tool.
#
# Original honeypi expects raw "*.fastq.gz" names and looks for its own trim_galore
# outputs as "<name>_val_1.fq.gz". Reads already trimmed (*_val_1.fq.gz / *_val_2.fq.gz)
# make it crash, so this script creates a folder of symlinks named *.fastq.gz plus a
# matching readpairs list. NB: honeypi will trim these reads a second time.
#
# Usage: make_honeypi_input.sh <trimmed_reads_dir> <output_dir> [readpairs_list_name]
#        (default list name: <trimmed_reads_dir name>_readpairslist.txt)
set -euo pipefail

if [[ $# -lt 2 ]]; then
    sed -n '2,/^set -e/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'
    exit 1
fi

src=$(realpath "$1")
dst=$2
list_name=${3:-$(basename "$src")_readpairslist.txt}

[[ -d $src ]] || { echo "ERROR: '$src' is not a directory" >&2; exit 1; }

shopt -s nullglob
r1_files=("$src"/*_val_1.fq.gz)
if [[ ${#r1_files[@]} -eq 0 ]]; then
    echo "No Trim Galore reads (*_val_1.fq.gz) found in $src - nothing to do." >&2
    exit 2
fi

mkdir -p "$dst"
dst=$(realpath "$dst")
list="$dst/$list_name"
{
    echo '# Lines beginning with "#" is ignored. '
    echo -e '# SampleID\tFilename for forward reads\tFilename for reverse reads'
} > "$list"

n=0
for r1 in "${r1_files[@]}"; do
    stem=$(basename "$r1" _val_1.fq.gz)
    r2="$src/${stem/_R1/_R2}_val_2.fq.gz"
    if [[ ! -e $r2 ]]; then
        echo "WARNING: no R2 for $(basename "$r1") (expected $(basename "$r2")) - skipped" >&2
        continue
    fi
    f1="$stem.fastq.gz"
    f2="${stem/_R1/_R2}.fastq.gz"
    ln -sf "$r1" "$dst/$f1"
    ln -sf "$r2" "$dst/$f2"

    # Sample ID: text before "_S<digits>_L<digits>" (Illumina), else before "_R1"
    sid=$(sed -E 's/_S[0-9]+_L[0-9]+.*$//; t; s/_R1.*$//' <<< "$stem")
    sid=${sid//_/-}   # honeypi sample IDs must not contain underscores
    echo -e "$sid\t$f1\t$f2" >> "$list"
    n=$((n + 1))
done

dups=$(grep -v '^#' "$list" | cut -f1 | sort | uniq -d | head -3 | tr '\n' ' ')
[[ -z $dups ]] || echo "WARNING: duplicate sample IDs in list: $dups" >&2

echo "Linked $n sample(s) into $dst"
echo "Readpairs list: $list"
echo "Run: honeypi -i $dst -o <outdir> --amplicontype ITS2 -l $list"
