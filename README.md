# nf-honeypi

Nextflow + Singularity pipeline for ITS2 (and ITS1) amplicon-based pollen/plant identification.
Implements the [honeypi](https://github.com/hsgweon/honeypi) workflow by H. Soon Gweon as a
reproducible, containerised Nextflow pipeline.

## Pipeline steps

```
samplesheet (CSV)
      |
      v
TRIM_GALORE  (per sample, Trim Galore v0.6.10)
      |
      v
DADA2        (joint denoising of all samples, outputs ASVs.fasta + count table)
      |
      v
ITSX         (ITS2/ITS1 extraction + length filter, replaces vsearch step)
      |
      v
CONSOLIDATE  (combines ITSxed + non-ITSxed sequences, pure Python)
      |
      v
DOWNLOAD_DB  (downloads Gweon RDP training data if not cached)
RDP_CLASSIFIER (classifies ASVs + reformats taxonomy, replaces BioPython/pandas)
      |
      v
FILTER_ASV_TABLE  (removes ASVs absent from consolidated FASTA, pure Python)
      |
      v
MERGE_DUPLICATES  (sums counts for identical taxonomy strings, pure Python)
      |
      v
results/counts/ASVs_counts_merged.txt
```

## Quick start

```bash
# 1. Clone
git clone https://github.com/ukceh-molecularecology/nf-honeypi
cd nf-honeypi

# 2. Edit assets/samplesheet.csv  (no underscores in sample IDs)

# 3. Run
nextflow run main.nf \
  -profile singularity \
  --input samplesheet.csv \
  --outdir results
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | *(required)* | Path to samplesheet CSV |
| `--outdir` | `results` | Output directory |
| `--its_region` | `ITS2` | ITS region to extract (`ITS1` or `ITS2`) |
| `--rdp_db_dir` | `null` | Pre-extracted RDP training data directory (skips download) |
| `--rdp_db_cache` | `<project>/db_cache` | Where to cache downloaded databases |
| `--dada2_trunc_len_r1` | `230` | R1 truncation length for DADA2 |
| `--dada2_trunc_len_r2` | `190` | R2 truncation length for DADA2 |
| `--dada2_max_ee` | `2,2` | Max expected errors for R1,R2 |
| `--its_min_len` | `100` | Minimum sequence length after ITS extraction (bp) |
| `--its_max_len` | `700` | Maximum sequence length after ITS extraction (bp) |
| `--rdp_confidence` | `0.5` | RDP confidence threshold for taxonomy assignment |
| `--r_lib_cache` | `~/.nextflow_r_libs` | R package cache (shared across runs) |

## Samplesheet format

```csv
sample,fastq_1,fastq_2
sample-A,/data/A_R1.fastq.gz,/data/A_R2.fastq.gz
sample-B,/data/B_R1.fastq.gz,/data/B_R2.fastq.gz
```

**Note:** sample IDs must not contain underscores (`_`). Use hyphens (`-`) instead.

## Outputs

| Path | Description |
|------|-------------|
| `results/counts/ASVs_counts_merged.txt` | Final taxonomy-merged count table |
| `results/counts/ASVs_counts_filtered.txt` | Per-ASV count table (pre-merge) |
| `results/rdp_classifier/ASVs_taxonomy.txt` | Taxonomy assignments |
| `results/dada2/ASVs.fasta` | Raw DADA2 ASVs |
| `results/consolidate/ASVs.fasta` | Consolidated ASVs (post-ITSx) |
| `results/dada2/dada2_stats.txt` | Read tracking through DADA2 |
| `results/itsx/itsx_out.summary.txt` | ITSx extraction summary |

## Requirements

- Nextflow >= 23.04
- Singularity (or Docker)
- Internet access for first run (database download), or `--rdp_db_dir` pointing to a local copy

## Credits

Original honeypi tool: H. Soon Gweon (https://github.com/hsgweon/honeypi)
Nextflow implementation: UKCEH Molecular Ecology
