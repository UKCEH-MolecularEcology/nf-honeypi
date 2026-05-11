# Changelog

## v2.0.0

Complete rewrite as a Nextflow DSL2 pipeline with Singularity containers.

Based on the original honeypi tool by H. Soon Gweon:
https://github.com/hsgweon/honeypi

### Changes from v1

- All steps run in containers — no local tool installation required
- DADA2 truncLen and maxEE are configurable parameters (not hardcoded)
- Sample ID validation at pipeline entry (no underscores) eliminates a class of crashes in the original RDP output parser
- `filterbyname.sh` (bbmap) replaced with pure Python in CONSOLIDATE
- `vsearch` length filter replaced with pure Python in ITSX
- BioPython + pandas in MERGE_DUPLICATES replaced with pure Python
- `progressbar2` + `requests` download replaced with stdlib `urllib`
- RDP database cached with Nextflow `storeDir` — downloaded once, reused across runs
- `--rdp_db_dir` parameter accepts a pre-existing local database directory
- Per-run read tracking written to `dada2_stats.txt`
