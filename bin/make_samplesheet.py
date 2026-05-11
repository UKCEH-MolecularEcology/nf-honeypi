#!/usr/bin/env python3
"""
make_samplesheet.py — build a honeypi samplesheet from a directory of FASTQ files.

Usage:
    python3 make_samplesheet.py /path/to/fastq_dir [--out samplesheet.csv] [--fix-ids]

Pairing rules (in order of preference):
    1. *_R1*.fastq.gz / *_R2*.fastq.gz
    2. *_R1*.fq.gz   / *_R2*.fq.gz
    3. *_1.fastq.gz  / *_2.fastq.gz
    4. *_1.fq.gz     / *_2.fq.gz

Sample ID = file basename with the R1 suffix stripped, then any trailing _ removed.
Underscores in sample IDs are converted to hyphens unless --keep-underscores is given
(honeypi requires hyphen-only IDs).
"""

import argparse
import re
import sys
from pathlib import Path


_R1_PATTERNS = [
    (r'(.+)_R1(_\S*)?\.fastq\.gz$', '_R2'),
    (r'(.+)_R1(_\S*)?\.fq\.gz$',    '_R2'),
    (r'(.+)_1\.fastq\.gz$',          '_2'),
    (r'(.+)_1\.fq\.gz$',             '_2'),
]


def find_pairs(directory: Path):
    pairs = {}
    for f in sorted(directory.iterdir()):
        if not f.is_file():
            continue
        name = f.name
        for pattern, r2_token in _R1_PATTERNS:
            m = re.match(pattern, name, re.IGNORECASE)
            if m:
                sample_raw = m.group(1)
                # Find the R2 partner
                r2_candidates = []
                for f2 in directory.iterdir():
                    n2 = f2.name
                    # Replace R1 with R2 in the original filename
                    r2_name = re.sub(r'(?i)_R1', r2_token, name, count=1)
                    r2_name2 = re.sub(r'(?i)_1\.', '_2.', name, count=1)
                    if n2 in (r2_name, r2_name2):
                        r2_candidates.append(f2)
                if r2_candidates:
                    pairs[sample_raw] = (f, r2_candidates[0])
                break
    return pairs


def sanitise_id(sample_id: str, fix: bool) -> str:
    sid = sample_id
    # Strip trailing underscores/dashes
    sid = sid.strip('_-')
    if fix:
        sid = sid.replace('_', '-')
    return sid


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('directory', help='Directory containing paired FASTQ files')
    p.add_argument('--out', default='samplesheet.csv', help='Output CSV path (default: samplesheet.csv)')
    p.add_argument('--fix-ids', action='store_true', default=True,
                   help='Replace underscores with hyphens in sample IDs (default: on)')
    p.add_argument('--keep-underscores', action='store_true',
                   help='Keep underscores in sample IDs (overrides --fix-ids)')
    args = p.parse_args()

    directory = Path(args.directory).resolve()
    if not directory.is_dir():
        sys.exit(f"ERROR: Not a directory: {directory}")

    fix = args.fix_ids and not args.keep_underscores
    pairs = find_pairs(directory)

    if not pairs:
        sys.exit(
            f"ERROR: No paired FASTQ files found in {directory}.\n"
            "Expected naming: *_R1*.fastq.gz / *_R2*.fastq.gz (or _1/_2 variants)."
        )

    out_path = Path(args.out)
    warn_underscores = []
    rows = []
    for sample_raw, (r1, r2) in sorted(pairs.items()):
        sid = sanitise_id(sample_raw, fix)
        if '_' in sid:
            warn_underscores.append(sid)
        rows.append((sid, r1, r2))

    with open(out_path, 'w') as fh:
        fh.write('sample,fastq_1,fastq_2\n')
        for sid, r1, r2 in rows:
            fh.write(f'{sid},{r1},{r2}\n')

    print(f"Written {len(rows)} samples to {out_path}")

    if warn_underscores:
        print(f"\nWARNING: {len(warn_underscores)} sample IDs still contain underscores "
              f"(use --fix-ids to auto-convert):")
        for sid in warn_underscores:
            print(f"  {sid}")
        print("honeypi will reject these IDs. Rename files or pass --fix-ids.")


if __name__ == '__main__':
    main()
