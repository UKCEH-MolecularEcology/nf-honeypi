#!/usr/bin/env python3
"""Merge ASVs with identical sequences - a pure-python port of native honeypi_mergeDuplicateASV.

  merge_duplicate_asvs.py --fasta ASVs.fasta --counts ASVs_counts_filtered.txt \
                          --taxonomy assigned_taxonomy.txt --outdir DIR [--strip-sample-regex REGEX]

Native behaviour reproduced exactly:
  * ASVs with the same sequence are collapsed, counts summed, first ASV ID (FASTA order) kept
  * rows are sorted by sequence (pandas groupby), the FASTA follows the same order
  * ASVs_counts.txt has NO label for the ID column in its header; ASVs_taxonomy.txt has no header
Extra (not native): ASVs_counts_merged.txt = counts summed per identical taxonomy string.
"""
import argparse, os, re, sys
from collections import OrderedDict

ap = argparse.ArgumentParser()
ap.add_argument("--fasta", required=True)
ap.add_argument("--counts", required=True)
ap.add_argument("--taxonomy", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--sample-ids", default=None,
                help="file with the original sample IDs (one per line); used to undo R's make.names() "
                     "mangling in the DADA2 count table header (X-prefix, '-' -> '.')")
ap.add_argument("--strip-sample-regex", default=None,
                help="regex removed from sample names (e.g. '-S[0-9]+-L[0-9]+$')")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

# FASTA: header (whole line after '>') -> sequence, in file order
uid2seq = OrderedDict()
hdr, buf = None, []
with open(a.fasta) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None:
                uid2seq[hdr] = "".join(buf)
            hdr, buf = line[1:].strip(), []
        elif hdr is not None:
            buf.append(line.strip())
if hdr is not None:
    uid2seq[hdr] = "".join(buf)

seq2first = {}
for uid, seq in uid2seq.items():
    seq2first.setdefault(seq, uid)
n_dup = len(uid2seq) - len(seq2first)
print(f"Number of ASVs with a duplicated sequence: {n_dup}", flush=True)

# counts table (first header cell is empty or 'ASV_ID'; header has one cell fewer if unlabelled)
with open(a.counts) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    ncol = None
    rows = []
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) == len(head) + 1:          # unlabelled header (native dada2 style)
            samples = head
        else:
            samples = head[1:]
        rows.append((p[0], [int(float(x)) for x in p[1:]]))
def make_names(x):                            # R's make.names() for the characters seen in sample IDs
    x = re.sub(r"[^A-Za-z0-9._]", ".", x)
    return x if re.match(r"^([A-Za-z]|\.(?![0-9]))", x) else "X" + x

if a.sample_ids:
    ids = [l.strip() for l in open(a.sample_ids) if l.strip()]
    lookup = {make_names(i): i for i in ids}
    known = set(ids)
    unresolved = [s for s in samples if s not in lookup and s not in known]
    if unresolved:
        sys.exit("ERROR: cannot map sample columns to known IDs: " + ", ".join(unresolved[:5]))
    samples = [s if s in known else lookup[s] for s in samples]
if a.strip_sample_regex:
    samples = [re.sub(a.strip_sample_regex, "", s) for s in samples]
if len(set(samples)) != len(samples):
    sys.exit("ERROR: sample names are not unique after renaming: " +
             ", ".join(s for s in set(samples) if samples.count(s) > 1))

by_seq = {}
for uid, vals in rows:
    seq = uid2seq.get(uid)
    if seq is None:                           # pandas drops rows whose ID has no sequence
        continue
    if seq in by_seq:
        by_seq[seq] = [x + y for x, y in zip(by_seq[seq], vals)]
    else:
        by_seq[seq] = list(vals)

ordered = sorted(by_seq)                      # pandas groupby sorts keys
final_ids = [seq2first[s] for s in ordered]

with open(os.path.join(a.outdir, "ASVs_counts.txt"), "w") as fo:
    fo.write("\t".join(samples) + "\n")
    for s, uid in zip(ordered, final_ids):
        fo.write(uid + "\t" + "\t".join(map(str, by_seq[s])) + "\n")

# taxonomy: ID <TAB> taxonomy <TAB> score (no header), reindexed to the merged IDs
tax = {}
with open(a.taxonomy) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            tax[p[0]] = (p[1], p[2] if len(p) > 2 else "")
with open(os.path.join(a.outdir, "ASVs_taxonomy.txt"), "w") as fo:
    for uid in final_ids:
        t, sc = tax.get(uid, ("", ""))
        fo.write(f"{uid}\t{t}\t{sc}\n")

with open(os.path.join(a.outdir, "ASVs.fasta"), "w") as fo:
    for s, uid in zip(ordered, final_ids):
        fo.write(f">{uid}\n{s}\n")

# extra: counts per identical taxonomy string, most abundant first
bytax = OrderedDict()
for s, uid in zip(ordered, final_ids):
    t = tax.get(uid, ("k__unclassified", ""))[0] or "k__unclassified"
    bytax[t] = [x + y for x, y in zip(bytax[t], by_seq[s])] if t in bytax else list(by_seq[s])
with open(os.path.join(a.outdir, "ASVs_counts_merged.txt"), "w") as fo:
    fo.write("\t".join(["taxonomy"] + samples) + "\n")
    for t, v in sorted(bytax.items(), key=lambda kv: -sum(kv[1])):
        fo.write(t + "\t" + "\t".join(map(str, v)) + "\n")

print(f"Merged {len(rows)} ASVs -> {len(final_ids)} unique sequences, {len(bytax)} unique taxonomies", flush=True)
