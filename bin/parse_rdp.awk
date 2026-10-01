# Reformat RDP classifier output exactly like native honeypi_reformatAssignedTaxonomy.
#   usage: awk -v threshold=0.5 -f parse_rdp.awk assigned_taxonomy_prelim.txt
# Output (no header): ASV_ID <TAB> "k__X; p__Y; c__Z; ..." <TAB> confidence of the last kept rank
# Ranks are kept until the first one below the threshold; none kept -> "Unassignable" / 1.0
BEGIN {
    FS = "\t"; OFS = "\t"
    if (threshold == "") threshold = 0.5
    prefix["domain"]    = "k__"; prefix["kingdom"]   = "k__"
    prefix["phylum"]    = "p__"; prefix["subphylum"] = "subp__"
    prefix["class"]     = "c__"; prefix["subclass"]  = "subc__"
    prefix["order"]     = "o__"; prefix["family"]    = "f__"
    prefix["genus"]     = "g__"; prefix["species"]   = "s__"
}
NF == 0 { next }
{
    taxonomy = ""; kept = 0; last_conf = "1.0"; stop = 0
    # fields: 1 id, 2 orientation, 3-5 root (Root, rootrank, conf), then (name, rank, conf) triplets
    for (i = 6; i + 2 <= NF && !stop; i += 3) {
        name = $i; rank = $(i + 1); conf = $(i + 2)
        if (!(rank in prefix)) {
            print rank, "Error in RDP Classifier produced output: not a valid taxonomic level." > "/dev/stderr"
            exit 1
        }
        if (conf + 0 >= threshold + 0) {
            n = split(name, parts, "|"); name = parts[n]; gsub(/ /, "_", name)
            taxonomy = taxonomy (kept ? "; " : "") prefix[rank] name
            kept++
            # python str(float(conf)): "1" -> "1.0", "0.9400" -> "0.94"
            last_conf = conf
            if (last_conf ~ /^[0-9]+$/) last_conf = last_conf ".0"
            else { sub(/0+$/, "", last_conf); if (last_conf ~ /\.$/) last_conf = last_conf "0" }
        } else stop = 1
    }
    if (kept == 0) { taxonomy = "Unassignable"; last_conf = "1.0" }
    print $1, taxonomy, last_conf
}
