#!/usr/bin/env bash
set -euo pipefail

gtf=${1:?Usage: parseGtftorefFlat.sh genes.gtf[.gz] prefix}
prefix=${2:-genes}

gp="${prefix}.genePredExt"
meta="${prefix}.transcript_meta.tsv"

# 1. GTF -> genePredExt
# name  = transcript_id
# name2 = gene_name
zcat -f "$gtf" \
| gtfToGenePred \
    -genePredExt \
    -geneNameAsName2 \
    -ignoreGroupsWithoutExons \
    stdin \
    "$gp"

# 2. transcript-level standard refFlat
awk 'BEGIN{OFS="\t"} {
    print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10
}' "$gp" > "${prefix}.transcript.refFlat"

# transcript_id_to_name
gawk -F'\t' '
BEGIN { OFS="\t" }
$0 !~ /^#/ {
    if (match($9, /transcript_id "([^"]+)"/, tid) &&
        match($9, /transcript_name "([^"]+)"/, tname)) {
        id = tid[1]
        name = tname[1]
        if (!seen[id]++) {
            print id, name
        }
    }
}
' $gtf > $prefix.transcript_id_to_name.tsv

# 3. transcript feature だけから metadata と代表 transcript 用 score を作る
zcat -f "$gtf" \
| awk -F'\t' '
BEGIN{OFS="\t"}

function attr(s, key,    m) {
    if (match(s, "(^|;[[:space:]]*)" key "[[:space:]]+\"([^\"]*)\"", m)) return m[2]
    if (match(s, "(^|;[[:space:]]*)" key "=([^;]*)", m)) return m[2]
    return ""
}

function has_tag(s, tag,    pat) {
    pat = "(^|;[[:space:]]*)tag[[:space:]]+\"" tag "\""
    return (s ~ pat)
}

$0 !~ /^#/ && NF >= 9 && $3 == "transcript" {
    tid = attr($9, "transcript_id")
    if (tid == "") next

    gid   = attr($9, "gene_id")
    gname = attr($9, "gene_name")
    tname = attr($9, "transcript_name")

    gtype = attr($9, "gene_type")
    if (gtype == "") gtype = attr($9, "gene_biotype")

    ttype = attr($9, "transcript_type")
    if (ttype == "") ttype = attr($9, "transcript_biotype")

    score = 0

    # 代表 transcript の優先順位
    if (has_tag($9, "MANE_Select")) score += 100000000
    if (has_tag($9, "Ensembl_canonical")) score += 50000000
    if ($9 ~ /tag "appris_principal/) score += 10000000
    if (has_tag($9, "CCDS")) score += 5000000
    if (has_tag($9, "basic")) score += 1000000
    if (ttype == "protein_coding") score += 100000

    if (gid == "") gid = gname
    if (gname == "") gname = gid
    if (tname == "") tname = tid
    if (gtype == "") gtype = "NA"
    if (ttype == "") ttype = "NA"

    print tid,gid,gname,tname,gtype,ttype,score
}
' > "$meta"

# 4. gene_id_to_name.tsv
{
    echo -e "gene_id\tgene_name"
    awk -F'\t' 'BEGIN{OFS="\t"} $2 != "" {print $2,$3}' "$meta" \
    | LC_ALL=C sort -u
} > "${prefix}.gene_id_to_name.tsv"

# 5. representative transcript を選ぶ
# tie-breaker として exon length を足す
awk -F'\t' '
BEGIN{OFS="\t"}

ARGIND == 1 {
    gid[$1]   = $2
    gname[$1] = $3
    tname[$1] = $4
    gtype[$1] = $5
    ttype[$1] = $6
    baseScore[$1] = $7
    next
}

function exon_len(starts, ends,    s,e,n,i,len) {
    n = split(starts, s, ",")
    split(ends, e, ",")
    len = 0
    for (i = 1; i <= n; i++) {
        if (s[i] != "" && e[i] != "") len += e[i] - s[i]
    }
    return len
}

ARGIND == 2 {
    tid = $1
    if (!(tid in gid)) next

    gene = gid[tid]
    if (gene == "") next

    len = exon_len($9, $10)
    score = baseScore[tid] + len

    if (!(gene in best) || score > bestScore[gene]) {
        best[gene] = tid
        bestScore[gene] = score

        line[gene] = gname[tid] OFS gene OFS $2 OFS $3 OFS \
                     $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS \
                     gtype[tid] OFS ttype[tid] OFS tname[tid] OFS tid
    }
}

END {
    for (gene in line) print line[gene]
}
' "$meta" "$gp" \
| LC_ALL=C sort -k3,3 -k5,5n \
> "${prefix}.gene.extended.refFlat.tsv.tmp"

# 6. header 付き gene-level extended refFlat
{
    echo -e "geneName\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tgene_type\ttranscript_type\treference_transcript_name\treference_transcript_id"
    cat "${prefix}.gene.extended.refFlat.tsv.tmp"
} > "${prefix}.gene.extended.refFlat.tsv"

# 7. 11列版 gene refFlat-like
awk -F'\t' 'BEGIN{OFS="\t"} NR > 1 {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
}' "${prefix}.gene.extended.refFlat.tsv" > "${prefix}.gene.refFlat"

rm -f "${prefix}.gene.extended.refFlat.tsv.tmp"

echo "[done]"
echo "  ${prefix}.transcript.refFlat"
echo "  ${prefix}.gene.refFlat"
echo "  ${prefix}.gene.extended.refFlat.tsv"
echo "  ${prefix}.gene_id_to_name.tsv"


# 8. protein_coding gene only: gene-level refFlat-like files
awk -F'\t' 'BEGIN{OFS="\t"}
NR == 1 {
    print $0
    next
}
$12 == "protein_coding" {
    print $0
}
' "${prefix}.gene.extended.refFlat.tsv" \
> "${prefix}.protein_coding.gene.extended.refFlat.tsv"

awk -F'\t' 'BEGIN{OFS="\t"}
NR > 1 {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
}
' "${prefix}.protein_coding.gene.extended.refFlat.tsv" \
> "${prefix}.protein_coding.gene.refFlat"


# 9. protein_coding gene only: transcript-level refFlat-like files
# Definition:
#   include all transcripts whose parent gene_type is protein_coding.
gawk -F'\t' '
BEGIN{OFS="\t"}

ARGIND == 1 {
    tid = $1
    gid = $2
    gname = $3
    tname = $4
    gtype = $5
    ttype = $6

    if (gtype == "protein_coding") {
        keep[tid] = 1
        TNAME[tid] = tname
        GID[tid] = gid
        GNAME[tid] = gname
        GTYPE[tid] = gtype
        TTYPE[tid] = ttype
    }
    next
}

ARGIND == 2 {
    tid = $1
    if (!(tid in keep)) next

    geneName = $12
    if (geneName == "") geneName = GNAME[tid]
    if (geneName == "") geneName = GID[tid]

    print geneName,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10
}
' "$meta" "$gp" \
> "${prefix}.protein_coding.transcript.refFlat"


{
    echo -e "transcriptName\ttranscriptId\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tgene_type\ttranscript_type\tgene_name\tgene_id"

    gawk -F'\t' '
    BEGIN{OFS="\t"}

    ARGIND == 1 {
        tid = $1
        gid = $2
        gname = $3
        tname = $4
        gtype = $5
        ttype = $6

        if (gtype == "protein_coding") {
            keep[tid] = 1
            TNAME[tid] = tname
            GID[tid] = gid
            GNAME[tid] = gname
            GTYPE[tid] = gtype
            TTYPE[tid] = ttype
        }
        next
    }

    ARGIND == 2 {
        tid = $1
        if (!(tid in keep)) next

        transcriptName = TNAME[tid]
        if (transcriptName == "") transcriptName = tid

        geneName = GNAME[tid]
        if (geneName == "") geneName = GID[tid]

        print transcriptName,tid,$2,$3,$4,$5,$6,$7,$8,$9,$10,
              GTYPE[tid],TTYPE[tid],geneName,GID[tid]
    }
    ' "$meta" "$gp"
} > "${prefix}.protein_coding.transcript.extended.refFlat.tsv"

rm $gp
