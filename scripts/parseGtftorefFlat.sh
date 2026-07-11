#!/usr/bin/env bash
set -euo pipefail

gtf=${1:?Usage: parseGtftorefFlat.sh genes.gtf[.gz] prefix}
prefix=${2:-chr}

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

# transcript_id_to_name
# Works with both .gtf and .gtf.gz.
zcat -f "$gtf" \
| gawk -F'	' '
BEGIN { OFS="	" }
$0 !~ /^#/ && NF >= 9 {
    if (match($9, /transcript_id "([^"]+)"/, tid) &&
        match($9, /transcript_name "([^"]+)"/, tname)) {
        id = tid[1]
        name = tname[1]
        if (!seen[id]++) {
            print id, name
        }
    }
}
' > "${prefix}.transcript_id_to_name.tsv"

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
gawk -F'\t' '
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
| LC_ALL=C sort -t$'\t' -k3,3 -k5,5n \
> "${prefix}.gene.refFlat.tmp"


# 6. header 付き gene-level extended refFlat
{
    echo -e "geneName\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tgene_type\ttranscript_type\treference_transcript_name\treference_transcript_id"
    cat "${prefix}.gene.refFlat.tmp"
} > "${prefix}.gene.refFlat"

# 7. 11列版 gene refFlat-like
awk -F'\t' 'BEGIN{OFS="\t"} NR > 1 {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
}' "${prefix}.gene.refFlat" > "${prefix}.gene.simple.refFlat"

rm -f "${prefix}.gene.refFlat.tmp"

# 8. protein_coding gene only: gene-level refFlat-like files
awk -F'\t' 'BEGIN{OFS="\t"}
NR == 1 {
    print $0
    next
}
$12 == "protein_coding" {
    print $0
}
' "${prefix}.gene.refFlat" \
> "${prefix}.protein_coding.gene.refFlat"

awk -F'\t' 'BEGIN{OFS="\t"}
NR > 1 {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
}
' "${prefix}.protein_coding.gene.refFlat" \
> "${prefix}.protein_coding.gene.simple.refFlat"


# 9. transcript-level refFlat-like files
# Outputs:
#   ${prefix}.transcript.refFlat
#   ${prefix}.transcript.simple.refFlat
#   ${prefix}.protein_coding.transcript.refFlat
#   ${prefix}.protein_coding.transcript.simple.refFlat
#
# Definition of protein_coding transcript set:
#   include transcripts whose parent gene_type is protein_coding.
#   If you want transcript_type == protein_coding as well, change the condition below.
out_transcript_ext="${prefix}.transcript.refFlat"
out_transcript_ref="${prefix}.transcript.simple.refFlat"
out_pc_transcript_ext="${prefix}.protein_coding.transcript.refFlat"
out_pc_transcript_ref="${prefix}.protein_coding.transcript.simple.refFlat"

# Avoid appending to old files when the script is rerun.
: > "$out_transcript_ref"
: > "$out_pc_transcript_ref"

gawk -F'\t' \
    -v out_ext="$out_transcript_ext" \
    -v out_ref="$out_transcript_ref" \
    -v out_pc_ext="$out_pc_transcript_ext" \
    -v out_pc_ref="$out_pc_transcript_ref" '
BEGIN {
    OFS="\t"

    print "transcriptName","transcriptId","chrom","strand", \
          "txStart","txEnd","cdsStart","cdsEnd","exonCount", \
          "exonStarts","exonEnds","gene_type","transcript_type", \
          "gene_name","gene_id" > out_ext

    print "transcriptName","transcriptId","chrom","strand", \
          "txStart","txEnd","cdsStart","cdsEnd","exonCount", \
          "exonStarts","exonEnds","gene_type","transcript_type", \
          "gene_name","gene_id" > out_pc_ext
}

ARGIND == 1 {
    tid = $1
    GID[tid]   = $2
    GNAME[tid] = $3
    TNAME[tid] = $4
    GTYPE[tid] = $5
    TTYPE[tid] = $6
    next
}

ARGIND == 2 {
    tid = $1

    chrom      = $2
    strand     = $3
    txStart    = $4
    txEnd      = $5
    cdsStart   = $6
    cdsEnd     = $7
    exonCount  = $8
    exonStarts = $9
    exonEnds   = $10

    geneName = GNAME[tid]
    if (geneName == "") geneName = $12
    if (geneName == "") geneName = GID[tid]
    if (geneName == "") geneName = tid

    geneId = GID[tid]
    if (geneId == "") geneId = geneName

    transcriptName = TNAME[tid]
    if (transcriptName == "") transcriptName = tid

    geneType = GTYPE[tid]
    if (geneType == "") geneType = "NA"

    transcriptType = TTYPE[tid]
    if (transcriptType == "") transcriptType = "NA"

    # 11-column standard refFlat: no header.
    print geneName,tid,chrom,strand,txStart,txEnd,cdsStart,cdsEnd, \
          exonCount,exonStarts,exonEnds >> out_ref

    # 15-column extended transcript table: header is written above.
    print transcriptName,tid,chrom,strand,txStart,txEnd,cdsStart,cdsEnd, \
          exonCount,exonStarts,exonEnds,geneType,transcriptType, \
          geneName,geneId >> out_ext

    if (geneType == "protein_coding") {
        print geneName,tid,chrom,strand,txStart,txEnd,cdsStart,cdsEnd, \
              exonCount,exonStarts,exonEnds >> out_pc_ref

        print transcriptName,tid,chrom,strand,txStart,txEnd,cdsStart,cdsEnd, \
              exonCount,exonStarts,exonEnds,geneType,transcriptType, \
              geneName,geneId >> out_pc_ext
    }
}
' "$meta" "$gp"

# Basic column checks.
# Basic column checks. All output files are tab-separated.
awk -F'\t' '
NF != 11 {
    print "Error:", FILENAME, "line", NR, "has", NF, "tab-separated columns"
    exit 1
}
' "$out_transcript_ref"

awk -F'\t' '
NF != 11 {
    print "Error:", FILENAME, "line", NR, "has", NF, "tab-separated columns"
    exit 1
}
' "$out_pc_transcript_ref"

awk -F'\t' '
NR > 1 && NF != 15 {
    print "Error:", FILENAME, "line", NR, "has", NF, "tab-separated columns"
    exit 1
}
' "$out_transcript_ext"

awk -F'\t' '
NR > 1 && NF != 15 {
    print "Error:", FILENAME, "line", NR, "has", NF, "tab-separated columns"
    exit 1
}
' "$out_pc_transcript_ext"

echo "[done]"
echo "  ${prefix}.gene.simple.refFlat"
echo "  ${prefix}.gene.refFlat"
echo "  ${prefix}.protein_coding.gene.refFlat"
echo "  ${prefix}.protein_coding.gene.simple.refFlat"
echo "  ${prefix}.gene_id_to_name.tsv"
echo "  ${prefix}.transcript_id_to_name.tsv"
echo "  $meta"


rm -f "$gp"

echo "  ${prefix}.transcript.simple.refFlat"
echo "  ${prefix}.transcript.refFlat"
echo "  ${prefix}.protein_coding.transcript.refFlat"
echo "  ${prefix}.protein_coding.transcript.simple.refFlat"
