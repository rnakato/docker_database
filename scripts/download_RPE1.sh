#!/bin/bash
cmdname=`basename $0`
pwd=`pwd`
function usage()
{
    echo "Data from the Giunta Lab (https://github.com/GiuntaLab/RPE1)" 1>&2
    echo "" 1>&2
    echo "Usage:   $cmdname <hap1|hap2> <outputdir>" 1>&2
    echo "Example: $cmdname hap1 RPE1_hap1_data" 1>&2
}

build=$1
outputdir=$2

if [ $# -ne 2 ]; then
  usage
  exit 1
fi

ex(){
    echo $1
    eval $1
}

if test $build != "hap1" -a $build != "hap2"; then
    echo "Specify the correct build."
    usage
    exit 1
fi

echo "Start downloading. Selected genome build: RPE1 $build"


wget="wget --timestamping"

mkdir -p $outputdir && cd $_
mkdir -p chromosomes GCcontents gtf_chrUCSC

URL=https://nakatolab.iqb.u-tokyo.ac.jp/Datafolder_for_sharing/DockerDatabase/RPE1/$build/

# genome
ex "$wget $URL/genome.fa.gz"
ex "faToTwoBit genome.fa.gz genome.fa.2bit"
ex "unpigz genome.fa.gz"
ex "makegenometable.pl genome.fa > genometable.txt"
ex "splitmultifasta genome.fa --dir chromosomes"
ex "samtools faidx genome.fa"

# gtf
ex "$wget $URL/chr.gtf.gz"
ex "unpigz chr.gtf.gz"
ex "mv chr.gtf gtf_chrUCSC/"
ex "$wget $URL/chr.gene.refFlat.gz"
ex "unpigz chr.gene.refFlat.gz"
ex "mv chr.gene.refFlat gtf_chrUCSC/"

echo -en "estimate GC contents..."
for chr in $(seq 1 22) X #Y M
do
    fa=chromosomes/chr$chr.fa
    s="$s $fa"
    for bin in 100 1000 10000 25000 50000 100000 500000 1000000; do
        GCcount $fa $bin > GCcontents/chr$chr-bs$bin
    done
done
echo "done."

head=gtf_chrUCSC/chr
#gtf2refFlat -g $head.gtf > $head.transcript.refFlat
#gtf2refFlat -u -g $head.gtf > $head.gene.refFlat
cat $head.gene.refFlat       | awk 'BEGIN { OFS="\t" } {if($4=="+") {print $3, $5, $5, $1}  else {print $3, $6, $6, $1} }'  | uniq | grep -v chrom > $head.gene.TSS.bed
#cat $head.transcript.refFlat | awk 'BEGIN { OFS="\t" } {if($4=="+") {print $3, $5, $5, $14} else {print $3, $6, $6, $14} }' | uniq | grep -v chrom > $head.transcript.TSS.bed
cat $head.gene.refFlat       | awk 'BEGIN { OFS="\t" } {if($4=="+") {print $3, $6, $6, $1}  else {print $3, $5, $5, $1} }'  | uniq | grep -v chrom > $head.gene.TES.bed
#cat $head.transcript.refFlat | awk 'BEGIN { OFS="\t" } {if($4=="+") {print $3, $5, $5, $14} else {print $3, $5, $5, $14} }' | uniq | grep -v chrom > $head.transcript.TES.bed

mkdir -p gtf_chrUCSC/genedensity
ex "makegenedensity.pl genometable.txt $head.gene.refFlat 500000"
mv chr*-bs500000 gtf_chrUCSC/genedensity

# mappability
url=https://nakatolab.iqb.u-tokyo.ac.jp/Datafolder_for_sharing/DockerDatabase/mappability
for k in 28 36 50
do
    wget -q $url/RPE1_${build}_mappability_Mosaics_${k}mer.tar.bz2
    tar xvfj RPE1_${build}_mappability_Mosaics_${k}mer.tar.bz2 >& /dev/null
    rm RPE1_${build}_mappability_Mosaics_${k}mer.tar.bz2
done
