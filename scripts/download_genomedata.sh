#!/bin/bash -e
cmdname=`basename $0`
pwd=`pwd`
function usage()
{
    echo "$cmdname <build> <outputdir>" 1>&2
    echo "  build:" 1>&2
    echo "         human (GRCh38, GRCh37, T2T)" 1>&2
    echo "         mouse (GRCm39, GRCm38)" 1>&2
    echo "         rat (mRatBN7.2)" 1>&2
    echo "         fly (BDGP6)" 1>&2
    echo "         zebrafish (GRCz11)" 1>&2
    echo "         chicken (GRCg6a)" 1>&2
    echo "         African clawed frog (xenLae2)" 1>&2
    echo "         C. elegans (WBcel235)" 1>&2
    echo "         S. cerevisiae (R64-1-1)" 1>&2
    echo "         S. pombe (SPombe)" 1>&2
    echo "  Example:" 1>&2
    echo "         $cmdname GRCh38 Ensembl-GRCh38" 1>&2
}

ncore=4
while getopts p: option
do
    case ${option} in
	p)
            ncore=${OPTARG}
            ;;
	*)
	    usage
	    exit 1
	    ;;
    esac
done
shift $((OPTIND - 1))

build=$1
outprefix=$2

if [ $# -ne 2 ]; then
  usage
  exit 1
fi

ex(){
    echo $1
    eval $1
}

echo "Start downloading. Selected build: $build"

download_mappability(){
    build=$1
    if test $build = "T2T"; then
	label=T2T
    else
	label=Ensembl-${build}
    fi
    for k in 28 36 50
    do
	wget https://www.nakatolab.iqb.u-tokyo.ac.jp/DockerDatabase/mappability/${label}_mappability_Mosaics_${k}mer.tar.bz2
	tar xvfj ${label}_mappability_Mosaics_${k}mer.tar.bz2
	rm ${label}_mappability_Mosaics_${k}mer.tar.bz2
    done
}

mkdir -p $outprefix && cd $_
download_mappability $build

if test $build = "GRCh38"; then
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit -O genome_full.2bit"
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes -O genometable_full.txt"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.chr.gtf.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gff3/homo_sapiens/Homo_sapiens.GRCh38.106.chr.gff3.gz"
#    ex "wget --timestamping http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz"
#    ex "wget --timestamping http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gff3.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    wget --timestamping https://www.nakatolab.iqb.u-tokyo.ac.jp/DockerDatabase/RepeatMasker/hg38.txt.gz

    chrs="$(seq 1 22) X Y M"
elif test $build = "GRCh37"; then
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit -O genome_full.2bit"
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes -O genometable_full.txt"
    ex "wget --timestamping http://ftp.ensembl.org/pub/grch37/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/grch37/release-106/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gff3.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/grch37/release-106/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/grch37/release-106/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh37.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    wget --timestamping https://www.nakatolab.iqb.u-tokyo.ac.jp/DockerDatabase/RepeatMasker/hg19.txt.gz
    chrs="$(seq 1 22) X Y M"
elif test $build = "T2T"; then
    # https://genomeinformatics.github.io/CHM13v2/
    # https://github.com/marbl/chm13
    ex "wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -O genome.fa.gz"
    ex "wget --timestamping https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.cen_mask.bed"
    ex "unpigz -f genome.fa.gz"
    zcat /opt/T2Tdata/chm13v2.gtf.gz > chm13v2.gtf
    zcat /opt/T2Tdata/chm13v2.refFlat.gz > chm13v2.refFlat
    chrs="$(seq 1 22) X Y M"
elif test $build = "GRCm39"; then
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.2bit -O genome_full.2bit"
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes -O genometable_full.txt"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gtf/mus_musculus/Mus_musculus.GRCm39.106.chr.gtf.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gff3/mus_musculus/Mus_musculus.GRCm39.106.chr.gff3.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    chrs="$(seq 1 19) X Y M"
elif test $build = "GRCm38"; then
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit -O genome_full.2bit"
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes -O genometable_full.txt"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.chr.gtf.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-102/gff3/mus_musculus/Mus_musculus.GRCm38.102.chr.gff3.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    wget --timestamping https://www.nakatolab.iqb.u-tokyo.ac.jp/DockerDatabase/RepeatMasker/mm10.txt.gz

    chrs="$(seq 1 19) X Y M"
elif test $build = "mRatBN7.2"; then
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.2bit -O genome_full.2bit"
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.chrom.sizes -O genometable_full.txt"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.106.chr.gtf.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gff3/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.106.chr.gff3.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/rattus_norvegicus/ncrna/Rattus_norvegicus.mRatBN7.2.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    chrs="$(seq 1 20) X Y M"
elif test $build = "GRCz11"; then
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.2bit -O genome_full.2bit"
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.chrom.sizes -O genometable_full.txt"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gtf/danio_rerio/Danio_rerio.GRCz11.106.chr.gtf.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gff3/danio_rerio/Danio_rerio.GRCz11.106.chr.gff3.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/danio_rerio/ncrna/Danio_rerio.GRCz11.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    chrs="$(seq 1 25) M"
elif test $build = "GRCg6a"; then
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/bigZips/galGal6.2bit -O genome_full.2bit"
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/bigZips/galGal6.chrom.sizes -O genometable_full.txt"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gtf/gallus_gallus/Gallus_gallus.GRCg6a.106.chr.gtf.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gff3/gallus_gallus/Gallus_gallus.GRCg6a.106.chr.gff3.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/gallus_gallus/cdna/Gallus_gallus.GRCg6a.cdna.all.fa.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/gallus_gallus/ncrna/Gallus_gallus.GRCg6a.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    chrs="$(seq 1 28) $(seq 30 33) W Z M"
elif test $build = "xenLae2"; then
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/xenLae2/bigZips/xenLae2.2bit -O genome_full.2bit"
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/xenLae2/bigZips/xenLae2.chrom.sizes -O genometable_full.txt"
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/xenLae2/bigZips/refMrna.fa.gz"
    chrs="1L 1S 2L 2S 3L 3S 4L 4S 5L 5S 6L 6S 7L 7S 8L 8S 9_10L 9_10S M"
elif test $build = "BDGP6"; then
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.2bit -O genome_full.2bit"
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes -O genometable_full.txt"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.106.chr.gtf.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.106.chr.gff3.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.32.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    chrs="2L 2R 3L 3R 4 X Y M"
elif test $build = "WBcel235"; then
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.2bit -O genome_full.2bit"
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.chrom.sizes -O genometable_full.txt"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.106.gtf.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.106.gff3.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/caenorhabditis_elegans/ncrna/Caenorhabditis_elegans.WBcel235.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    chrs="I II III IV V X M"
elif test $build = "R64-1-1"; then
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.2bit -O genome_full.2bit"
    ex "wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes -O genometable_full.txt"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.106.gtf.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.106.gff3.gz"
    ex "wget --timestamping https://sgd-prod-upload.s3.amazonaws.com/S000212118/SGD_features.tab.20110827.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
    ex "wget --timestamping http://ftp.ensembl.org/pub/release-106/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz"
    ex "unpigz -f SGD_features.tab.20110827.gz *gtf.gz *gff3.gz"
    ex "cp /opt/OriDB/S.cerevisiae.Oridb.txt ."
chrs="I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI M"
elif test $build = "SPombe"; then
    ex "wget https://www.pombase.org/data/genome_sequence_and_features/genome_sequence/Schizosaccharomyces_pombe_all_chromosomes.fa.gz -O genome_full.fa.gz"
    ex "unpigz -f genome_full.fa.gz"
    ex "sed -i -e 's/>I/>chrI/g' genome_full.fa"
    ex "sed -i -e 's/mitochondrial/chrM/g' genome_full.fa"
    ex "makegenometable.pl genome_full.fa > genometable_full.txt"
    ex "faToTwoBit genome_full.fa genome_full.2bit"
    ex "wget --timestamping https://www.pombase.org/data/genome_sequence_and_features/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz"
    ex "wget https://www.pombase.org/data/genome_sequence_and_features/feature_sequences/cds.fa.gz"
    ex "unpigz -f *.gff3.gz"
    ex "sed -i -e 's/mitochondrial/chrM/g' Schizosaccharomyces_pombe_all_chromosomes.gff3"
    cat Schizosaccharomyces_pombe_all_chromosomes.gff3 | awk '{if ($1=="I" || $1=="II"|| $1=="III"){ print "chr"$0}}' > Schizosaccharomyces_pombe_all_chromosomes.gff3.temp
    mv Schizosaccharomyces_pombe_all_chromosomes.gff3.temp Schizosaccharomyces_pombe_all_chromosomes.gff3
    ex "gffread Schizosaccharomyces_pombe_all_chromosomes.gff3 -T -o Schizosaccharomyces_pombe_all_chromosomes.gtf"
    ex "cp /opt/OriDB/S.pombe.Oridb.txt ."
    chrs="I II III M"
else
    echo "Specify the correct build."
    usage
    exit 1
fi

mkdir -p chromosomes GCcontents

if test $build != "T2T"; then
    ex "twoBitToFa genome_full.2bit genome_full.fa"
    ex "splitmultifasta genome_full.fa --dir chromosomes"
    ex "samtools faidx genome_full.fa"
else
    ex "splitmultifasta genome.fa --dir chromosomes"
fi

for chr in $chrs
do
    fa=chromosomes/chr$chr.fa
    s="$s $fa"
    for bin in 100 1000 10000 25000 50000 100000 500000 1000000; do
	ex "GCcount $fa $bin > GCcontents/chr$chr-bs$bin"
    done
done

ex "cat $s > genome.fa"
ex "makegenometable.pl genome.fa > genometable.txt"
ex "faToTwoBit genome.fa genome.2bit"
ex "samtools faidx genome.fa"

if test $build != "SPombe" -a $build != "xenLae2" -a $build != "T2T"; then
    ex "zcat *.cdna.all.fa.gz *.ncrna.fa.gz > rna.fa"
fi

if test $build = "T2T"; then
    ex "mkdir -p genedensity"
    ex "makegenedensity.pl genometable.txt chm13v2.refFlat 500000"
    ex "mv chr*-bs500000 genedensity"

    ex "mkdir -p gtf_chrUCSC"
    ex "cp chm13v2.gtf gtf_chrUCSC/chr.gtf"
    exit
fi

gtf=`ls *.gtf`
gff=`ls *.gff3`

ex "mkdir -p gtf_original"
ex "mv $gtf gtf_original/chr.gtf"
ex "mv $gff gtf_original/chr.gff3"

convert_gtf_to_refFlat(){
    dir=$1
    ex "mkdir -p $dir"
    ex "extract_proteincoding $dir/chr.gtf > $dir/chr.proteincoding.gtf"
    for head in $dir/chr $dir/chr.proteincoding
    do
        ex "gtf2refFlat -g $head.gtf > $head.transcript.refFlat"
        ex "gtf2refFlat -u -g $head.gtf > $head.gene.refFlat"
        cat $head.gene.refFlat | awk 'BEGIN { OFS="\t" } {if($4=="+") {print $3, $5, $5, $1} else {print $3, $6, $6, $1} }' | uniq | grep -v chrom > $head.gene.TSS.bed
        cat $head.transcript.refFlat | awk 'BEGIN { OFS="\t" } {if($4=="+") {print $3, $5, $5, $14} else {print $3, $6, $6, $14} }' | uniq | grep -v chrom > $head.transcript.TSS.bed
        cat $head.gene.refFlat | awk 'BEGIN { OFS="\t" } {if($4=="+") {print $3, $6, $6, $1} else {print $3, $5, $5, $1} }' | uniq | grep -v chrom > $head.gene.TES.bed
        cat $head.transcript.refFlat | awk 'BEGIN { OFS="\t" } {if($4=="+") {print $3, $5, $5, $14} else {print $3, $5, $5, $14} }' | uniq | grep -v chrom > $head.transcript.TES.bed
    done

    ex "mkdir -p $dir/genedensity"
    ex "makegenedensity.pl genometable.txt $head.gene.refFlat 500000"
    ex "mv chr*-bs500000 $dir/genedensity"
}

if test $build = "SPombe"; then
    for dir in gtf_original; do convert_gtf_to_refFlat $dir; done
elif test $build = "xenLae2"; then
    echo "build $build does not contain gtf files."
else
    ex "mkdir -p gtf_chrUCSC"
    ex "convertchr_fromEns2UCSC gtf_original/chr.gtf > gtf_chrUCSC/chr.gtf"
    if test $build = "R64-1-1"; then
	ex "sed -i -e 's/Mito/M/g' gtf_chrUCSC/chr.gtf"
    fi
    for dir in gtf_original gtf_chrUCSC; do convert_gtf_to_refFlat $dir; done
fi
