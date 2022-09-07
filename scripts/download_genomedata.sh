#!/bin/bash -e
cmdname=`basename $0`
pwd=`pwd`
function usage()
{
    echo "$cmdname <build> <outputdir>" 1>&2
    echo "  build (Ensembl|UCSC, you can specify either):" 1>&2
    echo "         human (GRCh38|hg38, GRCh37|hg19, T2T)" 1>&2
    echo "         mouse (GRCm39|mm39, GRCm38|mm10)" 1>&2
    echo "         rat (mRatBN7.2|rn7)" 1>&2
    echo "         fly (BDGP6|dm6)" 1>&2
    echo "         zebrafish (GRCz11|danRer11)" 1>&2
    echo "         chicken (GRCg6a|galGal6)" 1>&2
    echo "         African clawed frog (Xenopus_tropicalis|xenLae2)" 1>&2
    echo "         C. elegans (WBcel235|ce11)" 1>&2
    echo "         S. cerevisiae (R64-1-1|sacCer3)" 1>&2
    echo "         S. pombe (SPombe)" 1>&2
    echo "  Example:" 1>&2
    echo "         $cmdname GRCh38 Ensembl-GRCh38" 1>&2
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

echo "Start downloading. Selected genome build: $build"

wget="wget -nv --timestamping"

download_mappability(){
    label=$1
    for k in 28 36 50
    do
        wget -q https://www.nakatolab.iqb.u-tokyo.ac.jp/DockerDatabase/mappability/${label}_mappability_Mosaics_${k}mer.tar.bz2
        tar xvfj ${label}_mappability_Mosaics_${k}mer.tar.bz2 >& /dev/null
        rm ${label}_mappability_Mosaics_${k}mer.tar.bz2
    done
}

download_genome2bit(){
    build=$1
    url=https://hgdownload.soe.ucsc.edu/goldenPath/$build/bigZips
    genome=$url/$build.2bit
    gt=$url/$build.chrom.sizes
    ex "wget -nv $genome -O genome_full.2bit"
    ex "wget -nv $gt    -O genometable_full.txt"
}

mkdir -p $outputdir && cd $_
Ensembl_version=106

if test $build = "GRCh38" -o $build = "hg38"; then
    download_genome2bit hg38
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gtf/homo_sapiens/Homo_sapiens.GRCh38.$Ensembl_version.chr.gtf.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gff3/homo_sapiens/Homo_sapiens.GRCh38.$Ensembl_version.chr.gff3.gz"
#    ex "$wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz"
#    ex "$wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gff3.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    wget -q https://www.nakatolab.iqb.u-tokyo.ac.jp/DockerDatabase/RepeatMasker/hg38.txt.gz -O RepeatMasker.txt.gz
    download_mappability Ensembl-GRCh38
    chrs="$(seq 1 22) X Y M"
elif test $build = "GRCh37" -o $build = "hg19"; then
    download_genome2bit hg19
    ex "$wget http://ftp.ensembl.org/pub/grch37/release-$Ensembl_version/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz"
    ex "$wget http://ftp.ensembl.org/pub/grch37/release-$Ensembl_version/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gff3.gz"
    ex "$wget http://ftp.ensembl.org/pub/grch37/release-$Ensembl_version/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz"
    ex "$wget http://ftp.ensembl.org/pub/grch37/release-$Ensembl_version/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh37.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    wget -q https://www.nakatolab.iqb.u-tokyo.ac.jp/DockerDatabase/RepeatMasker/hg19.txt.gz -O RepeatMasker.txt.gz
    download_mappability Ensembl-GRCh37
    chrs="$(seq 1 22) X Y M"
elif test $build = "T2T"; then
    # https://genomeinformatics.github.io/CHM13v2/
    # https://github.com/marbl/chm13
    ex "wget -nv https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -O genome.fa.gz"
    ex "$wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.cen_mask.bed"
    ex "unpigz -f genome.fa.gz"
    zcat /opt/T2Tdata/chm13v2.gtf.gz > chm13v2.gtf
    zcat /opt/T2Tdata/chm13v2.refFlat.gz > chm13v2.refFlat
    download_mappability T2T
    chrs="$(seq 1 22) X Y M"
elif test $build = "GRCm39" -o $build = "mm39"; then
    download_genome2bit mm39
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gtf/mus_musculus/Mus_musculus.GRCm39.$Ensembl_version.chr.gtf.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gff3/mus_musculus/Mus_musculus.GRCm39.$Ensembl_version.chr.gff3.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    download_mappability Ensembl-GRCm39
    chrs="$(seq 1 19) X Y M"
elif test $build = "GRCm38" -o $build = "mm10"; then
    download_genome2bit mm10
    ex "$wget http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.chr.gtf.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-102/gff3/mus_musculus/Mus_musculus.GRCm38.102.chr.gff3.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    wget -q https://www.nakatolab.iqb.u-tokyo.ac.jp/DockerDatabase/RepeatMasker/mm10.txt.gz -O RepeatMasker.txt.gz
    download_mappability Ensembl-GRCm38
    chrs="$(seq 1 19) X Y M"
elif test $build = "mRatBN7.2" -o $build = "rn7"; then
    download_genome2bit rn7
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.$Ensembl_version.chr.gtf.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gff3/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.$Ensembl_version.chr.gff3.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/rattus_norvegicus/ncrna/Rattus_norvegicus.mRatBN7.2.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    download_mappability Ensembl-mRatBN7.2
    chrs="$(seq 1 20) X Y M"
elif test $build = "GRCz11" -o $build = "danRer11"; then
    download_genome2bit danRer11
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gtf/danio_rerio/Danio_rerio.GRCz11.$Ensembl_version.chr.gtf.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gff3/danio_rerio/Danio_rerio.GRCz11.$Ensembl_version.chr.gff3.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/danio_rerio/ncrna/Danio_rerio.GRCz11.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    download_mappability Ensembl-GRCz11
    chrs="$(seq 1 25) M"
elif test $build = "GRCg6a" -o $build = "galGal6"; then
    download_genome2bit galGal6
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gtf/gallus_gallus/Gallus_gallus.GRCg6a.$Ensembl_version.chr.gtf.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gff3/gallus_gallus/Gallus_gallus.GRCg6a.$Ensembl_version.chr.gff3.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/gallus_gallus/cdna/Gallus_gallus.GRCg6a.cdna.all.fa.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/gallus_gallus/ncrna/Gallus_gallus.GRCg6a.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    download_mappability Ensembl-GRCg6a
    chrs="$(seq 1 28) $(seq 30 33) W Z M"
elif test $build = "Xenopus_tropicalis" -o $build = "xenLae2"; then
    download_genome2bit xenLae2
    ex "$wget https://hgdownload.soe.ucsc.edu/goldenPath/xenLae2/bigZips/refMrna.fa.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gtf/xenopus_tropicalis/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.$Ensembl_version.chr.gtf.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gff3/xenopus_tropicalis/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.$Ensembl_version.chr.gff3.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/xenopus_tropicalis/cdna/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.cdna.all.fa.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/xenopus_tropicalis/ncrna/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    download_mappability Ensembl-xenLae2
    chrs="1L 1S 2L 2S 3L 3S 4L 4S 5L 5S 6L 6S 7L 7S 8L 8S 9_10L 9_10S M"
elif test $build = "BDGP6" -o $build = "dm6"; then
    download_genome2bit dm6
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.$Ensembl_version.chr.gtf.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.$Ensembl_version.chr.gff3.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.32.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    download_mappability Ensembl-BDGP6
    chrs="2L 2R 3L 3R 4 X Y M"
elif test $build = "WBcel235" -o $build = "ce11"; then
    download_genome2bit ce11
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.$Ensembl_version.gtf.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.$Ensembl_version.gff3.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/caenorhabditis_elegans/ncrna/Caenorhabditis_elegans.WBcel235.ncrna.fa.gz"
    ex "unpigz -f *gtf.gz *gff3.gz"
    download_mappability Ensembl-WBcel235
    chrs="I II III IV V X M"
elif test $build = "R64-1-1" -o $build = "sacCer3"; then
    download_genome2bit sacCer3
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.$Ensembl_version.gtf.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.$Ensembl_version.gff3.gz"
    ex "$wget https://sgd-prod-upload.s3.amazonaws.com/S000212118/SGD_features.tab.20110827.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
    ex "$wget http://ftp.ensembl.org/pub/release-$Ensembl_version/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz"
    ex "unpigz -f SGD_features.tab.20110827.gz *gtf.gz *gff3.gz"
    ex "cp /opt/OriDB/S.cerevisiae.Oridb.txt ."
    download_mappability Ensembl-R64-1-1
    chrs="I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI M"
elif test $build = "SPombe"; then
    ex "wget -nv https://www.pombase.org/data/genome_sequence_and_features/genome_sequence/Schizosaccharomyces_pombe_all_chromosomes.fa.gz -O genome_full.fa.gz"
    ex "unpigz -f genome_full.fa.gz"
    ex "sed -i -e 's/>I/>chrI/g' genome_full.fa"
    ex "sed -i -e 's/mitochondrial/chrM/g' genome_full.fa"
    ex "makegenometable.pl genome_full.fa > genometable_full.txt"
    ex "faToTwoBit genome_full.fa genome_full.2bit"
    ex "$wget https://www.pombase.org/data/genome_sequence_and_features/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz"
    ex "$wget https://www.pombase.org/data/genome_sequence_and_features/feature_sequences/cds.fa.gz"
    ex "unpigz -f *.gff3.gz"
    ex "sed -i -e 's/mitochondrial/chrM/g' Schizosaccharomyces_pombe_all_chromosomes.gff3"
    cat Schizosaccharomyces_pombe_all_chromosomes.gff3 | awk '{if ($1=="I" || $1=="II"|| $1=="III"){ print "chr"$0}}' > Schizosaccharomyces_pombe_all_chromosomes.gff3.temp
    mv Schizosaccharomyces_pombe_all_chromosomes.gff3.temp Schizosaccharomyces_pombe_all_chromosomes.gff3
    ex "gffread Schizosaccharomyces_pombe_all_chromosomes.gff3 -T -o Schizosaccharomyces_pombe_all_chromosomes.gtf"
    ex "cp /opt/OriDB/S.pombe.Oridb.txt ."
    download_mappability SPombe
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

echo -en "estimate GC contents..."
for chr in $chrs
do
    fa=chromosomes/chr$chr.fa
    s="$s $fa"
    for bin in 100 1000 10000 25000 50000 100000 500000 1000000; do
        GCcount $fa $bin > GCcontents/chr$chr-bs$bin
    done
done
echo "done."

ex "cat $s > genome.fa"
ex "makegenometable.pl genome.fa > genometable.txt"
ex "faToTwoBit genome.fa genome.2bit"
ex "samtools faidx genome.fa"

if test $build != "SPombe" -a $build != "xenLae2" -a $build != "T2T"; then
    ex "zcat *.cdna.all.fa.gz *.ncrna.fa.gz > rna.fa"
fi

if test $build = "T2T"; then
    mkdir -p genedensity
    ex "makegenedensity.pl genometable.txt chm13v2.refFlat 500000"
    mv chr*-bs500000 genedensity

    ex "mkdir -p gtf_chrUCSC"
    ex "cp chm13v2.gtf gtf_chrUCSC/chr.gtf"
    head=gtf_chrUCSC/chr
    gtf2refFlat -g $head.gtf > $head.transcript.refFlat
    gtf2refFlat -u -g $head.gtf > $head.gene.refFlat
    cat $head.gene.refFlat       | awk 'BEGIN { OFS="\t" } {if($4=="+") {print $3, $5, $5, $1}  else {print $3, $6, $6, $1} }'  | uniq | grep -v chrom > $head.gene.TSS.bed
    cat $head.transcript.refFlat | awk 'BEGIN { OFS="\t" } {if($4=="+") {print $3, $5, $5, $14} else {print $3, $6, $6, $14} }' | uniq | grep -v chrom > $head.transcript.TSS.bed
    cat $head.gene.refFlat       | awk 'BEGIN { OFS="\t" } {if($4=="+") {print $3, $6, $6, $1}  else {print $3, $5, $5, $1} }'  | uniq | grep -v chrom > $head.gene.TES.bed
    cat $head.transcript.refFlat | awk 'BEGIN { OFS="\t" } {if($4=="+") {print $3, $5, $5, $14} else {print $3, $5, $5, $14} }' | uniq | grep -v chrom > $head.transcript.TES.bed
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

    mkdir -p $dir/genedensity
    ex "makegenedensity.pl genometable.txt $head.gene.refFlat 500000"
    mv chr*-bs500000 $dir/genedensity
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
