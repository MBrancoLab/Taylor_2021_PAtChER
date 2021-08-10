#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=48:0:0
#$ -l h_vmem=16G

module load use.own
module load samtools-1.12

prefix=P2102EP
rmasker=GRCh38/Repeatmasker_hg38_4.0.5.bed

##merge bam files
samtools merge $prefix.bam *unpaired.bam

##sort
samtools sort -o $prefix-sorted.bam $prefix.bam

##subset
samtools view -b -s 1.5 $prefix-sorted.bam > $prefix-half.bam
samtools view -b -s 2.5 $prefix-half.bam > $prefix-quarter.bam

##index
samtools index $prefix-sorted.bam
samtools index $prefix-half.bam
samtools index $prefix-quarter.bam

##get unique reads
samtools view -h $prefix-sorted.bam | grep '^@\|:u\|:s' - | samtools view -b - > $prefix-unique.bam
samtools view -h $prefix-half.bam | grep '^@\|:u\|:s' - | samtools view -b - > $prefix-halfU.bam
samtools view -h $prefix-quarter.bam | grep '^@\|:u\|:s' - | samtools view -b - > $prefix-quarterU.bam

##index unique
samtools index $prefix-unique.bam
samtools index $prefix-halfU.bam
samtools index $prefix-quarterU.bam

##get coverage
samtools bedcov -d 1 ~/genomes/$rmasker $prefix-sorted.bam > $prefix-allcov.txt
samtools bedcov -d 1 ~/genomes/$rmasker $prefix-unique.bam > $prefix-unicov.txt
samtools bedcov -d 1 ~/genomes/$rmasker $prefix-half.bam > $prefix-halfA.txt
samtools bedcov -d 1 ~/genomes/$rmasker $prefix-halfU.bam > $prefix-halfU.txt
samtools bedcov -d 1 ~/genomes/$rmasker $prefix-quarter.bam > $prefix-quarterA.txt
samtools bedcov -d 1 ~/genomes/$rmasker $prefix-quarterU.bam > $prefix-quarterU.txt

##make bigwig
module load python/3.6.3
module load localtools
bamCoverage --bam $prefix-sorted.bam -o $prefix-all.bw --binSize 200 --normalizeUsing CPM
bamCoverage --bam $prefix-unique.bam -o $prefix-uni.bw --binSize 200 --normalizeUsing CPM
