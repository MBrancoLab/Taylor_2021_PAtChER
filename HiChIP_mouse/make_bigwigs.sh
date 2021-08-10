#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=24:0:0
#$ -l h_vmem=16G
#$ -t 1-6

module load samtools

prefix=$(sed -n "${SGE_TASK_ID}p" prefix_list.txt)

##merge bam files
samtools merge $prefix.bam ${prefix}*unpaired.bam

##sort
samtools sort -o $prefix-sorted.bam $prefix.bam
samtools index $prefix-sorted.bam

##get unique reads
samtools view -h $prefix-sorted.bam | grep '^@\|:u\|:s' - | samtools view -b - > $prefix-unique.bam
samtools index $prefix-unique.bam

##make bigwig
module load python/3.6.3
module load use.own
module load localtools
bamCoverage --bam $prefix-sorted.bam -o $prefix-all.bw --binSize 200 --normalizeUsing CPM
bamCoverage --bam $prefix-unique.bam -o $prefix-uni.bw --binSize 200 --normalizeUsing CPM
