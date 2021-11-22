#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=12:0:0
#$ -l h_vmem=8G

module load python/3.8.5
module load samtools
module load use.own
module load htseq

base_name=MCF7

##merge and sort bam files
samtools merge ${base_name}_long.bam *-long_unp.bam
samtools merge ${base_name}_prox.bam *-prox_unp.bam
samtools sort -o ${base_name}_long-sort.bam ${base_name}_long.bam
samtools sort -o ${base_name}_prox-sort.bam ${base_name}_prox.bam
samtools index ${base_name}_long-sort.bam
samtools index ${base_name}_prox-sort.bam


##count reads within regions of interest
htseq-count -f bam -s no -a 0 -t Loop -i ID ${base_name}_long-sort.bam original_ROIs.gtf > ${base_name}-long_ori_counts.txt
htseq-count -f bam -s no -a 0 -t Loop -i ID ${base_name}_long-sort.bam longRange_dup_ROIs.gtf > ${base_name}-long_dup_counts.txt
htseq-count -f bam -s no -a 0 -t Loop -i ID ${base_name}_prox-sort.bam original2_ROIs.gtf > ${base_name}-prox_ori_counts.txt
htseq-count -f bam -s no -a 0 -t Loop -i ID ${base_name}_prox-sort.bam proximal_dup_ROIs.gtf > ${base_name}-prox_dup_counts.txt
