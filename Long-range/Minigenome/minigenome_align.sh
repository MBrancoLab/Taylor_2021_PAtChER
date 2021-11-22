#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=72:0:0
#$ -l h_vmem=8G
#$ -t 1-16

module load python/3.6.3
module load samtools
module load use.own
module load htseq

##data
R1=$(sed -n "${SGE_TASK_ID}p" fastq_list.txt)
R2=$(echo $R1 | sed s/_1.part/_2.part/)
base_name=$(basename $R1 | sed s/_1.part//)

##map to long-range minigenome
. ~/PAtChER/patcherenv/bin/activate
patcher -g longRange_minigenome.fa -r1 $R1 -r2 $R2 -d 20000 -o ${base_name}-long.bam -b

##map to proximal minigenome
patcher -g proximal_minigenome.fa -r1 $R1 -r2 $R2 -d 20000 -o ${base_name}-prox.bam -b

##unpair reads
for bam in ${base_name}-long.bam ${base_name}-prox.bam
do
	unpair -i $bam -o $(echo $bam | sed s/.bam/_unp.bam/)
done
deactivate
