#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=24:0:0
#$ -l h_vmem=24G
#$ -t 1-160


module load python/3.6.3

R1=$(sed -n "${SGE_TASK_ID}p" fastq_list.txt)
R2=$(echo $R1 | sed s/_R1.part/_R2.part/)
base=$(echo $R1 | sed s/R1.part//)

##align
. ~/PAtChER/patcherenv/bin/activate
patcher -g ~/genomes/mm10/mm10.fa -r1 $R1 -r2 $R2 -d 10000 -o $base.bam -b

##unpair reads
unpair -i $base.bam -o ${base}-unpaired.bam

