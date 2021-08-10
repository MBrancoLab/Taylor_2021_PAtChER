#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=24:0:0
#$ -l h_vmem=24G
#$ -t 1-186

module load trimgalore/0.6.5
module load python/3.6.3


##define files
R1=$(sed -n "${SGE_TASK_ID}p" fastq_list.txt)
R2=$(echo $R1 | sed s/L1_1/L1_2/)
base=$(echo $R1 | sed s/_L1_1.part//)


##hard-trim reads down to 75bp
trim_galore --hardtrim5 75 $R1
trim_galore --hardtrim5 75 $R2

##map with PAtChER
. ~/PAtChER/patcherenv/bin/activate
patcher -g ~/genomes/GRCh38/hg38.fa -r1 $R1.75bp_5prime.fq -r2 $R2.75bp_5prime.fq -d 20000 -o $base.bam -b

##unpair reads
unpair -i $base.bam -o $base-unpaired.bam
deactivate
