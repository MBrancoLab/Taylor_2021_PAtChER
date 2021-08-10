#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 4
#$ -l h_rt=12:0:0
#$ -l h_vmem=4G


module load bowtie2/2.4.1
module load samtools
module load R
module load use.own
module load hicup

hicup --config hicup.conf
