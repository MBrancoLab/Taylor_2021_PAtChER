#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G

module load python/3.6.3
module load use.own
module load localtools


##human
for te in L1HS_1kb L1PA2_1kb SVA_A_1kb SVA_F_1kb
do
	multiBigwigSummary BED-file -b hg38_PE75_20kb_accuracy.bw --BED $te.bed -o $te.npz --outRawCounts $te-meanAcc.txt
done

##mouse
for te in L1Tf_1kb L1A_1kb IAPEz-int_1kb MERVL-int_1kb
do
  	multiBigwigSummary BED-file -b mm10_PE75_10kb_accuracy.bw --BED $te.bed -o $te.npz --outRawCounts $te-meanAcc.txt
done
