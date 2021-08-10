#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G

module load python/3.6.3
module load use.own
module load localtools

computeMatrix scale-regions -S hg38_PE75_20kb_accuracy.bw -R L1HS_5kb.bed -m 6000 -a 2000 -b 2000 -bs 100 -o accuracy-l1.mat
plotProfile --matrixFile accuracy-l1.mat --outFileName accuracy-l1.png --outFileNameData accuracy-l1_trend.txt --averageType "mean"

