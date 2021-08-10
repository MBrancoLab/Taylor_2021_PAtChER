#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=4:0:0
#$ -l h_vmem=2G

module load python/3.6.3
module load use.own
module load localtools

computeMatrix scale-regions -S P2102EP_INPUT-all.bw P2102EP_INPUT-uni.bw -R L1HS_5kb.bed -m 6000 -a 2000 -b 2000 -bs 100 -o P2102Ep-l1.mat
plotProfile --matrixFile P2102Ep-l1.mat --outFileName P2102Ep-l1.png --outFileNameData P2102Ep-l1_trend.txt --averageType "mean"

computeMatrix scale-regions -S ESC-all.bw ESC-uni.bw -R MERVL-int_4kb.bed -m 6000 -a 2000 -b 2000 -bs 100 -o ESC-mervl.mat
plotProfile --matrixFile ESC-mervl.mat --outFileName ESC-mervl.png --outFileNameData ESC-mervl_trend.txt --averageType "mean"
