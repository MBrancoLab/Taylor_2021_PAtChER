#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G

module load python/3.6.3
module load use.own
module load localtools

te=MERVL-int_4kb

##bigwig files not in repository due to size
computeMatrix scale-regions -S AML12-K9me3_norm-all.bw AML12-K9me3_norm-uni.bw AML12-K9me3_norm-rand.bw \
	-R $te.bed \
	-m 5000 \
	-a 2000 \
	-b 2000 \
	-bs 100 \
	-o $te.mat

##Figure 5E
plotHeatmap --matrixFile $te.mat \
	--outFileName $te.svg \
	--outFileSortedRegions ${te}_order.txt \
	--sortUsing sum \
	--sortUsingSamples 1 \
	--heatmapHeight 12 \
	--heatmapWidth 3 \
	--whatToShow "heatmap and colorbar"
