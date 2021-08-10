#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G

module load python/3.6.3
module load use.own
module load localtools

te=RLTR10-int_1kb

##bigwig files not in repository due to size
computeMatrix scale-regions -S E14_meDIP-all_norm.bw E14_meDIP-uni_norm.bw \
	-R $te.bed \
	-m 2000 \
	-a 1000 \
	-b 1000 \
	-bs 100 \
	-o $te.mat

##Supp Fig 5D
plotHeatmap --matrixFile $te.mat \
	--outFileName $te.svg \
	--outFileSortedRegions ${te}_order.txt \
	--sortUsing sum \
	--sortUsingSamples 1 \
	--heatmapHeight 12 \
	--heatmapWidth 3 \
	--whatToShow "heatmap and colorbar"
