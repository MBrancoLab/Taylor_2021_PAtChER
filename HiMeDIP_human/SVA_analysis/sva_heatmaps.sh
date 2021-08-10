#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G

module load python/3.6.3
module load use.own
module load localtools

te=SVA_1kb

##bigwig files not in repository due to size
computeMatrix scale-regions -S Pool_4_MCF7_ratio-all.bw Pool4_2102EP_ratio-all.bw Pool_4_MCF7_ratio-uni.bw Pool4_2102EP_ratio-uni.bw \
	-R $te.bed \
	-m 1500 \
	-a 1000 \
	-b 1000 \
	-bs 100 \
	-o $te.mat

##Figure 5B
plotHeatmap --matrixFile $te.mat \
	--outFileName $te.svg \
	--outFileSortedRegions ${te}_order.txt \
	--sortUsing sum \
	--sortUsingSamples 1 \
	--heatmapHeight 12 \
	--heatmapWidth 3 \
	--whatToShow "heatmap and colorbar"
