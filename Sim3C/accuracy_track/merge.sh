#!/bin/sh
#$ -cwd
#$ -l h_vmem=4G
#$ -pe smp 1
#$ -l h_rt=1:0:0

module load python/3.6.3
module load use.own
module load localtools

#cp sim-1-unpaired.bedgraph ./sim3c-all.bg
#cp sim-1-unpaired_matched.bedgraph ./sim3c-correct.bg
#for i in {2..100}
#do
#	paste sim3c-all.bg sim-$i-unpaired.bedgraph | awk -v OFS='\t' '{print $1, $2, $3, $4+$8}' > temp.bg
#	mv temp.bg sim3c-all.bg
#	paste sim3c-correct.bg sim-$i-unpaired_matched.bedgraph | awk -v OFS='\t' '{print $1, $2, $3, $4+$8}' > temp.bg
#       mv temp.bg sim3c-correct.bg
#done

#sort -k1,1 -k2,2n sim3c-all.bg > sim3c-all-sorted.bg
#~/BrancoLab/Software/bedGraphToBigWig sim3c-all-sorted.bg ~/genomes/GRCh38/hg38.genome sim3c-all.bw

#sort -k1,1 -k2,2n sim3c-correct.bg > sim3c-corr-sorted.bg 
#~/BrancoLab/Software/bedGraphToBigWig sim3c-corr-sorted.bg ~/genomes/GRCh38/hg38.genome sim3c-correct.bw

bigwigCompare -b1 sim3c-correct.bw -b2 sim3c-all.bw --operation "ratio" -o sim3c-zero.bw --outFileFormat "bigwig" --pseudocount 0
