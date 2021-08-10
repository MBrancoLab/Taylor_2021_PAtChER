#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=4:0:0
#$ -l h_vmem=8G
#$ -t 1-2

module load python/3.6.3
module load use.own
module load localtools
module load macs2
module load bedtools


ip=$(sed -n "${SGE_TASK_ID}p" ip_list.txt)
input=$(sed -n "${SGE_TASK_ID}p" input_list.txt)
norm=$(echo $ip | sed s/-IP/_norm/)

##make bigwig
bigwigCompare -b1 $ip -b2 $input --operation "ratio" -o $norm --outFileFormat "bigwig" --pseudocount 0.01 --skipZeroOverZero

##convert to bedGraph
bdg=$(echo $norm | sed s/.bw/.bdg/)
~/BrancoLab/Software/bigWigToBedGraph $norm $bdg

##call peaks
peaks=$(echo $norm | sed s/.bw/_peaks.txt/)
macs2 bdgbroadcall -i $bdg -o $peaks -C 4 -c 10 -l 300

##intersect with Repeatmasker
sed -e '1d' $peaks | cut -f 1-3 - > temp-$SGE_TASK_ID.bed
overRM=$(echo $norm | sed s/.bw/_intersect_rm.txt/)
intersectBed -wao -f 0.25 -a temp-$SGE_TASK_ID.bed -b ~/genomes/mm10/Repeatmasker_mm10_4.0.5_noSimple.bed > $overRM
rm temp-$SGE_TASK_ID.bed

