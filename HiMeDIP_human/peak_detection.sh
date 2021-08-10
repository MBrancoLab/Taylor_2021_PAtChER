#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=4:0:0
#$ -l h_vmem=8G
#$ -t 1-4

module load python/3.6.3
module load use.own
module load localtools
module load macs2
module load bedtools


ip=$(sed -n "${SGE_TASK_ID}p" ip_list.txt)
input=$(sed -n "${SGE_TASK_ID}p" input_list.txt)
norm=$(echo $ip | sed s/_IP/_ratio/)

##make bigwig
bigwigCompare -b1 $ip -b2 $input --operation "ratio" -o $norm --outFileFormat "bigwig" --pseudocount 0.01 --skipZeroOverZero

##convert to bedGraph
bdg=$(echo $norm | sed s/.bw/.bdg/)
~/BrancoLab/Software/bigWigToBedGraph $norm $bdg

##call peaks
peaks=$(echo $norm | sed s/.bw/_peaks.txt/)
macs2 bdgbroadcall -i $bdg -o $peaks -C 2 -c 10


##define gtf files
if echo $ip | grep -q MCF7; then
	up=MCF7_ATLAS_5prime.gtf
	int=MCF7_ATLAS_internal.gtf
else
	up=2102Ep_ATLAS_5prime.gtf
        int=2102Ep_ATLAS_internal.gtf
fi

##intersect peaks with gtf files
overUp=$(echo $norm | sed s/.bw/_intersect_up.txt/)
overInt=$(echo $norm | sed s/.bw/_intersect_int.txt/)
sed -e '1d' $peaks | cut -f 1-3 > temp-$SGE_TASK_ID.bed
intersectBed -c -a $up -b temp-$SGE_TASK_ID.bed > $overUp
intersectBed -c -a $int -b temp-$SGE_TASK_ID.bed > $overInt


##intersect with Repeatmasker
overRM=$(echo $norm | sed s/.bw/_intersect_rm.txt/)
intersectBed -wao -f 0.25 -a temp-$SGE_TASK_ID.bed -b ~/genomes/GRCh38/Repeatmasker_hg38_4.0.5_noSimple.bed > $overRM
rm temp-$SGE_TASK_ID.bed

