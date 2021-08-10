#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=24:0:0
#$ -l h_vmem=24G
#$ -t 1-4

module load samtools
module load python/3.6.3
module load use.own
module load localtools

prefix=$(sed -n "${SGE_TASK_ID}p" prefix_list.txt)

##merge bam files
samtools merge $prefix.bam ${prefix}*unpaired.bam

##sort
samtools sort -o $prefix-sorted.bam $prefix.bam
samtools index $prefix-sorted.bam

##get unique reads
samtools view -h $prefix-sorted.bam | grep '^@\|:u\|:s' - | samtools view -b - > $prefix-unique.bam
samtools index $prefix-unique.bam

##make bigwig
bamCoverage --bam $prefix-sorted.bam -o $prefix-all.bw --binSize 100 --normalizeUsing CPM
bamCoverage --bam $prefix-unique.bam -o $prefix-uni.bw --binSize 100 --normalizeUsing CPM

##count reads within L1 regions of interest
module unload python/3.6.3
module load python/3.8.5
module load htseq
if echo $prefix | grep -q 2102EP; then
	l1gtf=2102Ep_ATLAS_internal.gtf
	upgtf=2102Ep_ATLAS_5prime.gtf
else
	l1gtf=MCF7_ATLAS_internal.gtf
	upgtf=MCF7_ATLAS_5prime.gtf
fi
for bam in $prefix-sorted.bam $prefix-unique.bam
do
	htseq-count -f bam -s no -a 0 -t L1_ROI -i ID $bam $l1gtf > $(echo $bam | sed s/.bam/_L1counts.txt/)
	htseq-count -f bam -s no -a 0 -t L1_ROI -i ID $bam $upgtf > $(echo $bam | sed s/.bam/_UPcounts.txt/)
done
