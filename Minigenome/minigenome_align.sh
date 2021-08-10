#!/bin/sh
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=72:0:0
#$ -l h_vmem=2G
#$ -t 1-8

module load python/3.6.3
module load samtools
module load use.own
module load htseq

##data
R1=$(sed -n "${SGE_TASK_ID}p" fastq_list.txt)
R2=$(echo $R1 | sed s/_R1.fastq.gz/_R2.fastq.gz/)
base_name=$(basename $R1 | sed s/_R1.fastq.gz//)

##map to mm10-based minigenome
. ~/PAtChER/patcherenv/bin/activate
patcher -g promoter_minigenome_ori.fa -r1 $R1 -r2 $R2 -d 10000 -o ${base_name}-nondup.bam -b &

##map to modified minigenome with duplications
patcher -g promoter_minigenome.fa -r1 $R1 -r2 $R2 -d 10000 -o ${base_name}-10kb.bam -b &
patcher -g promoter_minigenome.fa -r1 $R1 -r2 $R2 -d 20000 -o ${base_name}-20kb.bam -b &
patcher -g promoter_minigenome.fa -r1 $R1 -r2 $R2 -d 40000 -o ${base_name}-40kb.bam -b &
wait

##unpair reads
for bam in ${base_name}-nondup.bam ${base_name}-10kb.bam ${base_name}-20kb.bam ${base_name}-40kb.bam
do
	unpair -i $bam -o $(echo $bam | sed s/.bam/_unp.bam/)
done
deactivate

##split alignment type
samtools view -h ${base_name}-nondup_unp.bam | grep '^@\|:u\|:s' - > ${base_name}-nondup.sam
samtools view -h ${base_name}-10kb_unp.bam | grep '^@\|:u\|:s' - > ${base_name}-unique.sam
for bam in ${base_name}-10kb_unp.bam ${base_name}-20kb_unp.bam ${base_name}-40kb_unp.bam
do
	samtools view -h $bam | grep '^@\|:r' - > $(echo $bam | sed s/_unp.bam/R.sam/)
	samtools view -h $bam | grep '^@\|:p' - > $(echo $bam | sed s/_unp.bam/P.sam/)
done

##count reads within regions of interest
module unload python/3.6.3
module load python/3.8.5
htseq-count -f sam -s no -a 0 -t Promoter -i ID ${base_name}-nondup.sam nonduplicated_ROIs.gtf > ${base_name}-nondup_counts.txt
for sam in ${base_name}-unique.sam ${base_name}-10kbR.sam ${base_name}-20kbR.sam ${base_name}-40kbR.sam ${base_name}-10kbP.sam ${base_name}-20kbP.sam ${base_name}-40kbP.sam
do
	htseq-count -f sam -s no -a 0 -t Promoter -i ID $sam original_ROIs.gtf > $(echo $sam | sed s/.sam/_oriCounts.txt/)
	htseq-count -f sam -s no -a 0 -t Promoter -i ID $sam duplicated_ROIs.gtf > $(echo $sam | sed s/.sam/_dupCounts.txt/)
done
