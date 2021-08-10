#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=24:0:0
#$ -l h_vmem=24G
#$ -t 1-10


##generate reads
module load python/2.7.15
~/.local/bin/sim3C --seed ${SGE_TASK_ID} --dist uniform -n 10000000 -l 75 -e DpnII -m hic --simple-reads --efficiency 0.9 --anti-rate 0 --profile-name profile-${SGE_TASK_ID}.tsv hg38.fa sim-${SGE_TASK_ID}.fq

##split pairs
sed -ne '1~8{N;N;N;p}' sim-${SGE_TASK_ID}.fq > sim-${SGE_TASK_ID}-R1.fastq
sed -ne '5~8{N;N;N;p}' sim-${SGE_TASK_ID}.fq > sim-${SGE_TASK_ID}-R2.fastq

##make bed file
module unload python/2.7.15
module load python/3.6.3
python get_original_pos.py -r1 sim-${SGE_TASK_ID}-R1.fastq -r2 sim-${SGE_TASK_ID}-R2.fastq

##align
. ~/PAtChER/patcherenv/bin/activate
patcher -g hg38.fa -r1 sim-${SGE_TASK_ID}-R1.fastq -r2 sim-${SGE_TASK_ID}-R2.fastq -d 20000 -o sim-${SGE_TASK_ID}.bam -b

##unpair reads
unpair -i sim-${SGE_TASK_ID}.bam -o sim-${SGE_TASK_ID}-unpaired.bam

##split alignment type
module load samtools
samtools view -h sim-${SGE_TASK_ID}-unpaired.bam | grep '@\|:r' - > sim-${SGE_TASK_ID}-recov.sam
samtools view -h sim-${SGE_TASK_ID}-unpaired.bam | grep '@\|:p' - > sim-${SGE_TASK_ID}-prob.sam

##convert to bed
module load bedtools
bamToBed -i sim-${SGE_TASK_ID}-unpaired.bam > sim-${SGE_TASK_ID}-unpaired.bed
bamToBed -i sim-${SGE_TASK_ID}-recov.sam > sim-${SGE_TASK_ID}-recov.bed
bamToBed -i sim-${SGE_TASK_ID}-prob.sam > sim-${SGE_TASK_ID}-prob.bed

##intersect with L1HS annotation
intersectBed -u -a sim-${SGE_TASK_ID}-unpaired.bed -b L1HS.bed > sim-${SGE_TASK_ID}-L1HS.bed

##match reads
python match_reads.py -u sim-${SGE_TASK_ID}-unpaired.bed -r1 sim-${SGE_TASK_ID}-R1.bed -r2 sim-${SGE_TASK_ID}-R2.bed
python match_reads.py -u sim-${SGE_TASK_ID}-recov.bed -r1 sim-${SGE_TASK_ID}-R1.bed -r2 sim-${SGE_TASK_ID}-R2.bed
python match_reads.py -u sim-${SGE_TASK_ID}-prob.bed -r1 sim-${SGE_TASK_ID}-R1.bed -r2 sim-${SGE_TASK_ID}-R2.bed
python match_reads.py -u sim-${SGE_TASK_ID}-L1HS.bed -r1 sim-${SGE_TASK_ID}-R1.bed -r2 sim-${SGE_TASK_ID}-R2.bed
