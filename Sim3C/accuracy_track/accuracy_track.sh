#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=24:0:0
#$ -l h_vmem=24G
#$ -t 1-100


##generate reads
module load python/2.7.15
~/.local/bin/sim3C --seed ${SGE_TASK_ID} --dist uniform -n 10000000 -l 75 -e DpnII -m hic --simple-reads --efficiency 0.9 --anti-rate 0 --profile-name profile-${SGE_TASK_ID}.tsv hg38.fa sim-${SGE_TASK_ID}.fq

##split pairs
sed -ne '1~8{N;N;N;p}' sim-${SGE_TASK_ID}.fq > sim-${SGE_TASK_ID}-R1.fastq
sed -ne '5~8{N;N;N;p}' sim-${SGE_TASK_ID}.fq > sim-${SGE_TASK_ID}-R2.fastq
rm sim-${SGE_TASK_ID}.fq

##make bed file
module unload python/2.7.15
module load python/3.6.3
python get_original_pos.py -r1 sim-${SGE_TASK_ID}-R1.fastq -r2 sim-${SGE_TASK_ID}-R2.fastq

##align
. ~/PAtChER/patcherenv/bin/activate
patcher -g hg38.fa -r1 sim-${SGE_TASK_ID}-R1.fastq -r2 sim-${SGE_TASK_ID}-R2.fastq -d 20000 -o sim-${SGE_TASK_ID}.bam -b

##unpair reads
unpair -i sim-${SGE_TASK_ID}.bam -o sim-${SGE_TASK_ID}-unpaired.bam
rm sim-${SGE_TASK_ID}.bam

##convert to bed
module load bedtools
bamToBed -i sim-${SGE_TASK_ID}-unpaired.bam > sim-${SGE_TASK_ID}-unpaired.bed

##match reads
python match_reads.py -u sim-${SGE_TASK_ID}-unpaired.bed -r1 sim-${SGE_TASK_ID}-R1.bed -r2 sim-${SGE_TASK_ID}-R2.bed -wm

##get coverage
coverageBed -a probes.bed -b sim-${SGE_TASK_ID}-unpaired.bed | cut -f 1-4 - > sim-${SGE_TASK_ID}-unpaired.bedgraph
coverageBed -a probes.bed -b sim-${SGE_TASK_ID}-unpaired_matched.bed | cut -f 1-4 - > sim-${SGE_TASK_ID}-unpaired_matched.bedgraph

