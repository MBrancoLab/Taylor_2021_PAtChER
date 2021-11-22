################################################################
# Branco lab
# Calculating chromatin loops
# Olivia Grant
###############################################################

################################################################################
# using juicer https://github.com/aidenlab/juicer
# this code depends on a particular file structure, please see above github for details
################################################################################

# first, create restriction enzyme proposed cutting sites
cd ~/genome/
bwa index hg38.fa
python /usr/local/juicer/misc/generate_site_positions.py DpnII hg38 hg38.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}'  hg38_DpnII.txt > hg38.chrom.sizes
awk 'BEGIN{OFS="\t"}{print $1, $NF}'  hg38_DpnII.txt
gawk 'BEGIN{OFS="\t"}{print $1, $NF}'  hg38_DpnII.txt
less hg38.chrom.sizes


################################################################################
# first working with MCF7 cell line
################################################################################

################################################################################
# make directory for back up files
# make sure the files are named as such "filename_R1.fq"
cd /home/og16379/juicer_dir/juicer/Miguel_hg38/
gunzip MCF7_EKDL190133930-2a-AK4825_HVT3NDSXX_R1.fq.gz
gunzip MCF7_EKDL190133930-2a-AK4825_HVT3NDSXX_R2.fq.gq

#call preprocessing pipeline
/home/og16379/juicer_dir/juicer/CPU/juicer.sh -t 30 -z ~/juicer_dir/juicer/genome/hg38.fa  -p ~/juicer_dir/juicer/genome/hg38.chrom.sizes -y ~/juicer_dir/juicer/genome/hg38_DpnII.txt -d /storage/projects/ZabetLab/livAnalysis/juicerFiles

#generate the hic files
java -Xmx10g -jar ~/juicer_dir/juicer/scripts/common/juicer_tools.jar pre -q 30 -r 10000/storage/projects/ZabetLab/livAnalysis/juicerFiles/aligned/merged_nodups.txt /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic ~/juicer_dir/juicer/genome/hg38.chrom.sizes

# calculate significant interactions
# make compartments directory
# you can put this in a shell script, but to follow along the code easily
# i share the code here

mkdir compartments
#!/bin/bash
java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p  KR  /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr1 chr1 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr1.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr2 chr2 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr2.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr3 chr3 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr3.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr4 chr4 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr4.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr5 chr5 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr5.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr6 chr6 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr6.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr7 chr7 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr7.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr8 chr8 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr8.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr9 chr9 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr9.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr10 chr10 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr10.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr11 chr11 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr11.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr12 chr12 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr12.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr13 chr13 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr13.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr14 chr14 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr14.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr15 chr15 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr15.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr16 chr16 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr16.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr17 chr17 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr17.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr18 chr18 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr18.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr19 chr19 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr19.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr20 chr20 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr20.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr21 chr21 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr21.txt


java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chr22 chr22 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chr22.txt


java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chrX chrX BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chrX.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
chrY chrY BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel_chrY.txt

#gunzip all the files
cd compartments
gzip *.txt

################################################################################
# then, generate loops
# here, we do use a shell script. see github for details of this script
################################################################################

qsub ./storage/projects/ZabetLab/livAnalysis/scripts/brancoLoopsMCF7.sh

java -Xmx10g -jar ~/juicer_dir/juicer/scripts/common/juicer_tools.jar apa -r 5000,2000,1000,500 -u \
/storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic  \
/storage/projects/ZabetLab/livAnalysis/juicerFiles/branco_hiccups_loops_10Kb/merged_loops.bedpe \
/storage/projects/ZabetLab/livAnalysis/juicerFiles/1_10000_hiccups_loops_10Kb_Branco

################################################################################
################################################################################
# second, working with P2102EP cell line
################################################################################
################################################################################


#  i moved everything out of the dir, to ensure juicer worked nicely, cleaning them out completely
# next, add the next cell line in relevant directories and gunzip as we did above
gunzip P2102EP_EKDL190133930-2a-AK5227_HVT3NDSXX_L1_R1.fq.gz
gunzip P2102EP_EKDL190133930-2a-AK5227_HVT3NDSXX_L1_R2.fq.gz

#call preprocessing pipeline
/home/og16379/juicer_dir/juicer/CPU/juicer.sh -t 30 -z ~/juicer_dir/juicer/genome/hg38.fa  -p ~/juicer_dir/juicer/genome/hg38.chrom.sizes -y ~/juicer_dir/juicer/genome/hg38_DpnII.txt -d /storage/projects/ZabetLab/livAnalysis/juicerFiles

#generate the hic files
java -Xmx10g -jar ~/juicer_dir/juicer/scripts/common/juicer_tools.jar pre -q 30 -r 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/aligned/merged_nodups.txt /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic ~/juicer_dir/juicer/genome/hg38.chrom.sizes

################################################################################
# then generate significant interactions for each chr
# we are doing this for the cell line P2102EP
################################################################################
# make a new directory to store the files
mkdir compartments

#!/bin/bash
java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p  KR  /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr1 chr1 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr1.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr2 chr2 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr2.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr3 chr3 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr3.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr4 chr4 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr4.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr5 chr5 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr5.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr6 chr6 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr6.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr7 chr7 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr7.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr8 chr8 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr8.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr9 chr9 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr9.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr10 chr10 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr10.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr11 chr11 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr11.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr12 chr12 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr12.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr13 chr13 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr13.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr14 chr14 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr14.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr15 chr15 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr15.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr16 chr16 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr16.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr17 chr17 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr17.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr18 chr18 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr18.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr19 chr19 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr19.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr20 chr20 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr20.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr21 chr21 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr21.txt


java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chr22 chr22 BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chr22.txt


java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chrX chrX BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chrX.txt

java -Xmx10g -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar dump oe \
-p KR /storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
chrY chrY BP 10000 /storage/projects/ZabetLab/livAnalysis/juicerFiles/compartments/miguel2_chrY.txt

#gunzip all the files
cd compartments
gzip *.txt

################################################################################
# then, generate loops for P2102EP cell line
################################################################################

# run shell script to calculate the loops
qsub ./storage/projects/ZabetLab/livAnalysis/scripts/brancoLoopsP2102EP.sh

java -Xmx10g -jar ~/juicer_dir/juicer/scripts/common/juicer_tools.jar apa -r 5000,2000,1000,500 -u \
/storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoP2102EP_HiC.hic  \
/storage/projects/ZabetLab/livAnalysis/juicerFiles/branco_hiccups_loops_10Kb_P2102EP/merged_loops.bedpe \
/storage/projects/ZabetLab/livAnalysis/juicerFiles/2_10000_hiccups_loops_10Kb_Branco


################################################################################
# complete
################################################################################
