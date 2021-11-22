#!/bin/sh
#$ -cwd
#$ -l h_vmem=2G
#$ -l h_rt=12:0:0

module load java

sed 1d chromatin_loops.bedpe | while read line
do
	chr1=$(echo $line | awk '{print $1}')
	start1=$(echo $line | awk '{print $2}')
	end1=$(echo $line | awk '{print $3}')
	res=$(($end1-$start1))

	##counts for proximal region
	java -jar juicer_tools.1.7.6_jcuda.0.8.jar dump observed NONE brancoMCF7_HiC.hic ${chr1}:${start1}:${end1} ${chr1}:${start1}:${end1} BP $res temp.txt
	grep ${start1}$'\t'${start1} temp.txt >> proximal_counts.txt
done

rm temp.txt
