#!/bin/bash
# FILE: bramcoLoopsMCF7.sh
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q gpu.q
#$ -l gpu=1
source /usr/local/gpuallocation.sh

java -Xms512m -Xmx2048m -jar /home/og16379/juicer_dir/juicer/scripts/common/juicer_tools.jar hiccups \
-r 5000,10000,25000  -k KR -f .1,.1,.1 -p 4,2,1  -i 7,5,3 -t 0.02,1.5,1.75,2  -d 20000,20000,50000   \
/storage/projects/ZabetLab/livAnalysis/juicerFiles/brancoMCF7_HiC.hic \
/storage/projects/ZabetLab/livAnalysis/juicerFiles/branco_hiccups_loops_10Kb
