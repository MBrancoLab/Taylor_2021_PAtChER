#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=8:0:0
#$ -l h_vmem=16G

module load samtools
module load bedtools
module load use.own
module load python/3.6.3
module load localtools
module load macs2

dist=20kb

##merge
samtools merge HiC_Bonev_nondup.bam ESC*nondup_unp.bam
samtools merge HiChIP_Mumbach_nondup.bam mES_HiChIP_H3K27ac*nondup_unp.bam
samtools merge HiC_Bonev_$dist.bam ESC*${dist}_unp.bam
samtools merge HiChIP_Mumbach_$dist.bam mES_HiChIP_H3K27ac*${dist}_unp.bam

##simulate non-enriched reads from duplicated region
for d in {3000..15000..1000} 20000 30000 50000
do
	##reads that would be correctly assigned to duplicated region
	grep "distSet \"${d}\"" original_ROIs.gtf > temp_ori.gtf
	bedtools intersect -a HiC_Bonev_$dist.bam -b temp_ori.gtf > temp_${d}_ori.bam
	samtools view -H temp_${d}_ori.bam > temp_${d}_correct.sam
       samtools view temp_${d}_ori.bam | awk -v d=$d -F"\t" '{OFS=FS}{$4=$4+d+3000 ; print}' >> temp_${d}_correct.sam
	rm temp_ori.gtf

	##reads that would be misassigned to original region
	grep "distSet \"${d}\"" duplicated_ROIs.gtf > temp_dup.gtf
	bedtools intersect -a HiC_Bonev_$dist.bam -b temp_dup.gtf > temp_${d}_dup.bam
	samtools view -H temp_${d}_dup.bam > temp_${d}_misassigned.sam
	samtools view temp_${d}_dup.bam | awk -v d=$d -F"\t" '{OFS=FS}{$4=$4-d-3000 ; print}' >> temp_${d}_misassigned.sam
	rm temp_dup.gtf
done

##add simulated reads to the rest
samtools merge HiC_Bonev_corrected.bam HiC_Bonev_$dist.bam temp_*_correct.sam temp_*_misassigned.sam
samtools merge HiChIP_Mumbach_corrected.bam HiChIP_Mumbach_$dist.bam temp_*_correct.sam temp_*_misassigned.sam#rm temp_*_ori.bam temp_*_dup.bam temp_*_correct.sam temp_*_misassigned.sam


##make bigwig tracks
for bam in HiC_Bonev_nondup.bam HiChIP_Mumbach_nondup.bam HiC_Bonev_corrected.bam HiChIP_Mumbach_corrected.bam
do

	##sort
	sorted=$(echo $bam | sed s/.bam/-sorted.bam/)
	samtools sort -o $sorted $bam
	samtools index $sorted

	##get unique reads
	unique=$(echo $sorted | sed s/sorted/unique/)
	samtools view -h $sorted | grep '^@\|:u\|:s' - | samtools view -b - > $unique
	samtools index $unique

	##make bigwig
	bwall=$(echo $bam | sed s/.bam/-all.bw/)
	bwuni=$(echo $unique | sed s/-unique.bam/-uni.bw/)
	bamCoverage --bam $sorted -o $bwall --binSize 100 --normalizeUsing CPM
	bamCoverage --bam $unique -o $bwuni --binSize 100 --normalizeUsing CPM
done

##normalise
bigwigCompare -b1 HiChIP_Mumbach_nondup-uni.bw -b2 HiC_Bonev_nondup-uni.bw --operation "ratio" -o H3K27ac-nondup_norm.bw --outFileFormat "bigwig" --pseudocount 0.01  --skipZeroOverZero
bigwigCompare -b1 HiChIP_Mumbach_corrected-all.bw -b2 HiC_Bonev_corrected-all.bw --operation "ratio" -o H3K27ac-all_norm.bw --outFileFormat "bigwig" --pseudocount 0.01  --skipZeroOverZero
bigwigCompare -b1 HiChIP_Mumbach_corrected-uni.bw -b2 HiC_Bonev_corrected-uni.bw --operation "ratio" -o H3K27ac-uni_norm.bw --outFileFormat "bigwig" --pseudocount 0.01  --skipZeroOverZero

##convert to bedGraph
for bw in H3K27ac-nondup_norm.bw H3K27ac-all_norm.bw H3K27ac-uni_norm.bw
do
	bdg=$(echo $bw | sed s/bw/bdg/)
	~/BrancoLab/Software/bigWigToBedGraph $bw $bdg
done

##call peaks and intersect with ROI annotations
for t in {1..10}
do
	macs2 bdgbroadcall -i H3K27ac-nondup_norm.bdg -o nondup_peaks_${t}.txt -C $t -c 20
	sed -e '1d' nondup_peaks_${t}.txt | cut -f 1-3 - > temp.bed
	bedtools intersect -c -a nonduplicated_ROIs.gtf -b temp.bed | cut -f 9,10 - > nondup_temp.txt

	macs2 bdgbroadcall -i H3K27ac-all_norm.bdg -o ori_peaks_${t}.txt -C $t -c 20
	sed -e '1d' ori_peaks_${t}.txt | cut -f 1-3 - > temp.bed
        bedtools intersect -c -a original_ROIs.gtf -b temp.bed | cut -f 10 - | paste nondup_temp.txt - > patcher_temp.txt
        bedtools intersect -c -a duplicated_ROIs.gtf -b temp.bed | cut -f 10 - | paste patcher_temp.txt - > patcher_intersection_${t}.txt

	macs2 bdgbroadcall -i H3K27ac-uni_norm.bdg -o ori_peaks_${t}.txt -C $t -c 20
        sed -e '1d' ori_peaks_${t}.txt | cut    -f 1-3 - > temp.bed
        bedtools intersect -c -a original_ROIs.gtf -b temp.bed  | cut -f 10 - | paste nondup_temp.txt - > unique_temp.txt
        bedtools intersect -c -a duplicated_ROIs.gtf -b temp.bed  | cut -f 10 - | paste unique_temp.txt - > unique_intersection_${t}.txt

	rm temp.bed nondup_temp.txt patcher_temp.txt unique_temp.txt
done
