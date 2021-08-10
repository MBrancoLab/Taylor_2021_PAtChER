# Taylor_2021_PAtChER

Scripts and files from Taylor et al 2021 bioRxiv preprint


## Contact frequency

Hi-C data were mapped with HiCUP (hicup_align.sh). The histogram of contact frequencies in Figure 1B was generated by HiC_histogram.R.


## Coverage

*Primary processing*

Hi-C data were aligned with PAtChER (patcher_align.sh). The get_coverage.sh script then: a) subsets the mapped data to test the effect of read depth (Supp. Fig. 1D), b) extracts uniquely mapped reads to be analysed in parallel, c) extracts the coverage for all RepeatMasker entries, and d) produces bigwig files that were uploaded to the WashU Epigenome Browser to generate Figure 1F.

*Polymorphic TEs*

Data from Nellaker et al (PMID: 22703977) was used to remove reference B6 insertions that are absent in all sequenced 129 genomes. The poly_TE.R script generates a list of LINE and IAP elements that are absent from the 129 strain. After liftover to mm10, this file (refTE_not129_mm10.bed) was used to remove these insertions from the RepeatMasker annotation (with bedtools intersect).

*Coverage distribution*

The coverage.R script takes the RepeatMasker coverage information and plots the coverage distribution for different groups of repeats depending on their age (Figure 1D) or family (Supplementary Figure 1).

*Coverage across TEs*

The te_trends.sh script uses the bigwig files to generate mean read depth profiles across selected TE families. Figure 1E was plotted by te_trends.R.


## Sim3C

*Primary processing*

A subset of real Hi-C and Sim3C data were compared in terms of interaction distances in the sim3C_distribution.R script, which yields Supplementary Figure 2A. To assess mapping accuracy, the sim3c.sh script: a) generates Sim3C data, b) records the original position of the reads in a bed file using get_original_pos.py, c) maps the data using PAtChER, d) separates alignment types, e) matches aligned reads to simulated reads using match_reads.py (which can also be applied to a subset of reads intersecting TE annotations).

*Overall accuracy*

The accuracy_analysis.R script uses the outputs from match_reads.py (compiled in accuracy_mm10.txt and accuracy_hg38.txt) to generate genome-wide accuracy plot (Figure 2A), as well as TE family-specific plots (Figure 2B).

*Element-specific accuracy*

Genome-wide accuracy tracks were generated using a more extensive Sim3C dataset (accuracy_track.sh). These were used for Figure 2C and Supplementary Figure 2C. The accuracy tracks were also used by mean_accuracy.sh to extract a mean accuracy value per element of selected TE families. The accuracy distributions in Supplementary Figure 2B were then generated by accuracy_analysis.R. The average accuracy profiles in Figure 2D and Supplementary Figure 2D were also extracted from the accuracy tracks using te_trends.sh, and plotted using te_trends.R.


## Minigenome

*Generating the minigenomes*

Based on selected regions surrounding gene promoters (gene_promoters.bed), make_minigenome.R extracts the respective sequences, generates the duplicate regions are variable distances, and merges all sequences. Both the mm10 and modified versions of the minigenomes are written, as well as GTF annotations of the duplicated regions.

*Mapping accuracy*

Hi-C or HiChIP data were aligned using minigenome.sh, which then counts the reads at the regions of interest defined above, split by read type. In accuracy_analysis.R, these counts are processed to generate the accuracy and read recovery plots in Figure 2F and Supplementary Figure 2E.

*H3K27ac HiChIP analysis*

For enrichment and peak analyses, read counts that would have come from the artificially duplicate region are estimated, assuming no enrichment (see paper methods for more details). In the enrichment analysis this is done within enrichment_analysis.sh, which then generates the plots in Figures 3C,D. In the peak analysis this is done in peak_detection.sh, which then generates normalised bigWig tracks (used in Figure 3B) and performs peak detection using MACS2. Peaks overlapping the regions of interest within the minigenome are analysed in peak_analysis.R to generate the ROC curves in Supplementary Figure 3, and the plot in Figure 3E.


## HiMeDIP human

*Primary processing*

LINE-1 DNA methylation data from targeted bisulphite sequencing were parsed into GTF files (for different regions of interest) using make_gtf.R. After alignment of HiMeDIP data using patcher_align.sh, reads were counted at the bisulphite-covered regions with count_reads.sh.

*Enrichment analysis*

Read counts from above were further processed with enrichment_analysis.R, wherein counts are normalised and coupled to the bisulphite data. The same script then generates plots correlating enrichment with bisulphite data (Figure 4B and Supplementary Figure 4A), as well as plots correlating the difference in DNA methylation measured by the two techniques (Figure 4C and Supplementary Figure 4B).

*Peak analysis*

Normalised bigWig files were generated in peak_detection.sh (used in Figures 4E and 5C, and Supplementary Figure 5B), which then performs peak detection using MACS2, and intersects them with the LINE-1 regions, as well as with RepeatMasker. These peak intersections were then analysed in peak_analysis.R to generate the plots in Figures 4D and 5A.

*SVA analysis*

The SVA heatmaps in Figure 5B were generated with sva_heatmaps.sh. To correlate SVA methylation with gene expression, gene TSSs were extracted with get_TSS.R, and the closest SVA element determined with bedtools (generating TSS_closest_SVA.bed). HiMeDIP read counts at SVAs, and RNA-seq reads counts at genes were generated with htseq-count. SVA_analysis.R then compiles these data to generate the scatter plot in Supplementary Figure 5A.


## HiChIP mouse

*Primary processing*

All data were aligned with PAtChER using patcher_align.sh and furher processed with make_bigwigs.sh to extract unique reads and generate bigWig files.

*ESC meDIP*

Normalised bigWig files were generated in peak_detection.sh, which then performs peak detection using MACS2, and intersects them with RepeatMasker. These peak intersections were then analysed in peak_analysis.R to generate the plot in Supplementary Figure 5C. The RLTR10-int heatmaps in Supplementary Figure 5D were generated with rltr10_heatmaps.sh.

*AML12 H3K9me3*

Normalised bigWig files were generated in peak_detection.sh (used in Figure 5F), which then performs peak detection using MACS2, and intersects them with RepeatMasker. These peak intersections were then analysed in peak_analysis.R to generate the plot in Figure 5D. The MERVL-int heatmaps in Figure 5E were generated with mervl_heatmaps.sh.

*ESC H3K27ac*

The MERVL-int heatmaps in Figure 5G were generated with mervl_heatmaps.sh, which also performs k-means clustering of the elements. The nearest genes to elements in these clusters were extracted using bedtools, and RNA-seq reads within genes counted using htseq-counts. These data were processed in MERVL_context.R to generate the plots in Figures 5H and 5I. The same script takes an annotation of spatial A and B compartments to generate the plot in Figure 5J.
