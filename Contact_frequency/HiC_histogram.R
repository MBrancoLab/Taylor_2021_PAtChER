##Plot contact frequency, as assessed from HiCUP mapping of Hi-C data (hicup_align.sh script)
##Histogram values were produced in Seqmonk


library('tidyverse')


h = read_tsv('HiC_histogram.txt', col_names=c('lower','upper','count'), skip=1)

mid = rowMeans(h[,1:2])/1000
nf = h$Count/h$Count[1]

quartz(w=3,h=3)
h %>% mutate(mid=(lower+upper)/2000, nf=count/h$count[1]) %>%
ggplot(aes(mid, nf)) +
	theme_classic() +
	xlab('Distance (kb)') + ylab('Normalised frequency') +
	geom_point() +
	geom_line()  #Figure 1B
	