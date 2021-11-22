##Analyse HiMeDIP peaks intersected with RepeatMasker (peak_detection.sh script)


library(tidyverse)


##track-specific peaks intersected with repeatmasker (w/o simple repeats)

cnames = c('chr.peak', 'start.peak', 'end.peak',
	'chr.rep', 'start.rep', 'end.rep', 'repname', 'score', 'strand', 'overlap')
all = read_tsv('E14_meDIP-all_norm_intersect_rm.txt.gz', col_names=cnames) %>% add_column(type='all')
uni = read_tsv('E14_meDIP-uni_norm_intersect_rm.txt.gz', col_names=cnames) %>% add_column(type='uni')
data = rbind(all, uni) %>% mutate(is.rep=repname!='.', coord=paste(chr.peak,start.peak))


##number of peaks

#n.peaks = data %>% group_by(type, is.rep) %>% summarise(n=length(unique(coord)))

#quartz(w=4,h=4)
#ggplot(n.peaks, aes(type, n, fill=is.rep)) +
#	geom_bar(stat='identity', colour='black')


##numbers of peaks per family

fam.peaks = data %>%	 filter(repname!='.') %>%
	group_by(type, repname) %>% summarise(n=length(unique(coord))) %>%
	pivot_wider(names_from=type, values_from=n, values_fill=0)

quartz(w=3.5,h=3.5)
ggplot(fam.peaks, aes(log10(uni+1), log2(all+1)-log2(uni+1))) +
	geom_point(size=0.8, colour='lightblue') +
	theme_classic() +
	geom_hline(yintercept=0,linetype='dashed') +
	xlab('log10 unique peaks') +
	ylab('log2 fold change')  #Supp Fig 6C