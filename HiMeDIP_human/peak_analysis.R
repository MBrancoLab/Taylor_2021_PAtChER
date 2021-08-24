##Analyse HiMeDIP peaks intersected with either L1 gtf files or RepeatMasker (peak_detection.sh script)


library(tidyverse)


##read intersections with L1 regions

data = tibble(chr=character(), source=character(), feature=character(),
	start=integer(), end=integer(), met=numeric(), strand=character(),
	empty=character(), roi=character(), count=integer(),
	region=character(), cells=character(), type=character())

for (region in c('int','up')) {
	for (cells in c('Pool_4_MCF7','Pool4_2102EP')) {
		for (type in c('all','uni')) {
			sub = read_tsv(paste('peak_intersections/',cells,'_ratio-',type,'_intersect_',region,'.txt', sep=''),
				col_names=colnames(data)[1:10])
			sub = add_column(sub, region=region, cells=cells, type=type)
			data = rbind(data, sub)
		}
	}
}


##summarise L1 data

sum.data = data %>% mutate(met_group=cut(met, c(0,25,75,100), include.lowest=TRUE)) %>%
	group_by(cells, type, region, met_group) %>%
	summarise(total=sum(count>0), perc=sum(count>0)/length(count)*100)


##plot L1 data

plot.peaks = function(cell, reads, y.lim=c(0,70)) {
	sum.data %>% filter(cells==cell, type==reads) %>%
	ggplot(aes(met_group, perc, fill=factor(region,levels=c('up','int')))) +
		geom_bar(stat='identity', position=position_dodge(), colour='black') +
		theme_classic() +
		ylim(y.lim[1], y.lim[2]) +
		xlab('% methylation') +
		ylab('% with peaks') +
		labs(fill='') +
		scale_fill_discrete(labels=c("5' end",'Internal'))
}

quartz(w=3.5,h=4)
plot.peaks(cell='Pool_4_MCF7', reads='uni')  #Figure 4D

quartz(w=3.5,h=4)
plot.peaks(cell='Pool_4_MCF7', reads='all')  #Figure 4D


##read intersections with Repeatmasker (no simple repeats)
cnames = c('chr.peak', 'start.peak', 'end.peak',
	'chr.rep', 'start.rep', 'end.rep', 'repname', 'score', 'strand', 'overlap')
all = read_tsv('peak_intersections/Pool_4_MCF7_ratio-all_intersect_rm.txt.gz', col_names=cnames) %>%
	add_column(type='all')
uni = read_tsv('peak_intersections/Pool_4_MCF7_ratio-uni_intersect_rm.txt.gz', col_names=cnames) %>%
	add_column(type='uni')
data2 = rbind(all, uni) %>% mutate(is.rep=repname!='.', coord=paste(chr.peak,start.peak))


##total number of peaks
#n.peaks = data2 %>% group_by(type, is.rep) %>% summarise(n=length(unique(coord)))

#quartz(w=4,h=4)
#ggplot(n.peaks, aes(type, n, fill=is.rep)) +
#	geom_bar(stat='identity', colour='black')


##numbers of peaks per family
fam.peaks = data2 %>% filter(repname!='.') %>%
	group_by(type, repname) %>% summarise(n=length(unique(coord))) %>%
	pivot_wider(names_from=type, values_from=n, values_fill=0)

quartz(w=3.5,h=3.5)
ggplot(fam.peaks, aes(log10(uni+1), log2(all+1)-log2(uni+1))) +
	geom_point(size=0.8, colour='lightblue') +
	theme_classic() +
	ylim(-1.4,1.4) +
	geom_hline(yintercept=0,linetype='dashed') +
	xlab('log10 unique peaks') +
	ylab('log2 fold change')  #Figure 5A

