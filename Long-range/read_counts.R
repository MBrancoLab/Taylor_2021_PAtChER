##Process read counts from juicer dump
##Compares reads involved in chromatin loops with those in proximal interactions



library(tidyverse)


##read data

cl = read.delim('chromatin_loops_MCF7.bedpe', as.is=T)
pc = read.delim('proximal_counts_MCF7.txt',header=F)

#cl = read.delim('chromatin_loops_2102Ep.bedpe', as.is=T)
#pc = read.delim('proximal_counts_2102Ep.txt',header=F)


##make tibble

hic = tibble(chr = rep(cl$chr1,2),
			start = rep(cl$x1,2),
			end = rep(cl$x2,2),
			counts = c(cl$observed, pc$V3),
			type = rep(c('loops','proximal'),each=nrow(cl)))


##plot raw read counts (Supp Fig 3A)

quartz(w=4,h=3)
mutate(hic, size=(end-start)/1000) %>%
ggplot(aes(x=factor(size), y=counts, fill=type)) +
	geom_boxplot(outlier.colour=NA) +
	theme_classic() +
	ylim(0,2500) +
	xlab('Bin size (kb)') +
	ylab('Read count')


##plot ratio between counts (Supp Fig 3B)

quartz(w=4, h=3)
mutate(cl, ratio=log2(pc$V3/observed)) %>%
ggplot(aes(x=ratio)) +
	geom_histogram(fill='grey', colour='black') +
	theme_classic() +
	ylab('Frequency') +
	scale_x_continuous('Proximal/loop read count ratio',breaks=seq(2,8,2),labels=2^(seq(2,8,2)))

