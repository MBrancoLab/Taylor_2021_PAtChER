##Process read counts at original and duplicated positions (minigenome_align.sh script)
##Counts are separated into unique, R and P. The latter two are PAtChER-aligned read types:
##R - Region-based alignment; one hit within the region defined by Ds
##P - Proximity-based alignment; more than one hit within Ds, nearest is used.

##Count data used here to assess mapping accuracy and read recovery at original location


##read count data

dataset = 'Bonev'

if (dataset=='Bonev') files = list.files('./counts', pattern='^ESC', full.names=T)
if (dataset=='Mumbach') files = list.files('./counts', pattern='^mES_HiChIP_H3K27ac', full.names=T)

sp.names = strsplit(files, split='-')
group = factor(unlist(lapply(sp.names, function(x) x[2])))

counts = data.frame(roi = read.delim(files[1], header=F)$V1) #get ROI names
for (g in levels(group)) {
	g.files = files[group==g]
	g.count = numeric(nrow(counts))
	for (f in g.files) g.count = g.count + read.delim(f, header=F)$V2
	counts = cbind(counts, g.count)
}
colnames(counts)[-1] = gsub('\\.txt','',levels(group))
counts = counts[!grepl('^__',counts$roi),] #remove extra lines from htseq-count


##get distances between original and duplicated ROIs

gtf = read.delim('original_ROIs.gtf', header=F, as.is=T)
attrib = strsplit(gtf$V9, split='[; ]')
gtf.id = unlist(lapply(attrib, function(x) x[5]))
gtf.dist = unlist(lapply(attrib, function(x) x[2]))

dist = as.numeric(gtf.dist[match(counts$roi, gtf.id)])/1000


##normalised counts per distance

norm.counts = function(prefix, include=c('unique','R','P')) {
	ori = dup = numeric(nrow(counts))
	if ('unique' %in% include) {
		ori = ori+counts[,colnames(counts)=='unique_oriCounts']
		dup = dup+counts[,colnames(counts)=='unique_dupCounts']
	}
	if ('R' %in% include) {
		ori = ori+counts[,colnames(counts)==paste(prefix,'R_oriCounts',sep='')]
		dup = dup+counts[,colnames(counts)==paste(prefix,'R_dupCounts',sep='')]
	}
	if ('P' %in% include) {
		ori = ori+counts[,colnames(counts)==paste(prefix,'P_oriCounts',sep='')]
		dup = dup+counts[,colnames(counts)==paste(prefix,'P_dupCounts',sep='')]
	}
	ori = ori/counts$nondup_counts*100
	dup = dup/counts$nondup_counts*100
	out = data.frame(roi=counts$roi, ori, dup, dist)
	return(out)
}


##tidy data to plot 

library(tidyverse)

sel.dist = c(5, 10, 20, 30, 50)

pdata = rbind(norm.counts('10kb', include='unique') %>% add_column(Ds='0kb'),
			norm.counts('10kb') %>% add_column(Ds='10kb'),
			norm.counts('20kb') %>% add_column(Ds='20kb'),
			norm.counts('40kb') %>% add_column(Ds='40kb')) %>%
			filter(dist %in% sel.dist) %>%
			mutate(acc = ori/(dup+ori)*100)

##summarise

sdata = pdata %>% group_by(dist,Ds) %>%
	summarise(mrec=mean(ori), sdrec=sd(ori), macc=mean(acc, na.rm=T), sdacc=sd(acc, na.rm=T))


##plot recovery (20kb only)

quartz(w=4,h=3)
sdata %>% filter(Ds=='0kb' | Ds=='20kb') %>%
ggplot(aes(factor(dist), mrec, fill=Ds)) +
	theme_classic() +
	geom_bar(stat='identity', position=position_dodge(), colour='black', width=0.8) +
	geom_errorbar(aes(ymin=mrec, ymax=mrec+sdrec), width=0.2, position=position_dodge(0.8)) +
	scale_fill_manual(values=c('grey','orange')) +
	xlab('Inter-repeat distance (kb)') +
	ylab('% Recovered reads') +
	labs(fill='')  #Figure 2F


##plot accuracy (20kb only)

quartz(w=4,h=3)
sdata %>% filter(Ds=='0kb' | Ds=='20kb') %>%
ggplot(aes(factor(dist), macc, fill=Ds)) +
	theme_classic() +
	geom_bar(stat='identity', position=position_dodge(), colour='black', width=0.8) +
	geom_errorbar(aes(ymin=macc, ymax=macc+sdacc), width=0.2, position=position_dodge(0.8)) +
	scale_fill_manual(values=c('grey','orange')) +
	xlab('Inter-repeat distance (kb)') +
	ylab('Accuracy (%)') +
	labs(fill='')  #Figure 2F
	

##plot recovery (all)

quartz(w=5,h=3)
sdata %>%
ggplot(aes(factor(dist), mrec, fill=Ds)) +
	theme_classic() +
	geom_bar(stat='identity', position=position_dodge(), colour='black', width=0.8) +
	geom_errorbar(aes(ymin=mrec, ymax=mrec+sdrec), width=0.2, position=position_dodge(0.8)) +
	scale_fill_manual(values=c('grey','yellow2','orange','red')) +
	xlab('Inter-repeat distance (kb)') +
	ylab('% Recovered reads') +
	labs(fill='')  #Supp Fig 2E


##plot accuracy (all)

quartz(w=5,h=3)
sdata %>%
ggplot(aes(factor(dist), macc, fill=Ds)) +
	theme_classic() +
	geom_bar(stat='identity', position=position_dodge(), colour='black', width=0.8) +
	geom_errorbar(aes(ymin=macc, ymax=macc+sdacc), width=0.2, position=position_dodge(0.8)) +
	scale_fill_manual(values=c('grey','yellow2','orange','red')) +
	xlab('Inter-repeat distance (kb)') +
	ylab('Accuracy (%)') +
	labs(fill='')  #Supp Fig 2E
	
