##Process read counts at original and duplicated positions (minigenome_align.sh script)
##Counts are separated into unique, R and P. The latter two are PAtChER-aligned read types:
##R - Region-based alignment; one hit within the region defined by Ds
##P - Proximity-based alignment; more than one hit within Ds, nearest is used.

##Count data used here to measure enrichment in H3K27ac HiChIP data
##Read counts that would've come from the artificially duplicate region are estimated assuming no enrichment


##read count data

parse.files = function(files) {
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
	
	return(counts)
}

input = parse.files(list.files('./counts', pattern='^ESC', full.names=T))
input[,-1] = input[,-1]/1.3  #normalise by total read count
ip = parse.files(list.files('./counts', pattern='^mES_HiChIP_H3K27ac', full.names=T))


##get distances between original and duplicated ROIs

gtf = read.delim('original_ROIs.gtf', header=F, as.is=T)
attrib = strsplit(gtf$V9, split='[; ]')
gtf.id = unlist(lapply(attrib, function(x) x[5]))
gtf.dist = unlist(lapply(attrib, function(x) x[2]))

dist = as.numeric(gtf.dist[match(ip$roi, gtf.id)])/1000


##calculate enrichment

enrichment = function(prefix) {	
	ori.cols = grepl(prefix,colnames(ip)) & grepl('ori',colnames(ip))
	dup.cols = grepl(prefix,colnames(ip)) & grepl('dup',colnames(ip))
	
	ori.ip = ip$unique_oriCounts +  #unique counts
		rowSums(ip[,ori.cols]) +  #RP counts
		rowSums(input[,dup.cols])  #misassigned reads that would've come from duplicated region
	
	ori.input = input$unique_oriCounts +
		rowSums(input[,ori.cols]) +
		rowSums(input[,dup.cols])
	
	dup.ip = ip$unique_dupCounts +  #unique counts
		rowSums(ip[,dup.cols]) +  #RP counts
		rowSums(input[,ori.cols])  #correct reads that would've come from duplicated region
	
	dup.input = input$unique_dupCounts +
		rowSums(input[,dup.cols]) +
		rowSums(input[,ori.cols])	
	
	out = data.frame(roi = ip$roi,
		ori = log2(ori.ip+1) - log2(ori.input+1),
		dup = log2(dup.ip+1) - log2(dup.input+1),
		nondup = log2(ip$nondup_counts+1) - log2(input$nondup_counts+1),
		dist)
	return(out)
}


##tidy data 

library(tidyverse)

enr20 = enrichment('20kb')
enr20t = tibble(roi = rep(enr20$roi,2),
	type = rep(c('ori','dup'), each=nrow(enr20)),
	dup = c(enr20$ori, enr20$dup),
	nondup = rep(enr20$nondup,2),
	dist = rep(enr20$dist,2))


##plots

quartz(w=4, h=3.5)
enr20t %>% filter(dist>=10) %>%
ggplot(aes(nondup,dup,colour=type)) +
	theme_classic() +
	geom_point(size=0.8) +
	geom_smooth(method=lm)  #Figure 3C

quartz(w=4, h=3.5)
mutate(enr20t, dclass=cut(dist, c(0,10,20,50))) %>%
filter(nondup>=1) %>%
ggplot(aes(dclass,dup,fill=type)) +
	theme_classic() +
	geom_boxplot(outlier.shape=NA) +
	ylim(-0.2,2.5)  #Figure 3D

