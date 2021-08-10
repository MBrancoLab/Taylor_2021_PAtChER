##Process read counts from HiMeDIP (count_reads.sh script)
##Compare with bisulphite data summarised in GTF files


library(tidyverse)


##file parsing function

parse.files = function(gtf.file, counts.file, read.type='all') {
	
	gtf = read_tsv(gtf.file, col_names=c('roi.chr','source','feature','roi.start',
		'roi.end','met','strand','frame','l1.info')) %>%
		select(-source, -feature, -frame) %>%
		mutate(id=unlist(lapply(strsplit(l1.info,'\"'), function(x) x[4])),
			met.cat=cut(met,c(0,25,75,100),include.lowest=TRUE))
	
	counts = read_tsv(counts.file)
	total = summarise_at(counts, vars(-id), sum)
	out = inner_join(gtf, counts, by='id')
	
	if (read.type=='all') {
		out = mutate(out, ip.rpm = ip.a.counts/total$ip.a.counts*1e6,
			input.rpm = input.a.counts/total$input.a.counts*1e6)
	}
	if (read.type=='unique') {
		out = mutate(out, ip.rpm = ip.u.counts/total$ip.u.counts*1e6,
			input.rpm = input.u.counts/total$input.u.counts*1e6)
	}
	if (read.type=='recov') {
		out = mutate(out, ip.rpm = (ip.a.counts-ip.u.counts)/total$ip.a.counts*1e6,
			input.rpm = (input.a.counts-input.u.counts)/total$input.a.counts*1e6)
	}
	
	out	= mutate(out, norm = log2(ip.rpm+0.01) - log2(input.rpm+0.01))
	return(out)
}


##enrichment plot function

enr.plot = function(data, min.rpm=0.01) {
	data %>% filter(input.rpm>min.rpm) %>%
	ggplot(aes(met.cat, norm, fill=type)) +
		theme_classic() +
		geom_boxplot(outlier.shape=NA, position=position_dodge(0.8), width=0.7) +
		ylim(-3,4) +
		xlab('% methylation') +
		ylab('log2 enrichment') +
		labs(fill='') +
		scale_fill_discrete(labels=c('Internal',"5' end"))
}


##differential methylation plot function

diff.plot = function(data, min.rpm=0.01, roi='L1', colour='#F8766D') {
	data %>% filter(rpm.ep>min.rpm & rpm.mcf>min.rpm, type==roi) %>%
	ggplot(aes(bs.diff, medip.diff)) +
		theme_classic() +
		xlim(-100, 100) +
		ylim(-4, 6) +
		geom_point(size=0.5, colour=colour) +
		xlab('ATLAS-BS difference') +
		ylab('Hi-MeDIP difference') +
		geom_smooth(method=lm, colour='black')
}



##parse all read counts (Figure 4B/C)

ep.l1 = parse.files(gtf.file = '2102Ep_ATLAS_internal.gtf',
	counts.file = 'counts/2102Ep_L1counts.txt',
	read.type = 'all')
ep.up = parse.files(gtf.file = '2102Ep_ATLAS_5prime.gtf',
	counts.file = 'counts/2102Ep_UPcounts.txt',
	read.type = 'all')
ep = rbind(ep.l1, ep.up) %>% add_column(type = rep(c('L1','UP'), c(nrow(ep.l1), nrow(ep.up))))

mcf.l1 = parse.files(gtf.file = 'MCF7_ATLAS_internal.gtf',
	counts.file = 'counts/MCF7_L1counts.txt',
	read.type = 'all')
mcf.up = parse.files(gtf.file = 'MCF7_ATLAS_5prime.gtf',
	counts.file = 'counts/MCF7_UPcounts.txt',
	read.type = 'all')
mcf = rbind(mcf.l1, mcf.up) %>% add_column(type = rep(c('L1','UP'), c(nrow(mcf.l1), nrow(mcf.up))))


##parse unique/non-unique counts depending on region (Supp Fig 4A/B)

#ep.l1 = parse.files(gtf.file = '2102Ep_ATLAS_internal.gtf',
#	counts.file = 'counts/2102Ep_L1counts.txt',
#	read.type = 'recov')
#ep.up = parse.files(gtf.file = '2102Ep_ATLAS_5prime.gtf',
#	counts.file = 'counts/2102Ep_UPcounts.txt',
#	read.type = 'unique')
#ep = rbind(ep.l1, ep.up) %>% add_column(type = rep(c('L1','UP'), c(nrow(ep.l1), nrow(ep.up))))

#mcf.l1 = parse.files(gtf.file = 'MCF7_ATLAS_internal.gtf',
#	counts.file = 'counts/MCF7_L1counts.txt',
#	read.type = 'recov')
#mcf.up = parse.files(gtf.file = 'MCF7_ATLAS_5prime.gtf',
#	counts.file = 'counts/MCF7_UPcounts.txt',
#	read.type = 'unique')
#mcf = rbind(mcf.l1, mcf.up) %>% add_column(type = rep(c('L1','UP'), c(nrow(mcf.l1), nrow(mcf.up))))


##plot enrichment

quartz(w=3.8,h=3)
enr.plot(ep)  #Figure 4B / Supp Fig 4A

quartz(w=3.8,h=3)
enr.plot(mcf)  #Figure 4B / Supp Fig 4A


##plot differential methylation

all = inner_join(ep %>% select(l1.info,met.ep=met,rpm.ep=input.rpm,norm.ep=norm,type),
	mcf %>% select(l1.info,met.mcf=met,rpm.mcf=input.rpm,norm.mcf=norm,type),
	by=c('l1.info','type')) %>%
	mutate(bs.diff=met.mcf-met.ep, medip.diff=norm.mcf-norm.ep)

quartz(w=3,h=3)
diff.plot(all)  #Figure 4C / Supp Fig 4B

quartz(w=3,h=3)
diff.plot(all, roi='UP', colour='#00BFC4')  #Figure 4C / Supp Fig 4B







