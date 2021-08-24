##Plot average read depth across selected TE families
##Source data produced by deeptools (te_trends.sh script)


library(tidyverse)

##functions

parse_file = function(file) {
	df = read.delim(file, row.names=1, skip=1, header=F)[,-1]
	tibble(bins=rep(as.numeric(df[1,]),2),
		rpm=c(as.numeric(df[2,]),as.numeric(df[3,])),
		type=rep(c('all','uni'),each=ncol(df)))
}

plot_trend = function(data) {
	mutate(data, log.rpm=log2(rpm)) %>%
	ggplot() +
		theme_classic() +
		ylab('log2 mean RPM') +
		labs(color='') +
		geom_line(aes(bins, log.rpm, color=type), size=1) +
		scale_colour_manual(values=c('orange','grey')) +
		scale_x_continuous(name='', breaks=c(0,20,80,100), labels=c('-2kb','start','end','+2kb'))
}

get_dpn = function(seq_file) {
	seq = scan(seq_file, character(), sep='\n', skip=1)
	dpn = gregexpr('GATC', paste(seq,collapse=''), ignore.case=TRUE)
	dpn.tib = tibble(norm=(dpn[[1]] + 2000) / (sum(nchar(seq))+4000) * 100, val=-4)
	return(dpn.tib)
}


##plots

mervl = parse_file('ESC-mervl_trend.txt')
mervl.dpn = get_dpn('MERVL_DF0003918.fa')
quartz(w=4.5,h=2.5)
pm = plot_trend(mervl)  #Figure 1E
pm + geom_segment(data=mervl.dpn, aes(x=norm,y=val,xend=norm,yend=val-0.5),col='red')

l1hs = parse_file('P2102Ep-l1_trend.txt')
l1hs.dpn = get_dpn('L1HS_L19092.fa')
quartz(w=4.5,h=2.5)
pl = plot_trend(l1hs)  #Figure 1E
pl + geom_segment(data=l1hs.dpn, aes(x=norm,y=val,xend=norm,yend=val-0.5),col='red')