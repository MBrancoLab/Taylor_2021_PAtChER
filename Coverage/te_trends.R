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

plot_trend = function(data, dist=5) {
	mutate(data, log.rpm=log2(rpm)) %>%
	ggplot(aes(bins, log.rpm, color=type)) +
		theme_classic() +
		ylab('log2 mean RPM') +
		labs(color='') +
		geom_line(size=1) +
		scale_colour_manual(values=c('orange','grey')) +
		scale_x_continuous(name='', breaks=c(0,20,80,100), labels=c('-2kb','start','end','+2kb'))
}


##plots

mervl = parse_file('ESC-mervl_trend.txt')
quartz(w=4.5,h=2.5)
plot_trend(mervl)  #Figure 1E

l1hs = parse_file('P2102Ep-l1_trend.txt')
quartz(w=4.5,h=2.5)
plot_trend(l1hs)  #Figure 1E

