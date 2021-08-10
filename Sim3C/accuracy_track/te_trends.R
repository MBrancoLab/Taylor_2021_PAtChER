##Plot mean accuracy across selected TE families
##Data extracted from bigWig accuracy track (te_trends.sh script)
##bigWig accuracy track not included in repository due to size


library(tidyverse)

##functions

parse_file = function(file) {
	df = read.delim(file, row.names=1, skip=1, header=F)[,-1]
	tibble(bins=as.numeric(df[1,]),
		acc=as.numeric(df[2,])*100)
}

plot_trend = function(data) {
	ggplot(data, aes(bins, acc)) +
		theme_classic() +
		ylab('Accuracy (%)') +
		labs(color='') +
		geom_line(colour='#4A90E2', size=1) +
		ylim(0,100) +
		scale_x_continuous(name='', breaks=c(0,20,80,100), labels=c('-2kb','start','end','+2kb'))
}


##plot L1HS

l1hs = parse_file('l1hs_trend.txt')
quartz(w=4,h=2.8)
plot_trend(l1hs) + geom_vline(xintercept=c(20,80), linetype='dashed')  #Figure 2D



##plot mouse TEs

l1a = parse_file('l1a_trend.txt')
iap = parse_file('iapez_trend.txt')
mervl = parse_file('mervl_trend.txt')

quartz(w=4,h=2.8)
plot_trend(l1a) +
	geom_line(data=iap, mapping=aes(bins,acc), size=1, colour=grey(0.6)) +
	geom_line(data=mervl, mapping=aes(bins,acc), size=1, colour='red') +
	geom_vline(xintercept=c(20,80), linetype='dashed')  #Supp Fig 2D
