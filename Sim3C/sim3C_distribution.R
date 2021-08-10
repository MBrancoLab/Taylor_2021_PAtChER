##Compare real Hi-C contact frequencies (mapped by HiCUP) with Sim3C simulated data


library(tidyverse)


##real Hi-C mapped data (mESCs)

r1 = read_tsv('real-R1.bed.gz', col_names=c('chr1','start1','end1','id','score1','strand1')) %>%
	mutate(id = substr(id,1,nchar(id)-2))
r2 = read_tsv('real-R2.bed.gz', col_names=c('chr2','start2','end2','id','score2','strand2'))%>%
	mutate(id = substr(id,1,nchar(id)-2))

real = inner_join(r1, r2, by='id') %>%
	filter(chr1==chr2) %>%
	mutate(dist = abs(start1-start2))


##simulated data

s1 = read_tsv('sim-R1.bed.gz', col_names=c('chr1','start1','end1','id','score1','strand1'))
s2 = read_tsv('sim-R2.bed.gz', col_names=c('chr2','start2','end2','id','score2','strand2'))

sim = inner_join(s1, s2, by='id') %>%
	filter(chr1==chr2, grepl(':3C',id)) %>%
	mutate(dist = abs(start1-start2))

sim.sub = sim[1:nrow(real),] #match number of data points between the two sets


##plot

quartz(w=3.5, h=3.5)
real %>% filter(dist<100000) %>%
ggplot(aes(x=dist)) +
	theme_classic() +
	xlab('Distance (bp)') +
	ylab('Frequency') +
	stat_bin(breaks=seq(0,100000,5000), colour='black', fill='white') +
	stat_bin(data=filter(sim.sub,dist<100000), mapping=aes(x=dist),
			breaks=seq(0,100000,5000), geom='line', colour='blue', size=1)  #Supp Fig 2A
