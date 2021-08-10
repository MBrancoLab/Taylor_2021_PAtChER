##Analyse genomic context of MERVL elements clustered by deeptools (mervl_heatmaps.sh script)
##Information on neares gene from closestBed (MERVL-int_4kb_closestGene.txt file)
##RNA-seq gene counts generated using htseq-counts


library(tidyverse)


##MERVL-int file

mervl = read_tsv('MERVL-int_4kb_closestGene.txt', col_names=FALSE) %>%
	select(X1,X2,X3,X6,X13,X22,X23) %>%
	mutate(gene_id=unlist(lapply(strsplit(X22,split='"'), function(x) x[2])))
colnames(mervl) = c('chr','start','end','strand','cluster','gene','dist','gene_id')


##add RNA-seq data

rna = read_tsv('E14-geneCounts.txt', col_names=c('gene_id','counts')) %>%
	mutate(rna.rpm = log2(1e6*(counts+1)/(sum(counts)+1))) %>%
	select(-counts)

mervl = inner_join(mervl, rna, by='gene_id')


##plot distance

sum.dist = mervl %>% mutate(dcat = cut(dist, c(-1,0,1e4,1e5,max(dist)))) %>%
	group_by(cluster, dcat) %>%
	summarise(n=length(dcat))

quartz(w=3.5,h=3.5)
sum.dist %>% filter(cluster!='cluster_4') %>%
ggplot(aes(cluster, n, fill=dcat)) +
	theme_classic() +
	xlab('') +
	ylab('Fraction') +
	labs(fill='Distance') +
	scale_fill_manual(values=colorRampPalette(c('lightblue','navy'))(4),
		labels=c('overlapping','<10kb','10-100kb','>100kb')) +
	geom_bar(stat='identity', position=position_fill(reverse=TRUE), col='black', width=0.7)  #Figure 5H


##plot expression

quartz(w=3,h=3.5)
mervl %>% filter(dist<100000,cluster!='cluster_4') %>%
ggplot(aes(cluster,rna.rpm)) + 
	theme_classic() +
	geom_violin(fill='lightblue') +
	stat_summary(fun=median, geom='point') +
	xlab('') +
	ylab('log2 RPM')  #Figure 5I


##overlap with A/B compartment annotation

mervl %>% select(chr,start,end,strand,cluster) %>%
write_tsv(file='temp.bed',col_names=FALSE)

system('~/Documents/bedtools2/bin/intersectBed -wao -a temp.bed -b AB_compartments.txt > comp.bed')
comp = read_tsv('comp.bed', col_names=FALSE)
unlink(c('temp.bed','comp.bed'))


##plot A/B comparment data

sum.comp = comp %>% group_by(X5,X9) %>%
	summarise(n=length(X9))
colnames(sum.comp) = c('cluster','comp','n')

quartz(w=3.5,h=3.5)
sum.comp %>% filter(cluster!='cluster_4') %>%
ggplot(aes(cluster, n, fill=comp)) +
	theme_classic() +
	xlab('') +
	ylab('Fraction') +
	labs(fill='Compartment') +
	scale_fill_manual(values=c('grey','green','red'),
		labels=c('unnanotated','A','B')) +
	geom_bar(stat='identity', position=position_fill(), col='black', width=0.7)   #Figure 5J
