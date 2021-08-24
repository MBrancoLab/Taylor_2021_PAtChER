##PAtChER mapping accuracy from alignment of Sim3C data
##Accuracy data tables generated using matched_reads.py (see sim3c.sh script)
##Data on mean accuracy per element generated from bigWig accuracy track (mean_accuracy.sh script)
##bigWig accuracy track not included in repository due to size


library('tidyverse')

hs = read_tsv('accuracy_hg38.txt') %>% mutate(acc=matched/total*100)
mm = read_tsv('accuracy_mm10.txt') %>% mutate(acc=matched/total*100)


##barplot for overall accuracy

bplot = function(data) {
	data %>% filter(type %in% c('all','recov')) %>%
	ggplot(aes(type, acc, fill=dist)) +
		theme_classic() +
		xlab('') + ylab('Accuracy') +
		geom_bar(stat='identity', width=0.8, col='black', position=position_dodge()) +
		scale_fill_manual(values=hcl.colors(3)) +
		geom_hline(yintercept=data$acc[data$type=='random'], linetype='dashed')
}

quartz(w=3.5,h=3)
bplot(hs)  #Figure 2A

quartz(w=3.5,h=3)
bplot(mm)  #Figure 2A



##reads vs accuracy for selected families

fam.plot = function(data, sel.te) {	
	data %>% filter(type %in% sel.te) %>%
		group_by(type) %>% mutate(norm=total/min(total)) %>%
	
	ggplot(aes(norm, acc, color=type)) +
		theme_classic() +
		theme(panel.grid.major=element_line(size = 0.3,linetype = 'solid',colour = "grey")) +
		geom_point(aes(shape=dist)) +
		ylim(0,100) +
		xlab('Fold change in reads') +
		ylab('Accuracy') +
		geom_line()
}

quartz(w=3.5,h=3)
fam.plot(hs, sel.te=c('l1hs','l1pa2','sva_a','sva_f','hervk_int'))   #Figure 2B

quartz(w=3.5,h=3)
fam.plot(mm, sel.te=c('l1tf','l1a','etnerv-int','iapez-int','mervl-int'))  #Figure 2B



##accuracy per element

hs.el = rbind(read_tsv('accuracy_families/L1HS_1kb-meanAcc.txt') %>% add_column(te='L1HS'),
			read_tsv('accuracy_families/L1PA2_1kb-meanAcc.txt') %>% add_column(te='L1PA2'),
			read_tsv('accuracy_families/SVA_A_1kb-meanAcc.txt') %>% add_column(te='SVA_A'),
			read_tsv('accuracy_families/SVA_F_1kb-meanAcc.txt') %>% add_column(te='SVA_F'))
colnames(hs.el) = c('chr','start','end','acc','te')

mm.el = rbind(read_tsv('accuracy_families/L1Tf_1kb-meanAcc.txt') %>% add_column(te='L1Tf'),
			read_tsv('accuracy_families/L1A_1kb-meanAcc.txt') %>% add_column(te='L1A'),
			read_tsv('accuracy_families/IAPEz-int_1kb-meanAcc.txt') %>% add_column(te='IAPEz-int'),
			read_tsv('accuracy_families/MERVL-int_1kb-meanAcc.txt') %>% add_column(te='MERVL-int'))
colnames(mm.el) = c('chr','start','end','acc','te')


element.plot = function(data) {
	ggplot(data, aes(te,acc*100)) +
		theme_classic() +
		geom_violin(fill='orange', scale='width', width=0.8) +
		xlab('') +
		ylab('Accuracy (%)') +
		stat_summary(fun='median', geom='point', position=position_dodge(0.8))
}


quartz(w=3,h=3.5)
element.plot(hs.el)  #Supp Fig 2B

quartz(w=3,h=3.5)
element.plot(mm.el)  #Supp Fig 2B

