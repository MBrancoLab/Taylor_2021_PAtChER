##Process H3K27ac peaks detected with MACS2 (peak_detection.sh)
##macs2 bdgbroadcall was run with different values of -C 


library(tidyverse)


##read intersections between ROIs and peaks

data = tibble(roi=character(), nondup=integer(), ori=integer(), dup=integer(), cutoff=integer(), type=character())

for (t in 1:10) { #values for -C
	all = read_tsv(paste('roi_peaks/patcher_intersection_',t,'.txt', sep=''), col_names=colnames(data)[1:4])
	all = add_column(all, cutoff=rep(t,nrow(all)), type=rep('all',nrow(all)))
	uni = read_tsv(paste('roi_peaks/unique_intersection_',t,'.txt', sep=''), col_names=colnames(data)[1:4])
	uni = add_column(uni, cutoff=rep(t,nrow(uni)), type=rep('uni',nrow(uni)))
	data = rbind(data, all, uni)
}

data = data %>% mutate(correct=(nondup>0 & ori>0), incorrect=(nondup>0 & dup>0),
	dist=as.integer(unlist(lapply(strsplit(roi,'\\"'), function(x) x[2]))))


##make ROCs separated by inter-repeat distance

roc.data = data %>% mutate(dfac = cut(dist, c(0,15000,50000))) %>%
	group_by(cutoff, dfac, type) %>%
	summarise(recov=sum(correct)/sum(nondup>0)*100, fdr=sum(incorrect)/length(incorrect)*100)

quartz(w=5, h=4)
roc.data %>% filter(type=='all') %>%
	ggplot(aes(fdr,recov,colour=dfac)) +
	theme_bw() +
	xlab('False positive rate') +
	ylab('True positive rate') +
	geom_point() +
	scale_colour_manual(values=c('orange','red'), name='Gap', labels=c('≤15kb','>15kb'))  #Supp Fig 3

quartz(w=5, h=4)
roc.data %>% filter(type=='uni') %>%
	ggplot(aes(fdr,recov,colour=dfac)) +
	theme_bw() +
	xlab('False positive rate') +
	ylab('True positive rate') +
	geom_point() +
	scale_colour_manual(values=c('grey','black'),name='Gap', labels=c('≤15kb','>15kb'))  #Supp Fig 3


##calculate AUC

auc = function(roc) {
	fpr = c(roc$fdr, 0)
	tpr = c(roc$recov, 0)
	fpr.diff = fpr[1:(length(fpr)-1)] - fpr[2:length(fpr)]
	tpr.mean = (tpr[1:(length(tpr)-1)] + tpr[2:length(tpr)])/2
	total.area = sum(fpr.diff*tpr.mean)
	return(total.area/10000)
}

roc.data %>% filter(dfac=='(0,1.5e+04]', type=='all') %>% auc()
roc.data %>% filter(dfac=='(1.5e+04,5e+04]', type=='all') %>% auc()
roc.data %>% filter(dfac=='(0,1.5e+04]', type=='uni') %>% auc()
roc.data %>% filter(dfac=='(1.5e+04,5e+04]', type=='uni') %>% auc()


##plot just -C 6 peak data

quartz(h=3.5,w=2.8)
roc.data %>% filter(cutoff==6) %>%
	ggplot(aes(dfac,recov,fill=type)) +
	theme_classic() +
	geom_bar(stat='identity', position=position_dodge(), colour='black') +
	scale_fill_manual(values=c('orange','grey')) +
	scale_x_discrete(labels=c('≤15kb','15kb')) +
	ylim(0,80) +
	ylab('True positive rate')  #Figure 3E

quartz(h=3.5,w=2.8)
roc.data %>% filter(cutoff==6) %>%
	ggplot(aes(dfac,fdr,fill=type)) +
	theme_classic() +
	geom_bar(stat='identity', position=position_dodge(), colour='black') +
	scale_fill_manual(values=c('orange','grey')) +
	scale_x_discrete(labels=c('≤15kb','15kb')) +
	ylim(0,80) +
	ylab('False positive rate')  #Figure 3E

