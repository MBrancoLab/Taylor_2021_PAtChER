 ##Analysis of coverage files produced by bedcov (get_coverage.sh script)
##Coverage files not included in the repository due to size


library("tidyverse")


##functions

get_data = function(all.file, uni.file) {
	cnames = c('chr','start','end','te','div','strand','count','cov')
	all = read_tsv(all.file, col_names=cnames) %>% mutate(perc=cov/(end-start)*100) %>% add_column(type='all') 
	uni = read_tsv(uni.file, col_names=cnames) %>% mutate(perc=cov/(end-start)*100) %>% add_column(type='uni') 

	return(list(all, uni))
}

plotcov = function(data, lmin=1000, div.breaks=c(0,3,6)) {
	sub = lapply(data, function(x) x %>% filter(end-start>lmin))
	
	rbind(sub[[1]], sub[[2]]) %>%
	mutate(div_fac = cut(div,c(div.breaks,max(div)), include.lowest=T)) %>%	
	ggplot(aes(div_fac, perc, fill=factor(type,levels=c('uni','all')))) +
		theme_classic() +
		xlab('% divergence') + ylab('% coverage') + labs(fill='') +
		scale_fill_manual(values=c('grey','orange')) +
		geom_violin(scale='width', width=0.8) +
		stat_summary(fun='median', geom='point', position=position_dodge(0.8))
}

plotte = function(data, sel.tes, lmin=1000) {
	sub = lapply(data, function(x) x %>% filter(te %in% sel.tes, end-start>lmin))

	rbind(sub[[1]], sub[[2]]) %>%	
	ggplot(aes(te, perc, fill=factor(type,levels=c('uni','all')))) +
		theme_classic() +
		xlab('') + ylab('% coverage') + labs(fill='') +
		scale_fill_manual(values=c('grey','orange')) +
		geom_violin(scale='width', width=0.8) +
		stat_summary(fun='median', geom='point', position=position_dodge(0.8))
}

scatter = function(data, sel.tes, lmin=1000) {
	sub = data[[1]] %>% add_column(uni=data[[2]]$perc) %>%
		filter(te %in% sel.tes, end-start>lmin)
	sub$dcol = densCols(sub$uni, sub$perc,colramp=hcl.colors)

	ggplot(sub, aes(uni,perc)) +
		theme_classic() +
		xlab('Unique coverage') + ylab('PAtChER coverage') +
		geom_point(size=0.3,aes(colour=dcol)) +
		scale_color_identity() +
		geom_abline(linetype='dashed')
}


##human

hs = get_data('P2102EP-allcov.txt','P2102EP-unicov.txt')

quartz(w=3.5,h=3)
plotcov(hs)  #Figure 1D

quartz(w=3.5,h=3)
plotte(hs, c('L1HS','L1PA2','SVA_A'))  #Supp Fig 1A

quartz(w=3,h=3)
scatter(hs, c('L1HS','L1PA2'))  #Supp Fig 1A


##mouse

mus = get_data('ESC-allcov.txt','ESC-unicov.txt')

quartz(w=3.5,h=3)
plotcov(mus)  #Figure 1D

quartz(w=3.5,h=3)
plotte(mus, c('L1MdTf_III','IAPEz-int','MERVL-int'))  #Supp Fig 1B

quartz(w=3,h=3)
scatter(mus, c('L1MdTf_I','L1MdTf_II','L1MdTf_III','L1MdA_I','L1MdA_II',
	'L1MdA_III','L1MdA_IV','L1MdA_V','L1MdA_VI','L1MdA_VII','L1MdGf_I','L1MdGf_II'))   #Supp Fig 1B


##files wih IAP or L1Md elements excluding those that are absent in 129 strain

ref129 = get_data('ESC_129refTE-allcov.txt','ESC_129refTE-unicov.txt')

quartz(w=3,h=3)
plotte(ref129, c('IAPEz-int','L1MdTf_III'))  #Supp Fig 1C


##effect of sequencing depth

full = hs
half = get_data('P2102EP-halfA.txt','P2102EP-halfU.txt')
quarter = get_data('P2102EP-quarterA.txt','P2102EP-quarterU.txt')

#full = mus
#half = get_data('ESC-halfA.txt','ESC-halfU.txt')
#quarter = get_data('ESC-quarterA.txt','ESC-quarterU.txt')

young.long = function(x,d) x %>% filter(end-start>1000, div<3) %>% add_column(depth=d)
full.sub = lapply(full, function(x) young.long(x,d='full'))
half.sub = lapply(half, function(x) young.long(x,d='half'))
quarter.sub = lapply(quarter, function(x) young.long(x,d='quarter'))
depth = rbind(quarter.sub[[1]],quarter.sub[[2]],half.sub[[1]],half.sub[[2]],full.sub[[1]],full.sub[[2]])

quartz(w=3.5,h=3)
ggplot(depth, aes(factor(depth,levels=c('quarter','half','full')), perc, fill=factor(type,levels=c('uni','all')))) +
	theme_classic() +
	xlab('') + ylab('% coverage') + labs(fill='') +
	scale_fill_manual(values=c('grey','orange')) +
	geom_violin(scale='width', width=0.8) +
	stat_summary(fun='median', geom='point', position=position_dodge(w=0.8))  #Supp Fig 1D
	