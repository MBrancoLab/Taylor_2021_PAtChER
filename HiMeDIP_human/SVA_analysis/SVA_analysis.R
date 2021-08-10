##Process counts from HiMeDIP at SVAs and RNA-seq at genes (using htseq-counts)
##Merge the two datasets based on proximity of gene TSSs to SVAs (TSS_closest_SVA.bed, generated with bedtools)


library(tidyverse)


##file parsing function

parse.files = function(counts.file, closest.tss.gtf, rna.counts.file) {
		
	##normalise meDIP read counts
	counts = read_tsv(counts.file)
	total = summarise_at(counts, vars(-id), sum)
	norm = counts %>%
		mutate(ip.a.rpm = ip.a.counts/total$ip.a.counts*1e6,
		input.a.rpm = input.a.counts/total$input.a.counts*1e6,
		ip.u.rpm = ip.u.counts/total$ip.u.counts*1e6,
		input.u.rpm = input.u.counts/total$input.u.counts*1e6) %>%
		mutate(norm.all = log2(ip.a.rpm+0.01) - log2(input.a.rpm+0.01),
			norm.uni = log2(ip.u.rpm+0.01) - log2(input.u.rpm+0.01)) %>%
		select(-ip.a.counts, -input.a.counts, -ip.u.counts, -input.u.counts)
	
	##join TSS info and meDIP counts
	tss = read_tsv(closest.tss.gtf, col_names=FALSE) %>%
		mutate(id=unlist(lapply(strsplit(X15,'\"'), function(x) x[4])),
				sva.length = X11-X10) %>%
		select(X4, X5, id, sva.length, X16)
	colnames(tss) = c('gene_id','gene_name','id','sva.length','dist')
	tss.dip = inner_join(tss, norm, by='id')

	##add RNA-seq counts
	rna = read_tsv(rna.counts.file, col_names=c('gene_id','counts')) %>%
		mutate(rna.rpm = log2(1e6*(counts+1)/(sum(counts)+1))) %>%
		select(-counts)
	out = inner_join(tss.dip, rna, by='gene_id')

	return(out)
}


##parse files

ep = parse.files(counts.file = '2102Ep-SVAcounts.txt',
	closest.tss.gtf = 'TSS_closest_SVA.bed',
	rna.counts.file = '2102Ep-geneCounts.txt')

mcf = parse.files(counts.file = 'MCF7-SVAcounts.txt',
	closest.tss.gtf = 'TSS_closest_SVA.bed',
	rna.counts.file = 'MCF7-geneCounts.txt')


##get differential methylation and expression

all = inner_join(ep %>% select(gene_id, gene_name, id, sva.length, dist,
								rpm.a.ep=input.a.rpm, rpm.u.ep=input.u.rpm,
								norm.a.ep=norm.all, norm.u.ep=norm.uni,
								rna.ep=rna.rpm),
				mcf %>% select(gene_id, gene_name, id, sva.length, dist,
								rpm.a.mcf=input.a.rpm, rpm.u.mcf=input.u.rpm,
								norm.a.mcf=norm.all, norm.u.mcf=norm.uni,
								rna.mcf=rna.rpm),
	by=c('gene_id','gene_name','id','sva.length','dist')) %>%
	mutate(met.a.diff = norm.a.ep-norm.a.mcf,
			met.u.diff = norm.u.ep-norm.u.mcf,
			rna.diff = rna.ep-rna.mcf,
			sva = substr(id,1,5))


##select specific subgroups

sva100 = all %>% filter(rna.ep>0 | rna.mcf>0,
					dist<100000)
					
#def100 = all %>% filter(rna.ep>0 | rna.mcf>0,
#					dist<100000,
#					sva %in% c('SVA_D', 'SVA_E', 'SVA_F'))

#abc100 = all %>% filter(rna.ep>0 | rna.mcf>0,
#					dist<100000,
#					sva %in% c('SVA_A', 'SVA_B', 'SVA+C'))


##plot

t = 0.02

quartz(w=3.5, h=3.5)
sva100 %>% filter(sva.length>1000, rpm.a.mcf>t, rpm.a.ep>t) %>%
ggplot(aes(met.a.diff, rna.diff)) +
	theme_classic() +
	geom_point(size=0.5) +
	xlab('log2 meDIP FC') +
	ylab('log2 RNA FC') +
	geom_hline(yintercept=0, colour='blue', linetype='dashed') +
	geom_vline(xintercept=0, colour='blue', linetype='dashed')  #Supp Fig 5A




	