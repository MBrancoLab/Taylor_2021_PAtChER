##Generate a minigenome of selected gene promoters, with duplications at set distances
##Generate annotation tracks for the minigenomes
##Requires HOMER to extract sequences from mm10 reference genome


##get sequences
##regions are 203kb (3kb promoter ±100kb either side)

homer.path = '~/Documents/homer/'
system(paste(homer.path, 'bin/homerTools extract gene_promoters.bed ',
 homer.path, 'data/genomes/mm10/ -fa > gene_promoters.fa', sep=''))


##read in sequences

fa = scan('gene_promoters.fa',character())
genome = as.list(fa[seq(2,length(fa),2)])
names(genome) = paste('prom',1:(length(fa)/2),sep='')


##set distances between duplicate sequences

d = c(seq(3000,15000,1000),20000,30000,50000)
dset = rep(d,each=60)


##make genome with duplicates
##duplicated regions of interest (ROI) are 3kb

dup.gen = list()
for (i in 1:length(dset)) {
	seq = c(substr(genome[[i]],1,103000+dset[i]), #to the end of the ROI + spacer
	 substr(genome[[i]],100001,103000), #add duplicated ROI
	 substr(genome[[i]],103001+dset[i],nchar(genome[[i]]))) #add the remaining sequence
	dup.gen[[i]] = paste(seq,collapse='')
}

dup.chrom = paste(unlist(dup.gen),collapse='')
write('>chrD','promoter_minigenome.fa')
write(dup.chrom,'promoter_minigenome.fa',append=T)


##make non-duplicated minigenome

nondup.chrom = paste(unlist(genome),collapse="")
write('>chrO','promoter_minigenome_ori.fa')
write(nondup.chrom,'promoter_minigenome_ori.fa',append=T)


##make annotation files

make.gtf = function(chr_name, start, end, distance, filename) {
	attribute = paste('distSet "',distance,
	 '"; ID "prom',1:length(start),
	 '"',sep='')
	
	gtf = data.frame(seqname=rep(chr_name,length(start)),
	 source=rep('Minigenome',length(start)),
	 feature=rep('Promoter',length(start)),
	 as.integer(start),
	 as.integer(end),
	 score=rep('.',length(start)),
	 strand=rep('.',length(start)),
	 frame=rep('.',length(start)),
	 attribute)
	 
	write.table(gtf,filename,sep='\t',quote=F,col.names=F,row.names=F)
}

ori.start = seq(100001,206000*length(dup.gen),206000)
ori.end = ori.start+2999
make.gtf('chrD', ori.start, ori.end, dset, filename='original_ROIs.gtf')

dup.start = ori.end+dset+1
dup.end = dup.start+2999
make.gtf('chrD', dup.start, dup.end, dset, filename='duplicated_ROIs.gtf')

nondup.start = seq(100001,203000*length(dup.gen),203000)
nondup.end = nondup.start+2999
make.gtf('chrO', nondup.start, nondup.end, dset, filename='nonduplicated_ROIs.gtf')
