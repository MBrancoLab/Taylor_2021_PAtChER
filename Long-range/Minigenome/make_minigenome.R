##Generate a minigenome of selected chromatin loops, with duplications at either the other end of the loop or 50kb away
##Generate annotation tracks for the minigenomes
##Requires HOMER to extract sequences from hg38 reference genome


##chromatin loops

cl = read.delim('chromatin_loops.bedpe', as.is=T)


##select 10kb interacting regions > 200kb away

filt = cl[cl$x2-cl$x1==10000 & cl$y1-cl$x1>200000,]


##randomly sample 200 regions

set.seed(143)
sel = filt[sample(1:nrow(filt),200),]


##make bed files with 200kb regions around loop anchors

bed1 = data.frame(chr=paste('chr',sel$chr1,sep=''), start=sel$x1-95000, end=sel$x2+95000)
bed2 = data.frame(chr=paste('chr',sel$chr2,sep=''), start=sel$y1-95000, end=sel$y2+95000)
write.table(bed1, 'loop_end1.bed', sep='\t', quote=F, row.names=F, col.names=F)
write.table(bed2, 'loop_end2.bed', sep='\t', quote=F, row.names=F, col.names=F)


##get sequences

homer.path = '~/Documents/homer/'

system(paste(homer.path, 'bin/homerTools extract loop_end1.bed ',
 homer.path, 'data/genomes/hg38/ -fa > loop_end1.fa', sep=''))
system(paste(homer.path, 'bin/homerTools extract loop_end2.bed ',
 homer.path, 'data/genomes/hg38/ -fa > loop_end2.fa', sep=''))


##read in sequences

read.seq = function(fa.file) {
	fa = scan(fa.file,character())
	seq = as.list(fa[seq(2,length(fa),2)])
	names(seq) = paste('loop',1:(length(fa)/2),sep='')
	return(seq)
}

end1 = read.seq('loop_end1.fa')
end2 = read.seq('loop_end2.fa')


##make genome with regions from one end duplicated on the other end
##duplicated regions of interest (ROI) are 3kb

long.gen = list()
for (i in 1:length(end1)) {
	seq = c(end1[[i]], #end1 is left intact
	 substr(end2[[i]],1,100000), #first half of end2
	 substr(end1[[i]],98501,101500), #add duplicated ROI from end1
	 substr(end2[[i]],100001,200000)) #second half of end2
	long.gen[[i]] = paste(seq,collapse='')
}

long.chrom = paste(unlist(long.gen),collapse='')
write('>chrL','longRange_minigenome.fa')
write(long.chrom,'longRange_minigenome.fa',append=T)


##make genome with duplicated regions within 50kb of end1

prox.gen = list()
for (i in 1:length(end1)) {
	seq = c(substr(end1[[i]],1,150000), #up to the point of insertion
	 substr(end1[[i]],98501,101500), #add duplicated ROI
	 substr(end1[[i]],150001,200000), #rest of end1 sequence
	 end2[[i]]) #add end2, to be similar to long-range minigenome
	prox.gen[[i]] = paste(seq,collapse='')
}

prox.chrom = paste(unlist(prox.gen),collapse='')
write('>chrP','proximal_minigenome.fa')
write(prox.chrom,'proximal_minigenome.fa',append=T)


##make annotation files

make.gtf = function(chr_name, start, end, filename) {
	attribute = paste('ID "loop',1:length(start),
	 '"',sep='')
	
	gtf = data.frame(seqname=rep(chr_name,length(start)),
	 source=rep('Minigenome',length(start)),
	 feature=rep('Loop',length(start)),
	 as.integer(start),
	 as.integer(end),
	 score=rep('.',length(start)),
	 strand=rep('.',length(start)),
	 frame=rep('.',length(start)),
	 attribute)
	 
	write.table(gtf,filename,sep='\t',quote=F,col.names=F,row.names=F)
}

ori.start = seq(98501,403000*length(long.gen),403000)
ori.end = ori.start+2999
make.gtf('chrL', ori.start, ori.end, filename='original_ROIs.gtf')

lr.start = seq(300001,403000*length(long.gen),403000)
lr.end = lr.start+2999
make.gtf('chrL', lr.start, lr.end, filename='longRange_dup_ROIs.gtf')

ori2.start = seq(98501,403000*length(prox.gen),403000)
ori2.end = ori2.start+2999
make.gtf('chrP', ori.start, ori.end, filename='original2_ROIs.gtf')

pr.start = seq(150001,403000*length(prox.gen),403000)
pr.end = pr.start+2999
make.gtf('chrP', pr.start, pr.end, filename='proximal_dup_ROIs.gtf')



