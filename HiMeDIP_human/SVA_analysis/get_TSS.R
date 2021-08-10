##extract TSS from gene annotation GTF
##output bed format

gene = read.delim('gencode.v38.genesOnly.gtf.gz', as.is=T, header=F)
gene.info = strsplit(gene$V9, split='; ')
gene.id = unlist(lapply(gene.info, function(x) {
		id = x[grep('gene_id',x)]
		return(substr(id, 9, nchar(id)))
	}))
gene.name = unlist(lapply(gene.info, function(x) {
		nam = x[grep('gene_name',x)]
		return(substr(nam, 11, nchar(nam)))
	}))

rev = gene$V7=='-'
tss = gene$V4
tss[rev] = gene$V5[rev]

tss.bed = data.frame(chr=gene$V1, start=tss, end=tss+1, gene.id, gene.name, strand=gene$V7)
write.table(tss.bed, 'TSS.bed', sep='\t', quote=F, row.names=F, col.names=F)
