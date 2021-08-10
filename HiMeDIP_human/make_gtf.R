##Convert processed bisulphite ATLAS-seq data to GTF with average methylation per region
##Files for 5' and internal L1 regions generated


make.gtf = function(cells, roi.start, roi.end, out.file) {

	##intersect with list of full-length (>5kb) L1 elements	
	system(paste('~/Documents/bedtools2/bin/intersectBed -wa -wb -a ',
		cells, '.all_L1.flanks_and_internal.hg38.bedGraph -b FL_L1_hg38.bed > ',
		'intersect_output.bed', sep=''))	
	ann = read.delim('intersect_output.bed', as.is=T, header=F)
	id = factor(paste(ann$V5,ann$V6,ann$V7,sep='_'))
	unlink('intersect_output.bed')
	
	##make list of methylation values per L1	
	l1 = list()
	for (i in 1:nlevels(id)) {
		l1[[i]] = ann[id==levels(id)[i],]
	}
	names(l1) = levels(id)
		
	##select region of interest within L1
	l1.int = lapply(l1, function(x) {
		if (x$V10[1]=='+') {
			left = x$V6[1] + roi.start
			right = x$V6[1] + roi.end
		} else {
			left = x$V7[1] - roi.end
			right = x$V7[1] - roi.start
		}
		x$roi.left = left
		x$roi.right = right
		return (x[x$V2>=left & x$V3<=right,])
	})
	
	##filter for minimum number of CpGs
	n.cpg = unlist(lapply(l1.int,nrow))
	l1.cpg = l1.int[n.cpg>=3]
	
	##get average methylation
	met = unlist(lapply(l1.cpg,function(x) mean(x$V4)))
	
	##generate gtf file
	attribute = paste('Family "', unlist(lapply(l1.cpg, function(x) x$V8[1])),
		 '"; ID "', names(l1.cpg),
		 '"', sep='')
	gtf = data.frame(seqname=unlist(lapply(l1.cpg, function(x) x$V1[1])),
		source=rep('ATLAS_BSseq',length(l1.cpg)),
		feature=rep('L1_ROI',length(l1.cpg)),
		start=unlist(lapply(l1.cpg, function(x) x$roi.left[1])),
		end=unlist(lapply(l1.cpg, function(x) x$roi.right[1])),
		score=met,
		strand=unlist(lapply(l1.cpg, function(x) x$V10[1])),
		frame=rep('.',length(l1.cpg)),
		attribute)
	 
	##write
	write.table(gtf, out.file,
		sep='\t', quote=F, col.names=F, row.names=F)
}



##files for internal L1 region

make.gtf('MCF7', roi.start=75, roi.end=210, out.file='MCF7_ATLAS_internal.gtf')
make.gtf('2102Ep', roi.start=75, roi.end=210, out.file='2102Ep_ATLAS_internal.gtf')


##files for 5' L1 region

make.gtf('MCF7', roi.start=0, roi.end=135, out.file='MCF7_ATLAS_5prime.gtf')
make.gtf('2102Ep', roi.start=0, roi.end=135, out.file='2102Ep_ATLAS_5prime.gtf')

