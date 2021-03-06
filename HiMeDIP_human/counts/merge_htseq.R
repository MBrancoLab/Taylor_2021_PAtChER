
##Merges htseq-count files, to simplify
##Original htseq-count files not included in repository


merge.files = function(ip.all.file, input.all.file, ip.uni.file, input.uni.file) {
	ip.all = read.delim(ip.all.file, header=FALSE, as.is=TRUE)
	input.all = read.delim(input.all.file, header=FALSE, as.is=TRUE)
	ip.uni = read.delim(ip.uni.file, header=FALSE, as.is=TRUE)
	input.uni = read.delim(input.uni.file, header=FALSE, as.is=TRUE)
	
	merged = data.frame(id = ip.all$V1,
		ip.a.counts = ip.all$V2,
		input.a.counts = input.all$V2,
		ip.u.counts = ip.uni$V2,
		input.u.counts = input.uni$V2)
	
	return(merged)
}


##internal L1 counts

ep.l1 = merge.files(ip.all.file='Pool4_2102EP_IP-sorted_L1counts.txt',
	input.all.file='P2102EP_INPUT-sorted_L1counts.txt',
	ip.uni.file='Pool4_2102EP_IP-unique_L1counts.txt',
	input.uni.file='P2102EP_INPUT-unique_L1counts.txt')
write.table(ep.l1, '2102Ep_L1counts.txt', sep='\t', quote=FALSE, row.names=FALSE)

mcf.l1 = merge.files(ip.all.file='Pool_4_MCF7_IP-sorted_L1counts.txt',
	input.all.file='Pool_4_MCF7INPUT2-sorted_L1counts.txt',
	ip.uni.file='Pool_4_MCF7_IP-unique_L1counts.txt',
	input.uni.file='Pool_4_MCF7INPUT2-unique_L1counts.txt')
write.table(mcf.l1, 'MCF7_L1counts.txt', sep='\t', quote=FALSE, row.names=FALSE)


##upstream L1 counts

ep.up = merge.files(ip.all.file='Pool4_2102EP_IP-sorted_UPcounts.txt',
	input.all.file='P2102EP_INPUT-sorted_UPcounts.txt',
	ip.uni.file='Pool4_2102EP_IP-unique_UPcounts.txt',
	input.uni.file='P2102EP_INPUT-unique_UPcounts.txt')
write.table(ep.up, '2102Ep_UPcounts.txt', sep='\t', quote=FALSE, row.names=FALSE)

mcf.up = merge.files(ip.all.file='Pool_4_MCF7_IP-sorted_UPcounts.txt',
	input.all.file='Pool_4_MCF7INPUT2-sorted_UPcounts.txt',
	ip.uni.file='Pool_4_MCF7_IP-unique_UPcounts.txt',
	input.uni.file='Pool_4_MCF7INPUT2-unique_UPcounts.txt')
write.table(mcf.up, 'MCF7_UPcounts.txt', sep='\t', quote=FALSE, row.names=FALSE)
