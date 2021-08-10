##Select reference B6 insertions that are absent in all sequenced 129 genomes
##Source data from Nellaker et al (PMID: 22703977)


library(tidyverse)

##read data
te = rbind(read_tsv('LINE.shared.tab'), read_tsv('IAP-I.shared.tab'))

##get reference insertions
ref = te %>% filter(C57B6_ref==1)

##get insertions absent in 129
del.129 = ref %>% filter(`129P2`=='DEL' | `129S1`=='DEL' | `129S5`=='DEL') %>%
	select(chrm, start, stop)

##write output
write_tsv(del.129, 'refTE_not129_mm9.bed', col_names=FALSE)