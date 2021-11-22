##Process read counts at original and duplicated positions for regions involved in chromatin loops
##Duplicated regions were placed either at the end of loops ("long") or 50kb away ("prox")



library(tidyverse)


##read count data

long = inner_join(read_tsv('MCF7-long_ori_counts.txt', col_names=c('loop','ori')),
				read_tsv('MCF7-long_dup_counts.txt', col_names=c('loop','dup')),
				by='loop') %>% 
		add_column(type='long')

prox = inner_join(read_tsv('MCF7-prox_ori_counts.txt', col_names=c('loop','ori')),
				read_tsv('MCF7-prox_dup_counts.txt', col_names=c('loop','dup')),
				by='loop') %>% 
		add_column(type='prox')

data = rbind(long, prox) %>% mutate(acc = ori/(ori+dup)*100)


##plot accuracy (Supp Fig 3D)

quartz(w=2.5,h=3)
data %>%
ggplot(aes(x=type, y=acc)) +
	theme_classic() +
	geom_boxplot(outlier.shape=NA, fill='orange',width=0.5) +
	ylab('Accuracy (%)') +
	labs(fill='')
	

##stats

wilcox.test(data$acc[data$type=='long'], data$acc[data$type=='prox'])