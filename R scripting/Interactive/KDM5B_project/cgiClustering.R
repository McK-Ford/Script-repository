library(tidyverse)
mytheme = theme(
  panel.background=element_rect(fill="white"),
  text=element_text(color="black",face="bold",family="sans"),
  axis.text=element_text(color="black"),
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
  axis.ticks=element_line(color="black"),
  plot.margin=unit(c(0.25, 0.25, .25, 0.25),"cm"),
  plot.title=element_text(vjust=2),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  axis.title.y=element_text(vjust=0),
  axis.title.x=element_text(vjust=0),
  panel.border = element_blank(),
  axis.line=element_line()
)

#get clusters by w or wout CGI:
annoDF5_wCGI$CGI_present = rep("CGI", dim(annoDF5_wCGI)[1])
annoDF5_woutCGI$CGI_present = rep("No CGI", dim(annoDF5_woutCGI)[1])
annoDF5_union = rbind(annoDF5_wCGI, annoDF5_woutCGI)

annoDF_sum = annoDF5_union %>% add_count(V15) %>% group_by(V15, V28, CGI_present) %>% mutate(count_anno_by_clust=n()) %>% mutate(per_id=count_anno_by_clust/n) %>% summarise(percentage=first(per_id)) #then should summarize into 3 column table.
plot1 = ggplot(data=annoDF_sum, aes(y=percentage, x=V15)) + mytheme + geom_col(aes(fill=V28)) + xlab("cluster") +labs(fill="annotation") + facet_wrap(~CGI_present)

annoDF4_wCGI$CGI_present = rep("CGI", dim(annoDF4_wCGI)[1])
annoDF4_woutCGI$CGI_present = rep("No CGI", dim(annoDF4_woutCGI)[1])
annoDF4_union = rbind(annoDF4_wCGI, annoDF4_woutCGI)

introns4_union <- annoDF4_union %>%
  filter(grepl('Intron(.*)', V16))
not_introns4_union  <- annoDF4_union %>%
  filter(!grepl('Intron(.*)', V16))
exons4_union <- not_introns4_union %>%
  filter(grepl('Exon(.*)', V16))
not_intron_or_exon4_union <- not_introns4_union %>%
  filter(!grepl('Exon(.*)', V16))
first_intron4_union  <- introns4_union %>%
  filter(grepl('(.*)intron 1(.*)', V16))
not_first_intron4_union <- introns4_union %>%
  filter(!grepl('(.*)intron 1(.*)', V16))
first_exon4_union <- exons4_union %>%
  filter(grepl('(.*)exon 1(.*)', V16))
not_first_exon4_union <- exons4_union %>%
  filter(!grepl('(.*)exon 1(.*)', V16))
not_intron_or_exon4_union$simple_anno = not_intron_or_exon4_union$V16
first_intron4_union$simple_anno = rep("1st Intron", dim(first_intron4_union)[1])
first_exon4_union$simple_anno = rep("1st Exon", dim(first_exon4_union)[1])
not_first_intron4_union$simple_anno = rep("Other Intron", dim(not_first_intron4_union)[1])
not_first_exon4_union$simple_anno = rep("Other Exon", dim(not_first_exon4_union)[1])

annoDF4_union= rbind(not_intron_or_exon4_union, not_first_exon4_union, not_first_intron4_union, first_exon4_union, first_intron4_union)

annoDF_sum4 = annoDF4_union %>% add_count(V15) %>% group_by(V15, simple_anno, CGI_present) %>% mutate(count_anno_by_clust=n()) %>% mutate(per_id=count_anno_by_clust/n) %>% summarise(percentage=min(per_id)) #then should summarize into 3 column table.
plot2 = ggplot(data=annoDF_sum4, aes(y=percentage, x=V15)) + mytheme + geom_col(aes(fill=simple_anno)) + xlab("cluster") +labs(fill="annotation") + facet_wrap(~CGI_present)

pdf("ClusterCGIinfo.pdf")
plot1
plot2
dev.off()
