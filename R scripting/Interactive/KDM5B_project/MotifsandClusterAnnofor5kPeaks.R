#edited 12/28 for annotation purposes

#BiocManager::install("ChIPseeker")
library(ChIPseeker)
library(BSgenome)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)

library(EnsDb.Hsapiens.v75)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
#BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

totalgrange = makeGRangesFromDataFrame(point25percent_peaks_4clust, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE) #this annotation program requires the peaks to be in granges format.

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(totalgrange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")

annoGraph=plotAnnoBar(peakAnno)
annoGraph+scale_fill_manual(values=c("Promoter"="Red", "5' UTR"="Orange","3' UTR"="Yellow","1st Exon"="Green", "Other Exon"="Blue", "1st Intron"="Purple", "Other Intron"="Orchid", "Downstream (<=300)"="Brown", "Distal Intergenic"="Grey"))

annoDF = as.data.frame(peakAnno)

totalgrange = makeGRangesFromDataFrame(point25percent_peaks_5clust, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE) #this annotation program requires the peaks to be in granges format.

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(totalgrange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")

annoGraph=plotAnnoBar(peakAnno)
annoGraph+scale_fill_manual(values=c("Promoter"="Red", "5' UTR"="Orange","3' UTR"="Yellow","1st Exon"="Green", "Other Exon"="Blue", "1st Intron"="Purple", "Other Intron"="Orchid", "Downstream (<=300)"="Brown", "Distal Intergenic"="Grey"))

annoDF5 = as.data.frame(peakAnno)

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

library(tidyverse)
annoDF_sum = annoDF %>% add_count(V13) %>% group_by(V13, annotation) %>% mutate(count_anno_by_clust=n()) #no that doesn't work bc the intron exon labels are messy, but it works for everything except them so I just need to make them less messy

introns <- annoDF %>%
  filter(grepl('Intron(.*)', annotation))
not_introns  <- annoDF %>%
  filter(!grepl('Intron(.*)', annotation))
exons <- not_introns %>%
  filter(grepl('Exon(.*)', annotation))
not_intron_or_exon <- not_introns %>%
  filter(!grepl('Exon(.*)', annotation))
first_intron  <- introns %>%
  filter(grepl('(.*)intron 1(.*)', annotation))
not_first_intron <- introns %>%
  filter(!grepl('(.*)intron 1(.*)', annotation))
first_exon  <- exons %>%
  filter(grepl('(.*)exon 1(.*)', annotation))
not_first_exon <- exons %>%
  filter(!grepl('(.*)exon 1(.*)', annotation))
not_intron_or_exon$simple_anno = not_intron_or_exon$annotation
first_intron$simple_anno = rep("1st Intron", dim(first_intron)[1])
first_exon$simple_anno = rep("1st Exon", dim(first_exon)[1])
not_first_intron$simple_anno = rep("Other Intron", dim(not_first_intron)[1])
not_first_exon$simple_anno = rep("Other Exon", dim(not_first_exon)[1])

annoDF= rbind(not_intron_or_exon, not_first_exon, not_first_intron, first_exon, first_intron)
annoDF_sum = annoDF %>% add_count(V13) %>% group_by(V13, simple_anno) %>% mutate(count_anno_by_clust=n()) %>% mutate(per_id=count_anno_by_clust/n) %>% summarise(percentage=first(per_id)) #then should summarize into 3 column table.
ggplot(data=annoDF_sum, aes(y=percentage, x=simple_anno)) + mytheme + geom_col(aes(fill=V13), position=position_dodge()) + xlab("annotation") +labs(fill="cluster")

introns <- annoDF5 %>%
  filter(grepl('Intron(.*)', annotation))
not_introns  <- annoDF5 %>%
  filter(!grepl('Intron(.*)', annotation))
exons <- not_introns %>%
  filter(grepl('Exon(.*)', annotation))
not_intron_or_exon <- not_introns %>%
  filter(!grepl('Exon(.*)', annotation))
first_intron  <- introns %>%
  filter(grepl('(.*)intron 1(.*)', annotation))
not_first_intron <- introns %>%
  filter(!grepl('(.*)intron 1(.*)', annotation))
first_exon  <- exons %>%
  filter(grepl('(.*)exon 1(.*)', annotation))
not_first_exon <- exons %>%
  filter(!grepl('(.*)exon 1(.*)', annotation))
not_intron_or_exon$simple_anno = not_intron_or_exon$annotation
first_intron$simple_anno = rep("1st Intron", dim(first_intron)[1])
first_exon$simple_anno = rep("1st Exon", dim(first_exon)[1])
not_first_intron$simple_anno = rep("Other Intron", dim(not_first_intron)[1])
not_first_exon$simple_anno = rep("Other Exon", dim(not_first_exon)[1])

annoDF5= rbind(not_intron_or_exon, not_first_exon, not_first_intron, first_exon, first_intron)
annoDF5_sum = annoDF5 %>% add_count(V13) %>% group_by(V13, simple_anno) %>% mutate(count_anno_by_clust=n()) %>% mutate(per_id=count_anno_by_clust/n) %>% summarise(percentage=first(per_id)) #then should summarize into 3 column table.
ggplot(data=annoDF5_sum, aes(y=percentage, x=simple_anno)) + mytheme + geom_col(aes(fill=V13), position=position_dodge()) + xlab("annotation") +labs(fill="cluster")

#write table
write.table(annoDF, "annoDF4.txt", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(annoDF5, "annoDF5.txt", quote=FALSE, col.names=FALSE, row.names=FALSE)
