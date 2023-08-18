#edited 12/28 for annotation purposes

clusters5 <- read.delim("~/clusters5.bed", header=FALSE, comment.char="#")
clusters4 <- read.delim("~/clusters4.bed", header=FALSE, comment.char="#")

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

totalgrange = makeGRangesFromDataFrame(clusters4, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE) #this annotation program requires the peaks to be in granges format.

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(totalgrange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")

annoGraph=plotAnnoBar(peakAnno)
annoGraph+scale_fill_manual(values=c("Promoter"="Red", "5' UTR"="Orange","3' UTR"="Yellow","1st Exon"="Green", "Other Exon"="Blue", "1st Intron"="Purple", "Other Intron"="Orchid", "Downstream (<=300)"="Brown", "Distal Intergenic"="Grey"))

C5_1 = clusters5 %>% filter(V13=="cluster_1")
C5_2 = clusters5 %>% filter(V13=="cluster_2")
C5_3 = clusters5 %>% filter(V13=="cluster_3")
C5_4 = clusters5 %>% filter(V13=="cluster_4")
C5_5 = clusters5 %>% filter(V13=="cluster_5")

C4_1 = clusters4 %>% filter(V13=="cluster_1")
C4_2 = clusters4 %>% filter(V13=="cluster_2")
C4_3 = clusters4 %>% filter(V13=="cluster_3")
C4_4 = clusters4 %>% filter(V13=="cluster_4")
library(GenomicRanges)
C4_1grange = makeGRangesFromDataFrame(C4_1, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE) #this annotation program requires the peaks to be in granges format.

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(C4_1grange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")

annoGraph=plotAnnoBar(peakAnno)
annoGraph+scale_fill_manual(values=c("Promoter"="Red", "5' UTR"="Orange","3' UTR"="Yellow","1st Exon"="Green", "Other Exon"="Blue", "1st Intron"="Purple", "Other Intron"="Orchid", "Downstream (<=300)"="Brown", "Distal Intergenic"="Grey"))

C4_2grange = makeGRangesFromDataFrame(C4_2, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE)
peakAnno2 <- annotatePeak(C4_2grange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
annoGraph=plotAnnoBar(peakAnno2)
annoGraph+scale_fill_manual(values=c("Promoter"="Red", "5' UTR"="Orange","3' UTR"="Yellow","1st Exon"="Green", "Other Exon"="Blue", "1st Intron"="Purple", "Other Intron"="Orchid", "Downstream (<=300)"="Brown", "Distal Intergenic"="Grey"))

C4_4grange = makeGRangesFromDataFrame(C4_4, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE)
peakAnno4 <- annotatePeak(C4_4grange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
annoGraph=plotAnnoBar(peakAnno4)
annoGraph+scale_fill_manual(values=c("Promoter"="Red", "5' UTR"="Orange","3' UTR"="Yellow","1st Exon"="Green", "Other Exon"="Blue", "1st Intron"="Purple", "Other Intron"="Orchid", "Downstream (<=300)"="Brown", "Distal Intergenic"="Grey"))

C4_3grange = makeGRangesFromDataFrame(C4_3, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE)
peakAnno3 <- annotatePeak(C4_3grange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
annoGraph=plotAnnoBar(peakAnno3)
annoGraph+scale_fill_manual(values=c("Promoter"="Red", "5' UTR"="Orange","3' UTR"="Yellow","1st Exon"="Green", "Other Exon"="Blue", "1st Intron"="Purple", "Other Intron"="Orchid", "Downstream (<=300)"="Brown", "Distal Intergenic"="Grey"))

C5_1grange = makeGRangesFromDataFrame(C5p2_1, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE)
peakAnnoC51 <- annotatePeak(C5_1grange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
annoGraph=plotAnnoBar(peakAnnoC51)
annoGraph+scale_fill_manual(values=c("Promoter"="Red", "5' UTR"="Orange","3' UTR"="Yellow","1st Exon"="Green", "Other Exon"="Blue", "1st Intron"="Purple", "Other Intron"="Orchid", "Downstream (<=300)"="Brown", "Distal Intergenic"="Grey"))

C5_2grange = makeGRangesFromDataFrame(C5p2_2, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE)
peakAnnoC52 <- annotatePeak(C5_2grange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
annoGraph=plotAnnoBar(peakAnnoC52)
annoGraph+scale_fill_manual(values=c("Promoter"="Red", "5' UTR"="Orange","3' UTR"="Yellow","1st Exon"="Green", "Other Exon"="Blue", "1st Intron"="Purple", "Other Intron"="Orchid", "Downstream (<=300)"="Brown", "Distal Intergenic"="Grey"))

C5_3grange = makeGRangesFromDataFrame(C5p2_3, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE)
peakAnnoC53 <- annotatePeak(C5_3grange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
annoGraph=plotAnnoBar(peakAnnoC53)
annoGraph+scale_fill_manual(values=c("Promoter"="Red", "5' UTR"="Orange","3' UTR"="Yellow","1st Exon"="Green", "Other Exon"="Blue", "1st Intron"="Purple", "Other Intron"="Orchid", "Downstream (<=300)"="Brown", "Distal Intergenic"="Grey"))

C5_4grange = makeGRangesFromDataFrame(C5p2_4, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE)
peakAnnoC54 <- annotatePeak(C5_3grange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
annoGraph=plotAnnoBar(peakAnnoC54)
annoGraph+scale_fill_manual(values=c("Promoter"="Red", "5' UTR"="Orange","3' UTR"="Yellow","1st Exon"="Green", "Other Exon"="Blue", "1st Intron"="Purple", "Other Intron"="Orchid", "Downstream (<=300)"="Brown", "Distal Intergenic"="Grey"))

C5_5grange = makeGRangesFromDataFrame(C5p2_5, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE)
peakAnnoC55 <- annotatePeak(C5_5grange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
annoGraph=plotAnnoBar(peakAnnoC55)
annoGraph+scale_fill_manual(values=c("Promoter"="Red", "5' UTR"="Orange","3' UTR"="Yellow","1st Exon"="Green", "Other Exon"="Blue", "1st Intron"="Purple", "Other Intron"="Orchid", "Downstream (<=300)"="Brown", "Distal Intergenic"="Grey"))

PA1 = as.data.frame(peakAnno)
PA2 = as.data.frame(peakAnno2)
PA3 = as.data.frame(peakAnno3)
PA4 = as.data.frame(peakAnno4)
write.table(annotated, "annotated.txt", row.names=FALSE, col.names=FALSE, sep="\t")
#graphs are under 12/20 in notebook
#############################################333
annotated <- read.delim("~/annotated.txt", header=FALSE)
BiocManager::install("JASPAR2020")
#BiocManager::install("motifmatchr")
#BiocManager::install("TFBSTools")

library(JASPAR2020)
library(motifmatchr)
library(TFBSTools)

opts <- list()
opts[["species"]] <- "Homo sapiens"
PFMatrixList <- getMatrixSet(JASPAR2020, opts) #this retrieves the JASPAR motifs, 117 TF associated motifs. #there we go, that makes a huge difference. Jaspar 2020 is much bigger than jasper 2014, who would expect that?

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

annotated <- read.delim("~/annotated.txt", header=FALSE)
annotated$uniq_id <- paste0(annotated$V2, annotated$V26) #because the gene name is technically not unique, ie there are multiple KDM5B peaks associated with a single gene.
peaks <- GRanges(seqnames=annotated$V1, ranges=IRanges(start=annotated$V2, end=annotated$V3))
peaks$uniq_ID = annotated$uniq_id
motifmatchs <- matchMotifs(PFMatrixList, peaks, genome="hg38")
allmotifs_matrix = motifMatches(motifmatchs)
allmotifs_matrix_df = as.data.frame(as.matrix(allmotifs_matrix)) #as.matrix is needed to convert it into a normal matrix, then it can be turned in a data frame.

IDs = ID(PFMatrixList)
names = name(PFMatrixList)
#should be in order
colnames(allmotifs_matrix_df) <- names
annotated_w_TFs = cbind(annotated, allmotifs_matrix_df) #order has been preserved so a simple binding works.
#r converts true to 1 and false to 0, therefore sum can calculate number true/false. 
sum(annotated_w_TFs$SPIB)
#6639
#What would be the high throughput version of this?
library(tidyverse)
total = annotated_w_TFs %>% pivot_longer(29:661, names_to="TF", values_to="is_TF_present") %>% group_by(TF) %>% summarize("total_true"=sum(is_TF_present))
counts = annotated_w_TFs %>% group_by(V15) %>% summarize(n=n()) #use this to label and to get percentages, normalize bc we'd expect there to be less 
cluster1 = annotated_w_TFs %>% filter(V15=="cluster_1") %>% pivot_longer(29:661, names_to="TF", values_to="is_TF_present") %>% group_by(TF) %>% summarize("cluster1_true_of1470"=sum(is_TF_present))
cluster2 = annotated_w_TFs %>% filter(V15=="cluster_2") %>% pivot_longer(29:661, names_to="TF", values_to="is_TF_present") %>% group_by(TF) %>% summarize("cluster2_true_of1454"=sum(is_TF_present))
cluster3 = annotated_w_TFs %>% filter(V15=="cluster_3") %>% pivot_longer(29:661, names_to="TF", values_to="is_TF_present") %>% group_by(TF) %>% summarize("cluster3_true_of5570"=sum(is_TF_present))
cluster4 = annotated_w_TFs %>% filter(V15=="cluster_4") %>% pivot_longer(29:661, names_to="TF", values_to="is_TF_present") %>% group_by(TF) %>% summarize("cluster4_true_of9074"=sum(is_TF_present))

TFs_summary = cbind(total, cluster1[2], cluster2[2], cluster3[2], cluster4[2])
TFs_summary$total_percent_true <- TFs_summary$total_true/17568
TFs_summary$cluster1_percent_true <- TFs_summary$cluster1_true_of1470/1470
TFs_summary$cluster2_percent_true <- TFs_summary$cluster2_true_of1454/1454
TFs_summary$cluster3_percent_true <- TFs_summary$cluster3_true_of5570/5570
TFs_summary$cluster4_percent_true <- TFs_summary$cluster4_true_of9074/9074

#okay now to represent w/ bar graph? would be a lot of columns though
#not sure how to represent it more clearly though?
TFs_summary_for_graphing <- TFs_summary %>% pivot_longer(7:11, names_to="cluster_info", values_to="percentage_true")
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
ggplot(data=TFs_summary_for_graphing, aes(y=reorder(TF, percentage_true), x=percentage_true)) + mytheme + geom_col(aes(fill=cluster_info), position=position_dodge())

TFs_summary_for_graphing2 = TFs_summary_for_graphing %>% ungroup() %>% group_by(TF) %>% mutate(variability=sd(percentage_true))
ggplot(data=TFs_summary_for_graphing2, aes(y=reorder(TF, variability), x=percentage_true)) + mytheme + geom_col(aes(fill=cluster_info), position=position_dodge())

TFs_summary
#get some interesting subsets:
#what's the average difference between clusters 1 and 3?
summary(TFs_summary$cluster1_percent_true - TFs_summary$cluster3_percent_true)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.04097  0.01572  0.03440  0.03608  0.05348  0.13413 
#what about some of the others?
summary(TFs_summary$cluster1_percent_true - TFs_summary$cluster2_percent_true)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.160910 -0.082631 -0.043188 -0.035940 -0.004599  0.323438 
summary(TFs_summary$cluster1_percent_true - TFs_summary$cluster4_percent_true)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.04303  0.01184  0.02817  0.04007  0.04996  0.40591 
summary(TFs_summary$cluster2_percent_true - TFs_summary$cluster3_percent_true)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.33649  0.02371  0.08294  0.07202  0.12798  0.27060 
summary(TFs_summary$cluster3_percent_true - TFs_summary$cluster4_percent_true)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.099959 -0.036908 -0.014247  0.003981  0.021702  0.418963 

#so some filtering...
#Groups I want: 
#1. 'cluster-independent' TFs. These genes would be pretty similar between the clusters, and probably would be ignored for further analysis
#2. 'probably just promoter' TFs vs 'dif in 1 and 3' TFs. Might be interesting.
#Figure out further from there?
TFs_summary_for_graphing2 = TFs_summary_for_graphing %>% ungroup() %>% group_by(TF) %>% mutate(variability=sd(percentage_true))
TFs2 = TFs_summary_for_graphing2 %>% pivot_wider(names_from = cluster_info, values_from = percentage_true)
summary(TFs2$variability)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.003936 0.027943 0.039563 0.043418 0.053775 0.189434 
#lets remove the half that are most similar
cluster_independent <- TFs2 %>% filter(variability < 0.045)
cluster_dependent <- TFs2 %>% filter(variability > 0.045)

#what's the average difference between clusters 1 and 3?
summary(cluster_dependent$cluster1_percent_true - cluster_dependent$cluster3_percent_true)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.04097  0.02517  0.05272  0.04996  0.07344  0.13413
#what about some of the others?
summary(cluster_dependent$cluster1_percent_true - cluster_dependent$cluster2_percent_true)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.16091 -0.10679 -0.08815 -0.05081 -0.04366  0.32344
summary(cluster_dependent$cluster1_percent_true - cluster_dependent$cluster4_percent_true)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.043031  0.001688  0.028900  0.053317  0.071287  0.405907 
summary(cluster_dependent$cluster2_percent_true - cluster_dependent$cluster3_percent_true)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.3365  0.1162  0.1361  0.1008  0.1661  0.2706 
summary(cluster_dependent$cluster3_percent_true - cluster_dependent$cluster4_percent_true)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.099959 -0.056543 -0.037772  0.003361  0.003650  0.418963 

#probably promoter driven associations?
prob_pro <- cluster_dependent %>% filter((cluster1_percent_true - cluster3_percent_true)<0.05)
summary(prob_pro$cluster1_percent_true - prob_pro$cluster2_percent_true)
summary(prob_pro$cluster2_percent_true - prob_pro$cluster4_percent_true)

dif_in_1n3 <- cluster_dependent %>% filter((cluster1_percent_true - cluster3_percent_true)>0.05)
summary(dif_in_1n3$cluster1_percent_true - dif_in_1n3$cluster2_percent_true)
summary(dif_in_1n3$cluster2_percent_true - dif_in_1n3$cluster4_percent_true)
#cluster independent, prob pro,dif_in_1n3
prob_pro_long <- prob_pro %>% pivot_longer(8:12, names_to="cluster_info", values_to="percentage_true")
ggplot(data=prob_pro_long, aes(y=TF, x=percentage_true)) + mytheme + geom_col(aes(fill=cluster_info), position=position_dodge())
dif_in_1n3_long <- dif_in_1n3 %>% pivot_longer(8:12, names_to="cluster_info", values_to="percentage_true")
ggplot(data=dif_in_1n3_long, aes(y=TF, x=percentage_true)) + mytheme + geom_col(aes(fill=cluster_info), position=position_dodge())

cluster_independent$subgroup <- rep("cluster_independent", dim(cluster_independent)[1])
cluster1_predominant <- TFs2 %>% filter(cluster1_percent_true>(cluster2_percent_true+variability))
cluster2_predominant <- TFs2 %>% filter(cluster2_percent_true>(cluster1_percent_true+variability))
similar_except_in_4 <- TFs2 %>% filter(cluster2_percent_true < (cluster1_percent_true+variability) & cluster1_percent_true < (cluster2_percent_true+variability))

cluster1_predominant$subgroup <- rep("cluster1_predominant", dim(cluster1_predominant)[1])
cluster2_predominant$subgroup <- rep("cluster2_predominant", dim(cluster2_predominant)[1])
similar_except_in_4$subgroup <- rep("cluster_independent", dim(similar_except_in_4)[1])
classed_anno <- rbind(cluster1_predominant, cluster2_predominant, similar_except_in_4)
write.table(
  classed_anno,
  file="classed_anno.txt",
  row.names=FALSE,
  quote=FALSE,
  sep="\t"
)
