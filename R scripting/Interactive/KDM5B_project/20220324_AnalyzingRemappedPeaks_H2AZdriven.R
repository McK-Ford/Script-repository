#####
library(tidyverse)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
H3K4me3_driven_intersection <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/H3K4me3_driven_intersection.txt",
  header=FALSE, na.strings='.')
#####
H3K4me3_driven_intersection$V7[is.na(H3K4me3_driven_intersection$V7)] = "placeholder"
Nointersect_H3K4me3 = H3K4me3_driven_intersection %>% filter(V7=="placeholder")
H3K4me3_single_peak_per_peak = H3K4me3_driven_intersection %>%
  filter(V7!="placeholder") %>%
  group_by(V6, V7) %>% mutate(max_col = max(V11)) %>%
  ungroup() %>% filter(V11==max_col) %>% dplyr::select(-c("max_col"))
H3K4me3_single_peak_per_peak = rbind(H3K4me3_single_peak_per_peak,
                                     Nointersect_H3K4me3)
collapsed_H3K4me3 <- data.frame(
  do.call(paste, c(H3K4me3_single_peak_per_peak[1:6], sep=",")),
  H3K4me3_single_peak_per_peak$V7,
  do.call(paste, c(H3K4me3_single_peak_per_peak[8:13], sep=",")))
names(collapsed_H3K4me3) <- c("col_H3K4me3", "bed_id", "col_info")
collapsed_H3K4me3$bed_id <- str_remove_all(
  collapsed_H3K4me3$bed_id,
  "ZS1_19_MCF7_EV_|ZS1_19_MCF7_0025EV_|\\.stringent\\.bed")
wide_collapsed_H3K4me3 <- collapsed_H3K4me3 %>%
  pivot_wider(
    names_from = bed_id, values_from = col_info,
    names_prefix = "samp_")

######
#now I need to arrange these all identically
wide_collapsed_H3K4me3 = wide_collapsed_H3K4me3[c(1,2,3,5)]
names(wide_collapsed_H3K4me3) = c("H3K4me3", "H2AZ", "KDM5B", "H3K27me3")
H3K4me3_stats = wide_collapsed_H3K4me3

H3K4me3_stats$KDM5B_status = H3K4me3_stats$KDM5B
H3K4me3_stats$KDM5B_status[is.na(H3K4me3_stats$KDM5B)]=FALSE
H3K4me3_stats$KDM5B_status[!is.na(H3K4me3_stats$KDM5B)]=TRUE

H3K4me3_stats$H3K27me3_status = H3K4me3_stats$H3K27me3
H3K4me3_stats$H3K27me3_status[is.na(H3K4me3_stats$H3K27me3)] = FALSE
H3K4me3_stats$H3K27me3_status[!is.na(H3K4me3_stats$H3K27me3)] = TRUE

H3K4me3_stats$H2AZ_status = H3K4me3_stats$H2AZ
H3K4me3_stats$H2AZ_status[is.na(H3K4me3_stats$H2AZ)] = FALSE
H3K4me3_stats$H2AZ_status[!is.na(H3K4me3_stats$H2AZ)] = TRUE

#####
H3K4me3_stats = separate(
  H3K4me3_stats, col=KDM5B,
  into=paste0("KDM5B_", c("chr", "start", "end", "tot_score", "max_score",
                          "max_reg")), sep=",", remove=TRUE, convert=TRUE)
H3K4me3_stats = separate(
  H3K4me3_stats, col=H3K27me3,
  into=paste0("H3K27me3_", c("chr", "start", "end","tot_score", "max_score",
                             "max_reg")), sep=",", remove=TRUE, convert=TRUE)
H3K4me3_stats = separate(
  H3K4me3_stats, col=H2AZ,
  into=paste0("H2AZ_", c("chr", "start", "end","tot_score", "max_score",
                         "max_reg")), sep=",", remove=TRUE, convert=TRUE)

H3K4me3_stats = separate(
  H3K4me3_stats, col=H3K4me3,
  into=paste0("H3K4me3_", c("chr", "start", "end","tot_score", "max_score",
                         "max_reg")), sep=",", remove=TRUE, convert=TRUE)
######
# Now I have my table (which I need to write to file).
#Steps to repeat from KDM5B peaks stats include:
#Annotate. Can do that with Ben's pipeline instead this time, and just define
#regions based on that info.
#get peakdifs
#summary of peakdifs
#additionally summary of categories
#don't forget CpG islands
#also I need to graph heatmaps at peaks not just at genes
H3K4me3_stats$peak_id <- paste0("peak_", seq(1:dim(H3K4me3_stats)[1]))
save_state = H3K4me3_stats #now if I mess up I can load it from this

intersection_grange = makeGRangesFromDataFrame(
  H3K4me3_stats,
  start.field = "H3K4me3_start",
  end.field = "H3K4me3_end",
  seqnames.field = "H3K4me3_chr",
  keep.extra.columns = TRUE)

tx = TxDb.Hsapiens.UCSC.hg38.knownGene
hs = org.Hs.eg.db
source("seqTools.annotBedtoTxDbGene.R")

annot = annotBedtoTxDbGene(intersection_grange, tx, prefix = "hg38", org = hs)
intersection_anno = as.data.frame(annot)

write.table(
  H3K4me3_stats,
  file="H3K4me3_stats.txt",
  row.names=FALSE,
  quote=FALSE,
  sep="\t"
)
write.table(
  intersection_anno,
  file="H3K4me3_anno.txt",
  row.names=FALSE,
  quote=FALSE,
  sep="\t"
)

intersection_anno <- read.delim("~/intersection_anno.txt")
intersection_anno_clean <- intersection_anno %>%
  dplyr::select(-c(5,9,15,21, 31:36, 38, 43, 45, 47))
promoter_subset <- intersection_anno_clean %>%
  filter(abs(hg38.tsDist)<=2000) #that's most of the genes as well? How are they doing the TSS distance?
#okay we'll call this anno 1.
intersection_anno_1 = intersection_anno_clean
#We could try chipseeker or homer instead...
txdb = tx
peakAnno <- annotatePeak(intersection_grange,
                         tssRegion=c(-2000, 2000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
intersection_anno_2 = as.data.frame(peakAnno)
peakAnno2 <- annotatePeak(intersection_grange,
                         tssRegion=c(-2000, 2000),
                         TxDb=txdb, annoDb="org.Hs.eg.db", overlap="all")
intersection_anno_3 = as.data.frame(peakAnno2)

#we'll start by looking at anno2.
get_simplified_anno <- function(tabs){
  genic <- as.data.frame(tabs) %>%
    filter(grepl('Intron|Exon', annotation))
  not_introns  <- as.data.frame(tabs) %>%
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
  annoDF= rbind(not_intron_or_exon, not_first_exon,
                not_first_intron, first_exon, first_intron)
  return(annoDF)
}

genic <- as.data.frame(intersection_anno_2) %>%
  filter(grepl('Intron|Exon|intron|exon|UTR', annotation))
non_genic <- as.data.frame(intersection_anno_2) %>%
  filter(!grepl('Intron|Exon|intron|exon|UTR', annotation))
promoter <- as.data.frame(intersection_anno_2) %>%
  filter(grepl('Promoter', annotation))
#nope okay this is basically only annotating things as promoters, I'm guessing
#the intergenic stuff isn't registering bc it's so much less...

#just do H2AZ or KDM5B instead...

H2AZ_driven_intersection <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/H2AZ_driven_intersection.txt",
  header=FALSE, na.strings='.')
#####
H2AZ_driven_intersection$V7[is.na(H2AZ_driven_intersection$V7)] = "placeholder"
Nointersect_H2AZ = H2AZ_driven_intersection %>% filter(V7=="placeholder")
H2AZ_single_peak_per_peak = H2AZ_driven_intersection %>%
  filter(V7!="placeholder") %>%
  group_by(V6, V7) %>% mutate(max_col = max(V11)) %>%
  ungroup() %>% filter(V11==max_col) %>% dplyr::select(-c("max_col"))
H2AZ_single_peak_per_peak = rbind(H2AZ_single_peak_per_peak,
                                     Nointersect_H2AZ)
collapsed_H2AZ <- data.frame(
  do.call(paste, c(H2AZ_single_peak_per_peak[1:6], sep=",")),
  H2AZ_single_peak_per_peak$V7,
  do.call(paste, c(H2AZ_single_peak_per_peak[8:13], sep=",")))
names(collapsed_H2AZ) <- c("col_H2AZ", "bed_id", "col_info")
collapsed_H2AZ$bed_id <- str_remove_all(
  collapsed_H2AZ$bed_id,
  "ZS1_19_MCF7_EV_|ZS1_19_MCF7_0025EV_|\\.stringent\\.bed")
wide_collapsed_H2AZ <- collapsed_H2AZ %>%
  pivot_wider(
    names_from = bed_id, values_from = col_info,
    names_prefix = "samp_")

######
#now I need to arrange these all identically
wide_collapsed_H2AZ = wide_collapsed_H2AZ[c(1,2,3,4)]
names(wide_collapsed_H2AZ) = c("H2AZ", "H3K4me3", "KDM5B", "H3K27me3")
H2AZ_stats = wide_collapsed_H2AZ

H2AZ_stats$KDM5B_status = H2AZ_stats$KDM5B
H2AZ_stats$KDM5B_status[is.na(H2AZ_stats$KDM5B)]=FALSE
H2AZ_stats$KDM5B_status[!is.na(H2AZ_stats$KDM5B)]=TRUE

H2AZ_stats$H3K27me3_status = H2AZ_stats$H3K27me3
H2AZ_stats$H3K27me3_status[is.na(H2AZ_stats$H3K27me3)] = FALSE
H2AZ_stats$H3K27me3_status[!is.na(H2AZ_stats$H3K27me3)] = TRUE

H2AZ_stats$H3K4me3_status = H2AZ_stats$H3K4me3
H2AZ_stats$H3K4me3_status[is.na(H2AZ_stats$H3K4me3)] = FALSE
H2AZ_stats$H3K4me3_status[!is.na(H2AZ_stats$H3K4me3)] = TRUE

#####
H2AZ_stats = separate(
  H2AZ_stats, col=KDM5B,
  into=paste0("KDM5B_", c("chr", "start", "end", "tot_score", "max_score",
                          "max_reg")), sep=",", remove=TRUE, convert=TRUE)
H2AZ_stats = separate(
  H2AZ_stats, col=H3K27me3,
  into=paste0("H3K27me3_", c("chr", "start", "end","tot_score", "max_score",
                             "max_reg")), sep=",", remove=TRUE, convert=TRUE)
H2AZ_stats = separate(
  H2AZ_stats, col=H2AZ,
  into=paste0("H2AZ_", c("chr", "start", "end","tot_score", "max_score",
                         "max_reg")), sep=",", remove=TRUE, convert=TRUE)

H2AZ_stats = separate(
  H2AZ_stats, col=H3K4me3,
  into=paste0("H3K4me3_", c("chr", "start", "end","tot_score", "max_score",
                            "max_reg")), sep=",", remove=TRUE, convert=TRUE)
######
# Now I have my table (which I need to write to file).
#Steps to repeat from KDM5B peaks stats include:
#Annotate. Can do that with Ben's pipeline instead this time, and just define
#regions based on that info.
#get peakdifs
#summary of peakdifs
#additionally summary of categories
#don't forget CpG islands
#also I need to graph heatmaps at peaks not just at genes
H2AZ_stats$peak_id <- paste0("peak_", seq(1:dim(H2AZ_stats)[1]))
save_state = H2AZ_stats #now if I mess up I can load it from this

intersection_grange = makeGRangesFromDataFrame(
  H2AZ_stats,
  start.field = "H2AZ_start",
  end.field = "H2AZ_end",
  seqnames.field = "H2AZ_chr",
  keep.extra.columns = TRUE)

tx = TxDb.Hsapiens.UCSC.hg38.knownGene
hs = org.Hs.eg.db
source("seqTools.annotBedtoTxDbGene.R")

annot = annotBedtoTxDbGene(intersection_grange, tx, prefix = "hg38", org = hs)
intersection_anno = as.data.frame(annot)

intersection_anno_clean <- intersection_anno %>%
  dplyr::select(-c(5,9,15,21, 31:36, 38, 43, 45, 47))
promoter_subset <- intersection_anno_clean %>%
  filter(abs(hg38.tsDist)<=2000) #that's most of the genes as well? How are they doing the TSS distance?
#okay we'll call this anno 1.
intersection_anno_1 = intersection_anno_clean
#We could try chipseeker or homer instead...
txdb = tx
peakAnno <- annotatePeak(intersection_grange,
                         tssRegion=c(-2000, 2000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
intersection_anno_2 = as.data.frame(peakAnno)
peakAnno2 <- annotatePeak(intersection_grange,
                          tssRegion=c(-2000, 2000),
                          TxDb=txdb, annoDb="org.Hs.eg.db", overlap="all")
intersection_anno_3 = as.data.frame(peakAnno2)

gene_body <- as.data.frame(intersection_anno_2) %>%
  filter(grepl('Intron|Exon|intron|exon|UTR', annotation))
promoter <- as.data.frame(intersection_anno_2) %>%
  filter(grepl('Promoter', annotation))
intergenic <-  as.data.frame(intersection_anno_2) %>%
  filter(grepl('Downstream|Distal', annotation))
#5986 in gene body, 5990 intergenic, 12861 promoter
gene_body$simple_anno = rep("gene body", dim(gene_body)[1])
promoter$simple_anno = rep("promoter", dim(promoter)[1])
intergenic$simple_anno = rep("intergenic", dim(intergenic)[1])
H2AZ_simpanno = rbind(gene_body, promoter, intergenic)

mytheme = theme(
  panel.background=element_rect(fill="white"),
  text=element_text(color="black",face="bold",family="sans"),
  axis.text=element_text(color="black"),
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

###########
H2AZ_simpanno$H3K27me3_width = H2AZ_simpanno$H3K27me3_end - H2AZ_simpanno$H3K27me3_start
H2AZ_simpanno$H3K4me3_width = H2AZ_simpanno$H3K4me3_end - H2AZ_simpanno$H3K4me3_start
H2AZ_simpanno$KDM5B_width = H2AZ_simpanno$KDM5B_end - H2AZ_simpanno$KDM5B_start

H2AZ = str_split(H2AZ_simpanno$H2AZ_max_reg, "[:punct:]")
H2AZ_simpanno_H2AZ_reg = as.data.frame(do.call(rbind, H2AZ))
H2AZ_simpanno_H2AZ_reg$H2AZ_max_center = (
  as.numeric(H2AZ_simpanno_H2AZ_reg$V3) - as.numeric(
    H2AZ_simpanno_H2AZ_reg$V2))/2 + as.numeric(H2AZ_simpanno_H2AZ_reg$V2)

H3K4me3 = str_split(H2AZ_simpanno$H3K4me3_max_reg, "[:punct:]") #this lets me get the max region.
H2AZ_simpanno_H3K4me3_reg = as.data.frame(do.call(rbind, H3K4me3))
H2AZ_simpanno_H3K4me3_reg$H3K4me3_max_center = (
  as.numeric(H2AZ_simpanno_H3K4me3_reg$V3) - as.numeric(
    H2AZ_simpanno_H3K4me3_reg$V2))/2 + as.numeric(H2AZ_simpanno_H3K4me3_reg$V2)

KDM5B = str_split(H2AZ_simpanno$KDM5B_max_reg, "[:punct:]") #this lets me get the max region.
H2AZ_simpanno_KDM5B_reg = as.data.frame(do.call(rbind, KDM5B))
H2AZ_simpanno_KDM5B_reg$KDM5B_max_center = (
  as.numeric(H2AZ_simpanno_KDM5B_reg$V3) - as.numeric(
    H2AZ_simpanno_KDM5B_reg$V2))/2 + as.numeric(H2AZ_simpanno_KDM5B_reg$V2)

Peaks_w_peakdif <- cbind(
  H2AZ_simpanno,
  H2AZ_simpanno_H2AZ_reg,
  H2AZ_simpanno_H3K4me3_reg,
  H2AZ_simpanno_KDM5B_reg)

Peaks_w_peakdif$KDM5B_H2AZ_dif = Peaks_w_peakdif$KDM5B_max_center - Peaks_w_peakdif$H2AZ_max_center
Peaks_w_peakdif$KDM5B_H3K4me3_dif = Peaks_w_peakdif$KDM5B_max_center - Peaks_w_peakdif$H3K4me3_max_center
Peaks_w_peakdif$H2AZ_H3K4me3_dif = Peaks_w_peakdif$H2AZ_max_center - Peaks_w_peakdif$H3K4me3_max_center

names(Peaks_w_peakdif) = c(
  names(
    Peaks_w_peakdif)[1:46],
  "H2AZ_max_chrom", "H2AZ_max_start", "H2AZ_max_end",
  names(Peaks_w_peakdif)[50],
  "H3K4me3_max_chrom", "H3K4me3_max_start", "H3K4me3_max_end",
  names(Peaks_w_peakdif)[54],
  "KDM5B_max_chrom", "KDM5B_max_start", "KDM5B_max_end",
  names(Peaks_w_peakdif)[58:61])
Peaks_w_peakdif$strandedanno = paste0(Peaks_w_peakdif$simple_anno, Peaks_w_peakdif$geneStrand)

## ref https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/
#### plotting ####
library(rstatix)
Peaks_w_peakdif$strandedanno=as.factor(Peaks_w_peakdif$strandedanno)
Peaks_w_peakdif$simple_anno=as.factor(Peaks_w_peakdif$simple_anno)

filtered_for_analytics = Peaks_w_peakdif %>% filter(H3K4me3_status==TRUE)
filtered_for_analytics_KDM5B = Peaks_w_peakdif %>% filter(KDM5B_status==TRUE)
filtered_for_analytics_KDM5B_H3K4me3 = Peaks_w_peakdif %>% filter(H3K4me3_status==TRUE & KDM5B_status==TRUE)
krusk_test_1 <- filtered_for_analytics %>% kruskal_test(H2AZ_H3K4me3_dif~strandedanno)
pwc <- filtered_for_analytics %>% 
  dunn_test(H2AZ_H3K4me3_dif~strandedanno, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "strandedanno")
library(ggpubr)
H2_H3_1g = ggplot(data = filtered_for_analytics, aes(x=strandedanno, y=H2AZ_H3K4me3_dif)) +
  geom_boxplot() + 
  scale_y_continuous(limits=c(-2500,2500), n.breaks=11) +
  mytheme + stat_pvalue_manual(
    pwc, hide.ns = TRUE, y.position=c(2300, 2400, 2500)) +
  labs(
    subtitle = get_test_label(
      krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc)) +
  geom_hline(yintercept = 0, linetype="dotted")

krusk_test_1 <- filtered_for_analytics_KDM5B %>% kruskal_test(KDM5B_H2AZ_dif~strandedanno)
pwc <- filtered_for_analytics_KDM5B %>% 
  dunn_test(KDM5B_H2AZ_dif~strandedanno, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "strandedanno")
K_H2_1g = ggplot(data = filtered_for_analytics_KDM5B, aes(x=strandedanno, y=KDM5B_H2AZ_dif)) +
  geom_boxplot() + 
  scale_y_continuous(limits=c(-2500,2500), n.breaks=11) +
  mytheme + stat_pvalue_manual(
    pwc, hide.ns = TRUE, y.position=c(2500)) + labs(
      subtitle = get_test_label(
        krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc)) +
  geom_hline(yintercept = 0, linetype="dotted")

krusk_test_1 <- filtered_for_analytics_KDM5B_H3K4me3 %>% kruskal_test(KDM5B_H3K4me3_dif~strandedanno)
pwc <- filtered_for_analytics_KDM5B_H3K4me3 %>% 
  dunn_test(KDM5B_H3K4me3_dif~strandedanno, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "strandedanno")
K_H3_1g =ggplot(data = filtered_for_analytics_KDM5B_H3K4me3, aes(x=strandedanno, y=KDM5B_H3K4me3_dif)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-2500,2500), n.breaks=11) +
  mytheme + stat_pvalue_manual(
    pwc, hide.ns = TRUE, y.position=c(2500)) + 
  labs(subtitle = get_test_label(
    krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc)) +
  geom_hline(yintercept = 0, linetype="dotted")

#### abs version ####
krusk_test_1 <- filtered_for_analytics %>% kruskal_test(abs_HZ_H3~simple_anno)
filtered_for_analytics$abs_HZ_H3 <- abs(filtered_for_analytics$H2AZ_H3K4me3_dif)
pwc <- filtered_for_analytics %>% 
  dunn_test(abs_HZ_H3~simple_anno, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "simple_anno")
H2_H3_2g = ggplot(data = filtered_for_analytics, aes(x=simple_anno, y=abs_HZ_H3)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(0,2500), n.breaks=6) +
  mytheme + stat_pvalue_manual(
    pwc, hide.ns = TRUE, y.position=c(2400, 2500)) +
  labs(subtitle = get_test_label(
    krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc))

krusk_test_1 <- filtered_for_analytics_KDM5B %>% kruskal_test(abs_KD_HZ~simple_anno)
filtered_for_analytics_KDM5B$abs_KD_HZ <- abs(filtered_for_analytics_KDM5B$KDM5B_H2AZ_dif)
pwc <- filtered_for_analytics_KDM5B %>% 
  dunn_test(abs_KD_HZ~simple_anno, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "simple_anno")
K_H2_2g = ggplot(data = filtered_for_analytics_KDM5B, aes(x=simple_anno, y=abs_KD_HZ)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(0,2500), n.breaks=6) +
  mytheme + stat_pvalue_manual(
    pwc, hide.ns = TRUE) +
  labs(subtitle = get_test_label(
    krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc))

krusk_test_1 <- filtered_for_analytics_KDM5B_H3K4me3 %>% kruskal_test(abs_KD_H3~simple_anno)
filtered_for_analytics_KDM5B_H3K4me3$abs_KD_H3 <- abs(filtered_for_analytics_KDM5B_H3K4me3$KDM5B_H3K4me3_dif)
pwc <- filtered_for_analytics_KDM5B_H3K4me3 %>% 
  dunn_test(abs_KD_H3~simple_anno, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "simple_anno")
K_H3_2g = ggplot(data = filtered_for_analytics_KDM5B_H3K4me3, aes(x=simple_anno, y=abs_KD_H3)) +
  geom_boxplot() + scale_y_continuous(limits=c(0,2500), n.breaks=6) +
  mytheme + stat_pvalue_manual(pwc, hide.ns = TRUE) + labs(subtitle = get_test_label(
    krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc))

Peaks_w_peakdif$geneStrand[Peaks_w_peakdif$geneStrand==1]=
  "+"
Peaks_w_peakdif$geneStrand[Peaks_w_peakdif$geneStrand==2]=
  "-"
write.table(
  Peaks_w_peakdif,
  file="H2AZ_peaks_stats.txt",
  row.names=FALSE,
  quote=FALSE,
  sep="\t"
)

pdf(file="H2AZ_peaks_stats_graphs.pdf")
H2_H3_1g
H2_H3_2g
K_H2_1g
K_H2_2g
K_H3_1g
K_H3_2g
dev.off()

#test the next - is this dif for KDM5B true/false? Let's start with a facet wrap
H2_H3_KDM5Bg = ggplot(data = filtered_for_analytics, aes(x=simple_anno, y=abs_HZ_H3)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(0,2500), n.breaks=6) +
  mytheme +
  facet_wrap(~KDM5B_status)
H2_H3_KDM5Bg2 = ggplot(data = filtered_for_analytics, aes(x=strandedanno, y=H2AZ_H3K4me3_dif)) +
  geom_boxplot() + 
  scale_y_continuous(limits=c(-2500,2500), n.breaks=11) +
  mytheme + 
  geom_hline(yintercept = 0, linetype="dotted") +
  facet_wrap(~KDM5B_status)
KDM5B_fill = ggplot(data = Peaks_w_peakdif, aes(x=simple_anno, fill=KDM5B_status)) +
  geom_bar(position="fill") +
  mytheme
H3K4me3_fill = ggplot(data = Peaks_w_peakdif, aes(x=simple_anno, fill=H3K4me3_status)) +
  geom_bar(position="fill") +
  mytheme
H3K27me3_fill = ggplot(data = Peaks_w_peakdif, aes(x=simple_anno, fill=H3K27me3_status)) +
  geom_bar(position="fill") +
  mytheme

pdf(file="H2AZ_peaks_stats_bar_graphs.pdf")
KDM5B_fill
H3K4me3_fill
H3K27me3_fill
dev.off()


H2AZ_peaks_stats <- read.delim("~/H2AZ_peaks_stats.txt")
promoter_bed_for_graphing = H2AZ_peaks_stats %>%
  filter(simple_anno=="promoter") %>%
  select(c(1:3,30,6,36))
gene_body_bed_for_graphing = H2AZ_peaks_stats %>%
  filter(simple_anno=="gene body") %>%
  select(c(1:3,30,6,36))
intergenic_bed_for_graphing = H2AZ_peaks_stats %>%
  filter(simple_anno=="intergenic") %>%
  select(c(1:3,30,6,36))
intergenic_bed_for_graphing$geneStrand = "+"
gene_body_bed_for_graphing$geneStrand[gene_body_bed_for_graphing$geneStrand==1]=
  "+"
gene_body_bed_for_graphing$geneStrand[gene_body_bed_for_graphing$geneStrand==2]=
  "-"
promoter_bed_for_graphing$geneStrand[promoter_bed_for_graphing$geneStrand==1]=
  "+"
promoter_bed_for_graphing$geneStrand[promoter_bed_for_graphing$geneStrand==2]=
  "-"

write.table(
  gene_body_bed_for_graphing,
  file="gene_body_H2AZ.txt",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)
write.table(
  promoter_bed_for_graphing,
  file="promoter_H2AZ.txt",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)
write.table(
  intergenic_bed_for_graphing,
  file="intergenic_H2AZ.txt",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)

##For TSS_TES graphs, start with promoter only bc gene body likely to smear and disappear.
T1 = H2AZ_peaks_stats %>% filter(simple_anno=="promoter") %>% select(c(1,33,34, 30, 6, 36))
T1$geneStrand[T1$geneStrand==1]=
  "+"
T1$geneStrand[T1$geneStrand==2]=
  "-"
write.table(
  T1,
  file="limited_promoters.txt",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)

