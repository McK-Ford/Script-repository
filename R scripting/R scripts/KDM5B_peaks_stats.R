#### Load packages ####
library(tidyverse)
library(ChIPseeker)
library(BSgenome)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(ggpubr)
library(rstatix)
#### import data ####
point25percent_peaks_4clust <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/analysis/KDM5B 5K peaks/point25percent_peaks_4clust.bed",
  header=FALSE, comment.char="#"
  )
KDM5B_H2AZnH3K4me3_intersection <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/analysis/KDM5B 5K peaks/h3k4me3_h2az_offset/KDM5B_H2AZnH3K4me3_intersection.bed",
  header=FALSE
  )
KDM5B_noH2AZnoH3K4me3_intersection <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/analysis/KDM5B 5K peaks/h3k4me3_h2az_offset/KDM5B_noH2AZnoH3K4me3_intersection.bed",
  header=FALSE
)
KDM5B_H2AZnoH3K4me3_intersection <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/analysis/KDM5B 5K peaks/h3k4me3_h2az_offset/KDM5B_H2AZnoH3K4me3_intersection.bed",
  header=FALSE
)
KDM5B_H3K4me3noH2AZ_intersection <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/analysis/KDM5B 5K peaks/h3k4me3_h2az_offset/KDM5B_H3K4me3noH2AZ_intersection.bed",
  header=FALSE
)

#### label+join data ####
KDM5B_H2AZnH3K4me3_intersection$H2AZ_status <- rep(
  "True", dim(KDM5B_H2AZnH3K4me3_intersection)[1])
KDM5B_H2AZnoH3K4me3_intersection$H2AZ_status <- rep(
  "True", dim(KDM5B_H2AZnoH3K4me3_intersection)[1])
KDM5B_H3K4me3noH2AZ_intersection$H2AZ_status <- rep(
  "False", dim(KDM5B_H3K4me3noH2AZ_intersection)[1])
KDM5B_noH2AZnoH3K4me3_intersection$H2AZ_status <- rep(
  "False", dim(KDM5B_noH2AZnoH3K4me3_intersection)[1])

KDM5B_H2AZnH3K4me3_intersection$H3K4me3_status <- rep(
  "True", dim(KDM5B_H2AZnH3K4me3_intersection)[1])
KDM5B_H2AZnoH3K4me3_intersection$H3K4me3_status <- rep(
  "False", dim(KDM5B_H2AZnoH3K4me3_intersection)[1])
KDM5B_H3K4me3noH2AZ_intersection$H3K4me3_status <- rep(
  "True", dim(KDM5B_H3K4me3noH2AZ_intersection)[1])
KDM5B_noH2AZnoH3K4me3_intersection$H3K4me3_status <- rep(
  "False", dim(KDM5B_noH2AZnoH3K4me3_intersection)[1])

#### join the tables together ####
KDM5B_H2AZnoH3K4me3_intersection_trim1 <- KDM5B_H2AZnoH3K4me3_intersection %>%
  select(-c(5:10,14))
KDM5B_H2AZnH3K4me3_intersection_trim1 <- KDM5B_H2AZnH3K4me3_intersection %>%
  select(-c(5:10,14,20))
KDM5B_noH2AZnoH3K4me3_intersection_trim1 <- KDM5B_noH2AZnoH3K4me3_intersection %>%
  select(-c(5:10))
KDM5B_H3K4me3noH2AZ_intersection_trim1 <- KDM5B_H3K4me3noH2AZ_intersection %>%
  select(-c(5:10,14))

#bc some tables don't have H2AZ/H3K4me3 peaks, need to fill empty spaces w/ NAs.
NA_l1 <- rep(NA, dim(KDM5B_H2AZnoH3K4me3_intersection_trim1)[1])
NA_t1 <- cbind(NA_l1, NA_l1, NA_l1, NA_l1, NA_l1)
KDM5B_H2AZnoH3K4me3_intersection_trim2 <- cbind(
  KDM5B_H2AZnoH3K4me3_intersection_trim1[, c(1:12)],
  NA_t1, KDM5B_H2AZnoH3K4me3_intersection_trim1[, c(13:14)])

NA_l2 <- rep(NA, dim(KDM5B_H3K4me3noH2AZ_intersection_trim1)[1])
NA_t2 <- cbind(NA_l2, NA_l2, NA_l2, NA_l2, NA_l2)
KDM5B_H3K4me3noH2AZ_intersection_trim2 <- cbind(
  KDM5B_H3K4me3noH2AZ_intersection_trim1[, c(1:12)],
  NA_t2, KDM5B_H3K4me3noH2AZ_intersection_trim1[, c(13:14)])

NA_l3 <- rep(NA, dim(KDM5B_noH2AZnoH3K4me3_intersection_trim1)[1])
NA_t3 <- cbind(NA_l3, NA_l3, NA_l3, NA_l3,
               NA_l3, NA_l3, NA_l3, NA_l3, NA_l3, NA_l3)
KDM5B_noH2AZnoH3K4me3_intersection_trim2 <- cbind(
  KDM5B_noH2AZnoH3K4me3_intersection_trim1[, c(1:7)],
  NA_t3, KDM5B_noH2AZnoH3K4me3_intersection_trim1[, c(8:9)])

col_names <- c("Chrom","KDM5B_s","KDM5B_e", "KDM5B_maxreg",
               "KDM5B_totscore","KDM5B_maxscore",
               "clusters4", "H2AZ_s", "H2AZ_e",
               "H2AZ_totscore","H2AZ_maxscore", "H2AZ_maxreg",
               "H3K4me3_s","H3K4me3_e","H3K4me3_totscore",
               "H3K4me3_maxscore", "H3K4me3maxreg", "H2AZ_status",
               "H3K4me3_status")
names(KDM5B_H2AZnH3K4me3_intersection_trim1)<-col_names
names(KDM5B_H2AZnoH3K4me3_intersection_trim2)<-col_names
names(KDM5B_H3K4me3noH2AZ_intersection_trim2)<-col_names
names(KDM5B_noH2AZnoH3K4me3_intersection_trim2)<-col_names

KDM5B_peaks_tab <- rbind(
  KDM5B_noH2AZnoH3K4me3_intersection_trim2,
  KDM5B_H3K4me3noH2AZ_intersection_trim2,
  KDM5B_H2AZnoH3K4me3_intersection_trim2,
  KDM5B_H2AZnH3K4me3_intersection_trim1)

#### annotate ####
clusters4_lt <- point25percent_peaks_4clust$V13
names(clusters4_lt) <- point25percent_peaks_4clust$V2
KDM5B_peaks_tab$clusters4 <- unname(
  (clusters4_lt[as.character(KDM5B_peaks_tab$KDM5B_s)]))
#okay 1, before profiling lets see what this looks like. 
summary_tab = KDM5B_peaks_tab %>%
  group_by(clusters4,H2AZ_status,H3K4me3_status) %>% summarize(num=n())

totalgrange = makeGRangesFromDataFrame(
  KDM5B_peaks_tab,
  start.field = "KDM5B_s",
  end.field = "KDM5B_e",
  seqnames.field = "Chrom",
  keep.extra.columns = TRUE)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(totalgrange,
                         tssRegion=c(-2000, 2000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
annoDF = as.data.frame(peakAnno)

get_simplified_anno <- function(tabs){
  introns <- as.data.frame(tabs) %>%
    filter(grepl('Intron(.*)', annotation))
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

anno_df_w_simple_col = get_simplified_anno(annoDF)

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

anno_df_w_simple_col$KDM5B_width = anno_df_w_simple_col$end - anno_df_w_simple_col$start
anno_df_w_simple_col$H3K4me3_width = anno_df_w_simple_col$H3K4me3_e - anno_df_w_simple_col$H3K4me3_s
anno_df_w_simple_col$H2AZ_width = anno_df_w_simple_col$H2AZ_e - anno_df_w_simple_col$H2AZ_s

KDM5B_peaks_stats = anno_df_w_simple_col #so technically what happened here is 
#this is one of the spots where I wrote the table to file and read it back in, 
#but that makes the clean code confusing.
H2AZ = str_split(KDM5B_peaks_stats$H2AZ_maxreg, "[:punct:]")
KDM5B_peaks_H2AZ_reg = as.data.frame(do.call(rbind, H2AZ))
KDM5B_peaks_H2AZ_reg$H2AZ_max_center = (
  as.numeric(KDM5B_peaks_H2AZ_reg$V3) - as.numeric(
    KDM5B_peaks_H2AZ_reg$V2))/2 + as.numeric(KDM5B_peaks_H2AZ_reg$V2)

H3K4me3 = str_split(KDM5B_peaks_stats$H3K4me3maxreg, "[:punct:]") #this lets me get the max region.
KDM5B_peaks_H3K4me3_reg = as.data.frame(do.call(rbind, H3K4me3))
KDM5B_peaks_H3K4me3_reg$H3K4me3_max_center = (
  as.numeric(KDM5B_peaks_H3K4me3_reg$V3) - as.numeric(
    KDM5B_peaks_H3K4me3_reg$V2))/2 + as.numeric(KDM5B_peaks_H3K4me3_reg$V2)

KDM5B = str_split(KDM5B_peaks_stats$KDM5B_maxreg, "[:punct:]") #this lets me get the max region.
KDM5B_peaks_KDM5B_reg = as.data.frame(do.call(rbind, KDM5B))
KDM5B_peaks_KDM5B_reg$KDM5B_max_center = (
  as.numeric(KDM5B_peaks_KDM5B_reg$V3) - as.numeric(
    KDM5B_peaks_KDM5B_reg$V2))/2 + as.numeric(KDM5B_peaks_KDM5B_reg$V2)

Peaks_w_peakdif <- cbind(
  KDM5B_peaks_stats,
  KDM5B_peaks_H2AZ_reg,
  KDM5B_peaks_H3K4me3_reg,
  KDM5B_peaks_KDM5B_reg)

Peaks_w_peakdif$KDM5B_H2AZ_dif = Peaks_w_peakdif$KDM5B_max_center - Peaks_w_peakdif$H2AZ_max_center
Peaks_w_peakdif$KDM5B_H3K4me3_dif = Peaks_w_peakdif$KDM5B_max_center - Peaks_w_peakdif$H3K4me3_max_center
Peaks_w_peakdif$H2AZ_H3K4me3_dif = Peaks_w_peakdif$H2AZ_max_center - Peaks_w_peakdif$H3K4me3_max_center

names(Peaks_w_peakdif) = c(
  names(
    Peaks_w_peakdif)[1:37],
  "H2AZ_max_chrom", "H2AZ_max_start", "H2AZ_max_end",
  names(Peaks_w_peakdif)[41],
  "H3K4me3_max_chrom", "H3K4me3_max_start", "H3K4me3_max_end",
  names(Peaks_w_peakdif)[45],
  "KDM5B_max_chrom", "KDM5B_max_start", "KDM5B_max_end",
  names(Peaks_w_peakdif)[49:52])

Promoter_peaks <- Peaks_w_peakdif %>%
  filter(grepl('Promoter(.*)', simple_anno))
Non_promoter_peaks <- Peaks_w_peakdif %>%
  filter(!grepl('Promoter(.*)', simple_anno))

summary_check = Promoter_peaks %>% group_by(clusters4,geneStrand) %>% summarize(num=n())

Promoter_strand_plus <- Promoter_peaks %>%
  filter(geneStrand=="1")
Promoter_strand_minus <- Promoter_peaks %>%
  filter(geneStrand=="2")

summary(Promoter_peaks$H2AZ_H3K4me3_dif)
#Min.   1st Qu.   Median   Mean   3rd Qu.  Max.   NA's 
#-292440.0 -530.0   0.0  -228.3   470.0   88750.0   929 
summary(Promoter_peaks$KDM5B_H3K4me3_dif)
#Min.     1st Qu.   Median  Mean   3rd Qu.  Max.      NA's 
#-13090.00  -380.00  -25.00 -21.66  330.00  15800.00  929 
summary(Promoter_peaks$KDM5B_H2AZ_dif)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-150605    -390      20     259     405  295040     199  

summary(Non_promoter_peaks$H2AZ_H3K4me3_dif)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-242795    -230       0   -2138     300   32850    1328 
summary(Non_promoter_peaks$KDM5B_H3K4me3_dif)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-4805.0  -315.0   -40.0   209.4   350.0 11815.0    1328
summary(Non_promoter_peaks$KDM5B_H2AZ_dif)
#Min.    1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-159490    -775   90    9164    4408    280960     453 

summary(Promoter_strand_plus$H2AZ_H3K4me3_dif)
#Min.      1st Qu.  Median    Mean  3rd Qu.   Max.     NA's 
#-292440.00  -790.00  -210.00  -670.34  26.25 14820.00 457 
summary(Promoter_strand_plus$KDM5B_H3K4me3_dif)
#Min.      1st Qu.  Median  Mean 3rd Qu.   Max.      NA's 
#-150605.00 -500.00 -110.00 -36.04  301.25 250725.00  97 
summary(Promoter_strand_plus$KDM5B_H2AZ_dif)
#Min.    1st Qu.  Median Mean  3rd Qu.  Max.     NA's 
#-71845.0 -260.0  130.0  545.4   492.5 295040.0  102

summary(Promoter_strand_minus$H2AZ_H3K4me3_dif)
#Min.     1st Qu.  Median  Mean   3rd Qu.   Max.   NA's 
#-258060.0   -110.0  190.0 236.4   780.0   88750.0   472 
summary(Promoter_strand_minus$KDM5B_H3K4me3_dif)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -8435 -240     115     101     460   15800     472 
summary(Promoter_strand_minus$KDM5B_H2AZ_dif)
#Min.    1st Qu.  Median    Mean  3rd Qu.   Max.      NA's 
#-150605.00  -500.00 -110.00 -36.04 301.25 250725.00 97 

Non_promoter_peaks$peak_category = rep("non_promoter", dim(Non_promoter_peaks)[1])
Promoter_strand_minus$peak_category = rep("promoter_minus", dim(Promoter_strand_minus)[1])
Promoter_strand_plus$peak_category = rep("promoter_plus", dim(Promoter_strand_plus)[1])
merged = rbind(Promoter_strand_plus, Promoter_strand_minus, Non_promoter_peaks)

## ref https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/
#### plotting ####
KDM5B_peaks_stats = merged
KDM5B_peaks_stats$peak_category=as.factor(KDM5B_peaks_stats$peak_category)
krusk_test_1 = KDM5B_peaks_stats %>% kruskal_test(
  H2AZ_H3K4me3_dif~peak_category)
krusk_test_1_eff = KDM5B_peaks_stats %>% kruskal_effsize(
  H2AZ_H3K4me3_dif~peak_category)
pwc <- KDM5B_peaks_stats %>% 
  dunn_test(H2AZ_H3K4me3_dif~peak_category, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "peak_category")
H2_H3_1g = ggplot(data = KDM5B_peaks_stats, aes(x=peak_category, y=H2AZ_H3K4me3_dif)) +
  geom_boxplot() + 
  scale_y_continuous(limits=c(-2500,2500), n.breaks=11) +
  mytheme + stat_pvalue_manual(
    pwc, hide.ns = TRUE, y.position=c(2300, 2400, 2500)) +
  labs(
    subtitle = get_test_label(
      krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc))

krusk_test_1 = KDM5B_peaks_stats %>% kruskal_test(KDM5B_H2AZ_dif~peak_category)
krusk_test_1_eff = KDM5B_peaks_stats %>% kruskal_effsize(KDM5B_H2AZ_dif~peak_category)
pwc <- KDM5B_peaks_stats %>% 
  dunn_test(KDM5B_H2AZ_dif~peak_category, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "peak_category")
K_H2_1g = ggplot(data = KDM5B_peaks_stats, aes(x=peak_category, y=KDM5B_H2AZ_dif)) +
  geom_boxplot() + 
  scale_y_continuous(limits=c(-2500,2500), n.breaks=11) +
  mytheme + stat_pvalue_manual(
    pwc, hide.ns = TRUE, y.position=c(2400, 2500)) + labs(
      subtitle = get_test_label(
        krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc))

krusk_test_1 = KDM5B_peaks_stats %>% kruskal_test(KDM5B_H3K4me3_dif~peak_category)
krusk_test_1_eff = KDM5B_peaks_stats %>% kruskal_effsize(KDM5B_H3K4me3_dif~peak_category)
pwc <- KDM5B_peaks_stats %>% 
  dunn_test(KDM5B_H3K4me3_dif~peak_category, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "peak_category")
K_H3_1g =ggplot(data = KDM5B_peaks_stats, aes(x=peak_category, y=KDM5B_H3K4me3_dif)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-2500,2500), n.breaks=11) +
  mytheme + stat_pvalue_manual(
    pwc, hide.ns = TRUE, y.position=c(2400, 2500)) + 
  labs(subtitle = get_test_label(
    krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc))

#### abs version ####
KDM5B_peaks_stats$peak_category=as.factor(KDM5B_peaks_stats$peak_category)
KDM5B_peaks_stats$abs_HZ_H3 <- abs(KDM5B_peaks_stats$H2AZ_H3K4me3_dif)
krusk_test_1 = KDM5B_peaks_stats %>% kruskal_test(abs_HZ_H3~peak_category)
krusk_test_1_eff = KDM5B_peaks_stats %>% kruskal_effsize(abs_HZ_H3~peak_category)
pwc <- KDM5B_peaks_stats %>% 
  dunn_test(abs_HZ_H3~peak_category, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "peak_category")
H2_H3_2g = ggplot(data = KDM5B_peaks_stats, aes(x=peak_category, y=abs_HZ_H3)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(0,2500), n.breaks=6) +
  mytheme + stat_pvalue_manual(
    pwc, hide.ns = TRUE, y.position=c(2400, 2500)) +
  labs(subtitle = get_test_label(
    krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc))

KDM5B_peaks_stats$abs_KD_HZ <- abs(KDM5B_peaks_stats$KDM5B_H2AZ_dif)
pwc <- KDM5B_peaks_stats %>% 
  dunn_test(abs_KD_HZ~peak_category, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "peak_category")
K_H2_2g = ggplot(data = KDM5B_peaks_stats, aes(x=peak_category, y=abs_KD_HZ)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(0,2500), n.breaks=6) +
  mytheme + stat_pvalue_manual(
    pwc, hide.ns = TRUE, y.position=c(2400, 2500)) +
  labs(subtitle = get_test_label(
    krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc))

KDM5B_peaks_stats$abs_KD_H3 <- abs(KDM5B_peaks_stats$KDM5B_H3K4me3_dif)
pwc <- KDM5B_peaks_stats %>% 
  dunn_test(abs_KD_H3~peak_category, p.adjust.method = "bonferroni") 
pwc <- pwc %>% add_xy_position(x = "peak_category")
K_H3_2g = ggplot(data = KDM5B_peaks_stats, aes(x=peak_category, y=abs_KD_H3)) +
  geom_boxplot() + scale_y_continuous(limits=c(0,2500), n.breaks=6) +
  mytheme + stat_pvalue_manual(pwc, hide.ns = TRUE) + labs(subtitle = get_test_label(
    krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc))


#########
ggplot(data=KDM5B_peaks_stats) + geom_bar(aes(x=clusters4, fill=simple_anno)) + mytheme

mini_tab_4_graphing = KDM5B_peaks_stats[c(9,34)]
tmp1 = mini_tab_4_graphing
tmp1$clusters4 = rep("total", dim(tmp1)[1])
mini_tab_4_graphing = rbind(mini_tab_4_graphing, tmp1)
stat1g = ggplot(data=mini_tab_4_graphing) + geom_bar(aes(x=clusters4, fill=simple_anno), position="fill") + mytheme + scale_y_continuous(expand=c(0,0))

mini_tab_4_graphing2 = KDM5B_peaks_stats %>% filter(abs(distanceToTSS)<2000) %>% select(c(9,27))
tmp1 = mini_tab_4_graphing2
tmp1$clusters4 = rep("total", dim(tmp1)[1])
mini_tab_4_graphing2 = rbind(mini_tab_4_graphing2, tmp1)
stat2g= ggplot(data=mini_tab_4_graphing2) + geom_bar(aes(x=clusters4, fill=as.factor(geneStrand)), position="fill") + mytheme + scale_y_continuous(expand=c(0,0))

ggplot(data=KDM5B_peaks_stats) + geom_bar(aes(x=clusters2, fill=H3K4me3_status), position="dodge") + mytheme + facet_wrap(~H2AZ_status)
ggplot(data=KDM5B_peaks_stats) + geom_bar(aes(x=clusters4, fill=H3K4me3_status)) + mytheme + facet_wrap(~H2AZ_status)

###################################
KDM5B_peaks_stats$clusters4=as.factor(KDM5B_peaks_stats$clusters4)
krusk_test_1 = KDM5B_peaks_stats %>% kruskal_test(H2AZ_H3K4me3_dif~clusters4)
krusk_test_1_eff = KDM5B_peaks_stats %>% kruskal_effsize(H2AZ_H3K4me3_dif~clusters4)
pwc2 <- KDM5B_peaks_stats %>% 
  dunn_test(H2AZ_H3K4me3_dif~clusters4, p.adjust.method = "bonferroni")
pwc2 <- pwc2 %>% add_xy_position(x = "clusters4")
H2_H3_3g = ggplot(data = KDM5B_peaks_stats, aes(x=clusters4, y=H2AZ_H3K4me3_dif)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-2500,2500), n.breaks=11) +
  mytheme + stat_pvalue_manual(
    pwc2, hide.ns = TRUE, y.position=c(2100, 2200, 2300, 2400, 2500)) +
  labs(subtitle = get_test_label(krusk_test_1, detailed = TRUE),
       caption = get_pwc_label(pwc2))

KDM5B_peaks_stats$clusters4=as.factor(KDM5B_peaks_stats$clusters4)
krusk_test_1 = KDM5B_peaks_stats %>% kruskal_test(KDM5B_H2AZ_dif~clusters4)
krusk_test_1_eff = KDM5B_peaks_stats %>% kruskal_effsize(KDM5B_H2AZ_dif~clusters4)
pwc2 <- KDM5B_peaks_stats %>% 
  dunn_test(KDM5B_H2AZ_dif~clusters4, p.adjust.method = "bonferroni")
pwc2 <- pwc2 %>% add_xy_position(x = "clusters4")
K_H2_3g = ggplot(data = KDM5B_peaks_stats, aes(x=clusters4, y=KDM5B_H2AZ_dif)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-2500,2500), n.breaks=11) +
  mytheme + stat_pvalue_manual(pwc2, hide.ns = TRUE) +
  labs(subtitle = get_test_label(
    krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc2))

KDM5B_peaks_stats$clusters4=as.factor(KDM5B_peaks_stats$clusters4)
krusk_test_1 = KDM5B_peaks_stats %>% kruskal_test(KDM5B_H3K4me3_dif~clusters4)
krusk_test_1_eff = KDM5B_peaks_stats %>% kruskal_effsize(KDM5B_H3K4me3_dif~clusters4)
pwc2 <- KDM5B_peaks_stats %>% 
  dunn_test(KDM5B_H3K4me3_dif~clusters4, p.adjust.method = "bonferroni")
pwc2 <- pwc2 %>% add_xy_position(x = "clusters4")
K_H3_3g = ggplot(data = KDM5B_peaks_stats, aes(x=clusters4, y=KDM5B_H3K4me3_dif)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-2500,2500), n.breaks=11) +
  mytheme + stat_pvalue_manual(
    pwc2, hide.ns = TRUE, y.position=c(2100, 2200, 2300, 2400, 2500)) +
  labs(subtitle = get_test_label(
    krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc2))


#abs version!!!!
krusk_test_1 = KDM5B_peaks_stats %>% kruskal_test(abs(H2AZ_H3K4me3_dif)~clusters4)
krusk_test_1_eff = KDM5B_peaks_stats %>% kruskal_effsize(abs(H2AZ_H3K4me3_dif)~clusters4)
KDM5B_peaks_stats$abs_HZ_H3 <- abs(KDM5B_peaks_stats$H2AZ_H3K4me3_dif)
pwc2 <- KDM5B_peaks_stats %>% 
  dunn_test(abs_HZ_H3~clusters4, p.adjust.method = "bonferroni")
pwc2 <- pwc2 %>% add_xy_position(x = "clusters4")
H2_H3_4g = ggplot(data = KDM5B_peaks_stats, aes(x=clusters4, y=abs_HZ_H3)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(0,2500), n.breaks=6) +
  mytheme + stat_pvalue_manual(
    pwc2, hide.ns = TRUE, y.position=c(2300, 2400, 2500)) +
  labs(subtitle = get_test_label(
    krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc2))

KDM5B_peaks_stats$abs_KD_HZ <- abs(KDM5B_peaks_stats$KDM5B_H2AZ_dif)
krusk_test_1 = KDM5B_peaks_stats %>% kruskal_test(abs_KD_HZ~clusters4)
krusk_test_1_eff = KDM5B_peaks_stats %>% kruskal_effsize(abs_KD_HZ~clusters4)
pwc2 <- KDM5B_peaks_stats %>% 
  dunn_test(abs_KD_HZ~clusters4, p.adjust.method = "bonferroni")
pwc2 <- pwc2 %>% add_xy_position(x = "clusters4")
K_H2_4g = ggplot(data = KDM5B_peaks_stats, aes(x=clusters4, y=abs_KD_HZ)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(0,2500), n.breaks=6) +
  mytheme + stat_pvalue_manual(
    pwc2, hide.ns = TRUE, y.position=c(2300, 2400, 2500)) +
  labs(subtitle = get_test_label(
    krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc2))


pwc2 <- KDM5B_peaks_stats %>% 
  dunn_test(abs_KD_H3~clusters4, p.adjust.method = "bonferroni")
pwc2 <- pwc2 %>% add_xy_position(x = "clusters4")
K_H3_4g = ggplot(data = KDM5B_peaks_stats, aes(x=clusters4, y=abs_KD_H3)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(0,2500), n.breaks=6) +
  mytheme + stat_pvalue_manual(
    pwc2, hide.ns = TRUE, y.position=c(2200, 2300, 2400, 2500)) +
  labs(subtitle = get_test_label(
    krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc2))





###########################
#I think I'd just directly intersected KDM5B peak stats with cgis.
KDM5_cpg <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/analysis/KDM5B 5K peaks/h3k4me3_h2az_offset/KDM5_cpg.txt", header=FALSE)
names(KDM5_cpg) <- c(names(KDM5B_peaks_stats[1:50]), "cpg_chrom", "cpg_start", "cpg_end", "cpg_id")
KDM5_cpg$cpg_status <- rep("true", dim(KDM5_cpg)[1])
##
KDM5_nocpg <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/analysis/KDM5B 5K peaks/h3k4me3_h2az_offset/KDM5_nocpg.txt", header=FALSE)
NA_l <- rep(NA, dim(KDM5_nocpg)[1])
NA_t <- cbind(NA_l, NA_l, NA_l, NA_l)
KDM5_nocpg <- cbind(KDM5_nocpg, NA_t)
names(KDM5_nocpg) <- c(names(KDM5B_peaks_stats[1:50]), "cpg_chrom", "cpg_start", "cpg_end", "cpg_id")
KDM5_nocpg$cpg_status <- rep("false", dim(KDM5_nocpg)[1])
KDM5B_peaks_stats <- rbind(KDM5_cpg, KDM5_nocpg)
#######
KDM5B_peaks_stats$cpg_status=as.factor(KDM5B_peaks_stats$cpg_status)
krusk_test_1 = KDM5B_peaks_stats %>% kruskal_test(H2AZ_H3K4me3_dif~cpg_status)
krusk_test_1_eff = KDM5B_peaks_stats %>% kruskal_effsize(H2AZ_H3K4me3_dif~cpg_status)
pwc2 <- KDM5B_peaks_stats %>% 
  wilcox_test(H2AZ_H3K4me3_dif~cpg_status, p.adjust.method = "bonferroni")
pwc2 <- pwc2 %>% add_xy_position(x = "cpg_status")
ggplot(data = KDM5B_peaks_stats, aes(x=cpg_status, y=H2AZ_H3K4me3_dif)) + geom_boxplot() + ylim(-2000,2100) + mytheme + stat_pvalue_manual(pwc2, hide.ns = FALSE, y.position=c(2000)) + labs(subtitle = get_test_label(krusk_test_1, detailed = TRUE), caption = get_pwc_label(pwc2))
###
ggplot(data=KDM5B_peaks_stats) + geom_bar(aes(x=cpg_status, fill=simple_anno), position="dodge") + mytheme
ggplot(data=KDM5B_peaks_stats) + geom_bar(aes(x=clusters4, fill=cpg_status), position="dodge") + mytheme

mini_tab_4_graphing3 = KDM5B_peaks_stats %>% select(c(8,55))
tmp1 = mini_tab_4_graphing3
tmp1$clusters4 = rep("total", dim(tmp1)[1])
mini_tab_4_graphing3 = rbind(mini_tab_4_graphing3, tmp1)
ggplot(data=mini_tab_4_graphing3) + geom_bar(aes(x=clusters4, fill=cpg_status), position="fill") + mytheme + scale_y_continuous(expand=c(0,0))
ggplot(data=KDM5B_peaks_stats) + geom_bar(aes(x=cpg_status, fill=H3K4me3_status), position="dodge") + mytheme + facet_wrap(~H2AZ_status)


H2AZ_bed <- KDM5B_peaks_stats %>% select(c(1,39,40)) %>% filter(!is.na(H2AZ_max_start))
H3K4me3_bed <- KDM5B_peaks_stats %>% select(c(1,43,44)) %>% filter(!is.na(H3K4me3_max_start))
KDM5B_bed <- KDM5B_peaks_stats %>% select(c(1,47,48)) %>% filter(!is.na(KDM5B_max_start))
write.table(
  H2AZ_bed,
  file="H2AZ_maxreg.bed",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)
write.table(
  H3K4me3_bed,
  file="H3K4me3_maxreg.bed",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)
write.table(
  KDM5B_bed,
  file="KDM5B_maxreg.bed",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)
write.table(
  KDM5B_peaks_stats,
  file="KDM5B_peaks_stats.txt",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)

pdf(file="KDM5B_peaks_stats_graphs.pdf")
H2_H3_1g
H2_H3_2g
H2_H3_3g
H2_H3_4g
K_H2_1g
K_H2_2g
K_H2_3g
K_H2_4g
K_H3_1g
K_H3_2g
K_H3_3g
K_H3_4g
stat1g
stat2g
stat3g
dev.off()

