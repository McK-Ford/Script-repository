#refgene is refseq gene
setwd("~/")
source("C:/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/metaplotter.R")
enableJIT(3)
Hg38_refseq_genes <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/Regions/Hg38 genes CpG intersection/Hg38_refseq_genes.txt.gz", header=TRUE)
Hg38_refseq_genes <- Hg38_refseq_genes %>%
  filter(!grepl('chrY|([\\w_]+)alt|random|fix|v1', chrom)) %>%
  group_by(name2) %>%
  filter(txStart==min(txStart), txEnd==max(txEnd)) %>%
  distinct(txStart, txEnd, .keep_all = TRUE) %>%
  ungroup %>%
  select(-c(1,7:11, 14:16))
#massively cleaned up gene list
H2AZall_genescaled_plus = get_anchored_scores_scaled(bed=Hg38_refseq_genes, bw="ZS1_19_MCF7_EV_H2AZ.dedup.rpkm.bw", anch1="txStart", anch2="txEnd", strand="+", bin_num=100) #this is too large bc genes are big, isn't it? I'd probably be better off running this in deeptools and downloading the matrices.
H2AZall_genenotscaled_plus = get_anchored_scores(bed=Hg38_refseq_genes, bw="ZS1_19_MCF7_EV_H2AZ.dedup.rpkm.bw", anch="txStart", b=1000, a=3000, strand="+", bs=50) #okay I have this fixed.
merged_bins_sum <- H2AZall_genenotscaled_plus %>%
  drop_na() %>%
  group_by(Bins) %>%
  summarise(score = mean(score), .groups = "drop")
#write this as trimmed to use
Hg38_refseq_genes2 <- Hg38_refseq_genes[c(2,4,5,1,6,3,7)]
write_bed_file(Hg38_refseq_genes2, "Hg38_refseq_uniqued.txt")

H2AZ_matrix <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/H2AZ_IVs_noheader.txt", header=FALSE, comment.char="#")
H2AZ_matrix_reformat <- H2AZ_matrix %>%
  pivot_longer(cols=1:length(H2AZ_matrix),
               names_to = "Bins",
               names_prefix="V",
               values_to = "score",
               values_drop_na = TRUE) %>% #pivots table into longform with bins as one column and values as another
  mutate(Bins=(as.numeric(Bins)), #because extra columns interfere with the bin naming (bin 1 is imported as X3 for example)
         .keep="unused")
H2AZ_av = H2AZ_matrix_reformat %>%
  group_by(Bins) %>%
  summarise(avg_score = mean(score))

H3K4me3_matrix <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/H3K4me3_IVs_noheader.txt", header=FALSE, comment.char="#")
H3K4me3_matrix_reformat <- H3K4me3_matrix %>%
  pivot_longer(cols=1:length(H3K4me3_matrix),
               names_to = "Bins",
               names_prefix="V",
               values_to = "score",
               values_drop_na = TRUE) %>% #pivots table into longform with bins as one column and values as another
  mutate(Bins=(as.numeric(Bins)), #because extra columns interfere with the bin naming (bin 1 is imported as X3 for example)
         .keep="unused")
H3K4me3_av = H3K4me3_matrix_reformat %>%
  group_by(Bins) %>%
  summarise(avg_score = mean(score))

KDM5B_matrix <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/KDM5B_IVs_noheader.txt", header=FALSE, comment.char="#")
KDM5B_matrix_reformat <- KDM5B_matrix %>%
  pivot_longer(cols=1:length(KDM5B_matrix),
               names_to = "Bins",
               names_prefix="V",
               values_to = "score",
               values_drop_na = TRUE) %>% #pivots table into longform with bins as one column and values as another
  mutate(Bins=(as.numeric(Bins)), #because extra columns interfere with the bin naming (bin 1 is imported as X3 for example)
         .keep="unused")
max(KDM5B_av$avg_score) #6.42
max(H3K4me3_av$avg_score) #41
max(H2AZ_av$avg_score) #12.49
av_for_graphing = KDM5B_av
av_for_graphing$H3K4me3_avg_score = H3K4me3_av$avg_score/4.1
av_for_graphing$H2AZ_avg_score = H2AZ_av$avg_score/1.24
av_for_graphing$KDM5B_avg_score = av_for_graphing$avg_score/0.642

av_for_graphing2 <- av_for_graphing %>%
  pivot_longer(cols=3:5,
               names_to = "category",
               values_to = "score",
               values_drop_na = TRUE)

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

ggplot(data=av_for_graphing2) + geom_line(aes(x=Bins, y=score, color=category)) + mytheme
