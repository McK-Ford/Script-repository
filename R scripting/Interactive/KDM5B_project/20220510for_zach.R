#for zach
source("/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/NewMetaplotterJoshBased.R")
enableJIT(3)
##############
## plotting ##
bamDirectory1="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/bams/hg38_dedup/"
bamDirectory2="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/"
refseq_curated_longest_cgi_minus <- read.delim(paste0(bamDirectory2, "refseq_curated_longest_cgi_minus.hg38.txt"))
refseq_curated_longest_cgi_plus <- read.delim(paste0(bamDirectory2, "refseq_curated_longest_cgi_plus.hg38.txt"))

TSS_CpG_end_plus= refseq_curated_longest_cgi_plus %>%
  select(1,5,3,8,2,7) #using CpG start as score placeholder
TSS_CpG_end_plus$name=paste0(TSS_CpG_end_plus$name, "_", TSS_CpG_end_plus$cpg_e)
TSS_CpG_end_minus= refseq_curated_longest_cgi_minus %>% select(1,2,6,8,3,7) #cgi end placeholder
TSS_CpG_end_minus$name=paste0(TSS_CpG_end_minus$name, "_", TSS_CpG_end_minus$cpg_s)

name_vec = c("chrom", "start", "end", "name", "score", "strand")
colnames(TSS_CpG_end_minus) = name_vec
colnames(TSS_CpG_end_plus) = name_vec
TSS_CpG_end = rbind(TSS_CpG_end_minus, TSS_CpG_end_plus)

bam_list_paired = list(paste0(bamDirectory1, "ZS1_19_MCF7_EV_H3K4me3.dedup.bam"),
                       paste0(bamDirectory1, "ZS1_19_MCF7_EV_KDM5B.dedup.bam"),
                       paste0(bamDirectory1, "ZS1_19_MCF7_SH1_H3K4me3.dedup.bam"),
                       paste0(bamDirectory1, "ZS1_19_MCF7_SH2_H3K4me3.dedup.bam"),
                       paste0(bamDirectory1, "ZS1_19_MCF7_EV_H2AZ.dedup.bam")
)

mat_list = lapply(bam_list_paired, get_score_matrix, bed=TSS_CpG_end, b=-2000, a=2000,
                  method="single_stranded_anchored")
mat_list_reads_only = lapply(bam_list_paired, get_score_matrix, bed=TSS_CpG_end, b=-2000, a=2000,
                  method="single_stranded_anchored", ReadsOnly = TRUE)

metaplot_list <- lapply(mat_list, data.table::melt, measure.vars=c(2:401),
                        variable.name="Bins", value.name="Score")
merged_mat=cbind(metaplot_list[[1]], metaplot_list[[2]][[3]],
                 metaplot_list[[3]][[3]], metaplot_list[[4]][[3]], metaplot_list[[5]][[3]])
colnames(merged_mat)=c("gene_ID", "Bins", "H3K4me3_EV", "KDM5B_EV",
                       "H3K4me3_SH1", "H3K4me3_SH2", "H2AZ_EV"
)
merged_mat[,gene_ID:=NULL]
mean_merged_mat = merged_mat[,lapply(.SD, mean), by=Bins]
mean_merged_mat$Bins = as.numeric(mean_merged_mat$Bins)*10-2000
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

mat_melt = data.table::melt(mean_merged_mat, measure.vars=2:6,
                            variable.name="Antibody", value.name = "Score")

geom = ggplot(data=mat_melt, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H2AZ_EV"="goldenrod3", "H3K4me3_EV"="darkseagreen4",
                              "KDM5B_EV"="firebrick4",
                              "H3K4me3_SH1"="steelblue2", "H3K4me3_SH2"="steelblue4"
                              ))
geom

#just our H3K4me3
mat_melt2 = mat_melt %>% 
  filter(Antibody=="H3K4me3_EV"|Antibody=="H3K4me3_SH1"|Antibody=="H3K4me3_SH2")
geom1 = ggplot(data=mat_melt2, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H3K4me3_EV"="darkseagreen4",
                              "H3K4me3_SH1"="steelblue2", "H3K4me3_SH2"="steelblue4"))
geom1


mat_melt3 = mat_melt %>% filter(Antibody=="KDM5B_EV")
geom2 = ggplot(data=mat_melt3, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("KDM5B_EV"="firebrick4"))
geom2

mat_melt4 = mat_melt %>% filter(Antibody=="H2AZ_EV")
geom = ggplot(data=mat_melt4, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H2AZ_EV"="goldenrod3"))
geom


#################
metaplot_list2 <- lapply(mat_list_reads_only, data.table::melt, measure.vars=c(2:401),
                        variable.name="Bins", value.name="Score")
merged_mat2=cbind(metaplot_list2[[1]], metaplot_list2[[2]][[3]],
                 metaplot_list2[[3]][[3]], metaplot_list2[[4]][[3]], metaplot_list2[[5]][[3]])
colnames(merged_mat2)=c("gene_ID", "Bins", "H3K4me3_EV", "KDM5B_EV",
                       "H3K4me3_SH1", "H3K4me3_SH2", "H2AZ_EV"
)
merged_mat2[,gene_ID:=NULL]
mean_merged_mat2 = merged_mat2[,lapply(.SD, mean), by=Bins]
mean_merged_mat2$Bins = as.numeric(mean_merged_mat2$Bins)*10-2000
mat_melt5 = data.table::melt(mean_merged_mat2, measure.vars=2:6,
                            variable.name="Antibody", value.name = "Score")

geom = ggplot(data=mat_melt5, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H2AZ_EV"="goldenrod3", "H3K4me3_EV"="darkseagreen4",
                              "KDM5B_EV"="firebrick4",
                              "H3K4me3_SH1"="steelblue2", "H3K4me3_SH2"="steelblue4"
  ))
geom

#just our H3K4me3
mat_melt6 = mat_melt5 %>% 
  filter(Antibody=="H3K4me3_EV"|Antibody=="H3K4me3_SH1"|Antibody=="H3K4me3_SH2")
geom3 = ggplot(data=mat_melt6, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from TSS")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H3K4me3_EV"="darkseagreen4",
                              "H3K4me3_SH1"="steelblue2", "H3K4me3_SH2"="steelblue4"))
geom3


mat_melt7 = mat_melt5 %>% filter(Antibody=="KDM5B_EV")
geom4 = ggplot(data=mat_melt7, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from TSS")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("KDM5B_EV"="firebrick4"))
geom4

mat_melt8 = mat_melt5 %>% filter(Antibody=="H2AZ_EV")
geom = ggplot(data=mat_melt8, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H2AZ_EV"="goldenrod3"))
geom

pdf(file="zach_graphs.pdf")
geom3
geom4
dev.off()


########################################################
KDM5B_peaks_stats <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/analysis/pre_fix_bamcoverage_error/KDM5B 5K peaks/KDM5B_peaks_stats.txt"
  )
plus_stats = KDM5B_peaks_stats %>% filter(geneStrand==1)
minus_stats = KDM5B_peaks_stats %>% filter(geneStrand==2)
plus_stats_col_trim = plus_stats %>% mutate(p1=geneStart+1) %>%
  select(c(1,24,51,32,7))
plus_stats_col_trim$strand = "+"

minus_stats_col_trim = minus_stats %>% mutate(p1=geneEnd-1) %>%
  select(c(1,51,25,32,7))
minus_stats_col_trim$strand = "-"

name_vec = c("chrom",  "start", "end",    "name",   "score",  "strand")

colnames(plus_stats_col_trim) = name_vec
colnames(minus_stats_col_trim) = name_vec

KDM5B_TSSes = rbind(plus_stats_col_trim, minus_stats_col_trim)
write.table(KDM5B_TSSes, file="KDM5B_TSSes.bed", quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
