#unscaled
#Running NewMetaplotterJoshBased.R
source("/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/NewMetaplotterJoshBased.R")
#does it load the packages?
##############
## plotting ##
##############
bamDirectory1="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/bams/hg38_dedup/"
bamDirectory2="/Users/kayle/Box/Vertinolab/McKayla Ford/"
refseq_curated_longest_cgi_minus <- read.delim("~/refseq_curated_longest_cgi_minus.hg38.txt")
refseq_curated_longest_cgi_plus <- read.delim("~/refseq_curated_longest_cgi_plus.hg38.txt")

TSS_CpG_end_plus= refseq_curated_longest_cgi_plus %>%
  select(1,5,3,8,2,7) #using CpG start as score placeholder
TSS_CpG_end_plus$name=paste0(TSS_CpG_end_plus$name, "_", TSS_CpG_end_plus$cpg_e)
TSS_CpG_end_minus= refseq_curated_longest_cgi_minus %>% select(1,2,6,8,3,7) #cgi end placeholder
TSS_CpG_end_minus$name=paste0(TSS_CpG_end_minus$name, "_", TSS_CpG_end_minus$cpg_s)

name_vec = c("chrom", "start", "end", "name", "score", "strand")
colnames(TSS_CpG_end_minus) = name_vec
colnames(TSS_CpG_end_plus) = name_vec
TSS_CpG_end = rbind(TSS_CpG_end_minus, TSS_CpG_end_plus)

bam_list_paired = list(paste0(bamDirectory1, "ZS1_19_MCF7_EV_H2AZ.dedup.bam"),
                       paste0(bamDirectory1, "ZS1_19_MCF7_EV_H3K4me3.dedup.bam"),
                       paste0(bamDirectory1, "ZS1_19_MCF7_EV_KDM5B.dedup.bam"),
                       paste0(bamDirectory1, "ZS1_19_MCF7_SH1_H3K4me3.dedup.bam"),
                       paste0(bamDirectory1, "ZS1_19_MCF7_SH2_H3K4me3.dedup.bam")
)
bam_list_notpaired = list(
  "/Users/kayle/Box/Vertinolab/McKayla Ford/MCF7.Polyak.G15.H3K4me3.WT.ChIPseq.dedup.bam",
      "/Users/kayle/Box/Vertinolab/McKayla Ford/MCF7.Polyak.G5.KDM5B.WT.ChIPseq.dedup.bam",
      "/Users/kayle/Box/Vertinolab/McKayla Ford/MCF7.Polyak.G15.H3K4me3.KDM5B_KO_N1.ChIPseq.dedup.bam",
      "/Users/kayle/Box/Vertinolab/McKayla Ford/MCF7.Polyak.G15.H3K4me3.KDM5B_KO_N2.ChIPseq.dedup.bam"
)
#wait was the polyak data paired end? I guess I'll see...
mat_list = lapply(bam_list_paired, get_score_matrix, bed=TSS_CpG_end, b=-1000, a=1000,
                  method="single_stranded_anchored")
mat_list2 = lapply(bam_list_notpaired, get_score_matrix, bed=TSS_CpG_end, b=-1000, a=1000,
                  method="single_stranded_anchored", pairedEnd=FALSE)
mat_list3 <- lapply(mat_list, data.table::melt, measure.vars=c(2:201),
                        variable.name="Bins", value.name="Score")
mat_list4 <- lapply(mat_list2, data.table::melt, measure.vars=c(2:201),
                        variable.name="Bins", value.name="Score")
merged_mat=cbind(mat_list3[[1]], mat_list3[[2]][[3]], mat_list3[[3]][[3]],
                 mat_list3[[4]][[3]], mat_list3[[5]][[3]], mat_list4[[1]][[3]],
                 mat_list4[[2]][[3]], mat_list4[[3]][[3]], mat_list4[[4]][[3]])
colnames(merged_mat)=c("gene_ID", "Bins", "H2AZ", "H3K4me3_EV", "KDM5B_EV",
                       "H3K4me3_SH1", "H3K4me3_SH2", "H3K4me3_pWT", "KDM5B_pWT",
                       "H3K4me3_pN1", "H3K4me3_pN2"
                       )
write.table(merged_mat, "merged_mat.txt")
merged_mat[,gene_ID:=NULL]
mean_merged_mat = merged_mat[,lapply(.SD, mean), by=Bins]
mean_merged_mat$Bins = as.numeric(mean_merged_mat$Bins)*10-1000
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

mat_melt = data.table::melt(mean_merged_mat, measure.vars=2:10,
                            variable.name="Antibody", value.name = "Score")

geom = ggplot(data=mat_melt, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H2AZ"="goldenrod3", "H3K4me3_EV"="darkseagreen4",
                        "KDM5B_EV"="firebrick4",
                     "H3K4me3_SH1"="steelblue2", "H3K4me3_SH2"="steelblue4",
                     "H3K4me3_pWT"="thistle4", "KDM5B_pWT"="tomato1",
                     "H3K4me3_pN1"="maroon2", "H3K4me3_pN2"="salmon3"))
geom

#just our H3K4me3
mat_melt3 = mat_melt %>% filter(Antibody=="H3K4me3_EV"|Antibody=="H3K4me3_SH1"|Antibody=="H3K4me3_SH2")
geom = ggplot(data=mat_melt3, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H3K4me3_EV"="darkseagreen4",
                              "H3K4me3_SH1"="steelblue2", "H3K4me3_SH2"="steelblue4"))
geom

mat_melt4 = mat_melt %>% filter(Antibody=="H3K4me3_pWT"|Antibody=="H3K4me3_pN1"|Antibody=="H3K4me3_pN2")
geom = ggplot(data=mat_melt4, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H3K4me3_pWT"="thistle4", "H3K4me3_pN1"="maroon2", "H3K4me3_pN2"="salmon3"))
geom

mat_melt5 = mat_melt %>% filter(Antibody=="KDM5B_pWT"|Antibody=="KDM5B_EV"|Antibody=="H2AZ")
geom = ggplot(data=mat_melt5, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H2AZ"="goldenrod3", "KDM5B_EV"="firebrick4",
                              "KDM5B_pWT"="tomato1"))
geom

#adding in the proseq data would be good
#possibly seeing if things look dif based on CpG-3end proseq counts - can classify things based on levels of transcription!
#don't forget the non-CGI genes would be interesting to look at! Need to decide how I'll determine them though.
# did a loj

test <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/test.txt", header=FALSE)
test_nointer = test %>% filter(V8==-1)
#that's 18375. Don't forget to deduplicate these...
test_inter_p = test %>% filter(V8!=-1 & V5=="+")
#that's 22.5k
test_inter_m = test %>% filter(V8!=-1 & V5=="-")
#that's 21.5k
test_inter_p = test_inter_p %>%
  group_by(V6) %>%
  mutate(testing_val = min(V8)) %>%
  filter(V8==testing_val) %>%
  mutate(testing_val2 = min(V2)) %>%
  filter(V2==testing_val2) %>%
  ungroup() %>%
  distinct(V6, .keep_all = TRUE)
#now its 8943
tru_nocgi_p = test_inter_p %>% filter(V8>V2) %>% select(c(1:3, 6, 4:5))
#okay 2022 true ones then, that happen to have gene body cgis
test_inter_m = test_inter_m %>%
  group_by(V6) %>%
  mutate(testing_val = max(V9)) %>%
  filter(V9==testing_val) %>%
  mutate(testing_val2 = max(V3)) %>%
  filter(V3==testing_val2) %>%
  ungroup() %>%
  distinct(V6, .keep_all = TRUE)
#now its 8624
tru_nocgi_m = test_inter_m %>% filter(V3>V9) %>% select(c(1:3, 6, 4:5))
#1893
#no inter cleanup
no_inter_p = test_nointer %>%
  filter(V5=="+") %>%
  group_by(V6) %>%
  mutate(testing_val = min(V2)) %>%
  filter(V2==testing_val) %>%
  ungroup() %>%
  distinct(V6, .keep_all = TRUE) %>% select(c(1:3, 6, 4:5))
#6566 more
no_inter_m = test_nointer %>%
  filter(V5=="-") %>%
  group_by(V6) %>%
  mutate(testing_val = max(V3)) %>%
  filter(V3==testing_val) %>%
  ungroup() %>%
  distinct(V6, .keep_all = TRUE) %>% select(c(1:3, 6, 4:5))
#6377 more
nointer_tab = rbind(tru_nocgi_p, tru_nocgi_m, no_inter_p, no_inter_m) #16858
nointer_tab <-nointer_tab %>% filter(!grepl("chrY|([\\w_]+)alt|random|fix|v1", V1)) %>% distinct(V6, .keep_all = TRUE) #15413
#########################
nointer_mat = lapply(bam_list_paired, get_score_matrix, bed=nointer_tab, b=-1000, a=1000,
                  method="single_stranded_anchored")
nointer_mat2 = lapply(bam_list_notpaired, get_score_matrix, bed=nointer_tab, b=-1000, a=1000,
                   method="single_stranded_anchored", pairedEnd=FALSE)
nointer_mat3 <- lapply(nointer_mat, data.table::melt, measure.vars=c(2:201),
                    variable.name="Bins", value.name="Score")
nointer_mat4 <- lapply(nointer_mat2, data.table::melt, measure.vars=c(2:201),
                    variable.name="Bins", value.name="Score")
merged_nointer_mat=cbind(nointer_mat3[[1]], nointer_mat3[[2]][[3]], nointer_mat3[[3]][[3]],
                 nointer_mat3[[4]][[3]], nointer_mat3[[5]][[3]], nointer_mat4[[1]][[3]],
                 nointer_mat4[[2]][[3]], nointer_mat4[[3]][[3]], nointer_mat4[[4]][[3]])
colnames(merged_nointer_mat)=c("gene_ID", "Bins", "H2AZ", "H3K4me3_EV", "KDM5B_EV",
                       "H3K4me3_SH1", "H3K4me3_SH2", "H3K4me3_pWT", "KDM5B_pWT",
                       "H3K4me3_pN1", "H3K4me3_pN2"
)

merged_nointer_mat[,gene_ID:=NULL]
mean_merged_nointer_mat = merged_nointer_mat[,lapply(.SD, mean), by=Bins]
mean_merged_nointer_mat$Bins = as.numeric(mean_merged_nointer_mat$Bins)*10-1000

mat_melt2 = data.table::melt(mean_merged_nointer_mat, measure.vars=2:10,
                            variable.name="Antibody", value.name = "Score")

geom = ggplot(data=mat_melt2, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H2AZ"="goldenrod3", "H3K4me3_EV"="darkseagreen4",
                              "KDM5B_EV"="firebrick4",
                              "H3K4me3_SH1"="steelblue2", "H3K4me3_SH2"="steelblue4",
                              "H3K4me3_pWT"="thistle4", "KDM5B_pWT"="tomato1",
                              "H3K4me3_pN1"="maroon2", "H3K4me3_pN2"="salmon3"))
geom

#just our H3K4me3
mat_melt3 = mat_melt2 %>% filter(Antibody=="H3K4me3_EV"|Antibody=="H3K4me3_SH1"|Antibody=="H3K4me3_SH2")
geom = ggplot(data=mat_melt3, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H3K4me3_EV"="darkseagreen4",
                              "H3K4me3_SH1"="steelblue2", "H3K4me3_SH2"="steelblue4"))
geom

mat_melt4 = mat_melt2 %>% filter(Antibody=="H3K4me3_pWT"|Antibody=="H3K4me3_pN1"|Antibody=="H3K4me3_pN2")
geom = ggplot(data=mat_melt4, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H3K4me3_pWT"="thistle4", "H3K4me3_pN1"="maroon2", "H3K4me3_pN2"="salmon3"))
geom

mat_melt5 = mat_melt2 %>% filter(Antibody=="KDM5B_pWT"|Antibody=="KDM5B_EV"|Antibody=="H2AZ")
geom = ggplot(data=mat_melt5, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H2AZ"="goldenrod3", "KDM5B_EV"="firebrick4",
                              "KDM5B_pWT"="tomato1"))
geom

write.table(merged_nointer_mat, "merged_nointer_mat.txt")
