#unscaled
#Running NewMetaplotterJoshBased.R
source("/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/NewMetaplotterJoshBased.R")
enableJIT(3)
##############
## plotting ##
##############
bamDirectory1="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/bams/hg38_dedup/"
bamDirectory2="/Users/kayle/Box/Vertinolab/McKayla Ford/"
refseq_curated_longest_cgi_minus <- read.delim(paste0(
  bamDirectory2,
  "Data/refseq_curated_longest_cgi_minus.hg38.txt"))
refseq_curated_longest_cgi_plus <- read.delim(paste0(
  bamDirectory2,
  "Data/refseq_curated_longest_cgi_plus.hg38.txt"))

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
                       paste0(bamDirectory1, "ZS1_19_MCF7_EV_KDM5B.dedup.bam"))
mat_list = lapply(bam_list_paired, get_score_matrix, bed=TSS_CpG_end, b=-1000, a=1000,
                  method="single_stranded_anchored")
mat_list3 <- lapply(mat_list, data.table::melt, measure.vars=c(2:201),
                    variable.name="Bins", value.name="Score")
merged_mat=cbind(mat_list3[[1]], mat_list3[[2]][[3]], mat_list3[[3]][[3]])
colnames(merged_mat)=c("gene_ID", "Bins", "H2AZ", "H3K4me3_EV", "KDM5B_EV"
)
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

mat_melt = data.table::melt(mean_merged_mat, measure.vars=2:4,
                            variable.name="Antibody", value.name = "Score")

geom = ggplot(data=mat_melt, aes(x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c(
    "H2AZ"="goldenrod3", "H3K4me3_EV"="darkseagreen4","KDM5B_EV"="firebrick4"
    ))
geom

mat_melt_H2AZ = mat_melt %>% filter(Antibody=="H2AZ")
geom = ggplot(data=mat_melt_H2AZ, aes(x=Bins, y=Score)) +
  mytheme +
  geom_line(size=1.5, color="cyan4") +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed")
geom

svg(file="h2az1.svg")
geom
dev.off()

mat_melt_h3k4me3 = mat_melt %>% filter(Antibody=="H3K4me3_EV")
geom = ggplot(data=mat_melt_h3k4me3, aes(x=Bins, y=Score)) +
  mytheme +
  geom_line(size=1.5, color="firebrick2") +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed")
geom

svg(file="h3k4me31.svg")
geom
dev.off()

mat_melt_KDM5B = mat_melt %>% filter(Antibody=="KDM5B_EV")
geom = ggplot(data=mat_melt_KDM5B, aes(x=Bins, y=Score)) +
  mytheme +
  geom_line(size=1.5, color="darkorchid4") +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed")
geom

svg(file="KDM5B1.svg")
geom
dev.off()

#scale it
summary(mean_merged_mat$H2AZ)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.3157  0.4789  0.6250  0.6324  0.7479  1.0411 
mean_merged_mat_scaled = mean_merged_mat
mean_merged_mat_scaled$H2AZ = mean_merged_mat_scaled$H2AZ/1.0411
summary(mean_merged_mat$H3K4me3_EV)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.219   2.565   3.297   3.563   4.884   5.774 
mean_merged_mat_scaled$H3K4me3_EV = mean_merged_mat_scaled$H3K4me3_EV/5.774
summary(mean_merged_mat$KDM5B_EV)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1848  0.2770  0.3747  0.4556  0.5450  1.0731 
mean_merged_mat_scaled$KDM5B_EV = mean_merged_mat_scaled$KDM5B_EV/1.0731

mat_melt_scaled = data.table::melt(mean_merged_mat_scaled, measure.vars=2:4,
                            variable.name="Antibody", value.name = "Score")

mytheme = theme(
  panel.background=element_rect(fill="white"),
  text=element_text(color="black",face="bold",family="sans", size = 18),
  axis.text=element_text(color="black"),
  axis.ticks=element_line(color="black"),
  plot.margin=unit(c(0.25, 0.25, .25, 0.25),"cm"),
  plot.title=element_text(vjust=2),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  axis.title.y=element_text(vjust=1),
  axis.title.x=element_text(vjust=0),
  panel.border = element_blank(),
  axis.line=element_line()
)

geom = ggplot(data=mat_melt_scaled, aes(x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1.5) +
  xlab("Distance from TSS") +
  ylab("Avg. Normalizaed Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(labels=c("H2AZ", "H3K4me3", "KDM5B"), values=c(
   "cyan4", "firebrick2","darkorchid2"
  )) +
  scale_x_continuous(n.breaks=9)
geom

svg(filename = "overlap.svg")
geom
dev.off()

nascentDirectory="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/"
bam_list_proseq = list(paste0(nascentDirectory, "MCF7_H9_plus.bam"),
                       paste0(nascentDirectory, "MCF7_H9_minus.bam"),
                       paste0(nascentDirectory, "MCF7_C11_plus.bam"),
                       paste0(nascentDirectory, "MCF7_C11_minus.bam"),
                       paste0(nascentDirectory, "MCF7_G11_plus.bam"),
                       paste0(nascentDirectory, "MCF7_G11_minus.bam"),
                       paste0(nascentDirectory, "MCF7_B7_plus.bam"),
                       paste0(nascentDirectory, "MCF7_B7_minus.bam"))
plus_query = lapply(bam_list_proseq, get_score_matrix, bed=TSS_CpG_end_plus, b=-1000, a=1000,
                    method="single_stranded_anchored", pairedEnd=FALSE)
minus_query = lapply(bam_list_proseq, get_score_matrix, bed=TSS_CpG_end_minus, b=-1000, a=1000,
                    method="single_stranded_anchored", pairedEnd=FALSE)

plus_melt <- lapply(plus_query, data.table::melt, measure.vars=c(2:201),
                    variable.name="Bins", value.name="Score")
minus_melt <- lapply(minus_query, data.table::melt, measure.vars=c(2:201),
                    variable.name="Bins", value.name="Score")

plus_sense <- cbind(plus_melt[[1]], plus_melt[[3]][[3]], plus_melt[[5]][[3]], plus_melt[[7]][[3]])
minus_sense <- cbind(minus_melt[[2]], minus_melt[[4]][[3]], minus_melt[[6]][[3]], minus_melt[[8]][[3]])
sense = rbind(plus_sense, minus_sense)

plus_antisense <- cbind(plus_melt[[2]], plus_melt[[4]][[3]], plus_melt[[6]][[3]], plus_melt[[8]][[3]])
minus_antisense <- cbind(minus_melt[[1]], minus_melt[[3]][[3]], minus_melt[[5]][[3]], minus_melt[[7]][[3]])
antisense = rbind(plus_antisense, minus_antisense)

merged_mat_ts=cbind(sense, antisense[,3:6])
colnames(merged_mat_ts)=c("gene_ID", "Bins", "H9_sense", "C11_sense", "G11_sense",
                       "B7_sense","H9_antisense", "C11_antisense", "G11_antisense",
                       "B7_antisense"
)
merged_mat_ts[,gene_ID:=NULL]
mean_merged_mat_ts = merged_mat_ts[,lapply(.SD, mean), by=Bins]
mean_merged_mat_ts$Bins = as.numeric(mean_merged_mat_ts$Bins)*10-1000
melt2 = data.table::melt(mean_merged_mat_ts, measure.vars=2:5,
                                variable.name="Subclone1", value.name = "Sense_TS")
melt3 = data.table::melt(melt2, measure.vars=2:5,
                         variable.name="Subclone2", value.name = "Antisense_TS")
melt3 = melt3[,-c(2,4)]
mean_melt3 = melt3[,lapply(.SD, mean), by=Bins]

mat_melt_ts = data.table::melt(mean_melt3, measure.vars=2:3,
                            variable.name="Antibody", value.name = "Score")

geom = ggplot(data=mat_melt_ts, aes(x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=2) +
  xlab("Distance from TSS") +
  ylab("Avg. Normalizaed Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(n.breaks=9)
geom


geom = ggplot(data=mean_merged_mat_ts, aes(x=Bins)) +
  mytheme +
  geom_line(size=2, aes(y=H9_antisense)) +
  geom_line(size=2, aes(y=H9_sense)) +
  xlab("Distance from TSS") +
  ylab("Avg. Normalizaed Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(n.breaks=9)
geom

H9_only = mean_merged_mat_ts[,c(1,2,6)]

mat_melt_ts = data.table::melt(H9_only, measure.vars=2:3,
                               variable.name="Direction", value.name = "Score")

geom = ggplot(data=mat_melt_ts, aes(x=Bins, y=Score, color=Direction)) +
  mytheme +
  geom_line(size=2) +
  xlab("Distance from TSS") +
  ylab("Avg. Normalizaed Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(n.breaks=9) +
  scale_color_manual(labels=c("Sense", "Antisense"), values = c("blue", "red"))
geom

svg(file="transcription.svg")
geom
dev.off()

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
nointer_mat3 <- lapply(nointer_mat, data.table::melt, measure.vars=c(2:201),
                       variable.name="Bins", value.name="Score")
merged_nointer_mat=cbind(nointer_mat3[[1]], nointer_mat3[[2]][[3]], nointer_mat3[[3]][[3]]
                         )
colnames(merged_nointer_mat)=c("gene_ID", "Bins", "H2AZ", "H3K4me3", "KDM5B")

merged_nointer_mat[,gene_ID:=NULL]
mean_merged_nointer_mat = merged_nointer_mat[,lapply(.SD, mean), by=Bins]
mean_merged_nointer_mat$Bins = as.numeric(mean_merged_nointer_mat$Bins)*10-1000

mat_melt2 = data.table::melt(mean_merged_nointer_mat, measure.vars=2:4,
                             variable.name="Antibody", value.name = "Score")



geom = ggplot(data=mat_melt2, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H2AZ"="goldenrod3", "H3K4me3_EV"="darkseagreen4",
                              "KDM5B_EV"="firebrick4"))
geom

####################
geom = ggplot(data=mean_merged_mat, aes(x=Bins, y=H3K4me3_EV)) +
  mytheme +
  geom_line(size=2, color="firebrick2") +
  xlab("Distance from TSS") +
  ylab("Avg. Normalizaed Tags: H3K4me3") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(n.breaks=9) + 
  ylim(0,6)
geom
svg(file="H3K4me31.svg")
geom
dev.off()
geom = ggplot(data=mean_merged_nointer_mat, aes(x=Bins, y=H3K4me3)) +
  mytheme +
  geom_line(size=2, color="firebrick2") +
  xlab("Distance from TSS") +
  ylab("Avg. Normalizaed Tags: H3K4me3") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(n.breaks=9) + 
  ylim(0,6)
geom
svg(file="H3K4me32.svg")
geom
dev.off()
geom = ggplot(data=mean_merged_mat, aes(x=Bins, y=H2AZ)) +
  mytheme +
  geom_line(size=2, color="cyan4") +
  xlab("Distance from TSS") +
  ylab("Avg. Normalizaed Tags: H2AZ") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(n.breaks=9) + 
  ylim(0,1.25)
geom
svg(file="H2AZ1.svg")
geom
dev.off()
geom = ggplot(data=mean_merged_nointer_mat, aes(x=Bins, y=H2AZ)) +
  mytheme +
  geom_line(size=2, color="cyan4") +
  xlab("Distance from TSS") +
  ylab("Avg. Normalizaed Tags: H2AZ") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(n.breaks=9) + 
  ylim(0,1.25)
geom
svg(file="H2AZ2.svg")
geom
dev.off()
geom = ggplot(data=mean_merged_mat, aes(x=Bins, y=KDM5B_EV)) +
  mytheme +
  geom_line(size=2, color="darkorchid2") +
  xlab("Distance from TSS") +
  ylab("Avg. Normalizaed Tags: KDM5B") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(n.breaks=9) + 
  ylim(0,1.25)
geom
svg(file="KDM5B1.svg")
geom
dev.off()
geom = ggplot(data=mean_merged_nointer_mat, aes(x=Bins, y=KDM5B)) +
  mytheme +
  geom_line(size=2, color="darkorchid2") +
  xlab("Distance from TSS") +
  ylab("Avg. Normalizaed Tags: KDM5B") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(n.breaks=9) + 
  ylim(0,1.25)
geom
svg(file="KDM5B2.svg")
geom
dev.off()
save.image(file="return.RData")
load(file="return.RData")
##################
merged_mat=cbind(mat_list3[[1]], mat_list3[[2]][[3]], mat_list3[[3]][[3]])
colnames(merged_mat)=c("gene_ID", "Bins", "H2AZ", "H3K4me3_EV", "KDM5B_EV"
)
merged_mat$Bins = as.numeric(merged_mat$Bins)*10-1000
mm = merged_mat %>%
  group_by(gene_ID) %>%
  mutate(tot_H3K4me3 = sum(H2AZ)) %>%
  ungroup()

mytheme2 = theme(
  text=element_text(color="black",face="bold",family="sans", size = 18),
  axis.text=element_text(color="black"),
  plot.title=element_text(vjust=2),
  axis.title.x=element_text(vjust=0),
  axis.text.y = element_blank()
)

h1 <- ggplot(data=mm, mapping=aes(
  x=Bins, y=reorder(gene_ID, tot_H3K4me3),fill=H3K4me3_EV)) +
  geom_raster() +
  theme_void() +
  scale_fill_gradient(name="H3K4me3", low = "white",
                      high = "firebrick4", oob=scales::squish, limits=c(0,6)) +
  mytheme2 +
  xlab("Distance from TSS") +
  ylab("")

h2 <- ggplot(data=mm, mapping=aes(
  x=Bins, y=reorder(gene_ID, tot_H3K4me3),fill=H2AZ)) + geom_raster() +
  theme_void() +
  scale_fill_gradient(name="H2AZ", low = "white",
                      high = "darkcyan", oob=scales::squish, limits=c(0,1.2)) +
  mytheme2 +
  xlab("Distance from TSS") +
  ylab("")

h3 <- ggplot(data=mm, mapping=aes(
  x=Bins, y=reorder(gene_ID, tot_H3K4me3),fill=KDM5B_EV)) + geom_raster() +
  theme_void() +
  scale_fill_gradient(name="KDM5B", low = "white",
                      high = "darkorchid4", oob=scales::squish, limits=c(0,1.2)) +
  mytheme2 +
  xlab("Distance from TSS") +
  ylab("")

svg(file="test.svg")
h1
dev.off()

svg(file="test2.svg")
h2
dev.off()

svg(file="test3.svg")
h3
dev.off()

merged_nointer_mat=cbind(nointer_mat3[[1]], nointer_mat3[[2]][[3]], nointer_mat3[[3]][[3]]
)
colnames(merged_nointer_mat)=c("gene_ID", "Bins", "H2AZ", "H3K4me3", "KDM5B")
merged_nointer_mat$Bins = as.numeric(merged_nointer_mat$Bins)*10-1000
mm2 = merged_nointer_mat %>%
  group_by(gene_ID) %>%
  mutate(tot_H3K4me3 = sum(H2AZ)) %>%
  ungroup()

h4 <- ggplot(data=mm2, mapping=aes(
  x=Bins, y=reorder(gene_ID, tot_H3K4me3),fill=H3K4me3)) +
  geom_raster() +
  theme_void() +
  scale_fill_gradient(name="H3K4me3", low = "white",
                      high = "firebrick4", oob=scales::squish, limits=c(0,6)) +
  mytheme2 +
  xlab("Distance from TSS") +
  ylab("")

h5 <- ggplot(data=mm2, mapping=aes(
  x=Bins, y=reorder(gene_ID, tot_H3K4me3),fill=H2AZ)) + geom_raster() +
  theme_void() +
  scale_fill_gradient(name="H2AZ", low = "white",
                      high = "darkcyan", oob=scales::squish, limits=c(0,1.2)) +
  mytheme2 +
  xlab("Distance from TSS") +
  ylab("")

h6 <- ggplot(data=mm2, mapping=aes(
  x=Bins, y=reorder(gene_ID, tot_H3K4me3),fill=KDM5B)) + geom_raster() +
  theme_void() +
  scale_fill_gradient(name="KDM5B", low = "white",
                      high = "darkorchid4", oob=scales::squish, limits=c(0,1.2)) +
  mytheme2 +
  xlab("Distance from TSS") +
  ylab("")

svg(file="test4.svg")
h4
dev.off()

svg(file="test5.svg")
h5
dev.off()

svg(file="test6.svg")
h6
dev.off()

pdf(file="test.pdf")
h1
h2
h3
dev.off()



bam_list_KO = list(
  paste0(bamDirectory1, "ZS1_19_MCF7_SH1_H3K4me3.dedup.bam"),
  paste0(bamDirectory1, "ZS1_19_MCF7_SH2_H3K4me3.dedup.bam"),
  paste0(bamDirectory1, "ZS1_19_MCF7_SH1_H2AZ.dedup.bam"),
  paste0(bamDirectory1, "ZS1_19_MCF7_SH2_H2AZ.dedup.bam")
)


mat_list_KO = lapply(bam_list_KO, get_score_matrix, bed=TSS_CpG_end, b=-1000, a=1000,
                  method="single_stranded_anchored")
mat_list_KO2 <- lapply(mat_list_KO, data.table::melt, measure.vars=c(2:201),
                    variable.name="Bins", value.name="Score")
merged_mat_KO=cbind(mat_list_KO2[[1]], mat_list_KO2[[2]][[3]], mat_list_KO2[[3]][[3]], mat_list_KO2[[4]][[3]])
colnames(merged_mat_KO)=c("gene_ID", "Bins", "H3K4me3_SH1", "H3K4me3_SH2", "H2AZ_SH1", "H2AZ_SH2"
)
merged_mat_KO[,gene_ID:=NULL]
mean_merged_mat_KO = merged_mat_KO[,lapply(.SD, mean), by=Bins]
mean_merged_mat_KO$Bins = as.numeric(mean_merged_mat_KO$Bins)*10-1000

mean_merged_mat_KO = cbind(mean_merged_mat_KO, mean_merged_mat[,4])
mean_merged_mat_KO$H3K4me3_SH1 = mean_merged_mat_KO$H3K4me3_SH1/5.774
mean_merged_mat_KO$H3K4me3_SH2 = mean_merged_mat_KO$H3K4me3_SH2/5.774
mean_merged_mat_KO$H2AZ_SH1 = mean_merged_mat_KO$H2AZ_SH1/1.0411
mean_merged_mat_KO$H2AZ_SH2 = mean_merged_mat_KO$H2AZ_SH2/1.0411
mat_melt_KO = data.table::melt(mean_merged_mat_KO, measure.vars=2:6,
                            variable.name="Antibody", value.name = "Score")

geom = ggplot(data=mat_melt_KO, aes(x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1.5) +
  xlab("Distance from TSS") +
  ylab("Avg. Normalized Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c(
    "H2AZ_SH1"="cyan1", "H2AZ_SH2"="lightblue4",
    "H3K4me3_SH1"="indianred1", "H3K4me3_SH2"="tomato4", "KDM5B_EV"="darkorchid4"
  )) +
  scale_y_continuous(breaks = c(.25, .5, .75, 1))
geom

svg(file="overlap2.svg")
geom
dev.off()
