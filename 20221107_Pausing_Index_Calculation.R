source("~/Script repository/My_Useful_Fns.R")

refseq <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/Regions/Untouched source data/Refseq_curated.hg38.txt.gz",
  header=TRUE
)
#remove columns bin (1), cds start - exon ends (7:11), more cds (14:16),
refseq2 <- refseq %>% select(-c(1, 7:11, 14:16)) #unless Ching-Hua used detect transcripts?
refseq3 <- refseq2[,c(2,4,5,7,6,3,1)]
colnames(refseq3) <- c("chrom", "start", "end", "symbol",
                       "score", "strand", "NM_ID")
refseq4p <- refseq3 %>%
  group_by(symbol) %>%
  filter(strand == "+") %>%
  mutate(start_site=min(start)) %>%
  filter(start==start_site) %>%
  select(-c(8)) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  ungroup
refseq4m <- refseq3 %>%
  group_by(symbol) %>%
  filter(strand == "-") %>%
  mutate(start_site=max(end)) %>%
  filter(end==start_site) %>%
  select(-c(8)) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  ungroup()

refseq4=rbind(refseq4p, refseq4m)
refseq4 <- refseq4 %>% filter(!grepl('chrY|([\\w_]+)alt|random|fix|v1', chrom))


T47 <- pI(bed=refseq4,
              bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/T47D_R1.BAM",
              pairedEnd=FALSE, pause_s=0, pause_e=150)
MCF7 <- pI(bed=refseq4,
               bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/MCF7_R1.BAM",
               pairedEnd = FALSE, pause_s=0, pause_e = 150)
MCF10A <- pI(bed=refseq4,
                 bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/MCF10A_R1.BAM",
                 pairedEnd = FALSE, pause_s=0, pause_e = 150)
SUM159 <- pI(bed=refseq4,
                 bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/SUM159_R1.BAM",
                 pairedEnd = FALSE, pause_s=0, pause_e = 150)
MB231 <- pI(bed=refseq4,
                bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/MB231_R1.BAM",
                pairedEnd = FALSE, pause_s=0, pause_e = 150)
MCF10A_pe <- pI(bed=refseq4,
                             bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/bamFiles/BAMdeDuped/MCF10A_deDuped.BMA",
                             pairedEnd = TRUE, pause_s=0, pause_e = 150) #much slower as expected
MCF10A_tgfb_se <- pI(bed=refseq4,
                    bam="../Box/Vertinolab/Ching-Hua Shih/project_treatmentTGFb/ProSeq/VER00065.PRO.2022.10.07.SE.NextSeq/proseq2.0/bamFiles/MCF10A_CTRL_A_R1_dedup_QC.sorted.BAM",
                    pairedEnd = FALSE, pause_s=0, pause_e=150)

#i am currently growing a list. I should fix that? smol so doesn't make big dif but still bad practice
pI_tab <- cbind(T47, MCF7[,10:14], SUM159[,10:14], MB231[,10:14], MCF10A[,10:14], MCF10A_pe[,10:14], MCF10A_tgfb_se[,10:14])
colnames(pI_tab) = c(colnames(pI_tab[,1:9]),
                     paste0("T47D",           "_", colnames(pI_tab[,10:14])),
                     paste0("MCF7",           "_", colnames(pI_tab[,10:14])),
                     paste0("SUM159",         "_", colnames(pI_tab[,10:14])),
                     paste0("MB231",          "_", colnames(pI_tab[,10:14])),
                     paste0("MCF10A",         "_", colnames(pI_tab[,10:14])),
                     paste0("MCF10A_pe",      "_", colnames(pI_tab[,10:14])),
                     paste0("MCF10A_tgfb_se", "_", colnames(pI_tab[,10:14])))
write.table(pI_tab, file="cell_lines_pI_3_no_filter.txt", quote=FALSE, sep="\t", row.names = FALSE)
#############################################################
pI_tab <- read.delim("~/cell_lines_pI_3_no_filter.txt", stringsAsFactors = FALSE)
pI_tab_has_sig = pI_tab %>%
  filter(T47D_body_counts != 0 | MCF7_body_counts != 0 | MCF10A_body_counts != 0 | SUM159_body_counts != 0 |
           MB231_body_counts != 0 | MCF10A_tgfb_se_body_counts != 0 | MCF10A_pe_body_counts != 0) 
#boxplot violinplot
#okay to make my life easier changing colnames, bc regex for these is a mess
namevec = c("pause_counts", "body_counts", "lengthnorm_pause", "lengthnorm_body", "pI")
colnames(pI_tab) = c(colnames(pI_tab[,1:9]),
                     paste0("T47D",           ".", namevec),
                     paste0("MCF7",           ".", namevec),
                     paste0("SUM159",         ".", namevec),
                     paste0("MB231",          ".", namevec),
                     paste0("MCF10A",         ".", namevec),
                     paste0("MCF10A_pe",      ".", namevec),
                     paste0("MCF10A_tgfb_se", ".", namevec))
pI_tab_longform = pI_tab %>%
  pivot_longer(cols=c(10:44), names_to = c("cell_line", ".value"), names_pattern = "(.+)\\.(.+)")

mytheme = theme(
  panel.background=element_blank(),
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
ggplot(data=pI_tab_longform) +
  geom_violin(aes(y=lengthnorm_body, x=cell_line)) +
  mytheme + 
  ylab("Body counts, length normalized") +
  xlab("") +
  scale_x_discrete(labels=c("T47D","MCF7","SUM159","MB231","MCF10A","MCF10A_pe","MCF10A_tgfb_se"))

ggplot(data=pI_tab_longform) +
  geom_point(aes(y=lengthnorm_body, x=pI, color=cell_line)) +
  mytheme + 
  ylab("Body counts, length normalized") +
  xlab("pI")
##
ggplot(data=pI_tab_longform) +
  geom_point(aes(y=lengthnorm_body, x=pI, color=cell_line)) +
  mytheme + 
  ylab("Body counts, length normalized") +
  xlab("pI") +
  ylim(0,0.5)
##
ggplot(data=pI_tab_longform) +
  geom_point(aes(y=lengthnorm_body, x=pI, color=cell_line)) +
  mytheme + 
  ylab("Body counts, length normalized") +
  xlab("pI") +
  ylim(0,0.01)
##
vlow_gb_counts <- pI_tab_longform %>% filter(pI>=1000 & body_counts>0)
ggplot(data=vlow_gb_counts) +
  geom_point(aes(y=lengthnorm_body, x=pI, color=lengthnorm_pause)) +
  mytheme + 
  ylab("Body counts, length normalized") +
  xlab("pI")
summary(vlow_gb_counts$pI)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1013    1192    1428    1613    1893    3910 
#obv this is already subsetted to the high pIs though
summary(vlow_gb_counts$lengthnorm_body)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000597 0.0001874 0.0003901 0.0009316 0.0010880 0.0056772
summary(vlow_gb_counts$lengthnorm_pause)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.06667 0.23667 0.58667 1.38204 1.96667 7.08667 
##
## What does it look like if I restrict it to needs more than 1 read every 1kb? 
vlow_gb_counts2 <- vlow_gb_counts %>% filter(lengthnorm_body>0.001) #leaves 20
ggplot(data=vlow_gb_counts2) +
  geom_point(aes(y=lengthnorm_body, x=pI, color=lengthnorm_pause)) +
  mytheme + 
  ylab("Body counts, length normalized") +
  xlab("pI")
summary(vlow_gb_counts2$lengthnorm_pause)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.267   2.253   3.250   3.602   4.618   7.087

pI_tab_longform$pI[pI_tab_longform$lengthnorm_body < 0.001] <- NA
#now there are going to be some genes we can remove entirely.
pI_tab_wideform = pI_tab_longform %>%
  pivot_wider(values_from=c(11:15), names_from = cell_line, names_sep=".", names_vary="slowest")
MCF10A_in_dif_seqs = pI_tab_wideform %>%
  select(c(1:9, 30:44))
MCF10A_in_dif_seqs_som_sig = MCF10A_in_dif_seqs %>%
  filter(!is.na(pI.MCF10A) | !is.na(pI.MCF10A_pe) | !is.na(pI.MCF10A_tgfb_se))
#talk about w/ paula
write.table(MCF10A_in_dif_seqs_som_sig, file="MCF10A_in_dif_seqs_pI.txt", quote=FALSE, sep="\t", row.names = FALSE)

project_cell_lines = pI_tab_wideform %>%
  select(c(1:34))
project_cell_lines_som_sig = project_cell_lines %>%
  filter(!is.na(pI.MCF10A) | !is.na(pI.MCF7) | !is.na(pI.MB231) | !is.na(pI.SUM159)| !is.na(pI.T47D))
#reurns 18k genes
#what about:
pI_tab_longform2 = project_cell_lines_som_sig %>%
  pivot_longer(cols=c(10:34), names_to = c(".value", "cell_line"), names_pattern = "(.+)\\.(.+)")


ggplot(data=pI_tab_longform2) +
  geom_violin(aes(y=pI, x=cell_line)) +
  mytheme + 
  ylab("pI") +
  xlab("")
write.table(project_cell_lines_som_sig, file="project_cell_lines_pI.txt", quote=FALSE, sep="\t", row.names = FALSE)
