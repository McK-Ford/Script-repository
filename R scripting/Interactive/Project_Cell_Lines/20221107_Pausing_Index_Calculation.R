source("~/Script repository/My_Useful_Fns.R")

gc <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/Regions/Untouched source data/Gencode39_knowngene.hg38.txt.gz",
  header=TRUE
)
#remove irrelevant columns
gc2 <- gc %>% dplyr::select(-c(7:17, 20, 22, 24:26)) #unless Ching-Hua used detect transcripts?
gc3 <- gc2[,c(1:3, 7, 5:6, 4, 8:10)]
colnames(gc3) <- c("chrom", "start", "end", "symbol",
                       "score", "strand", "ENST_ID", "ID3", "class", "type")
gc_genes_only <- gc3 %>%
  filter(type=="protein_coding" & !grepl('chrY|([\\w_]+)alt|random|fix|v1', chrom)) ###86K out of 266K
gc4p <- gc_genes_only %>%
  group_by(symbol) %>%
  filter(strand == "+") %>%
  mutate(start_site=min(start)) %>%
  filter(start==start_site) %>%
  dplyr::select(-c(11)) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  ungroup()
gc4m <- gc_genes_only %>%
  group_by(symbol) %>%
  filter(strand == "-") %>%
  mutate(start_site=max(end)) %>%
  filter(end==start_site) %>%
  dplyr::select(-c(11)) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  ungroup()

gc4=rbind(gc4p, gc4m)
gc4$symbol=paste0(gc4$symbol, "_", gc4$start)

T47 <- pI(bed=gc4,
              bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/T47D_R1.BAM",
              pairedEnd=FALSE, pause_s=0, pause_e=150)
MCF7 <- pI(bed=gc4,
               bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/MCF7_R1.BAM",
               pairedEnd = FALSE, pause_s=0, pause_e = 150)
MCF10A <- pI(bed=gc4,
                 bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/MCF10A_R1.BAM",
                 pairedEnd = FALSE, pause_s=0, pause_e = 150)
SUM159 <- pI(bed=gc4,
                 bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/SUM159_R1.BAM",
                 pairedEnd = FALSE, pause_s=0, pause_e = 150)
MB231 <- pI(bed=gc4,
                bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/MB231_R1.BAM",
                pairedEnd = FALSE, pause_s=0, pause_e = 150)

#i am currently growing a list. I should fix that? smol so doesn't make big dif but still bad practice
pI_tab <- cbind(T47, MCF7[,13:18], SUM159[,13:18], MB231[,13:18], MCF10A[,13:18])
colnames(pI_tab) = c(colnames(pI_tab[,1:12]),
                     paste0("T47D",           ".", colnames(pI_tab[,13:18])),
                     paste0("MCF7",           ".", colnames(pI_tab[,13:18])),
                     paste0("SUM159",         ".", colnames(pI_tab[,13:18])),
                     paste0("MB231",          ".", colnames(pI_tab[,13:18])),
                     paste0("MCF10A",         ".", colnames(pI_tab[,13:18])))
write.table(pI_tab, file="cell_lines_pI_3_no_filter.txt", quote=FALSE, sep="\t", row.names = FALSE)
#############################################################
pI_tab <- read.delim("~/cell_lines_pI_3_no_filter.txt", stringsAsFactors = FALSE)
#okay to make my life easier changing colnames, bc regex for these is a mess
pI_tab_longform = pI_tab %>%
  pivot_longer(cols=c(13:42), names_to = c("cell_line", ".value"), names_pattern = "(.+)\\.(.+)")

####
###
####

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
  geom_violin(aes(y=total_norm, x=cell_line)) +
  mytheme + 
  ylab("Normalized counts") +
  xlab("") +
  scale_x_discrete(labels=c("T47D","MCF7","SUM159","MB231","MCF10A")) #most are under 5

ggplot(data=pI_tab_longform) +
  geom_point(aes(y=total_norm, x=pI, color=cell_line)) +
  mytheme + 
  ylab("Body counts, length normalized") +
  xlab("pI")
##
summary(pI_tab_longform$pI)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.000   0.000   4.575     Inf  20.041     Inf    9362 
summary(pI_tab_longform$total_norm)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#   0.0000   0.0188   0.3597   1.4103   1.1316 555.9950
#okay how about this:
pI_tab_longform2 = pI_tab_longform %>%
  filter(!is.na(pI) & pI!="Inf")
summary(pI_tab_longform2$pI)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000    0.000    4.525   18.646   19.814 4002.022 
summary(pI_tab_longform2$total_norm)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.0003   0.0566   0.4755   1.5700   1.2642 555.9950  
ggplot(data=pI_tab_longform) +
  geom_violin(aes(y=log10(total_norm), x=cell_line)) +
  mytheme + 
  ylab("Normalized log counts") +
  xlab("") +
  scale_x_discrete(labels=c("T47D","MCF7","SUM159","MB231","MCF10A"))
ggplot(data=pI_tab_longform) +
  geom_violin(aes(y=log10(pI), x=cell_line)) +
  mytheme + 
  ylab("Normalized log pI") +
  xlab("") +
  scale_x_discrete(labels=c("T47D","MCF7","SUM159","MB231","MCF10A"))
pI_tab_longform2 = pI_tab_longform %>%
  filter(!is.na(pI) & pI!="Inf" & total_norm>=1) #puts us at 28841 rows, unknown num genes
ggplot(data=pI_tab_longform2) +
  geom_violin(aes(y=log10(pI), x=cell_line)) +
  mytheme + 
  ylab("Normalized log pI") +
  xlab("") +
  scale_x_discrete(labels=c("T47D","MCF7","SUM159","MB231","MCF10A"))
summary(pI_tab_longform2$pI)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000    1.873    7.540   16.314   18.858 2488.127
pI_tab_wideform = pI_tab_longform2 %>%
  pivot_wider(values_from=c(13:19), names_from = cell_line, names_sep=".",
              names_vary="slowest", values_fill = NA)
#8928 genes total
write.table(pI_tab_wideform, file="project_cell_lines_pI.txt", quote=FALSE, sep="\t", row.names = FALSE)
