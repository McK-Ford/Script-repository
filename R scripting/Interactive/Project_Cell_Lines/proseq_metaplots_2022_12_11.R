source("~/Script repository/My_Useful_Fns.R")
enableJIT(3)
tabf <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/PublicData/GSE118075_Matsuki_2019_MCF7_CAGE_netCAGE/getTSSes/tx_p.txt",
  header=FALSE
)
tabr <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/PublicData/GSE118075_Matsuki_2019_MCF7_CAGE_netCAGE/getTSSes/tx_m.txt",
  header=FALSE
)
colnames(tabf) = c("CGI_ID", "genechrom", "cgi_s", "cgi_e", "gene_s", "gene_e", "ensembl_ID", "strand",
                   "geneName", "unitprot_ID", "class", "transtype", "netcagetss", "tssscore")
colnames(tabr) = c("CGI_ID", "genechrom", "cgi_s", "cgi_e", "gene_s", "gene_e", "ensembl_ID", "strand",
                   "geneName", "unitprot_ID", "class", "transtype", "netcagetss", "tssscore")
tabf$scalestart=tabf$netcagetss
tabf$scaleend=tabf$cgi_e

tabf$tssstart=tabf$netcagetss
tabf$tssend=tabf$cgi_e

tabf$diststart=tabf$cgi_e
tabf$distend=tabf$gene_e

tabr$scalestart=tabr$cgi_s
tabr$scaleend=tabr$gene_e

tabr$tssstart=tabr$gene_s
tabr$tssend=tabr$netcagetss

tabr$diststart=tabr$gene_s
tabr$distend=tabr$cgi_s

tab=rbind(tabf, tabr)
tab$uniq_ID = paste0(tab$genechrom, "_", tab$netcagetss, "_", tab$geneName)

#scaled tab - netcagetss to cgi_e for plus, cgi_s to gene_e for minus.
#TSS tab - netcagetss to gene_e for plus, gene_s to netcagetss for minus.
#dist tab - cgi_e to gene_e for plus, gene_s to cgi_s for minus.
tab = tab %>%
  filter(gene_e-gene_s>=200)
tab_scaled = tab[,c(2, 15, 16, 21, 14, 8, 3:7, 9:12)]
tab_tss = tab[,c(2, 17, 18, 21, 14, 8, 3:7, 9:12)]
tab_cgi = tab[,c(2, 19, 20, 21, 14, 8, 3:7, 9:12)]
colnames(tab_scaled) <- c("chrom", "start", "end", "symbol",
                          "score", "strand", colnames(tab_scaled[7:15]))
colnames(tab_tss) <- c("chrom", "start", "end", "symbol",
                       "score", "strand", colnames(tab_tss[7:15]))
colnames(tab_cgi) <- c("chrom", "start", "end", "symbol",
                       "score", "strand", colnames(tab_cgi[7:15]))

val_tss = score_matrix(bed=tab_tss,
                       bam="C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/PublicData/GSE93229_Danko_2018_MCF7_PROseq/bamdeduped/MCF7_B7.deduped.bam",
                       pairedEnd=FALSE, b=-50, a=250, mode="sbp", revcomp=TRUE, method="single_anch", rnorm=FALSE)
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
setDT(val_tss)
mtss_long = data.table::melt(val_tss, measure.vars=c(16:45))
mtss_long$variable = as.numeric(as.character(mtss_long$variable))
mtss_avbybin_all = mtss_long %>%
  group_by(variable) %>%
  summarise(score=mean(value)) %>%
  ungroup()

tss1 = ggplot(data=mtss_avbybin_all, mapping=aes(
  x=variable, y=score)) + geom_line() +
  mytheme +
  scale_x_continuous(n.breaks = 6, name="Distance from TSS (bp)") +
  scale_y_continuous(name="Average read density") +
  ggtitle("GSE93229 (Danko) MCF7:B7 TSS-centered PRO-seq")
tss1
#plus version bed
#minus version bed
val_cgi = score_matrix(bed=tab_cgi,
                       bam="C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/PublicData/GSE93229_Danko_2018_MCF7_PROseq/bamdeduped/MCF7_B7.deduped.bam",
                       pairedEnd=FALSE, b=-250, a=250, mode="sbp", revcomp=TRUE, method="single_anch", rnorm=FALSE)
setDT(val_cgi)
mcgi_long = data.table::melt(val_cgi, measure.vars=c(16:65))
mcgi_long$variable = as.numeric(as.character(mcgi_long$variable))
mcgi_avbybin_all = mcgi_long %>%
  group_by(variable) %>%
  summarise(score=mean(value)) %>%
  ungroup()

cgi1 = ggplot(data=mcgi_avbybin_all, mapping=aes(
  x=variable, y=score)) + geom_line() +
  mytheme +
  scale_x_continuous(n.breaks = 6, name="Distance from 3' CGI (bp)") +
  scale_y_continuous(name="Average read density") +
  ggtitle("GSE93229 (Danko) MCF7:B7 CGI-centered PRO-seq")
cgi1
##############

val_tss2 = score_matrix(bed=tab_tss,
                       bam="C:/Users/kayle/Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/bamFiles/BAMdeDuped/MCF7_deDuped.BAM",
                       pairedEnd=FALSE, b=-50, a=250, mode="sbp", revcomp=TRUE, method="single_anch", rnorm=FALSE)
setDT(val_tss2)
mtss_long2 = data.table::melt(val_tss2, measure.vars=c(16:45))
mtss_long2$variable = as.numeric(as.character(mtss_long2$variable))
mtss_avbybin_all2 = mtss_long2 %>%
  group_by(variable) %>%
  summarise(score=mean(value)) %>%
  ungroup()

tss2 = ggplot(data=mtss_avbybin_all2, mapping=aes(
  x=variable, y=score)) + geom_line() +
  mytheme +
  scale_x_continuous(n.breaks = 6, name="Distance from TSS (bp)") +
  scale_y_continuous(name="Average read density") +
  ggtitle("VER0059 (Project Cell Lines) MCF7 TSS-centered PRO-seq")
tss2
#plus version bed
#minus version bed
val_cgi2 = score_matrix(bed=tab_cgi,
                       bam="C:/Users/kayle/Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/bamFiles/BAMdeDuped/MCF7_deDuped.BAM",
                       pairedEnd=FALSE, b=-250, a=250, mode="sbp", revcomp=TRUE, method="single_anch", rnorm=FALSE)
setDT(val_cgi2)
mcgi_long2 = data.table::melt(val_cgi2, measure.vars=c(16:65))
mcgi_long2$variable = as.numeric(as.character(mcgi_long2$variable))
mcgi_avbybin_all2 = mcgi_long2 %>%
  group_by(variable) %>%
  summarise(score=mean(value)) %>%
  ungroup()

cgi2 = ggplot(data=mcgi_avbybin_all2, mapping=aes(
  x=variable, y=score)) + geom_line() +
  mytheme +
  scale_x_continuous(n.breaks = 6, name="Distance from 3' CGI (bp)") +
  scale_y_continuous(name="Average read density") +
  ggtitle("VER0059 (Project Cell Lines) MCF7 CGI-centered PRO-seq")
cgi2


pdf(file="proseq_metaplots_2022_12_11.pdf")
tss1
cgi1
tss2
cgi2
print("Data is:")
print("For first two graphs, GSE93229. This is the Danko MCF7 'B7 subclone' single-end PRO-seq data.")
print("For second two graphs, VER0059. This is the project cell lines project, MCF7 line specifically, paired end library.")
dev.off()
