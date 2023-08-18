###Pre meeting attempt####
H2AZ_H3K4me3_closest <- read.delim("~/H2AZ_H3K4me3_closest.bed", header=FALSE)
summary(H2AZ_H3K4me3_closest$V13)

library(tidyverse)

H3K4me3_far_up <- H2AZ_H3K4me3_closest %>% filter(V13< -10000)
H3K4me3_far_down <- H2AZ_H3K4me3_closest %>% filter(V13>10000)

direct_olap <- H2AZ_H3K4me3_closest %>% filter(V13==0)
H3K4me3_close_up <- H2AZ_H3K4me3_closest %>% filter(V13>=-10000 & V13 < 0)
H3K4me3_close_down <- H2AZ_H3K4me3_closest %>% filter(V13<=10000 & V13 > 0)

direct_olap$H2AZ_center = round((direct_olap$V3-direct_olap$V2)/2)+direct_olap$V2
direct_olap$H3K4me3_center = round((direct_olap$V9-direct_olap$V8)/2)+direct_olap$V8
direct_olap$peak_dif = direct_olap$H3K4me3_center - direct_olap$H2AZ_center

summary(direct_olap$peak_dif)
#Min.   1st Qu.    Median 
#-124460.0    -450.0      10.0 
#Mean   3rd Qu.      Max. 
#19.1     465.0  168875.0 

library(ChIPseeker)
library(BSgenome)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(EnsDb.Hsapiens.v75)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

totalgrange = makeGRangesFromDataFrame(direct_olap, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(totalgrange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
annoDF = as.data.frame(peakAnno)

strand1 = annoDF %>% filter(geneStrand==1)
strand2 = annoDF %>% filter(geneStrand==2)

summary(strand1$peak_dif)
#Min. 1st Qu.  Median    Mean 
#-123530    -125     255     357 
#3rd Qu.    Max. 
#705  149130

summary(strand2$peak_dif)
#Min. 1st Qu.  Median    Mean 
#-124460    -735    -265    -339 
#3rd Qu.    Max. 
#135  168875 

save(H3K4me3_close_down, H3K4me3_close_up, H3K4me3_far_down, H3K4me3_far_up, strand1, strand2, file="H2AZ_H3K4me3_relation.RData")
load("H2AZ_H3K4me3_relation.RData")

library(tidyverse)
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
#What categories would I want for graphing? Maybe just lump all the non-overlapping together for now, as even for the closest ones may still not be associated ... might as well annotate tho, then I can assign them strands and also figure out what the annotation difs are.

library(ChIPseeker)
library(BSgenome)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(EnsDb.Hsapiens.v75)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

far = rbind(H3K4me3_far_down, H3K4me3_far_up)
totalgrange = makeGRangesFromDataFrame(far, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(totalgrange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
annoDF = as.data.frame(peakAnno)

strand1_far = annoDF %>% filter(geneStrand==1)
strand2_far = annoDF %>% filter(geneStrand==2)

close = rbind(H3K4me3_close_down, H3K4me3_close_up)
totalgrange = makeGRangesFromDataFrame(close, start.field = "V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = TRUE)
peakAnno <- annotatePeak(totalgrange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
annoDF = as.data.frame(peakAnno)

strand1_close = annoDF %>% filter(geneStrand==1)
strand2_close = annoDF %>% filter(geneStrand==2)

#okay now what? Well, what are the essential pieces of each table?

annoDF_sum = annoDF %>% add_count(V13) %>% group_by(V13, simple_anno) %>% mutate(count_anno_by_clust=n()) %>% mutate(per_id=count_anno_by_clust/n) %>% summarise(percentage=first(per_id)) #then should summarize into 3 column table.
ggplot(data=annoDF_sum, aes(y=percentage, x=simple_anno)) + mytheme + geom_col(aes(fill=V13), position=position_dodge()) + xlab("annotation") +labs(fill="cluster")

get_simplified_anno <- function(tabs){
    introns <- as.data.frame(tabs) %>%
      filter(grepl('Intron(.*)', "annotation"))
    not_introns  <- as.data.frame(tabs) %>%
      filter(!grepl('Intron(.*)', "annotation"))
    exons <- not_introns %>%
      filter(grepl('Exon(.*)', "annotation"))
    not_intron_or_exon <- not_introns %>%
      filter(!grepl('Exon(.*)', "annotation"))
    first_intron  <- introns %>%
      filter(grepl('(.*)intron 1(.*)', "annotation"))
    not_first_intron <- introns %>%
      filter(!grepl('(.*)intron 1(.*)', "annotation"))
    first_exon  <- exons %>%
      filter(grepl('(.*)exon 1(.*)', "annotation"))
    not_first_exon <- exons %>%
      filter(!grepl('(.*)exon 1(.*)', "annotation"))
    not_intron_or_exon$simple_anno = not_intron_or_exon$annotation
    first_intron$simple_anno = rep("1st Intron", dim(first_intron)[1])
    first_exon$simple_anno = rep("1st Exon", dim(first_exon)[1])
    not_first_intron$simple_anno = rep("Other Intron", dim(not_first_intron)[1])
    not_first_exon$simple_anno = rep("Other Exon", dim(not_first_exon)[1])
    annoDF= rbind(not_intron_or_exon, not_first_exon, not_first_intron, first_exon, first_intron)
  return(annoDF)
}

close2 = rbind(strand1_close, strand2_close)
simple_anno_close <- get_simplified_anno(tab=close2)
far2 = rbind(strand1_far, strand2_far)
simple_anno_far <- get_simplified_anno(tab=far2)

#what about with the actually overlapping ones? peak center is within 500 bp, 500-2000, >2000, before or after. Remember everything is relative to H2AZ. So plus strand negative dif means H3K4me3 is before H2AZ.
simple_anno_strand1 <- get_simplified_anno(tab=strand1)
simple_anno_strand2 <- get_simplified_anno(tab=strand2)
Peak_center_500_before_strand1 = strand1 %>% filter(peak_dif>=-500 & peak_dif<0) #826
Peak_center_500_after_strand1 = strand1 %>% filter(peak_dif<=500 & peak_dif >=0) #1342
Peak_center_500plus_before_strand1 = strand1 %>% filter(peak_dif< -500) #430
Peak_center_500plus_after_strand1 = strand1 %>% filter(peak_dif>500) #1365
###
Peak_center_500_after_strand2 = strand2 %>% filter(peak_dif>=-500 & peak_dif<0) #1202
Peak_center_500_before_strand2 = strand2 %>% filter(peak_dif<=500 & peak_dif >=0) #759
Peak_center_500plus_after_strand2 = strand2 %>% filter(peak_dif< -500) #1329
Peak_center_500plus_before_strand2 = strand2 %>% filter(peak_dif>500) #449
#probably need to revise my filters but this feels like a decent starting point, to capture the offset peaks. Maybe look at some from the categories to see? After all maybe for not offset it's actually best to do the centers as 250 off in either direction then call the rest offset? Or could use location of maximum peak density as in columns V6 and V12, but then I would need to figure out how to parse that. String split at the dash and colon? Then do the centers of that region instead of the centers of the peaks?

###
#Get the max signal area?
H2AZ_max_sig_list_str1 = str_split(strand1$V6, "[:punct:]")
H2AZ_max_sig_tab_str1 = as.data.frame(do.call(rbind, H2AZ_max_sig_list_str1))
H2AZ_max_sig_tab_str1$H2AZ_max_center = (as.numeric(H2AZ_max_sig_tab_str1$V3) - as.numeric(H2AZ_max_sig_tab_str1$V2))/2 + as.numeric(H2AZ_max_sig_tab_str1$V2)

H3K4me3_max_sig_list_str1 = str_split(strand1$V12, "[:punct:]")
H3K4me3_max_sig_tab_str1 = as.data.frame(do.call(rbind, H3K4me3_max_sig_list_str1))
H3K4me3_max_sig_tab_str1$H3K4me3_max_center = (as.numeric(H3K4me3_max_sig_tab_str1$V3) - as.numeric(H3K4me3_max_sig_tab_str1$V2))/2 + as.numeric(H3K4me3_max_sig_tab_str1$V2)

strand1_w_max_reg <- cbind(strand1, H2AZ_max_sig_tab_str1, H3K4me3_max_sig_tab_str1)
strand1_w_max_reg$peak_max_dif <- strand1_w_max_reg$H3K4me3_max_center - strand1_w_max_reg$H2AZ_max_center
##and strand 2:
H2AZ_max_sig_list_str2 = str_split(strand2$V6, "[:punct:]")
H2AZ_max_sig_tab_str2 = as.data.frame(do.call(rbind, H2AZ_max_sig_list_str2))
H2AZ_max_sig_tab_str2$H2AZ_max_center = (as.numeric(H2AZ_max_sig_tab_str2$V3) - as.numeric(H2AZ_max_sig_tab_str2$V2))/2 + as.numeric(H2AZ_max_sig_tab_str2$V2)

H3K4me3_max_sig_list_str2 = str_split(strand2$V12, "[:punct:]")
H3K4me3_max_sig_tab_str2 = as.data.frame(do.call(rbind, H3K4me3_max_sig_list_str2))
H3K4me3_max_sig_tab_str2$H3K4me3_max_center = (as.numeric(H3K4me3_max_sig_tab_str2$V3) - as.numeric(H3K4me3_max_sig_tab_str2$V2))/2 + as.numeric(H3K4me3_max_sig_tab_str2$V2)

strand2_w_max_reg <- cbind(strand2, H2AZ_max_sig_tab_str2, H3K4me3_max_sig_tab_str2)
strand2_w_max_reg$peak_max_dif <- strand2_w_max_reg$H3K4me3_max_center - strand2_w_max_reg$H2AZ_max_center
#so what do these look like?
summary(strand1_w_max_reg$peak_max_dif)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-32850.0    -50.0    200.0    374.8    685.0 242795.0 #similar to the other method but a bit more centralized
summary(strand2_w_max_reg$peak_max_dif)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-88750.0   -700.0   -190.0   -145.0     77.5 292440.0 

#these are easy to regenerate, just save the document

#### post meeting version ####
H2AZ_H3K4me3_closest <- read.delim("~/H2AZ_H3K4me3_closest.bed", header=FALSE)
library(tidyverse)
direct_olap <- H2AZ_H3K4me3_closest %>% filter(V13==0)

H2AZ_max_sig_list = str_split(direct_olap$V6, "[:punct:]")
H2AZ_max_sig_tab = as.data.frame(do.call(rbind, H2AZ_max_sig_list))
H2AZ_max_sig_tab$H2AZ_max_center = (as.numeric(H2AZ_max_sig_tab$V3) - as.numeric(H2AZ_max_sig_tab$V2))/2 + as.numeric(H2AZ_max_sig_tab$V2)

H3K4me3_max_sig_list = str_split(direct_olap$V12, "[:punct:]")
H3K4me3_max_sig_tab = as.data.frame(do.call(rbind, H3K4me3_max_sig_list))
H3K4me3_max_sig_tab$H3K4me3_max_center = (as.numeric(H3K4me3_max_sig_tab$V3) - as.numeric(H3K4me3_max_sig_tab$V2))/2 + as.numeric(H3K4me3_max_sig_tab$V2)

direct_olap_w_max_reg <- cbind(direct_olap, H2AZ_max_sig_tab, H3K4me3_max_sig_tab)
direct_olap_w_max_reg$peak_max_dif <- direct_olap_w_max_reg$H3K4me3_max_center - direct_olap_w_max_reg$H2AZ_max_center
summary(direct_olap_w_max_reg$peak_max_dif)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-88750.0   -435.0      0.0    122.5    455.0 292440.0 
summary(abs(direct_olap_w_max_reg$peak_max_dif))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.0    150.0    450.0    795.4    870.0 292440.0
#For now lets split it into abs quartiles then see what that looks like when graphed, how offset the various groups are.
names(direct_olap_w_max_reg) = c("H2AZ_chrom", "H2AZ_peak_start", "H2AZ_peak_end", "H2AZ_total_signal", "H2AZ_max_signal", "H2AZ_max_signal_region", "H3K4me3_chrom", "H3K4me3_peak_start", "H3K4me3_peak_end", "H3K4me3_total_signal", "H3K4me3_max_signal", "H3K4me3_max_signal_region", "peak_dif", "H2AZ_max_chrom", "H2AZ_max_start", "H2AZ_max_end", "H2AZ_max_center", "H3K4me3_max_chrom", "H3K4me3_max_start", "H3K4me3_max_end", "H3K4me3_max_center", "peak_max_dif")
q1 = direct_olap_w_max_reg %>% filter(abs(peak_max_dif)<=150)
q2 = direct_olap_w_max_reg %>% filter(abs(peak_max_dif)>150 & abs(peak_max_dif)<=450)
q3 = direct_olap_w_max_reg %>% filter(abs(peak_max_dif)>450 & abs(peak_max_dif)<=870)
q4 = direct_olap_w_max_reg %>% filter(abs(peak_max_dif)>870)
######### characterizing the overlaps ###########
#Characterize? Q4 may not be of use, based on visual perusal of high difs from earlier in the document, I suspect some of them come from what should have been two peaks being called as a single peak
merged_regions=direct_olap_w_max_reg %>%
  rowwise() %>%
  mutate(merge_start=min(H2AZ_peak_start, H3K4me3_peak_start), merge_end=max(H2AZ_peak_end, H3K4me3_peak_end))

merged_regions_reorder = merged_regions[,c(1,23,24,1,2:22)]

write_bed_file <- function(tab, b_name) {
  write.table(
    tab,
    file=b_name,
    row.names=FALSE,
    col.names=FALSE,
    quote=FALSE,
    sep="\t"
  )
  #You need an appropiately formatted (3 or 6 columns in correct order) table for this but these are good params for making a tab-separated genomic table of any form.
}

write_bed_file(merged_regions_reorder, "H2AZ_H3K4me3_olaps.txt")
##### how do the overlaps interact w cgis/KDM5B? ######
H2AZ_H3K4me3_wcgi <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/2021.11.17 ZS 1_19 Cut and Tag/h3k4me3_h2az_offset/H2AZ_H3K4me3_wcgi.txt", header=FALSE)
H2AZ_H3K4me3_woutcgi <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/H2AZ_H3K4me3_woutcgi.txt", header=FALSE)
H2AZ_H3K4me3_wkdm5b <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/H2AZ_H3K4me3_wkdm5b.txt", header=FALSE)
H2AZ_H3K4me3_woutkdm5b <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/H2AZ_H3K4me3_woutkdm5b.txt", header=FALSE)

library(tidyverse)

summary(H2AZ_H3K4me3_wcgi$V24)
#Min.  1st Qu.   Median     Mean  3rd Qu. Max.
#-88750.0   -460.0      0.0    237.3    480.0 292440.0 
summary(abs(H2AZ_H3K4me3_wcgi$V24))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0     155     470    1039     920  292440

summary(H2AZ_H3K4me3_woutcgi$V24)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-14820.00   -257.50      0.00     26.29    297.50   7610.00
summary(abs(H2AZ_H3K4me3_woutcgi$V24))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0   100.0   280.0   506.4   600.0 14820.0 

summary(H2AZ_H3K4me3_wkdm5b$V24)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-88750    -460       0    7594     595  292440 
summary(abs(H2AZ_H3K4me3_wkdm5b$V24))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0     160     530    9360    1070  292440 

summary(H2AZ_H3K4me3_woutkdm5b_w_max_reg$V24)
#Min.    1st Qu.     Median       Mean    3rd Qu. 
#-15040.000   -410.000      0.000      9.929    415.000 
#Max. 
#7610.000 
summary(abs(H2AZ_H3K4me3_woutkdm5b_w_max_reg$V24))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0   150.0   410.0   583.4   800.0 15040.0 
##### where are these overlaps located? How do peak centers dif in dif ones? #######
H2AZ_H3K4me3_wcgi$CGI_status = rep(TRUE, dim(H2AZ_H3K4me3_wcgi)[1])
H2AZ_H3K4me3_woutcgi$CGI_status = rep(FALSE, dim(H2AZ_H3K4me3_woutcgi)[1])
cgi_status_tab = rbind(H2AZ_H3K4me3_wcgi, H2AZ_H3K4me3_woutcgi)

H2AZ_H3K4me3_wkdm5b$KDM5B_status = rep(TRUE, dim(H2AZ_H3K4me3_wkdm5b)[1])
H2AZ_H3K4me3_woutkdm5b$KDM5B_status = rep(FALSE, dim(H2AZ_H3K4me3_woutkdm5b)[1])
kdm5b_status_tab = rbind(H2AZ_H3K4me3_wkdm5b, H2AZ_H3K4me3_woutkdm5b)

#CGI status is longer than KDM5B status though... Okay, original olaps are 7702.
merged_regions_reorder$uniq_id = paste0(merged_regions_reorder$H2AZ_peak_start, ".", merged_regions_reorder$H2AZ_peak_end, ".", merged_regions_reorder$H3K4me3_peak_start)
merged_regions_unique = merged_regions_reorder %>% distinct(uniq_id, .keep_all = TRUE) %>% select(-c(1))

cgi_status_tab$uniq_id = paste0(cgi_status_tab$V5, ".", cgi_status_tab$V6, ".", cgi_status_tab$V11)
cgi_unique = cgi_status_tab %>% distinct(uniq_id, .keep_all = TRUE)

cgi_lt = cgi_status_tab$CGI_status
names(cgi_lt) = cgi_status_tab$uniq_id
merged_regions_unique$CGI_status = unname((cgi_lt[as.character(merged_regions_unique$uniq_id)]))

kdm5b_status_tab$uniq_id = paste0(kdm5b_status_tab$V5, ".", kdm5b_status_tab$V6, ".", kdm5b_status_tab$V11)
kdm5b_unique = kdm5b_status_tab %>% distinct(uniq_id, .keep_all = TRUE)

kdm5b_lt = kdm5b_status_tab$KDM5B_status
names(kdm5b_lt) = kdm5b_status_tab$uniq_id
merged_regions_unique$KDM5B_status = unname((kdm5b_lt[as.character(merged_regions_unique$uniq_id)]))

write.table(
  merged_regions_unique,
  file="H2AZ_H3K4me3_with_cgi_kdm5b_stats.txt",
  row.names=FALSE,
  quote=FALSE,
  sep="\t"
)

#### Okay now annotate ####
library(ChIPseeker)
library(BSgenome)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(EnsDb.Hsapiens.v75)
library(clusterProfiler)

totalgrange = makeGRangesFromDataFrame(merged_regions_unique, start.field = "merge_start", end.field = "merge_end", seqnames.field = "H2AZ_chrom", keep.extra.columns = TRUE)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(totalgrange, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
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
  annoDF= rbind(not_intron_or_exon, not_first_exon, not_first_intron, first_exon, first_intron)
  return(annoDF)
}

anno_df_w_simple_col = get_simplified_anno(annoDF)
#Okay, now I have a table with a lot going on. Important categories include: KDM5B status, CGI status, annotation, and peak max diff. How do these relate? 1. Box plots could be a good idea.

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

ggplot(data = anno_df_w_simple_col) + geom_boxplot(aes(x=simple_anno, y=peak_max_dif)) + mytheme
ggplot(data = anno_df_w_simple_col) + geom_boxplot(aes(x=simple_anno, y=peak_max_dif)) + ylim(-500,500) + mytheme

ggplot(data = anno_df_w_simple_col) + geom_boxplot(aes(x=CGI_status, y=peak_max_dif)) + mytheme
ggplot(data = anno_df_w_simple_col) + geom_boxplot(aes(x=CGI_status, y=peak_max_dif)) + ylim(-5000,5000) + mytheme

ggplot(data = anno_df_w_simple_col) + geom_boxplot(aes(x=KDM5B_status, y=peak_max_dif)) + mytheme
ggplot(data = anno_df_w_simple_col) + geom_boxplot(aes(x=KDM5B_status, y=peak_max_dif)) + ylim(-5000,5000) + mytheme

ggplot(data = anno_df_w_simple_col) + geom_boxplot(aes(x=KDM5B_status, y=peak_max_dif, color=simple_anno)) + ylim(-5000,5000) + mytheme

ggplot(data=anno_df_w_simple_col) + geom_bar(aes(simple_anno)) + mytheme
ggplot(data=anno_df_w_simple_col) + geom_bar(aes(simple_anno)) + mytheme + ylim(0,150)
ggplot(data=anno_df_w_simple_col) + geom_bar(aes(CGI_status)) + mytheme
ggplot(data=anno_df_w_simple_col) + geom_bar(aes(KDM5B_status)) + mytheme

H2AZ_H3K4me3_annotated_homer <- read_delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/H2AZ_H3K4me3_annotated_homer.txt", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 1)
colnames(H2AZ_H3K4me3_annotated_homer) <- c("Peak_ID", "Chromosome", "Peak_start_pos", "Peak_end_pos", "Strand", "Peak_score", "FDR_Focus_Ratio_and_Size", "Annotation", "Detailed_Anno", "Dist_Refseq_TSS", "nTSS_nativeID", "nTSS_entrezID", "nTSS_unigeneID", "nTSS_Refseq_ID", "nTSS_ensemblID", "nTSS_geneSymbol", "nTSS_geneAlias", "nTSS_geneDescrip", "nTSS_geneType")

anno_df_w_simple_col$uniq_id2 <- paste0(anno_df_w_simple_col$start, ".", anno_df_w_simple_col$end)
H2AZ_H3K4me3_annotated_homer$uniq_id2 <- paste0((H2AZ_H3K4me3_annotated_homer$Peak_start_pos-1), ".", H2AZ_H3K4me3_annotated_homer$Peak_end_pos)
H2AZ_H3K4me3_annotated_homer_trimmed <- H2AZ_H3K4me3_annotated_homer %>% select(-c(1:7, 12:13))
anno_df_full <- merge.data.frame(anno_df_w_simple_col, H2AZ_H3K4me3_annotated_homer_trimmed, by.x="uniq_id2", by.y="uniq_id2")
anno_df_full <- anno_df_full %>% select(-c(1, 6, 11, 12, 17:19, 23))

#get simple homer anno
introns <- as.data.frame(anno_df_full) %>%
  filter(grepl('^intron(.*)', Annotation)) #
not_introns  <- as.data.frame(anno_df_full) %>%
  filter(!grepl('^intron(.*)', Annotation))
exons <- not_introns %>%
  filter(grepl('^exon(.*)', Annotation)) #
not_intron_or_exon <- not_introns %>%
  filter(!grepl('^exon(.*)', Annotation))
utr_3 <- not_intron_or_exon %>%
  filter(grepl("3[[:punct:]][[:space:]]UTR(.)*", Annotation)) #
not_utr_3 <- not_intron_or_exon %>%
  filter(!grepl("3[[:punct:]][[:space:]]UTR(.)*", Annotation))
utr_5 <- not_utr_3 %>%
  filter(grepl("5[[:punct:]][[:space:]]UTR(.)*", Annotation)) #
not_utr_5 <- not_utr_3 %>%
  filter(!grepl("5[[:punct:]][[:space:]]UTR(.)*", Annotation))
non_coding_genes <- not_utr_5 %>%
  filter(grepl("non[[:punct:]]coding(.)*", Annotation)) #
not_non_coding <- not_utr_5 %>%
  filter(!grepl("non[[:punct:]]coding(.)*", Annotation))
promoter <- not_non_coding %>%
  filter(grepl("promoter(.)*", Annotation)) #
not_promoter <- not_non_coding %>%
  filter(!grepl("promoter(.)*", Annotation))
TTS <- not_promoter %>%
  filter(grepl("TTS(.)*", Annotation)) #
notTTS <- not_promoter %>%
  filter(!grepl("TTS(.)*", Annotation)) #
#there's probably an easier way to do this but it's simpler to figure out the regular expressions this way.

introns$homer_simple_anno = rep("intron", dim(introns)[1])
exons$homer_simple_anno = rep("exon", dim(exons)[1])
utr_3$homer_simple_anno = rep("3_UTR", dim(utr_3)[1])
utr_5$homer_simple_anno = rep("5_UTR", dim(utr_5)[1])
non_coding_genes$homer_simple_anno = rep("non_coding", dim(non_coding_genes)[1])
promoter$homer_simple_anno = rep("promoter_tss", dim(promoter)[1])
TTS$homer_simple_anno = rep("TTS", dim(TTS)[1])
notTTS$homer_simple_anno = notTTS$Annotation

Annotated_olaps = rbind(introns, exons, utr_3, utr_5, non_coding_genes, promoter, TTS, notTTS)

write.table(
  Annotated_olaps,
  file="H2AZ_H3K4me3_with_cgi_kdm5b_stats.txt",
  row.names=FALSE,
  quote=FALSE,
  sep="\t"
)
#########################################################
t1 = Annotated_olaps %>% group_by(KDM5B_status, CGI_status, simple_anno) %>% summarize(num=n(), min_dif = min(peak_max_dif), q1 = quantile(peak_max_dif, 0.25), q2=median(peak_max_dif), q3=quantile(peak_max_dif, 0.75), max=max(peak_max_dif), spread=sd(peak_max_dif))

t2 = Annotated_olaps %>% group_by(KDM5B_status, CGI_status, homer_simple_anno) %>% summarize(num=n(), min_dif = min(peak_max_dif), q1 = quantile(peak_max_dif, 0.25), q2=median(peak_max_dif), q3=quantile(peak_max_dif, 0.75), max=max(peak_max_dif), spread=sd(peak_max_dif))

t3 = Annotated_olaps %>% group_by(KDM5B_status) %>% summarize(num=n(), min_dif = min(peak_max_dif), q1 = quantile(peak_max_dif, 0.25), q2=median(peak_max_dif), q3=quantile(peak_max_dif, 0.75), max=max(peak_max_dif), spread=sd(peak_max_dif))

t4 = Annotated_olaps %>% group_by(CGI_status) %>% summarize(num=n(), min_dif = min(peak_max_dif), q1 = quantile(peak_max_dif, 0.25), q2=median(peak_max_dif), q3=quantile(peak_max_dif, 0.75), max=max(peak_max_dif), spread=sd(peak_max_dif))

t5 = Annotated_olaps %>% group_by(simple_anno) %>% summarize(num=n(), min_dif = min(peak_max_dif), q1 = quantile(peak_max_dif, 0.25), q2=median(peak_max_dif), q3=quantile(peak_max_dif, 0.75), max=max(peak_max_dif), spread=sd(peak_max_dif))

ggplot(data = Annotated_olaps) + geom_boxplot(aes(x=KDM5B_status, y=peak_max_dif)) + ylim(-10000,10000) + mytheme
ggplot(data = Annotated_olaps) + geom_boxplot(aes(x=CGI_status, y=peak_max_dif)) + ylim(-10000,10000) + mytheme
ggplot(data = Annotated_olaps) + geom_boxplot(aes(x=simple_anno, y=peak_max_dif)) + ylim(-10000,10000) + mytheme
ggplot(data = Annotated_olaps) + geom_boxplot(aes(x=simple_anno, y=peak_max_dif)) + ylim(-2500,2500) + mytheme
ggplot(data = Annotated_olaps) + geom_boxplot(aes(x=homer_simple_anno, y=peak_max_dif)) + ylim(-7500,7500) + mytheme

H2AZ_H3K4me3_with_cgi_kdm5b_stats <- read.delim("~/H2AZ_H3K4me3_with_cgi_kdm5b_stats.txt")
library(tidyverse)
chisq.test(H2AZ_H3K4me3_with_cgi_kdm5b_stats$KDM5B_status, H2AZ_H3K4me3_with_cgi_kdm5b_stats$CGI_status)
#X-squared = 4.7234, df = 1, p-value = 0.02975
Annotated_olaps <- H2AZ_H3K4me3_with_cgi_kdm5b_stats

ggplot(data = Annotated_olaps) + geom_violin(aes(x=simple_anno, y=peak_max_dif)) + ylim(-5000,5000) + mytheme
ggplot(data = Annotated_olaps) + geom_violin(aes(x=simple_anno, y=peak_max_dif, color=KDM5B_status)) + ylim(-5000,5000) + mytheme
ggplot(data = Annotated_olaps) + geom_violin(aes(x=simple_anno, y=peak_max_dif, color=CGI_status)) + ylim(-5000,5000) + mytheme


ggplot(data = Annotated_olaps) + geom_violin(aes(x=CGI_status, y=peak_max_dif)) + ylim(-5000,5000) + mytheme
ggplot(data = Annotated_olaps) + geom_violin(aes(x=KDM5B_status, y=peak_max_dif)) + ylim(-5000,5000) + mytheme
ggplot(data = Annotated_olaps) + geom_violin(aes(x=homer_simple_anno, y=peak_max_dif)) + ylim(-5000,5000) + mytheme

ggplot(data = Annotated_olaps) + geom_violin(aes(x=homer_simple_anno, y=peak_max_dif, color=KDM5B_status)) + ylim(-5000,5000) + mytheme
ggplot(data = Annotated_olaps) + geom_violin(aes(x=homer_simple_anno, y=peak_max_dif, color=CGI_status)) + ylim(-5000,5000) + mytheme
