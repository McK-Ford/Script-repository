TSS_plus <- read.delim( #Reading in the CpG/Genes intersection file, from bedtools
  "plusmerge5.txt",
  header=FALSE,
  quote="")
TSS_minus <- read.delim( #Reading in the CpG/Genes intersection file, from bedtools
  "minusmerge5.txt",
  header=FALSE,
  quote="")

library(tidyverse)
#repair_last_col_m <- data.frame(unlist(str_split(TSS_minus$V34, " ", simplify = TRUE)))
#repair_last_col_p <- data.frame(unlist(str_split(TSS_plus$V34, " ", simplify = TRUE)))
#TSS_plus = cbind(TSS_plus, repair_last_col_p)
#TSS_minus = cbind(TSS_minus, repair_last_col_m)

TSS_minus_selected = TSS_minus %>% select(-c(2,5,9, 11:14, 18:21, 24, 26, 29:31))
TSS_plus_selected = TSS_plus %>% select(-c(3,5,9,11:14,18:21, 24, 26, 29:31))

colnames(TSS_minus_selected) = c("chrom", "tru_TSS", "TSS_score", "gene_s", "gene_e", "ensembl ID", "strand", "exon_sizes", "exon_starts", "UCSC_ID", "gene_name", "uniprot_ID", "tanscript_class", "transcript_type", "tag", "cgi_s", "cgi_e", "cgi_id")
colnames(TSS_plus_selected) = c("chrom", "tru_TSS", "TSS_score", "gene_s", "gene_e", "ensembl ID", "strand", "exon_sizes", "exon_starts", "UCSC_ID", "gene_name", "uniprot_ID", "tanscript_class", "transcript_type", "tag", "cgi_s", "cgi_e", "cgi_id")

#with the previous 'first in each group' problem.
#summary(as.numeric(TSS_minus_selected$tru_annotated_tss_dist))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1      11      42    2668     302  572336 
#summary(as.numeric(TSS_plus_selected$tru_annotated_tss_dist))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1      10      35    2313     217  686931

summary(abs(TSS_minus_selected$tru_TSS - TSS_minus_selected$gene_e))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.0      2.0     17.0   1200.0    105.8 372610.0 
#So there are still some way too high but not as many

summary(abs(TSS_plus_selected$tru_TSS - TSS_plus_selected$gene_s))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0       3      16    1133      86  686931 


TSS_plus_selected$gene_len <- TSS_plus_selected$gene_e - TSS_plus_selected$tru_TSS
TSS_minus_selected$gene_len <- TSS_minus_selected$tru_TSS - TSS_minus_selected$gene_s

TSS_plus_selected$cgi_len <- TSS_plus_selected$cgi_e - TSS_plus_selected$cgi_s
TSS_minus_selected$cgi_len <- TSS_minus_selected$cgi_e - TSS_minus_selected$cgi_s


TSS_plus_selected$p_gene_in_cgi <- (TSS_plus_selected$gene_len - pmax((TSS_plus_selected$cgi_s - TSS_plus_selected$tru_TSS),0) - pmax((TSS_plus_selected$gene_e-TSS_plus_selected$cgi_e),0))/TSS_plus_selected$gene_len
TSS_minus_selected$p_gene_in_cgi <- (TSS_minus_selected$gene_len - pmax((TSS_minus_selected$cgi_s - TSS_minus_selected$gene_s),0) - pmax((TSS_minus_selected$tru_TSS-TSS_minus_selected$cgi_e),0))/TSS_minus_selected$gene_len

TSS_plus_selected$cgi5_tru_tss_dist <- TSS_plus_selected$tru_TSS-TSS_plus_selected$cgi_s
TSS_plus_selected$tss_tru_cgi3_dist <- TSS_plus_selected$cgi_e-TSS_plus_selected$tru_TSS

TSS_minus_selected$cgi5_tss_dist <- TSS_minus_selected$cgi_e-TSS_minus_selected$tru_TSS
TSS_minus_selected$tss_cgi3_dist <- TSS_minus_selected$tru_TSS-TSS_minus_selected$cgi_s


#really all I need to pull is josh's, I can remake skew and proseq classes for this if it is easier.

#oh did i remeber to filter out chrom Y?

TSS_plus_selected_1name = TSS_plus_selected
library(data.table)
setDT(TSS_plus_selected_1name)
TSS_plus_selected_1name = TSS_plus_selected_1name[, .SD[which.max(TSS_score)], by = gene_name]

setDT(TSS_minus_selected)
TSS_minus_selected_1name = TSS_minus_selected[, .SD[which.max(TSS_score)], by = gene_name]

summary(abs(TSS_minus_selected_1name$tru_TSS - TSS_minus_selected_1name$gene_e))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0       3      20    1553     134  372610 

summary(abs(TSS_plus_selected_1name$tru_TSS - TSS_plus_selected_1name$gene_s))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0       3      20    1610     114  686931 

cgi_pause <- read.delim( #Reading in the CpG/Genes intersection file, from bedtools
  "cgi_pause_df_hg38_10_8_21.txt",
  quote="")

#Lookup #1, wendy classes
cgi_pause_uniq = cgi_pause %>% distinct(name2, .keep_all = TRUE)
wendylt = cgi_pause_uniq$class
names(wendylt) = cgi_pause_uniq$name2
TSS_plus_selected_1name$wendy_class = unname((wendylt[(TSS_plus_selected_1name$gene_name)]))
TSS_minus_selected_1name$wendy_class = unname((wendylt[(TSS_minus_selected_1name$gene_name)]))

for_quick_dt_maps_p <- TSS_plus_selected_1name %>% select(c(2,3,6))
for_quick_dt_maps_m <- TSS_minus_selected_1name %>% select(c(2,5,3))

write.table(for_quick_dt_maps_p, "test_p.bed", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(for_quick_dt_maps_m, "test_m.bed", quote=FALSE, row.names = FALSE, col.names = FALSE)
