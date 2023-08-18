###################
## Initial setup ##
###################
library(dplyr)
refseq__cpg_intersect <- read.delim("~/refseq_curated_genes_cpg_intersect.txt.gz", header=FALSE)
no_alt_chroms <- refseq__cpg_intersect %>%
  filter(!grepl('chrY|([\\w_]+)alt|random|fix|v1', V1)) %>%
  filter(!grepl('chrY|([\\w_]+)alt|random|fix|v1', V4))
names(no_alt_chroms) <- c("cpgchrom", "cgi_s", "cgi_e", "cpg_id",
                          "chrom", "gene_s", "gene_e",
                          "gene_name", "score", "strand")
no_alt_chroms = unique(no_alt_chroms)
##################################
## Processing plus strand genes ##
##################################
plus_strand_only <- no_alt_chroms %>%
  filter(strand=="+") #Filters to only plus strand genes
plus_TSS_intersect_only <- plus_strand_only %>%
  filter(cgi_s<gene_s) #As bedtools intersect intersects anywhere in gene, this is 'start of CpG island is
#before start of gene' which ensures the TSS is within the gene.
shortest_tss_in_cgi <- plus_TSS_intersect_only %>% group_by(cgi_s) %>% mutate(tss=max(gene_s)) %>% filter(tss==gene_s) %>% distinct(cgi_s, .keep_all=TRUE) #
#7550 genes
###################################
## Processing minus strand genes ##
###################################
#Majority of logic is same as with plus strand above - if not commented here, check analogous
#operation there.
minus_strand_only <- no_alt_chroms %>%
  filter(strand=="-")
minus_TSS_intersect_only <- minus_strand_only %>%
  filter(cgi_e>gene_e)
shortest_tss_in_cgi_m <- minus_TSS_intersect_only %>% group_by(cgi_s) %>% mutate(tss=min(gene_e)) %>% filter(tss==gene_e) %>% distinct(cgi_s, .keep_all = TRUE) #drops it to 7369

##############
## Metadata ##
##############
shortest_tss_in_cgi$gene_len <- shortest_tss_in_cgi$gene_e - shortest_tss_in_cgi$gene_s
shortest_tss_in_cgi_m$gene_len <- shortest_tss_in_cgi_m$gene_e - shortest_tss_in_cgi_m$gene_s
shortest_tss_in_cgi$cgi_len <- shortest_tss_in_cgi$cgi_e - shortest_tss_in_cgi$cgi_s
shortest_tss_in_cgi_m$cgi_len <- shortest_tss_in_cgi_m$cgi_e - shortest_tss_in_cgi_m$cgi_s
shortest_tss_in_cgi$cgi5_tss_dist <- shortest_tss_in_cgi$gene_s-shortest_tss_in_cgi$cgi_s
shortest_tss_in_cgi$tss_cgi3_dist <- shortest_tss_in_cgi$cgi_e-shortest_tss_in_cgi$gene_s
shortest_tss_in_cgi_m$cgi5_tss_dist <- shortest_tss_in_cgi_m$cgi_e-shortest_tss_in_cgi_m$gene_e
shortest_tss_in_cgi_m$tss_cgi3_dist <- shortest_tss_in_cgi_m$gene_e-shortest_tss_in_cgi_m$cgi_s
###########################
## Joining back together ##
###########################
refseq_tab_short <- rbind(shortest_tss_in_cgi, shortest_tss_in_cgi_m)
refseq_tab_longenough <- refseq_tab_short %>% filter(tss_cgi3_dist>=150) #that's actually a lot more genes than I expected, 13076.
write.table(
  refseq_tab_longenough,
  file="refseq_cgi_intersect_short_TSSes.hg38.txt",
  row.names=FALSE,
  col.names=TRUE,
  quote=FALSE,
  sep="\t"
)
