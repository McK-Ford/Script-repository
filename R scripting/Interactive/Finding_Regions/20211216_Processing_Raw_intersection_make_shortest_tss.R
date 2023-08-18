###################
## Initial setup ##
###################
library(dplyr)

ncbiRefSeq.sorted.txt <- read.delim("C:/Users/kayle/Downloads/ncbiRefSeq.sorted.txt.gz", header=FALSE)
ncbiRefSeq <- ncbiRefSeq.sorted.txt %>% select(-c(1,2,7,8,9,10,11,14,15,16)) %>% unique
no_alt_chroms <- ncbiRefSeq %>% filter(!grepl('chrY|([\\w_]+)alt|random|fix|v1', V3))

plus_strand_only <- no_alt_chroms %>%
  filter(V4=="+")
plus_strand_one_tes <- plus_strand_only %>% group_by(V13, V5) %>% summarise(chrom=first(V3), V6=max(V6), strand=first(V4), score=first(V12))
plus_strand_ordered <- cbind(plus_strand_one_tes$chrom, plus_strand_one_tes$V5, plus_strand_one_tes$V6, plus_strand_one_tes$V13, plus_strand_one_tes$score, plus_strand_one_tes$strand)

minus_strand_only <- no_alt_chroms %>%
  filter(V4=="-")
minus_strand_one_tes <- minus_strand_only %>% group_by(V13, V6) %>% summarise(chrom=first(V3), V5=min(V5), strand=first(V4), score=first(V12))
minus_strand_ordered <- cbind(minus_strand_one_tes$chrom, minus_strand_one_tes$V5, minus_strand_one_tes$V6, minus_strand_one_tes$V13, minus_strand_one_tes$score, minus_strand_one_tes$strand)

rejoined_for_intersection = rbind(plus_strand_ordered, minus_strand_ordered)
write.table(
  rejoined_for_intersection,
  file="refseq_igv_genes_for_intersection.hg38.txt",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)

######################################
#after intersection
refseq_igv_cpg_intersect.hg38 <- read.delim("~/R_work/refseq_igv_cpg_intersect.hg38.txt", header=FALSE)
#split by strand
plus = refseq_igv_cpg_intersect.hg38 %>% filter(V6=="+")
minus = refseq_igv_cpg_intersect.hg38 %>% filter(V6=="-")
plus_TSS_intersect_only <- plus %>%
  filter(V8<V2)
minus_TSS_intersect_only <- minus %>%
  filter(V9>V3)
shortest_tss_in_cgi_plus <- plus_TSS_intersect_only %>% group_by(V8) %>% mutate(tss=max(V2)) %>% filter(tss==V2) #
#8321 genes
shortest_tss_in_cgi_minus <- minus_TSS_intersect_only %>% group_by(V8) %>% mutate(tss=min(V3)) %>% filter(tss==V3)
#8167
##############
## Metadata ##
##############
shortest_tss_in_cgi_plus$gene_len <- shortest_tss_in_cgi_plus$V3 - shortest_tss_in_cgi_plus$V2
shortest_tss_in_cgi_minus$gene_len <- shortest_tss_in_cgi_minus$V3 - shortest_tss_in_cgi_minus$V2
shortest_tss_in_cgi_plus$cgi_len <- shortest_tss_in_cgi_plus$V9 - shortest_tss_in_cgi_plus$V8
shortest_tss_in_cgi_minus$cgi_len <- shortest_tss_in_cgi_minus$V9 - shortest_tss_in_cgi_minus$V8
shortest_tss_in_cgi_plus$cgi5_tss_dist <- shortest_tss_in_cgi_plus$V2-shortest_tss_in_cgi_plus$V8
shortest_tss_in_cgi_plus$tss_cgi3_dist <- shortest_tss_in_cgi_plus$V9-shortest_tss_in_cgi_plus$V2
shortest_tss_in_cgi_minus$cgi5_tss_dist <- shortest_tss_in_cgi_minus$V9-shortest_tss_in_cgi_minus$V3
shortest_tss_in_cgi_minus$tss_cgi3_dist <- shortest_tss_in_cgi_minus$V3-shortest_tss_in_cgi_minus$V8
###########################
## Joining back together ##
###########################
shortest_tss_in_cgi <- rbind(shortest_tss_in_cgi_minus, shortest_tss_in_cgi_plus)
shortest_tss_in_cgi_150 <- shortest_tss_in_cgi %>% filter(tss_cgi3_dist>=150)
dim(shortest_tss_in_cgi_150)
write.table(
  shortest_tss_in_cgi_150,
  file="filtered_shortTSS_cpgend_igv_refseq.txt",
  row.names=FALSE,
  col.names=TRUE,
  quote=FALSE,
  sep="\t"
)
