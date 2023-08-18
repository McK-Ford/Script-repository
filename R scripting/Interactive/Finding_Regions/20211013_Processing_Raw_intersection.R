###################
## Initial setup ##
###################
library(dplyr)
CpG_Genes_Intersection <- read.delim( #Reading in the CpG/Genes intersection file, from bedtools
  "gencode_cpg_intersection.txt",
  header=FALSE,
  quote="")
no_alt_chroms <- CpG_Genes_Intersection %>%
  filter(!grepl('chrY|([\\w_]+)alt|random|fix|v1', V1)) %>%
  filter(!grepl('chrY|([\\w_]+)alt|random|fix|v1', V4))
names(no_alt_chroms) <- c("cpgchrom", "cgi_s", "cgi_e", "cpg_id",
                          "chrom", "gene_s", "gene_e",
                          "ensembl_ID", "score", "strand",
                          "thickStart", "thickEnd", "reserved",
                          "blockCount", "blockSizes",
                          "chromStarts", "UCSC_gene_IDs",
                          "cdsStartStat", "cdsEndStat",
                          "exonFrames", "type", "geneName", 
                          "UniProt_ID", "gene_type", "class",
                          "source", "transtype", "tag", "level",
                          "tier")
##################################
## Processing plus strand genes ##
##################################
plus_strand_only <- no_alt_chroms %>%
  filter(strand=="+") #Filters to only plus strand genes
plus_TSS_intersect_only <- plus_strand_only %>%
  filter(cgi_s<gene_s) #As bedtools intersect intersects anywhere in gene, this is 'start of CpG island is
#before start of gene' which ensures the TSS is within the gene.
longest_tss_in_cgi <- plus_TSS_intersect_only %>% group_by(cgi_s) %>% mutate(tss=min(gene_s)) %>% filter(tss==gene_s) #drops it to 11067
plus_distinct <- longest_tss_in_cgi %>% select(-c(cpgchrom, thickStart, thickEnd, reserved, blockCount, blockSizes, cdsStartStat, cdsEndStat, type, source, tag, level, tier, gene_type))
only_coding <- plus_distinct %>% filter(class == "coding")
not_coding <- plus_distinct %>% filter(class != "coding")
###################################
## Processing minus strand genes ##
###################################
#Majority of logic is same as with plus strand above - if not commented here, check analogous
#operation there.
minus_strand_only <- no_alt_chroms %>%
  filter(strand=="-")
minus_TSS_intersect_only <- minus_strand_only %>%
  filter(cgi_e>gene_e)
longest_tss_in_cgi_m <- minus_TSS_intersect_only %>% group_by(cgi_s) %>% mutate(tss=max(gene_e)) %>% filter(tss==gene_e) #drops it to 10838
minus_distinct <- longest_tss_in_cgi_m %>% select(-c(cpgchrom, thickStart, thickEnd, reserved, blockCount, blockSizes, cdsStartStat, cdsEndStat, type, source, tag, level, tier, gene_type))
only_coding_m <- minus_distinct %>% filter(class == "coding")
not_coding_m <- minus_distinct %>% filter(class != "coding")
###########################
## Joining back together ##
###########################
merged_table <- rbind.data.frame(
  only_coding, only_coding_m
)
merged_table_noncoding <- rbind.data.frame(
  not_coding, not_coding_m
)
write.table(
  merged_table,
  file="Gencode_CGI_genes_intersection_table.hg38.txt",
  row.names=FALSE,
  col.names=TRUE,
  quote=FALSE,
  sep="\t"
)
write.table(
  merged_table_noncoding,
  file="Gencode_CGI_genes_noncoding_intersection_table.hg38.txt",
  row.names=FALSE,
  col.names=TRUE,
  quote=FALSE,
  sep="\t"
)

gencode_table <- read.delim( #Reading in the CpG/Genes intersection file, from bedtools
  "Gencode_CGI_genes_intersection_table.hg38.txt",
  quote="")
gencode_table$gene_len <- gencode_table$gene_e - gencode_table$gene_s
gencode_table$cgi_len <- gencode_table$cgi_e - gencode_table$cgi_s
gencode_table$p_gene_in_cgi <- (gencode_table$gene_len - pmax((gencode_table$cgi_s - gencode_table$gene_s),0) - pmax((gencode_table$gene_e-gencode_table$cgi_e),0))/gencode_table$gene_len

gencode_table_p <- gencode_table %>% filter(strand=="+")
gencode_table_p$cgi5_tss_dist <- gencode_table_p$gene_s-gencode_table_p$cgi_s
gencode_table_p$tss_cgi3_dist <- gencode_table_p$cgi_e-gencode_table_p$gene_s

gencode_table_m <- gencode_table %>% filter(strand=="-")
gencode_table_m$cgi5_tss_dist <- gencode_table_m$cgi_e-gencode_table_m$gene_e
gencode_table_m$tss_cgi3_dist <- gencode_table_m$gene_e-gencode_table_m$cgi_s

gencode_table <- rbind(gencode_table_m, gencode_table_p)

library(data.table)
setDT(gencode_table) 
tmp1 = gencode_table[, head(.SD, 1), by=.(gene_s, gene_e), ]
write.table(
  tmp1,
  file="Gencode_CGI_genes_intersection_table.hg38.txt",
  row.names=FALSE,
  col.names=TRUE,
  quote=FALSE,
  sep="\t"
)
