library(tidyverse)
library(data.table)

############################# works best with a unique id
cpgIsland_hg38 <- read.table(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/Regions/cpgIsland.hg38.txt")
cpgIsland_hg38$CGI_ID = paste0(
  cpgIsland_hg38$V1,
  cpgIsland_hg38$V2)
write.table(
  cpgIsland_hg38, "cpgIsland_hg38.txt", quote=FALSE,
  row.names=FALSE, col.names=FALSE, sep="\t"
)
###################################

Dir = "../Box/Vertinolab/McKayla Ford/Data/remappingStarts/getTSSes/getTSSes-selected/"

#Large files, iterate over chromosomes
chromlist = paste0("chr", seq(from=1, to=22))
chromlist = append(chromlist, "chrx")
out_list = list()
for (i in seq_along(chromlist)) {
  filename = paste0(
    Dir,
    "cpgs_starts_simple_plus_",
    chromlist[[i]],
    ".txt"
  )
  file = read.delim(
    filename,
    header=FALSE,
    quote=""
  )
  p_chr_dedup1 <- file %>%
    group_by(V2) %>%
    summarise(
      chrom=first(V1),
      start=first(V2),
      end=first(V3),
      uniq_ID=first(V4),
      score=sum(V8)
    )
  p_chr_dedup2 <- p_chr_dedup1 %>%
    ungroup() %>%
    group_by(uniq_ID) %>%
    mutate(score_max=max(score)) %>%
    filter(score==score_max) %>%
    summarize(
      chrom=first(chrom),
      start=min(start),
      end=first(end),
      uniq_ID=first(uniq_ID),
      score=first(score))
  out_list[[i]] = p_chr_dedup2
}
plus_vals_0 = do.call(rbind, out_list)

out_list = list()
for (i in seq_along(chromlist)) {
  filename = paste0(
    filename = paste0(
      Dir,
      "cpgs_starts_simple_minus_",
      chromlist[[i]],
      ".txt"
    )
  )
  file = read.delim(filename, header=FALSE, quote="")
  m_chr_dedup1 <- file %>%
    group_by(V2) %>%
    summarise(
      chrom=first(V1),
      end=first(V3),
      uniq_ID=first(V4),
      score=sum(V8)
    )
  m_chr_dedup2 <- m_chr_dedup1 %>%
    ungroup() %>%
    group_by(uniq_ID) %>%
    mutate(score_max=max(score)) %>%
    filter(score==score_max) %>%
    summarize(
      chrom=first(chrom),
      start=first(V2),
      end=max(end),
      score=first(score)
    )
  out_list[[i]] = m_chr_dedup2
}
minus_vals_0 = do.call(rbind, out_list)
##################################
### Read in raw CpG/genes intersection table
CpG_Genes_Intersection <- read.delim(
  "../Box/Vertinolab/McKayla Ford/Data/Regions/gencode_cpg_intersection.txt.gz",
  header=FALSE,
  quote="")
no_alt_chroms <- CpG_Genes_Intersection %>%
  filter(
    !grepl(
      'chrY|([\\w_]+)alt|random|fix|v1',
      V1
    )
  )
names(no_alt_chroms) <- c(
  "cpgchrom", "cgi_s", "cgi_e", "cpg_id",
  "chrom", "gene_s", "gene_e",
  "ensembl_ID", "score", "strand",
  "thickStart", "thickEnd", "reserved",
  "blockCount", "blockSizes",
  "chromStarts", "UCSC_gene_IDs",
  "cdsStartStat", "cdsEndStat",
  "exonFrames", "type", "geneName", 
  "UniProt_ID", "gene_type", "class",
  "source", "transtype", "tag", "level",
  "tier"
)


plus_strand_only <- no_alt_chroms %>%
  filter(strand=="+") #Filters to only plus strand genes
plus_TSS_intersect_only <- plus_strand_only %>%
  filter(cgi_s<gene_s & class!="pseudo")
minus_strand_only <- no_alt_chroms %>%
  filter(strand=="-")
minus_TSS_intersect_only <- minus_strand_only %>%
  filter(cgi_e>gene_e & class!="pseudo")
### Get unique CpG islands by strand and gene name.
save.image("~/tsscalling.RData")
#############################
plus_TSS_intersect_only$CGI_ID = paste0(
  plus_TSS_intersect_only$cpgchrom,
    plus_TSS_intersect_only$cgi_s)
plus_TSS_intersect_only$num_class[
  plus_TSS_intersect_only$class=="coding"
  ] = 0
plus_TSS_intersect_only$num_class[
  plus_TSS_intersect_only$class=="nonCoding"
] = 1
coding_preferred_p <- plus_TSS_intersect_only %>%
  ungroup() %>%
  group_by(CGI_ID) %>%
  mutate(trueclass=min(num_class)) %>%
  filter(trueclass==num_class)
deduped_CGIs_p <- coding_preferred_p %>%
  ungroup() %>%
  group_by(CGI_ID) %>%
  mutate(refstart=min(gene_s)) %>%
  filter(refstart==gene_s) %>%
  distinct(gene_s, .keep_all=TRUE) #what are the ones that don't intersect? CGIs with no genes coding or noncoding?
###################### what else, okay next olap.
minus_TSS_intersect_only$CGI_ID = paste0(
  minus_TSS_intersect_only$cpgchrom,
  minus_TSS_intersect_only$cgi_s)
minus_TSS_intersect_only$num_class[
  minus_TSS_intersect_only$class=="coding"
] = 0
minus_TSS_intersect_only$num_class[
  minus_TSS_intersect_only$class=="nonCoding"
] = 1
coding_preferred_m <- minus_TSS_intersect_only %>%
  ungroup() %>%
  group_by(CGI_ID) %>%
  mutate(trueclass=min(num_class)) %>%
  filter(trueclass==num_class)
deduped_CGIs_m <- coding_preferred_m %>%
  ungroup() %>%
  group_by(CGI_ID) %>%
  mutate(refstart=max(gene_s)) %>%
  filter(refstart==gene_s) %>%
  distinct(gene_s, .keep_all=TRUE)
plus_vals_0$CGI_ID = plus_vals_0$uniq_ID
minus_vals_0$CGI_ID = minus_vals_0$uniq_ID
p_merged = merge(deduped_CGIs_p, plus_vals_0, by="CGI_ID", all=TRUE)
m_merged = merge(deduped_CGIs_m, minus_vals_0, by="CGI_ID", all=TRUE)
#################################################################
p_genes_no_signal = p_merged %>% filter(is.na(score.y)) #1140
p_CGI_no_gene = p_merged %>% filter(is.na(gene_s)) #9456 ... thats a lot but I suppose the numbers do add up
p_genes_w_signal = p_merged %>% filter(!is.na(gene_s) & !is.na(score.y)) #8010
m_genes_no_signal = m_merged %>% filter(is.na(score.y)) #1067
m_CGI_no_gene = m_merged %>% filter(is.na(gene_s)) #9436
m_genes_w_signal = m_merged %>% filter(!is.na(gene_s) & !is.na(score.y)) #7852
#also worth keeping in mind, shared CGIs, bidirectional may not be antisense annotated but may have lots antisense signal

#okay what columns are actually necessary? Also do I want a 2-way version of the table?
# CGI ID, one chrom column (cpgchrom), cgi_s, cgi_e, gene_s, gene_e, ensembl_id, strand if we combine, gene_name,
#unit_prot_id, class, transtype, "start" and score.y.
p_merged_select = p_merged %>% select(c(1:4, 7:9, 11, 23, 24,26, 28, 37, 39))
m_merged_select = m_merged %>% select(c(1:4, 7:9, 11, 23, 24,26, 28, 38, 39))

colnames(p_merged_select) = c("CGI_ID", "genechrom", "cgi_s", "cgi_e", "gene_s", "gene_e", "ensembl_ID", "strand",
                              "geneName", "unitprot_ID", "class", "transtype", "netcagetss", "tssscore")
colnames(m_merged_select) = c("CGI_ID", "genechrom", "cgi_s", "cgi_e", "gene_s", "gene_e", "ensembl_ID", "strand",
                              "geneName", "unitprot_ID", "class", "transtype", "netcagetss", "tssscore")
p_merged_select$strand="+"
m_merged_select$strand="-"

Merged_table = rbind(p_merged_select, m_merged_select)
has_signal = Merged_table %>% filter(!is.na(tssscore))

wide_has_signal = has_signal %>% pivot_wider(names_from = strand, values_from = c(
  gene_s, gene_e, ensembl_ID, geneName, unitprot_ID, class, transtype, netcagetss, tssscore))
nogenes = wide_has_signal %>% filter(is.na(`gene_s_+`) & is.na(`gene_s_-`)) #15K, probably filter slightly based on strength
#then set aside as mystery transcribed CGIs.
bidir = wide_has_signal %>% filter(!is.na(`gene_s_+`) & !is.na(`gene_s_-`)) #3k
#though some are noncoding. Point is annotated bidir.
#honestly lets just assign some TSSes. Don't worry about antisense rn.

p_genes_no_signal = p_merged_select %>% filter(is.na(tssscore)) #1140
p_CGI_no_gene = p_merged_select %>% filter(is.na(gene_s)) #9456 ... thats a lot but I suppose the numbers do add up
p_genes_w_signal = p_merged_select %>% filter(!is.na(gene_s) & !is.na(tssscore)) #8010
m_genes_no_signal = m_merged_select %>% filter(is.na(tssscore)) #1067
m_CGI_no_gene = m_merged_select %>% filter(is.na(gene_s)) #9436
m_genes_w_signal = m_merged_select %>% filter(!is.na(gene_s) & !is.na(tssscore)) #7852

#genes with no signal get thrown out.
summary(p_CGI_no_gene$tssscore)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#1.00     1.00     3.00    36.18    15.00 14384.00 
summary(m_CGI_no_gene$tssscore) #lets filter out anything less than 3 reads for
#the enhancer ones
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#1.00     1.00     3.00    40.06    16.00 14030.00 
summary(p_genes_w_signal$tssscore)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#1.0     39.0    281.0   1004.4    891.5 846371.0 
summary(m_genes_w_signal$tssscore)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.0    36.0   279.0   826.2   856.0 84694.0 
#lets do 3 for all then.
p_CGI_no_gene_tab = p_CGI_no_gene %>% select(-c(2:12)) %>% filter(tssscore>=3)
m_CGI_no_gene_tab = m_CGI_no_gene %>% select(-c(2:12)) %>% filter(tssscore>=3)
p_CGI_withgenes_tab = p_genes_w_signal %>% filter(tssscore>=3)
m_CGI_withgenes_tab = m_genes_w_signal %>% filter(tssscore>=3)

write.table(
  p_CGI_no_gene_tab, "nonanno_tx_p.txt", quote=FALSE,
  row.names = FALSE, col.names = FALSE, sep = "\t"
)
write.table(
  m_CGI_no_gene_tab, "nonannno_tx_m.txt", quote=FALSE,
  row.names=FALSE, col.names=FALSE, sep = "\t"
)
write.table(
  p_CGI_withgenes_tab, "tx_p.txt", quote=FALSE,
  row.names = FALSE, col.names = FALSE, sep = "\t"
)
write.table(
  m_CGI_withgenes_tab, "tx_m.txt", quote=FALSE,
  row.names=FALSE, col.names=FALSE, sep = "\t"
)
