######################
## 1. Preprocessing ##
######################
library(tidyverse)
Refseq_curated <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/Regions/Untouched source data/Refseq_curated.hg38.txt.gz"
  )
#probably best to check a thing or two first, we don't actually have 90k genes.
Refseq_curated <- Refseq_curated %>% distinct(txStart, txEnd, .keep_all = TRUE) 
#puts us at 47043, much better starting point.
Refseq_curated2 <- Refseq_curated %>% select(c(3,5,6, 1, 4, 13)) %>% unique()
write.table(
  Refseq_curated2,
  file="Refseq_curated_tss4intersect.bed",
  col.names=FALSE,
  row.names=FALSE,
  quote=FALSE,
  sep="\t"
)
################
## 2. Postprocessing
library(tidyverse)
Refseq_curated_cgis_raw <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Refseq_curated_cgis_raw.txt",
  header=FALSE)
no_alt_chroms <- Refseq_curated_cgis_raw %>% filter(!grepl("chrY|([\\w_]+)alt|random|fix|v1", V1))
plus_strand_only <- no_alt_chroms %>%
  filter(V9=="+") %>%
  filter(V6>V2) %>%
  distinct(V6, .keep_all = TRUE) 
#Now we've eliminated variable TESes and made sure we only have CGIs that
#intersect the TSS. Now to trim by picking one TSS per. For now, lets just go
#with longest.
plus_one_tss <- plus_strand_only %>%
  group_by(V1, V2, V3, V4) %>%
  dplyr::summarize(geneS=min(V6),
                   geneE=dplyr::first(V7),
                   strand=dplyr::first(V9),
                   name=dplyr::first(V10)) %>%
  ungroup()
#I think I accidentally got rid of the NM accessions earlier. Oh well it will be fine.
#that's 7567 regions.
colnames(plus_one_tss)=c("chrom", "cpg_s", "cpg_e", "cpg_id",
                         "gene_s", "gene_e", "strand", "name")


minus_strand_only <- no_alt_chroms %>%
  filter(V9=="-") %>%
  filter(V3>V7)  %>%
  distinct(V7, .keep_all = TRUE)
minus_one_tss <- minus_strand_only %>%
  group_by(V1, V2, V3, V4) %>%
  dplyr::summarize(geneS=dplyr::first(V6),
                   geneE=max(V7),
                   strand=dplyr::first(V9),
                   name=dplyr::first(V10)) %>%
  ungroup()
colnames(minus_one_tss)=c("chrom", "cpg_s", "cpg_e", "cpg_id",
                         "gene_s", "gene_e", "strand", "name")


write.table(
  plus_one_tss,
  file="refseq_curated_longest_cgi_plus.hg38.txt",
  row.names=FALSE,
  col.names=TRUE,
  quote=FALSE,
  sep="\t"
)
write.table(
  minus_one_tss,
  file="refseq_curated_longest_cgi_minus.hg38.txt",
  row.names=FALSE,
  col.names=TRUE,
  quote=FALSE,
  sep="\t"
)
