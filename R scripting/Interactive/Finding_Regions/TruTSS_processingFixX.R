### Read in raw CpG/genes intersection table and do the basic processing of it. This is my standard filtering, as seen in "Processing raw intersection" file for commentary, it's just I don't save the intermediate step so I need to redo this part.
library(tidyverse)
CpG_Genes_Intersection <- read.delim(
  "gencode_cpg_intersection.txt",
  header=FALSE,
  quote="")
no_alt_chroms <- CpG_Genes_Intersection %>%
  filter(
    !grepl(
      'chrY|([\\w_]+)alt|random|fix|v1',
      V1
      )
    ) %>%
  filter(
    !grepl(
      'chrY|([\\w_]+)alt|random|fix|v1',
      V4
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
  filter(cgi_s<gene_s)
minus_strand_only <- no_alt_chroms %>%
  filter(strand=="-")
minus_TSS_intersect_only <- minus_strand_only %>%
  filter(cgi_e>gene_e)
### Get unique CpG islands by strand and gene name.
library(data.table)
setDT(plus_TSS_intersect_only)
dedup1_p = plus_TSS_intersect_only[
  , head(.SD, 1), c("cgi_s", "geneName")
  ]
#Figure out what islands are duplicated (which should have unique genes)
dupes_cgi_p = dedup1_p[
  duplicated(
    dedup1_p$cgi_s
    )|duplicated(
      dedup1_p$cgi_s,
      fromLast=TRUE
                               ),
  ]
#Figure out what genes are duplicated (which should have unique CGIs)
dupes_gene_p = dedup1_p[
  duplicated(
    dedup1_p$geneName
    )|duplicated(
      dedup1_p$geneName,
      fromLast=TRUE
      ),
  ]
#
setDT(minus_TSS_intersect_only)
dedup1_m = minus_TSS_intersect_only[
  , head(.SD, 1), c("cgi_s", "geneName")
  ]
dupes_cgi_m = dedup1_m[
  duplicated(
    dedup1_m$cgi_s
    )|duplicated(
      dedup1_m$cgi_s,
      fromLast=TRUE
      ),
  ]
dupes_gene_m = dedup1_m[
  duplicated(
    dedup1_m$geneName
    )|duplicated(
      dedup1_m$geneName,
      fromLast=TRUE
      ),
  ]
#so there are nearly 20k unique cgis that have TSSes (actually this isn't taking into consideration bidirectional, could be fewer...)
#1500 cpg islands (no, 1500 genes, 750 or less islands) with multiple genes contained, same strand. Filtered out for now, may figure out if there is anything interesting with them later.
summary(dupes_cgi_m$cgi_e - dupes_cgi_m$cgi_s)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#208.0   559.0   845.5  1378.3  1302.0 32637.0 
#so the CGI that contain multiple same strand genes are longer than normal. Maybe contain both start and end of some genes?
summary(dupes_cgi_p$cgi_e - dupes_cgi_p$cgi_s)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#202.0   578.5   860.0  1849.4  1433.0 27227.0 
#
### get deduplicated cpgs by CpGs
cgi_deduped_p = dedup1_p[
  !(
    duplicated(
      dedup1_p$cgi_s
      )|duplicated(
        dedup1_p$cgi_s,
        fromLast=TRUE
        )
    ),
  ]
cgi_deduped_m = dedup1_m[
  !(
    duplicated(dedup1_m$cgi_s
               )|duplicated(
                 dedup1_m$cgi_s,
                 fromLast=TRUE
                 )
    ),
  ]
#get a unique ID for labeling purposes
cgi_deduped_p$uniq_ID = paste0(
  cgi_deduped_p$geneName,
  cgi_deduped_p$cgi_s
  )
cgi_deduped_m$uniq_ID = paste0(
  cgi_deduped_m$geneName,
  cgi_deduped_m$cgi_s
  )
#get columns for beds and write the beds
cgi_deduped_p = cgi_deduped_p %>%
  select(
  chrom,
  cgi_s,
  cgi_e,
  uniq_ID,
  everything()
  )
cgi_deduped_m = cgi_deduped_m %>%
  select(
    chrom,
    cgi_s,
    cgi_e,
    uniq_ID,
    everything()
    )
write.table(
  cgi_deduped_p,
  "test_p.bed",
  quote=FALSE,
  row.names = FALSE,
  col.names = FALSE
  )
write.table(
  cgi_deduped_m,
  "test_m.bed",
  quote=FALSE,
  row.names = FALSE,
  col.names = FALSE
  )
################################
## After I process this in bluehive, I import the score files back into R.
################################
library(tidyverse)
library(data.table)

ptabX <- read.delim(
  "window_p_chrX_score_windowed.txt",
  header=FALSE,
  quote=""
  )
p_chrX_dedup1 <- ptabX %>%
  group_by(V2) %>%
  summarise(
    chrom=first(V1),
    end=first(V3),
    uniq_ID=first(V4),
    score=sum(V8)
    )
p_chrX_dedup2 <- p_chrX_dedup1 %>%
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
#Okay now that I've tested with plus X here's the loop for remaining chromosomes
chromlist = paste0("chr", seq(from=1, to=22))
out_list = list()
for (i in seq_along(chromlist)) {
  filename = paste0(
    "window_p_",
    chromlist[[i]],
    "_score_windowed.txt"
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
      end=first(V3),
      uniq_ID=first(V4),
      score=sum(V8)
      )
  p_chrX_dedup2 <- p_chrX_dedup1 %>%
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
  out_list[[i]] = p_chr_dedup2
}
plus_vals_0 = do.call(rbind, out_list)
plus_vals = rbind(plus_vals_0, p_chrX_dedup2, p_chr1_dedup2)

mtabX <- read.delim(
  "window_m_chrX_score_windowed.txt", header=FALSE, quote=""
  )
m_chrX_dedup1 <- mtabX %>%
  group_by(V2) %>%
  summarise(
    chrom=first(V1),
    end=first(V3),
    uniq_ID=first(V4),
    score=sum(V8)
    )
m_chrX_dedup2 <- m_chrX_dedup1 %>%
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

out_list = list()
for (i in seq_along(chromlist)) {
  filename = paste0(
    "window_m_", chromlist[[i]], "_score_windowed.txt"
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
minus_vals = rbind(minus_vals_0, m_chrX_dedup2)
minus_vals=unique(minus_vals)

#Also read in the unprocessed
test_p <- read.table(
  "~/Lab stuff/R_coding/backedup/test_p.bed", quote="", comment.char=""
  )
test_m <- read.table(
  "~/Lab stuff/R_coding/backedup/test_m.bed", quote="", comment.char=""
  )
test_p$uniq_ID = paste0(
  test_p$V4,
  test_p$V2
)
test_m$uniq_ID = paste0(
  test_m$V4,
  test_m$V2
)
test_p_x = test_p %>% filter(V1=="chrX")
test_m_x = test_m %>% filter(V1=="chrX")
#Get the processed values with their metaplots
p_with_verified = merge(test_p_x, p_chrX_dedup2, by="uniq_ID", all=TRUE)
m_with_verified = merge(test_m_x, m_chrX_dedup2, by="uniq_ID", all=TRUE)
#sort for has and doesn't has signal
pluschrx = p_with_verified %>% filter(V1=="chrX")
minuschrx =m_with_verified %>% filter(V1=="chrX") #alright good up to this point.
silentpx = pluschrx %>% filter(is.na(score))
silentmx = minuschrx %>% filter(is.na(score))
has_signal_px = pluschrx %>% filter(!is.na(score))
has_signal_mx = minuschrx %>% filter(!is.na(score))

silent_p = p_with_verified %>% filter(is.na(score))
silent_m = m_with_verified %>% filter(is.na(score))
has_signal_p = p_with_verified %>% filter(!is.na(score) & !is.na(V5))
has_signal_m = m_with_verified %>% filter(!is.na(score))
####################
col_trimmed_p = has_signal_px %>%
  select(-c(
    V5, V10, V12, V13, V14, V15, V20,
    V21, V22, V19, V24, V26, V29, V30, chrom
    ))
col_trimmed_p_silent = silentpx %>%
  select(-c(
    V5, V10, V12, V13, V14, V15, V20,
    V21, V22, V19, V24, V26, V29, V30, chrom
    ))
col_trimmed_m = has_signal_mx %>%
  select(-c(
    V5, V10, V12, V13, V14, V15, V20,
    V21, V22, V19, V24, V26, V29, V30, chrom
    ))
col_trimmed_m_silent = silentmx %>%
  select(-c(
    V5, V10, V12, V13, V14, V15, V20,
    V21, V22, V19, V24, V26, V29, V30, chrom
    ))
###########################3
colnames(col_trimmed_p) = c(
  "uniq_ID", "chrom", "cgi_s", "cgi_e", "gene_name", "cgi_ID",
  "gene_start", "gene_end", "ensembl_ID", "strand", "exon_len",
  "exon_starts", "name2", "name3", "gene_type", "transcript_type",
  "tags", "tss_s", "tss_e", "score" 
  )
colnames(col_trimmed_p_silent) =c(
  "uniq_ID", "chrom", "cgi_s", "cgi_e", "gene_name", "cgi_ID",
  "gene_start", "gene_end", "ensembl_ID", "strand", "exon_len",
  "exon_starts", "name2", "name3", "gene_type", "transcript_type",
  "tags", "tss_s", "tss_e", "score" 
  )
colnames(col_trimmed_m) =c(
  "uniq_ID", "chrom", "cgi_s", "cgi_e", "gene_name", "cgi_ID",
  "gene_start", "gene_end", "ensembl_ID", "strand", "exon_len",
  "exon_starts", "name2", "name3", "gene_type", "transcript_type",
  "tags", "tss_s", "tss_e", "score" 
  )
colnames(col_trimmed_m_silent) =c(
  "uniq_ID", "chrom", "cgi_s", "cgi_e", "gene_name", "cgi_ID",
  "gene_start", "gene_end", "ensembl_ID", "strand", "exon_len",
  "exon_starts", "name2", "name3", "gene_type", "transcript_type",
  "tags", "tss_s", "tss_e", "score" 
  )
#
col_trimmed_p$cgi_len = col_trimmed_p$cgi_e-col_trimmed_p$cgi_s
col_trimmed_p_silent$cgi_len = col_trimmed_p_silent$cgi_e-col_trimmed_p_silent$cgi_s
col_trimmed_m$cgi_len = col_trimmed_m$cgi_e-col_trimmed_m$cgi_s
col_trimmed_m_silent$cgi_len = col_trimmed_m_silent$cgi_e-col_trimmed_m_silent$cgi_s
#
col_trimmed_p$gene_len = col_trimmed_p$gene_e-col_trimmed_p$gene_s
col_trimmed_p_silent$gene_len = col_trimmed_p_silent$gene_e-col_trimmed_p_silent$gene_s
col_trimmed_m$gene_len = col_trimmed_m$gene_e-col_trimmed_m$gene_s
col_trimmed_m_silent$gene_len = col_trimmed_m_silent$gene_e-col_trimmed_m_silent$gene_s
#
col_trimmed_p$Updist = col_trimmed_p$tss_s-col_trimmed_p$cgi_s
col_trimmed_p_silent$Updist = col_trimmed_p_silent$tss_s-col_trimmed_p_silent$cgi_s
col_trimmed_m$Updist = col_trimmed_m$cgi_e-col_trimmed_m$tss_e
col_trimmed_m_silent$Updist = col_trimmed_m_silent$cgi_e-col_trimmed_m_silent$tss_e
#
col_trimmed_p$Downdist = col_trimmed_p$cgi_e-col_trimmed_p$tss_s
col_trimmed_p_silent$Downdist = col_trimmed_p_silent$cgi_e-col_trimmed_p_silent$tss_s
col_trimmed_m$Downdist = col_trimmed_m$tss_e-col_trimmed_m$cgi_s
col_trimmed_m_silent$Downdist = col_trimmed_m_silent$tss_e-col_trimmed_m_silent$cgi_s
#
#for graphing, I want the 
plus = rbind(col_trimmed_p, col_trimmed_p_silent)
minus = rbind(col_trimmed_m, col_trimmed_m_silent)
plus_filt = plus %>%
  filter(Downdist>=150) %>%
  arrange(-Downdist)
minus_filt = minus %>%
  filter(Downdist>=150) %>%
  arrange(-Downdist)

tstbed = cbind(
  minus_filt$chrom, minus_filt$cgi_s, minus_filt$tss_e,
  minus_filt$gene_name, minus_filt$score, minus_filt$strand
  )
write.table(
  tstbed, "minus_w_TSSes_testing.bed", quote=FALSE,
  row.names = FALSE, col.names = FALSE
  )

tstbed_p = cbind(
  plus_filt$chrom, plus_filt$tss_s, plus_filt$cgi_e,
  plus_filt$gene_name, plus_filt$score, plus_filt$strand
  )
write.table(
  tstbed_p, "plus_w_TSSes_testing.bed", quote=FALSE,
  row.names=FALSE, col.names=FALSE
  )

save(col_trimmed_p, col_trimmed_m, col_trimmed_p_silent,
     col_trimmed_m_silent, file="gettingTruTSSesX.RData")
