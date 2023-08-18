###################
## Initial setup ##
###################
library(dplyr)
library(stringr)

make_bed_file <- function(tab, b_name) {
  write.table(
    tab,
    file=b_name,
    row.names=FALSE,
    col.names=FALSE,
    quote=FALSE,
    sep="\t"
  )
}

CpG_Genes_Intersection <- read.delim( #Reading in the CpG/Genes intersection file, from bedtools
  "refseq_curated_genes_cpg_intersect.txt",
  header=FALSE,
  quote="")
no_alt_chroms <- CpG_Genes_Intersection %>%
  filter(!grepl('chrY|([\\w_]+)alt|random|fix|v1', V1))

get_info <- function(strand, tab) {
  stranded = tab %>% filter(V10==strand)
  if (strand=="+"){
      TSS_intersect_only = stranded %>%
      filter(V2<V6)
  }
  if (strand=="-"){
    TSS_intersect_only = stranded %>%
      filter(V3>V7)
  }
  uni = unique(TSS_intersect_only)
  dis = uni %>%
    group_by(V2) %>%
    summarize (cpg_end=max(V3), TSS=min(V6), TES=max(V7))
  dis$strand = rep_len(c(strand), dim(dis)[1])
  dis$score = rep_len(c(0), dim(dis)[1])
  dis = dis %>% rename("cpg_start"="V2")
  lt = distinct(uni, V2, .keep_all=TRUE)
  
  gene_names = lt$V8
  names(gene_names) = lt$V2
  dis$name = unname((gene_names[as.character(dis$cpg_start)]))
  chrs = lt$V1
  names(chrs) = lt$V2
  dis$chr = unname(chrs[as.character(dis$cpg_start)])
 
   rm_cpg_multis <- distinct(dis, cpg_end, .keep_all=TRUE)
  return (rm_cpg_multis)
}

get_stranded_TSS_Cpg_end_beds <- function(tab, prefix) {
  if (tab$strand[1] == "+") {
    t1 = tab %>% select(chr, TSS, cpg_end, name, score, strand)
    print(head(t1))
    make_bed_file(t1, (str_glue(prefix, "_", "plus.bed", sep="")))
  }
  if (tab$strand[1] == "-") {
    t1 = tab %>% select(chr, cpg_start, TES, name, score, strand)
    make_bed_file(t1, (str_glue(prefix, "_", "minus.bed", sep="")))
  }
  colnames(t1) = c('chr', 'start', 'end', 'name', 'score', 'strand')
  return (t1)
}

get_unidirectional <- function(tab1, tab2) {
 uni1 = anti_join(tab1, tab2, by="cpg_start")
 uni2 = anti_join(tab2, tab1, by="cpg_start")
 uni1_format = get_stranded_TSS_Cpg_end_beds(uni1, "unidir")
 uni2_format = get_stranded_TSS_Cpg_end_beds(uni2, "unidir")
 uni = rbind(uni1_format, uni2_format)
 make_bed_file(uni, "bidirectional islands/unidir.bed")
 }

get_bidirectional <- function(tab1, tab2) {
  bidir1 = semi_join(tab1, tab2, by="cpg_start")
  bidir2 = semi_join(tab2, tab1, by="cpg_start")
  bidir1_format = get_stranded_TSS_Cpg_end_beds(bidir1, "bidir")
  bidir2_format = get_stranded_TSS_Cpg_end_beds(bidir2, "bidir")
  bidir = rbind (bidir1_format, bidir2_format)
  make_bed_file(bidir, "bidirectional islands/bidir.bed")
  #close_far
  bidr1_sort = bidir1 %>% arrange(cpg_start)
  bidr2_sort = bidir2 %>% arrange(cpg_start)
  wide_bidir = cbind(bidr1_sort, bidr2_sort[3:5], bidr2_sort[7])
  colnames(wide_bidir) = c('cs', 'ce', 'TSS1', 'TES1', 'str1', 'scr', 'n1', 'chr',
                           'TSS2', 'TES2', 'str2', 'n2')
  wide_bidir$diff = wide_bidir$TSS1-wide_bidir$TES2
  bound_lower = quantile(wide_bidir$diff, 0.01)
  bound_upper = quantile(wide_bidir$diff, 0.99)
  rm_outliers = wide_bidir %>% filter (diff >= bound_lower & diff <= bound_upper)
  get_close = rm_outliers %>% filter (diff <=100)
  get_far = rm_outliers %>% filter (diff >=100)
  get_plusminus(get_close, "bidirectional islands/close_bidir")
  get_plusminus(get_far, "bidirectional islands/far_bidir") #probably don't need, most are far.
  }

get_plusminus <- function(tab, prefix) {
  tp = tab %>% select(chr, TSS1, ce, n1, scr, str1)
  make_bed_file(tp, (str_glue(prefix, "_", "plus.bed", sep="")))
  colnames(tp) = c('chr', 'start', 'end', 'name', 'score', 'strand')
  tm = tab %>% select(chr, cs, TES2, n2, scr, str2)
  make_bed_file(tm, (str_glue(prefix, "_", "minus.bed", sep="")))
  colnames(tm) = c('chr', 'start', 'end', 'name', 'score', 'strand')
  tabx <- rbind(tp, tm)
  make_bed_file(tabx, (str_glue(prefix, ".bed", sep="")))
}

plus <- get_info("+", no_alt_chroms)
minus <- get_info("-", no_alt_chroms)
get_unidirectional(plus, minus)
get_bidirectional(plus, minus)
