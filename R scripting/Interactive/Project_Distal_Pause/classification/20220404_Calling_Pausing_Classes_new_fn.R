#########################
## 1. Load in packages ##
#########################
library(Rsamtools)
library (tidyverse)
library(compiler)
library(data.table)
library(GenomicAlignments)
library(reshape2)
enableJIT(3)
##################
## 2. Functions ##
##################

get_score_for_classes <- function(bed, bam, b=NA, a=NA, bs=10, pairedEnd=TRUE){
  if (pairedEnd==TRUE) {
    bam_info=BamFile(bam,asMates=TRUE)
  } else {bam_info=BamFile(bam)}
  si=seqinfo(bam_info)
  read_counts= ifelse(
    pairedEnd==TRUE,
    countBam(
      bam, param = ScanBamParam(flag=scanBamFlag(isFirstMateRead = TRUE))
    )$records,
    countBam(
      bam)$records
  )
  chrs = list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
              "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
              "chr19", "chr20", "chr21", "chr22", "chrX")
  tmp_list=list()
  for (i in seq_along(chrs)) {
    bed_sub = bed %>% filter(bed[[1]]==chrs[[i]])
      ref_point = ifelse(bed_sub[[6]]=="+", bed_sub[[2]], bed_sub[[3]])
      hist_rows=bed_sub[[4]]
      bins = ifelse(
        bed_sub[6]=="+", list(seq(from=b, to=(a-bs), by=bs)),
        list(seq(from=a, to=(b+bs), by=-bs)))
      hist=round(do.call(rbind, bins)) #collapse bins vector into matrix
      bin_names = seq(from=b,to=a-bs,by=bs) #names need to match so can only use +
      hist=hist+ref_point
      dimnames(hist)=list(hist_rows, bin_names)
    hist[hist <= 0] = NA
    print(paste(chrs[[i]], Sys.time()))
    sbp = ScanBamParam(
      which = GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start = 0, end = si@seqlengths[si@seqnames == chrs[[i]]])))
    bam_aln = readGAlignments(bam, param = sbp)
    long_hist=reshape2::melt(hist, na.rm=TRUE) #longform lets us generate 'all
    #bin starts' and 'all bin ends' vectors for scoring.
      olap = countOverlaps(GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start = long_hist[[3]], end = long_hist[[3]]+ bs)), bam_aln, ignore.strand=TRUE)
    long_hist$value=olap
#    print(head(long_hist))
    hist=acast(long_hist, Var1~Var2) #get the hist back into wideform.
    tmp_list[[i]] = hist
  }
  hist = do.call(rbind, tmp_list) #collapse into one matrix
  hist = hist * 1e6 / read_counts #normalize
  hist=as.data.table(hist, keep.rownames=TRUE)
  return(hist)
}

#################
## 3. Data
bamDirectory="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/bams/hg38_dedup/"
refseq_curated_longest_cgi_minus <- read.delim("~/refseq_curated_longest_cgi_minus.hg38.txt")
refseq_curated_longest_cgi_plus <- read.delim("~/refseq_curated_longest_cgi_plus.hg38.txt")

#Lets say for first attempt at this again we'll do 0-150 surrounding TSS, and -100-50 surrounding 3' end, so
#we need 250 between start and end.

rclc_p = refseq_curated_longest_cgi_plus %>% filter(cpg_e - gene_s >=250) #start at 7567, puts us at 6110.
rclc_m = refseq_curated_longest_cgi_minus %>% filter(gene_e - cpg_s >=250) #start at 7389, puts us at 5982

TSS_CpG_end_plus= rclc_p %>%
  select(1,5,3,8,2,7) #using CpG start as score placeholder
CpG_end_plus= rclc_p %>% select(1,3,6,8,2,7)
TSS_CpG_end_plus$name=paste0(TSS_CpG_end_plus$name, "_", TSS_CpG_end_plus$cpg_e)
CpG_end_plus$name=paste0(CpG_end_plus$name, "_", CpG_end_plus$cpg_e)

TSS_CpG_end_minus= rclc_m %>% select(1,2,6,8,3,7) #cgi end placeholder
CpG_end_minus= rclc_m %>% select(1,5,2,8,3,7) 
TSS_CpG_end_minus$name=paste0(TSS_CpG_end_minus$name, "_", TSS_CpG_end_minus$cpg_s)
CpG_end_minus$name=paste0(CpG_end_minus$name, "_", CpG_end_minus$cpg_s)

nascentDirectory="/Users/kayle/Box/Vertinolab/McKayla Ford/"

plus_TSS <- get_score_for_classes(bed=TSS_CpG_end_plus, bam=paste0(
  nascentDirectory, "MCF7_H9_plus.bam"), b=0, a=150, bs=150, pairedEnd = FALSE) #lets see if this works
#nice thats super fast
plus_cgi <- get_score_for_classes(bed=CpG_end_plus, bam=paste0(
  nascentDirectory, "MCF7_H9_plus.bam"), b=-100, a=50, bs=150, pairedEnd = FALSE) 

#there is a small handful of weird ones in silent that have an exactly equal number of tags between the two,
#thats not 0. See LINC01560, DHX30, HHAT, FARS2 (MYB just looks messy in Danko data, clearly prox in 293T data
#, and it's not the only one) They look interesting in IGV too.

#local antisense signal? lets say we'll check IDK 100 bp upstream of pause region
plus_TSS_anti <- get_score_for_classes(bed=TSS_CpG_end_plus, bam=paste0(
  nascentDirectory, "MCF7_H9_minus.bam"), b=-100, a=50, bs=150, pairedEnd = FALSE)
plus_cgi_anti <- get_score_for_classes(bed=CpG_end_plus, bam=paste0(
  nascentDirectory, "MCF7_H9_minus.bam"), b=-200, a=-50, bs=150, pairedEnd = FALSE) #will this even work
#yep all I need is an anchor, it doesn't actually need to touch it.
#now we'll add extra bams in.
plus_TSS_C <- get_score_for_classes(bed=TSS_CpG_end_plus, bam=paste0(
  nascentDirectory, "MCF7_C11_plus.bam"), b=0, a=150, bs=150, pairedEnd = FALSE)
plus_cgi_C <- get_score_for_classes(bed=CpG_end_plus, bam=paste0(
  nascentDirectory, "MCF7_C11_plus.bam"), b=-100, a=50, bs=150, pairedEnd = FALSE) 
plus_TSS_G <- get_score_for_classes(bed=TSS_CpG_end_plus, bam=paste0(
  nascentDirectory, "MCF7_G11_plus.bam"), b=0, a=150, bs=150, pairedEnd = FALSE)
plus_cgi_G <- get_score_for_classes(bed=CpG_end_plus, bam=paste0(
  nascentDirectory, "MCF7_G11_plus.bam"), b=-100, a=50, bs=150, pairedEnd = FALSE) 
plus_TSS_B <- get_score_for_classes(bed=TSS_CpG_end_plus, bam=paste0(
  nascentDirectory, "MCF7_B7_plus.bam"), b=0, a=150, bs=150, pairedEnd = FALSE)
plus_cgi_B <- get_score_for_classes(bed=CpG_end_plus, bam=paste0(
  nascentDirectory, "MCF7_B7_plus.bam"), b=-100, a=50, bs=150, pairedEnd = FALSE) 
plus_sense_TSS=cbind(plus_TSS, plus_TSS_C[[2]], plus_TSS_B[[2]], plus_TSS_G[[2]])
plus_sense_cgi=cbind(plus_cgi, plus_cgi_C[[2]], plus_cgi_B[[2]], plus_cgi_G[[2]])
plus_sense_TSS[[6]] = plus_sense_TSS[[2]]+ plus_sense_TSS[[3]]+ plus_sense_TSS[[4]]+ plus_sense_TSS[[5]]
plus_sense_cgi[[6]] = plus_sense_cgi[[2]]+ plus_sense_cgi[[3]]+ plus_sense_cgi[[4]]+ plus_sense_cgi[[5]]


plus=data.frame(plus_sense_TSS[[1]], plus_sense_TSS[[6]], plus_sense_cgi[[6]])
colnames(plus)=c("gene_id", "tss_score", "cgi_end_score")
plus_dist= plus %>%
  filter(cgi_end_score > tss_score)
plus_prox= plus %>%
  filter(cgi_end_score < tss_score)
plus_sil=plus %>%
  filter(cgi_end_score == tss_score)

#RPL13?
  plus_dist_true = plus_dist %>% filter(tss_score*4>cgi_end_score & cgi_end_score>10)
  plus_prox_true = plus_prox %>% filter(tss_score>cgi_end_score & tss_score>10)
########################################################################333
minus_TSS <- get_score_for_classes(bed=TSS_CpG_end_minus, bam=paste0(
  nascentDirectory, "MCF7_H9_minus.bam"), b=0, a=150, bs=150, pairedEnd = FALSE) #lets see if this works
#nice thats super fast
minus_cgi <- get_score_for_classes(bed=CpG_end_minus, bam=paste0(
  nascentDirectory, "MCF7_H9_minus.bam"), b=-100, a=50, bs=150, pairedEnd = FALSE) 
