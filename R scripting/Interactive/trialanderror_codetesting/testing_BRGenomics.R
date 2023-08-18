###################################################
### code chunk number 1: load packages
###################################################
library(groHMM)
library(tidyverse)
library(data.table)
#################################################
### Pausing index
#################################################
T47 <- BamFile("../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/T47D_R1.BAM")
si=seqinfo(T47)
sbp = ScanBamParam(
  which = GRanges(seqnames = "chr1", ranges = IRanges(
    start = 0, end = si@seqlengths[si@seqnames == "chr1" ])))
#pretty sure G ranges is just completely ignoring the reverse complement argument
#not helpful guys.
bam_aln_m = readGAlignments(T47, param = sbp)
bam_aln_p = bam_aln_m[strand(bam_aln_m) == "-"]
vecp = rep_len("+", length(bam_aln_p))
rvecp = Rle(vecp)
strand(bam_aln_p) = rvecp

bam_aln_m2 = bam_aln_m[strand(bam_aln_m) == "+"]
vecp = rep_len("-", length(bam_aln_m2))
rvecp = Rle(vecp)
strand(bam_aln_m2) = rvecp

glst <- unlist(GRangesList(bam_aln_m2, bam_aln_p))

#############################
gncde <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/Regions/Untouched source data/Gencode39_knowngene.hg38.txt.gz",
  header=TRUE
  )
#remove columns score (5), thick start -chrom starts (7:12), cds start stat - type (14:17),
#gene name 2 (19), gene type(20), source(22), tag-tier (24-26),
gncde2 <- gncde %>% select(-c(5, 7:12, 14:17, 19, 20, 22, 24:26)) #unless Ching-Hua used detect transcripts?
gncde3p <- gncde2 %>%
  group_by(geneName) %>%
  mutate(start_site=min(chromStart)) %>%
  filter(chromStart==start_site & strand == "+")
gncde3m <- gncde2 %>%
  group_by(geneName) %>%
  mutate(start_site=max(chromEnd)) %>%
  filter(chromEnd==start_site & strand == "-")

gncde4=rbind(gncde3p, gncde3m) ##66409, but includes lots of pseudogenes and such.
#for testing:
gncde4.2 = gncde4 %>% filter(X.chrom=="chr1")

############################################
features <- GRanges(seqnames=gncde4.2$X.chrom, IRanges(gncde4.2$chromStart,gncde4.2$chromEnd), strand=gncde4.2$strand)
pi <- pausingIndex(features, glst)
pi

test <- cbind(gncde4.2, pi)

############################################
#happy with this now, now to run it higher thruput, easy to base off my metaplotter
#function, we'll say gencode processing thru gncde4 is outside of function

#CLEAN CODE
get_pI <- function(gen_tab, bam, debug=FALSE){
  bam_info=BamFile(bam)
  si=seqinfo(bam_info)
  chrs = list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
              "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
              "chr19", "chr20", "chr21", "chr22", "chrX")
  if (debug==TRUE){
    chrs = "chr1"
    print("Running in debug mode")
  }
  tmp_list=list()
  for (i in seq_along(chrs)) {
    print(paste(chrs[[i]], Sys.time()))
    gen_sub = gen_tab %>% filter(chrom==chrs[[i]])
    sbp = ScanBamParam(
      which = GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start = 0, end = si@seqlengths[si@seqnames == chrs[[i]]])))
    bam_aln = readGAlignments(bam_info, param = sbp)
    
    bam_aln_p = bam_aln[strand(bam_aln) == "-"]
    vecp = Rle(rep_len("+", length(bam_aln_p)))
    strand(bam_aln_p) = vecp
    
    bam_aln_m = bam_aln[strand(bam_aln) == "+"]
    vecp = Rle(rep_len("-", length(bam_aln_m)))
    strand(bam_aln_m) = vecp
    
    glst <- unlist(GRangesList(bam_aln_m, bam_aln_p))
    features <- GRanges(seqnames=gen_sub$chrom, IRanges(gen_sub$start,gen_sub$end), strand=gen_sub$strand)
    pi <- pausingIndex(features, glst)
    
    tmp_list[[i]] = cbind(gen_sub, pi)
  }
  pi_tab = do.call(rbind, tmp_list) #collapse into one matrix
  pi_tab=as.data.table(pi_tab, keep.rownames=TRUE)
  return(pi_tab)
}

#establishing gen tab
gncde <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/Regions/Untouched source data/Gencode39_knowngene.hg38.txt.gz",
  header=TRUE
)
gncde2 <- gncde %>% select(-c(5, 7:12, 14:17, 19, 20, 22, 24:26)) #unless Ching-Hua used detect transcripts?
gncde3p <- gncde2 %>%
  group_by(geneName) %>%
  mutate(start_site=min(chromStart)) %>%
  filter(chromStart==start_site & strand == "+") %>%
  select(-start_site)
gncde3m <- gncde2 %>%
  group_by(geneName) %>%
  mutate(start_site=max(chromEnd)) %>%
  filter(chromEnd==start_site & strand == "-") %>%
  select(-start_site)

gncde4=rbind(gncde3p, gncde3m) ##66409, but includes lots of pseudogenes and such.
colnames(gncde4) = c("chrom", "start", "end", "ID1", "strand",
                    "ID2", "symbol", "class", "type")
T47 <- get_pI(gen_tab=gncde4,
              bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/T47D_R1.BAM")
MCF7 <- get_pI(gen_tab=gncde4,
              bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/MCF7_R1.BAM")
MCF10A <- get_pI(gen_tab=gncde4,
              bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/MCF10A_R1.BAM")
SUM159 <- get_pI(gen_tab=gncde4,
              bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/SUM159_R1.BAM")
MB231 <- get_pI(gen_tab=gncde4,
              bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/MB231_R1.BAM")

pI_tab <- cbind(T47[,2:20], MCF7[,11:20], MCF10A[,11:20], SUM159[,11:20], MB231[,11:20])
colnames(pI_tab) = c(colnames(pI_tab[,1:9]),
                     paste0(colnames(pI_tab[,10:19]), "_", "T47D"),
                     paste0(colnames(pI_tab[,10:19]), "_", "MCF7"),
                     paste0(colnames(pI_tab[,10:19]), "_", "MCF10A"),
                     paste0(colnames(pI_tab[,10:19]), "_", "SUM159"),
                     paste0(colnames(pI_tab[,10:19]), "_", "MB231"))

write.table(pI_tab, file="cell_lines_pI_no_filter.txt", quote=FALSE, sep="\t", row.names = FALSE)
################################################### testing the new one
library(BRGenomics)

gncde <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/Regions/Untouched source data/Gencode39_knowngene.hg38.txt.gz",
  header=TRUE
)
gncde2 <- gncde %>% select(-c(5, 7:12, 14:17, 19, 20, 22, 24:26)) #unless Ching-Hua used detect transcripts?
gncde3p <- gncde2 %>%
  group_by(geneName) %>%
  mutate(start_site=min(chromStart)) %>%
  filter(chromStart==start_site & strand == "+") %>%
  select(-start_site)
gncde3m <- gncde2 %>%
  group_by(geneName) %>%
  mutate(start_site=max(chromEnd)) %>%
  filter(chromEnd==start_site & strand == "-") %>%
  select(-start_site)
gncde4=rbind(gncde3p, gncde3m) ##66409, but includes lots of pseudogenes and such.
colnames(gncde4) = c("chrom", "start", "end", "ID1", "strand",
                     "ID2", "symbol", "class", "type")
gncde5 = gncde4 %>% filter((end-start)>200) #breaks getpausingindices if don't filter bc genebodies takes out too small regions
features <- GRanges(seqnames=gncde5$chrom, IRanges(gncde5$start,gncde5$end), strand=gncde5$strand)
pr = promoters(features, upstream = 50, downstream = 200)
gb = genebodies(features, start=200, end=0)
t47d = import_bam("../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/T47D_R1.BAM",
                  revcomp = TRUE, trim.to="3p") #I feel like my method is faster and less memory-intensive.
  #Significanty faster, and my computer had to call garbage collection 2 times running this to clear memory.
  #maybe review their code and implement a modified version of this?
T47D = getPausingIndices(t47d, pr, gb)
#oh ffs they need the random etc chroms cleaned up.
t47d = tidyChromosomes(t47d) #at least thats fast
T47D = getPausingIndices(t47d, pr, gb) #oh it's fast bc the function that took a solid five minutes and two gc() calls failed
#to actually import the data. WTF. Maybe it overwrote it, some code doesn't play nice with rewriting its parent file.
# Let's just run a test of this my way getting granges w/ chr1.

T47 <- BamFile("../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/T47D_R1.BAM")
si=seqinfo(T47)
sbp = ScanBamParam(
  which = GRanges(seqnames = "chr1", ranges = IRanges(
    start = 0, end = si@seqlengths[si@seqnames == "chr1" ])))
#pretty sure G ranges is just completely ignoring the reverse complement argument
#not helpful guys.
bam_aln_m = readGAlignments(T47, param = sbp)
bam_aln_p = bam_aln_m[strand(bam_aln_m) == "-"]
vecp = rep_len("+", length(bam_aln_p))
rvecp = Rle(vecp)
strand(bam_aln_p) = rvecp

bam_aln_m2 = bam_aln_m[strand(bam_aln_m) == "+"]
vecp = rep_len("-", length(bam_aln_m2))
rvecp = Rle(vecp)
strand(bam_aln_m2) = rvecp

glst <- unlist(GRangesList(bam_aln_m2, bam_aln_p))
gr <- resize(glst, width = 1L, fix = "end") #based on their script there's no funny business with stranding and resize
#third times the charm?
gr2 = tidyChromosomes(gr)
T47D = getPausingIndices(gr2, pr, gb) #ffs you've got to be kidding me. oh i think i see
gr2$score = 1
T47D = getPausingIndices(gr2, pr, gb)
#yep I was right. have to add score column else doesn't work.
summary(T47D) #remember only chrom 1
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.00    0.00    1.74     Inf   11.46     Inf   55858 

#okay compare to groHMM take
pi <- pausingIndex(features, glst)
pi2 <- pausingIndex(features, gr)

summary(pi$Pause)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.0000  0.0000  0.0389  0.0000 50.3200 
summary(pi2$Pause)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.00000  0.00000  0.00000  0.03497  0.00000 48.94000 

#less nans if i do it this way, easier to compare
gncde5 = gncde4 %>% filter((end-start)>200 & chrom=="chr1")
features <- GRanges(seqnames=gncde5$chrom, IRanges(gncde5$start,gncde5$end), strand=gncde5$strand)
pr = promoters(features, upstream = 50, downstream = 200)
gb = genebodies(features, start=200, end=0)
T47D = getPausingIndices(gr2, pr, gb) #actual function is fast at least, the bams are the optimization problem
summary(T47D) #remember only chrom 1
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.000   1.743     Inf  11.461     Inf    1374 
#well these are pretty different, how does this even have inf values, I can see nan but not inf
T47D2 = T47D[ T47D != "Inf"]
summary(T47D2)
#okay why are these so different
test = cbind(pi, pi2, T47D, gncde5)
test2 = test[,-c(5,6,9,10,14:16,19,20)]
#you get INF in new method which I should be able to rework easily enough
#when the body counts are 0 (still not sure how GroHMM managed to calculate anything from body=0 tbh)
#so
test3 = test2 %>% filter(BodyCounts>0)
#ah I see, grohmm doesn't even really give the pause index it gives normalized pause counts and normalized body counts
#that's not made clear at all in the docs..
test3$pI1 = test3$Pause/test3$Body
test3$pI2 = test3$Pause.1/test3$Body.1

test3$p1vp2 = test3$pI2/test3$pI1 #so grohmm really only has small differences between bp 1 and full read
test3$nvp1 = test3$T47D/test3$pI1 #lots of dif here though, whole region vs very small region
######################## Alright, what is the best script I can make for BRGenomics PI########################################################
library(tidyverse)
library(data.table)
library(BRGenomics)
library(Rsamtools)
library(GenomicAlignments)

genes <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/Regions/Untouched source data/Gencode39_knowngene.hg38.txt.gz",
  header=TRUE
)
genes2 <- genes %>% select(-c(5, 7:12, 14:17, 19, 20, 22, 24:26)) #unless Ching-Hua used detect transcripts?
genes3p <- genes2 %>%
  group_by(geneName) %>%
  mutate(start_site=min(chromStart)) %>%
  filter(chromStart==start_site & strand == "+") %>%
  select(-start_site)
genes3m <- genes2 %>%
  group_by(geneName) %>%
  mutate(start_site=max(chromEnd)) %>%
  filter(chromEnd==start_site & strand == "-") %>%
  select(-start_site)
genes4=rbind(genes3p, genes3m) ##66409, but includes lots of pseudogenes and such.
colnames(genes4) = c("chrom", "start", "end", "ID1", "strand",
                     "ID2", "symbol", "class", "type")
genes4$uniq_id = paste0(genes4$start, genes4$symbol)

get_pI <- function(gen_tab, bam, pause_reg_before = 50,
                   pause_reg_after=200, debug=FALSE){
  #need table with chrom, start, end, strand, and uniq id.
  #also need a read1 only bam file.
  #requires dplyr and data.table
  print("establishing regions")
  size_filt_genes = gen_tab %>% filter((end-start)>pause_reg_after)
  print("referencing bam")
  bam_info=BamFile(bam)
  si=seqinfo(bam_info)
  chrs = list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
              "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
              "chr19", "chr20", "chr21", "chr22", "chrX")
  if (debug==TRUE){
    chrs = list("chr20", "chr21")
    print("Running in debug mode")
  }
  tmp_list=list()
  for (i in seq_along(chrs)) {
    print(paste(chrs[[i]], Sys.time()))
    print(paste("getting regions for chrom ", chrs[[i]]))
    chrom_filt_genes = size_filt_genes %>% filter(chrom==chrs[[i]])
    features = GRanges(seqnames = chrom_filt_genes$chrom,
                       IRanges(chrom_filt_genes$start, chrom_filt_genes$end),
                       strand=chrom_filt_genes$strand,
                       uniq_id=chrom_filt_genes$uniq_id)
    pr = promoters(features, upstream = pause_reg_before,
                   downstream = pause_reg_after)
    gb = genebodies(features, start=pause_reg_after, end=0)
    print("getting chrom-specific bam reads")
    sbp = ScanBamParam(
      which = GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start = 0, end = si@seqlengths[si@seqnames == chrs[[i]]])))
    bam_aln = readGAlignments(bam_info, param = sbp)
    print("Reverse complementing bams")
    bam_aln_p = bam_aln[strand(bam_aln) == "-"]
    vecp = Rle(rep_len("+", length(bam_aln_p)))
    strand(bam_aln_p) = vecp
    bam_aln_m = bam_aln[strand(bam_aln) == "+"]
    vecp = Rle(rep_len("-", length(bam_aln_m)))
    strand(bam_aln_m) = vecp
    gr <- unlist(GRangesList(bam_aln_m, bam_aln_p))
    print("resizing and tidying")
    gr2 <- resize(gr, width = 1L, fix = "end")
    gr3 = tidyChromosomes(gr2)
    gr3$score = 1
    pI = getPausingIndices(gr3, pr, gb)
    chrom_filt_genes$pI = pI
    tmp_list[[i]] = chrom_filt_genes
  }
  pi_tab = do.call(rbind, tmp_list) #collapse into one matrix
  pi_tab=as.data.table(pi_tab)
  return(pi_tab)
}


T47 <- get_pI(gen_tab=genes4,
              bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/T47D_R1.BAM")
MCF7 <- get_pI(gen_tab=genes4,
               bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/MCF7_R1.BAM")
MCF10A <- get_pI(gen_tab=genes4,
                 bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/MCF10A_R1.BAM")
SUM159 <- get_pI(gen_tab=genes4,
                 bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/SUM159_R1.BAM")
MB231 <- get_pI(gen_tab=genes4,
                bam="../Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/Split_bams/MB231_R1.BAM")

pI_tab <- cbind(T47, MCF7[,11], MCF10A[,11], SUM159[,11], MB231[,11])
colnames(pI_tab) = c(colnames(pI_tab[,1:10]), "T47D", "MCF7", "MCF10A", "SUM159", "MB231")

write.table(pI_tab, file="cell_lines_pI_2_no_filter.txt", quote=FALSE, sep="\t", row.names = FALSE)

###################################################
#### Processing
#first of all, what has no signal under any condition?
cell_lines_pI_no_filter <- read.delim("~/cell_lines_pI_2_no_filter.txt")
pI_tab = cell_lines_pI_2_no_filter
pI_tab_has_sig = pI_tab %>%
  filter(
    (T47D != "Inf" & T47D != "NaN") | (MCF7 != "Inf" & MCF7 != "NaN") |
    (MCF10A != "Inf" & MCF10A != "NaN") | (SUM159 != "Inf" & SUM159 != "NaN") |
      (MB231 != "Inf" & MB231 != "NaN")
  ) #drops size of df from 56k to 46k.
#still have some dups with identical TSSes too...
#Pick 1, prefer coding over noncoding/pseudo
#same way I did it for CGI TSSes
pI_tab_has_sig$num_class[
  pI_tab_has_sig$class=="coding"
] = 0
pI_tab_has_sig$num_class[
  pI_tab_has_sig$class=="nonCoding"
] = 1
pI_tab_has_sig$num_class[
  pI_tab_has_sig$class=="pseudo"
] = 2

pI_dedup <- pI_tab_has_sig %>%
  ungroup() %>%
  group_by(symbol) %>%
  mutate(trueclass=min(num_class)) %>%
  filter(trueclass==num_class) %>%
  ungroup() %>%
  distinct(symbol, .keep_all=TRUE) %>%
  select(-c(16,17))
#that's down to 42k.

lncRNA <- pI_dedup %>%
  filter(type=="lncRNA") #15k
coding <- pI_dedup %>%
  filter(type=="protein_coding") #15k, remaining 12k must be pseudogenes
###process coding first
summary(coding$T47D)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   0.000   4.414     Inf  14.548     Inf     778 
#1/0 = Inf
#0/0 = NaN
#0/1 = 0
#so as I thought, NaN and Inf are both useless to us.
#I feel like I should get rid of anything that's all only NaN, Inf, and basically 0...
#if I print out to 20 digits, first quartile is still 0. How small is 'basically 0?'
#well based on playing around with quantile, quantile(coding$T47D, 0.29879, na.rm=TRUE)
#the 29.879% is the first time I have any sig figs instead of printing 0, and that says 0.0006.
#probably pretty safe to say anything less than 0.0001 is effectively 0 for purpose of filtering.
pI_tab_has_sig = pI_tab %>%
  filter(
    (T47D != "Inf" & T47D != "NaN" & T47D > 0.0001) | (MCF7 != "Inf" & MCF7 != "NaN" & MCF7 > 0.0001) |
      (MCF10A != "Inf" & MCF10A != "NaN" & MCF10A > 0.0001) | (SUM159 != "Inf" & SUM159 != "NaN" & SUM159 > 0.0001) |
      (MB231 != "Inf" & MB231 != "NaN" & MB231 > 0.0001)
  )
#30k
pI_tab_has_sig$num_class[
  pI_tab_has_sig$class=="coding"
] = 0
pI_tab_has_sig$num_class[
  pI_tab_has_sig$class=="nonCoding"
] = 1
pI_tab_has_sig$num_class[
  pI_tab_has_sig$class=="pseudo"
] = 2

pI_dedup <- pI_tab_has_sig %>%
  ungroup() %>%
  group_by(symbol) %>%
  mutate(trueclass=min(num_class)) %>%
  filter(trueclass==num_class) %>%
  ungroup() %>%
  distinct(symbol, .keep_all=TRUE) %>%
  select(-c(16,17))
#that's down to 27k.

lncRNA <- pI_dedup %>%
  filter(type=="lncRNA") #8k
coding <- pI_dedup %>%
  filter(type=="protein_coding") #11.5k, remaining 7.5k must be pseudogenes
###process coding first
summary(coding$T47D)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   1.710   7.443     Inf  18.090     Inf     181 
#yeah a lot of the Na's were bc one cell line had maybe a 0.00000000000000000000001 signal somewhere,
#which isn't very meaningful.
#big downside of this one is it doesn't actually give us the num body reads in output...
#basically the things that are most trustworthy are the middle section.
