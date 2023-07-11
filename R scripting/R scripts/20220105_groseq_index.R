#make a groseq/proseq index

#########################
## 1. Load in packages ##
#########################
library(Rsamtools)
library (tidyverse)
library (rtracklayer)
library(compiler)
library(BSgenome) #make sure Hg38 is locally installed - makes functions much faster for only .7 gB.
library(data.table)
enableJIT(3)
#########################
## 2. Define functions ##
#########################
get_bw_by_chrom <- function(bw, chrom){
  partial_bw = import(bw, selection=GenomicSelection("hg38", chrom=chrom, colnames="score"))
  return(partial_bw)
}

get_bed_by_chrom_n_strand <- function(bed, chr, strd){
  partial_bed = bed %>% filter(chrom==chr) %>% filter(strand==strd)
  return(partial_bed)
}

get_anchored_scores <- function(bed, bw, anch, b, a, strand, bs=10){
  #  INPUT is a table (stored in txt or bed file, must contain a numerical coordinate column to use as anchors, a chromosome column (named in standard UCSC format ie chr1, chr2 etc), a strand column containing only + and -, and a gene name column. May need to rename columns as seen in the intersection table readin.), a bw file, how many bp you want before and after the anchor, the column of the table that your anchor is stored in, the method of combining bw bin scores (valid options are mean and sum), the size of the bins (default 10).
  #OUTPUT: A table with cgi and gene coordinate columns, chromosomes and gene names, bin numbers, and scores for each bin.
  #  1. Handle each chromosome separately (for speed and memory purposes), writing each processed chromosome to a list, after retrieving the bigwig data by chromosome using Rtracklayer's import bigwig, and retrieve the table data by chromosome and strand.
  chrs = list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
              "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
              "chr19", "chr20", "chr21", "chr22", "chrX")
  tmp_list = list()
  for (i in seq_along(chrs)) { #iterate over chroms
    #get our subsetted starting material
    partial_bw = get_bw_by_chrom(bw, chrs[[i]])
    partial_bed = get_bed_by_chrom_n_strand(bed, chrs[[i]], strand)
    #2. Retrieve desired start coordinates of windows based on our anchors. If plus strand, this is (anchor - the number of base pairs desired before anchor). If minus strand, this is (anchor - the number of base pairs desired after anchor) (since directionality is flipped.) 
    if (strand=="+") {
      anchored_window = partial_bed %>% mutate(
        win_start=(partial_bed[,anch]-b))
    }
    if (strand=="-") {
      anchored_window = partial_bed %>% mutate(
        win_start=(partial_bed[,anch]-a))
    }
    #3. Split window into bins. 1st bin for ea window is win_start to win_start + bs * bin_num (1 for first bin), 2nd is win_start + bs*(bin_num-1) to win_start + bs*(bin_num), etc. So bin coordinates can be determined from bin number. Get last_bin_start from window_width - bs, then generate sequence of all bin starts from 0 to last_bin_start in bs sized steps. Put in string to add to table (table doesn't play nicely w/ it in other formats), then break string apart into bin columns. Pivot into longform table. Calc end coordinates by adding bin width to start coordinates, -1 due to inclusive ranges.
    last_bin_start = b+a-bs
    start_bin_string = paste0(seq(from=0, to=last_bin_start, by=bs), collapse=",") 
    anch_win_w_bins = cbind(anchored_window, start_bin_string)
    num_bins = (a+b)/bs
    anch_win_w_bins_as_cols = separate(
      anch_win_w_bins, col = start_bin_string, into=paste(
        "X", seq(1:num_bins), sep=""
      ), sep=",", remove=TRUE, convert=TRUE
    )
    end_col = paste("X", num_bins, sep="")
    long_anch = pivot_longer(
      anch_win_w_bins_as_cols, cols=X1:end_col, names_to="Bins", names_prefix="X", values_to = "start_adjust"
    ) #bins are 5' to 3' - take that into account when handling minus down the line
    long_anch$new_start = long_anch$win_start+long_anch$start_adjust
    long_anch$new_end = long_anch$new_start + bs -1 
    #4. Convert bed file into GRanges, find overlaps using GRanges 'findOverlaps' command. Returns Hits object with Query and Subject Hits, (just fancy table w/ indices of all overlaps). Pull overlaps into table that you can bind together to have a bigwig score + regions table. 
    anchored_granges = GRanges (
      seqnames = long_anch$chrom,
      ranges=IRanges(
        start=long_anch$new_start,
        end=long_anch$new_end
      ),
      strand=long_anch$strand
    )
    #adding metadata to the granges. If I find a cleaner way to do this I will fix the function but for now this works.
    anchored_granges$gennames <- long_anch$gene_name
    anchored_granges$gen_5 <- long_anch$gene_l
    anchored_granges$gen_3 <- long_anch$gene_h
    anchored_granges$cgi_5 <- long_anch$cpg_l
    anchored_granges$cgi_3 <- long_anch$cpg_h
    anchored_granges$bins <- long_anch$Bins
    olaps = findOverlaps(anchored_granges, partial_bw)
    bw_df = data.frame(partial_bw[subjectHits(olaps)])
    regions_df = data.frame(anchored_granges[queryHits(olaps)])
    regions_w_raw_scores = cbind(bw_df, regions_df)
    names(regions_w_raw_scores) <- c(
      "chrom", "bs", "be", "bw", "star", "score", "chr2",
      "gs", "ge", "gw", "strand", "gn", "g5", "g3", "c5", "c3", "bins"
    )
    #5. bw bins don't match perfectly w/ desired window bins. May be larger (common for 0 score bins- many programs merge them together) or smaller, or just not line up cleanly. Therefore calculate an 'adjustor' (percentage of bw bin in window bin), then multiply score by adjustor, in data.table format for performance purposes. Then we get the sum of adjusted scores for each bin in each gene. Add the table to the list.
    regions_w_raw_scores$adjustor <- (
      regions_w_raw_scores$bw - pmax(
        (regions_w_raw_scores$gs-regions_w_raw_scores$bs),0
      ) - pmax(
        (regions_w_raw_scores$be-regions_w_raw_scores$ge),0)
    )/regions_w_raw_scores$bw
    regions_w_raw_scores$adjusted_score <- regions_w_raw_scores$score * regions_w_raw_scores$adjustor
    setDT(regions_w_raw_scores) 
    tmp1 = regions_w_raw_scores[, head(.SD, 1), by=.(bins, g5), .SDcols=c("chrom", "gs", "ge", "strand", "gn", "g3", "c5", "c3")]
    tmp2 = regions_w_raw_scores[, list((sum(adjusted_score)/sum(adjustor))), by=.(bins, g5)]
    summarized_regions_w_raw_scores = cbind(tmp1, tmp2)
    tmp_list[[i]] <- summarized_regions_w_raw_scores
  }
  olap_df <- do.call(rbind, tmp_list) #bind the list of tables together into a single table.
  olap_df$bins<-as.numeric(olap_df$bins)
  names(olap_df) = c("bins", "gene_5", "chrom", "bin_5", "bin_3", "strand", "gene_name", "gene_3", "cgi_5", "cgi_3", "bins2", "g5", "score")
  return(olap_df)
}

get_index <- function(tab){
  ntab = tab %>%
    group_by(gene_5) %>%
    summarise(gene_name=first(gene_name), strand=first(strand), gene_3=first(gene_3), cgi_5=first(cgi_5), cgi_3=first(cgi_3), score=mean(score)) #do I actually want this to be mean? Yes.
}

write_bed_file <- function(tab, b_name) {
  write.table(
    tab,
    file=b_name,
    row.names=FALSE,
    quote=FALSE,
    sep="\t"
  )
}

intersection_table <- read.delim( #Reading in the CpG/Genes intersection file, from bedtools
  "filtered_shortTSS_cpgend_igv_refseq.txt",
  quote="")
names(intersection_table) <- c("chrom", "gene_l", "gene_h", "gene_name", "score_filler", "strand", "chrom2", "cpg_l", "cpg_h", "cgi_id", "tss", "gene_len", "cgi_len", "cgi5_tss_dist", "tss_cgi3_dist")
#gene_name
uniq_id <- paste0(intersection_table$cpg_l, intersection_table$gene_name)
intersection_table$tru_name <- intersection_table$gene_name
intersection_table$gene_name <- uniq_id
##################
## 3. Get index ##
##################
###1. B7
plus_tss <- get_anchored_scores(bed=intersection_table, bw="MCF7_B7_fwd.bw", anch=2, b=0, a=100, strand="+")
minus_tss <- get_anchored_scores(bed=intersection_table, bw="MCF7_B7_rev.bw", anch=3, b=0, a=100, strand="-")
plus_cpg <- get_anchored_scores(bed=intersection_table, bw="MCF7_B7_fwd.bw", anch=9, b=50, a=50, strand="+")
minus_cpg <- get_anchored_scores(bed=intersection_table, bw="MCF7_B7_rev.bw", anch=8, b=50, a=50, strand="-")

plus_tss_index <- get_index(plus_tss)
plus_cpg_index <- get_index(plus_cpg)
plus_index <- plus_tss_index
plus_index$dist <- plus_cpg_index$score
plus_index_noolap <- plus_index %>% filter((cgi_3-gene_5)>=150)

minus_tss_index <- get_index(minus_tss)
minus_cpg_index <- get_index(minus_cpg)
minus_index <- minus_tss_index
minus_index$dist <- minus_cpg_index$score
minus_index_noolap <- minus_index %>% filter((gene_3-cgi_5)>=150)

#dist?
dim(plus_index_noolap %>% filter(dist>score)) #713 of 6701
dim(minus_index_noolap %>% filter(dist>score))#733 of 6486
dim(plus_index_noolap %>% filter(dist<score))#4721
dim(minus_index_noolap %>% filter(dist<score))#4521
dim(plus_index_noolap %>% filter(dist==score)) #1237
dim(minus_index_noolap %>% filter(dist==score))#1232

B7_plus_tab = plus_index_noolap
B7_minus_tab = minus_index_noolap

###2. H9
plus_tss <- get_anchored_scores(bed=intersection_table, bw="MCF7_H9_fwd.bw", anch=2, b=0, a=100, strand="+")
minus_tss <- get_anchored_scores(bed=intersection_table, bw="MCF7_H9_rev.bw", anch=3, b=0, a=100, strand="-")
plus_cpg <- get_anchored_scores(bed=intersection_table, bw="MCF7_H9_fwd.bw", anch=9, b=50, a=50, strand="+")
minus_cpg <- get_anchored_scores(bed=intersection_table, bw="MCF7_H9_rev.bw", anch=8, b=50, a=50, strand="-")

plus_tss_index <- get_index(plus_tss)
plus_cpg_index <- get_index(plus_cpg)
plus_index <- plus_tss_index
plus_index$dist <- plus_cpg_index$score
plus_index_noolap <- plus_index %>% filter((cgi_3-gene_5)>=150)

minus_tss_index <- get_index(minus_tss)
minus_cpg_index <- get_index(minus_cpg)
minus_index <- minus_tss_index
minus_index$dist <- minus_cpg_index$score
minus_index_noolap <- minus_index %>% filter((gene_3-cgi_5)>=150)

#dist?
dim(plus_index_noolap %>% filter(dist>score)) #651 of 6701
dim(minus_index_noolap %>% filter(dist>score))#653 of 6486
dim(plus_index_noolap %>% filter(dist<score))#4635
dim(minus_index_noolap %>% filter(dist<score))#4427
dim(plus_index_noolap %>% filter(dist==score)) #1415
dim(minus_index_noolap %>% filter(dist==score))#1406

H9_plus_tab = plus_index_noolap
H9_minus_tab = minus_index_noolap

###3. C11
plus_tss <- get_anchored_scores(bed=intersection_table, bw="MCF7_C11_fwd.bw", anch=2, b=0, a=100, strand="+")
minus_tss <- get_anchored_scores(bed=intersection_table, bw="MCF7_C11_rev.bw", anch=3, b=0, a=100, strand="-")
plus_cpg <- get_anchored_scores(bed=intersection_table, bw="MCF7_C11_fwd.bw", anch=9, b=50, a=50, strand="+")
minus_cpg <- get_anchored_scores(bed=intersection_table, bw="MCF7_C11_rev.bw", anch=8, b=50, a=50, strand="-")

plus_tss_index <- get_index(plus_tss)
plus_cpg_index <- get_index(plus_cpg)
plus_index <- plus_tss_index
plus_index$dist <- plus_cpg_index$score
plus_index_noolap <- plus_index %>% filter((cgi_3-gene_5)>=150)

minus_tss_index <- get_index(minus_tss)
minus_cpg_index <- get_index(minus_cpg)
minus_index <- minus_tss_index
minus_index$dist <- minus_cpg_index$score
minus_index_noolap <- minus_index %>% filter((gene_3-cgi_5)>=150)

#dist?
dim(plus_index_noolap %>% filter(dist>score)) #698 of 6701
dim(minus_index_noolap %>% filter(dist>score))#717 of 6486
dim(plus_index_noolap %>% filter(dist<score))#4673
dim(minus_index_noolap %>% filter(dist<score))#4462
dim(plus_index_noolap %>% filter(dist==score)) #1330
dim(minus_index_noolap %>% filter(dist==score))#1307

C11_plus_tab = plus_index_noolap
C11_minus_tab = minus_index_noolap

###4. G11
plus_tss <- get_anchored_scores(bed=intersection_table, bw="MCF7_G11_fwd.bw", anch=2, b=0, a=100, strand="+")
minus_tss <- get_anchored_scores(bed=intersection_table, bw="MCF7_G11_rev.bw", anch=3, b=0, a=100, strand="-")
plus_cpg <- get_anchored_scores(bed=intersection_table, bw="MCF7_G11_fwd.bw", anch=9, b=50, a=50, strand="+")
minus_cpg <- get_anchored_scores(bed=intersection_table, bw="MCF7_G11_rev.bw", anch=8, b=50, a=50, strand="-")

plus_tss_index <- get_index(plus_tss)
plus_cpg_index <- get_index(plus_cpg)
plus_index <- plus_tss_index
plus_index$dist <- plus_cpg_index$score
plus_index_noolap <- plus_index %>% filter((cgi_3-gene_5)>=150)

minus_tss_index <- get_index(minus_tss)
minus_cpg_index <- get_index(minus_cpg)
minus_index <- minus_tss_index
minus_index$dist <- minus_cpg_index$score
minus_index_noolap <- minus_index %>% filter((gene_3-cgi_5)>=150)

#dist?
dim(plus_index_noolap %>% filter(dist>score)) #635 of 6701
dim(minus_index_noolap %>% filter(dist>score))#611 of 6486
dim(plus_index_noolap %>% filter(dist<score))#4827
dim(minus_index_noolap %>% filter(dist<score))#4649
dim(plus_index_noolap %>% filter(dist==score)) #1239
dim(minus_index_noolap %>% filter(dist==score))#1226

G11_plus_tab = plus_index_noolap
G11_minus_tab = minus_index_noolap

##########################################
Int_tab_w_scores = intersection_table

B7 = rbind(B7_minus_tab, B7_plus_tab)
H9 = rbind(H9_minus_tab, H9_plus_tab)
C11 = rbind(C11_minus_tab, C11_plus_tab)
G11= rbind(G11_minus_tab, G11_plus_tab)

B7_score = B7$score
H9_score = H9$score
C11_score= C11$score
G11_score= G11$score

B7_dist = B7$dist
H9_dist = H9$dist
C11_dist= C11$dist
G11_dist= G11$dist

names(B7_score) = B7$gene_name
names(H9_score) = H9$gene_name
names(C11_score)= C11$gene_name
names(G11_score)= G11$gene_name

names(B7_dist) = B7$gene_name
names(H9_dist) = H9$gene_name
names(C11_dist)= C11$gene_name
names(G11_dist)= G11$gene_name

Int_tab_w_scores$B7_score = unname(B7_score[as.character(Int_tab_w_scores$gene_name)])
Int_tab_w_scores$H9_score = unname(H9_score[as.character(Int_tab_w_scores$gene_name)])
Int_tab_w_scores$C11_score= unname(C11_score[as.character(Int_tab_w_scores$gene_name)])
Int_tab_w_scores$G11_score= unname(G11_score[as.character(Int_tab_w_scores$gene_name)])

Int_tab_w_scores$B7_dist = unname(B7_dist[as.character(Int_tab_w_scores$gene_name)])
Int_tab_w_scores$H9_dist = unname(H9_dist[as.character(Int_tab_w_scores$gene_name)])
Int_tab_w_scores$C11_dist= unname(C11_dist[as.character(Int_tab_w_scores$gene_name)])
Int_tab_w_scores$G11_dist= unname(G11_dist[as.character(Int_tab_w_scores$gene_name)])
###############################################
Int_tab_w_scores$B7_score_dif = Int_tab_w_scores$B7_score - Int_tab_w_scores$B7_dist
Int_tab_w_scores$H9_score_dif = Int_tab_w_scores$H9_score - Int_tab_w_scores$H9_dist
Int_tab_w_scores$C11_score_dif = Int_tab_w_scores$C11_score - Int_tab_w_scores$C11_dist
Int_tab_w_scores$G11_score_dif = Int_tab_w_scores$G11_score - Int_tab_w_scores$G11_dist

No_na = Int_tab_w_scores %>% filter(!is.na(Int_tab_w_scores$B7_dist)) #catches all the other NAs too.
write.table(
  Int_tab_w_scores,
  file="danko_mcf7_subclones_proseq_short.txt",
  row.names=FALSE,
  quote=FALSE,
  sep="\t")

