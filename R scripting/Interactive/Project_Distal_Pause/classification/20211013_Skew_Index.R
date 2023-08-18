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
  #  INPUT is a table (stored in txt or bed file, must contain a numerical coordinate column to use as anchors, a chromosome column (named in standard UCSC format ie chr1, chr2 etc), a strand column containing only + and -, a gene name column, and preferably columns with the 5' and 3' ends of the genes and associated cpg islands, though if desired the function could likely be easily edited to remove this requirement. May need to rename columns as seen in the intersection table readin.), a bw file, how many bp you want before and after the anchor, the column of the table that your anchor is stored in, the method of combining bw bin scores (valid options are mean and sum), the size of the bins, the expected size of the bigwig bins.
  #OUTPUT: A table with cgi and gene coordinate columns, chromosomes and gene names, bin numbers, and scores for each bin.
  #  1. Because it makes the function less memory and time intensive to run on my personal machine (should be quicker on big things now too), we break everything down into small sub-problems. In this case, we handle each chromosome separately, writing each processed chromosome to a list (as expanding a list is significantly less memory and time intensive than expanding a dataframe). So we retrieve the bigwig data by chromosome using Rtracklayer's import bigwig, and retrieve the table data by chromosome and strand.
  
  chrs = list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
              "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
              "chr19", "chr20", "chr21", "chr22", "chrX")
  tmp_list = list()
  for (i in seq_along(chrs)) { #iterate over chroms
    #get our subsetted starting material
    partial_bw = get_bw_by_chrom(bw, chrs[[i]])
    partial_bed = get_bed_by_chrom_n_strand(bed, chrs[[i]], strand)
    #2. We retrieve the desired start coordinates of our windows based on our anchors. If plus strand, this is the anchor - the number of base pairs before the anchor desired. If minus strand, this is the anchor - the number of base pairs after the anchor desired (since directionality is flipped.) 
    if (strand=="+") {
      anchored_window = partial_bed %>% mutate(
        win_start=(partial_bed[,anch]-b))
    }
    if (strand=="-") {
      anchored_window = partial_bed %>% mutate(
        win_start=(partial_bed[,anch]-a))
    }
    #3. Then we split this window into bins. Logically speaking, if bins are for example 10 bp wide, the first bin for each window should be from window_start to window_start + 10, the next would be window_start+10 to window_start + 20, etc. Therefore, we retrieve the amount the start of the final bin for each region needs to be adjusted by (last_bin_start) by taking region width - bin size, then generate a sequence of bin starts from 0 to this last bin start by the bin size. We put this into a string so we can add it as a column to our windowed table (because the table wouldn't play nicely with it if I made it a vector), then break the string apart into however many columns we have bins. Then we can pivot it into longform data so if we have say 10 bins we have 10 rows of region data numbered bin 1-10 with associated start coordinate data. Finally, we can calculate the end coordinates by adding the bin width to the start coordinates (minus 1 bc Granges uses both ends inclusive ranges so if you don't subtract the ranges will be one bp bigger than desired.)
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
    #4. Convert the bed file into a GRanges object. Then find overlaps by using the GRanges 'findOverlaps' command. This will return a Hits object with Query and Subject Hits, which is really just a fancy table that gives you indices of all overlaps, so you can pull them out into a table that you can bind together to have a bigwig + regions table. 
    anchored_granges = GRanges (
      seqnames = long_anch$chrom,
      ranges=IRanges(
        start=long_anch$new_start,
        end=long_anch$new_end
      ),
      strand=long_anch$strand
    )
    #adding metadata to the granges. If I find a cleaner way to do this I will fix the function accordingly.
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
    #5. However, the bigwig bins will not overlap cleanly with the desired bins. They may be larger (common for 0 score bins- unfortunately many programs merge them together to save storage space) or smaller, or just not line up cleanly. Therefore we calculate an 'adjustor' which is meant to be the amount of the bw bin included in the desired bin (width of bw bin - how far bw bin start is before desired bin start (0 if negative) - how far bw bin end is after desired bin end (0 if negative) with the whole thing divided by the bw bin width. Example - bw bin width is 5. Bin of interest is lets say 100-110. Bw bins overlapping are 98-103, 103-108, and 108-113. Adjustors are therefore (5-2-0)/5 = 0.6, (5-0-0)/5=1, (5-0-3)/5 = 0.4. Then we multiply the score by the adjustor to get the adjusted score that takes into account what percentage of the bigwig bin actually overlaps the desired bin. Now, this is a big table that's difficult for the computer to work with, so we put it into data.table format. We retrieve the first row name, strand, chrom, coordinates data for every gene and desired bin, because that should be the same for every bigwig bin, and put those into a table (so instead of 3 rows of bin 1 of gene FUBAR, one for each bigwig bin, like we'd have above with our 100-110 example we only have one row of bin 1 of gene FUBAR.) Then, again based on the original table, we combine the adjusted scores for each bin of each gene. We bind these scores as a column back to the table with only one of each gene-bin combo and bind this to the list.
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
  olap_df <- do.call(rbind, tmp_list)
  olap_df$bins<-as.numeric(olap_df$bins)
  names(olap_df) = c("bins", "gene_5", "chrom", "bin_5", "bin_3", "strand", "gene_name", "gene_3", "cgi_5", "cgi_3", "bins2", "g5", "score")
  return(olap_df)
}

get_index <- function(tab){
  ntab = tab %>%
    group_by(gene_5) %>%
    summarise(gene_name=first(gene_name), strand=first(strand), gene_3=first(gene_3), cgi_5=first(cgi_5), cgi_3=first(cgi_3), score=mean(score))
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
  "Gencode_CGI_genes_intersection_table.hg38.txt",
  quote="")
#names(intersection_table) <- c("chrom", "cpg_l", "cpg_h", "cpg_id", "gene_l", "gene_h", "gene_name", "score_filler", "strand")
intersection_table <- intersection_table %>% dplyr::rename("cpg_l" = "cgi_s", "cpg_h" = "cgi_e", "gene_l" = "gene_s", "gene_h" = "gene_e", "score_filler" = "score")
#gene_name
uniq_id <- paste0(intersection_table$tss, intersection_table$ensembl_ID)
intersection_table$gene_name <- uniq_id
##################
## 3. Get index ##
##################
plus_tss <- get_anchored_scores(
  bed=intersection_table, bw="skew_10_plus.bw", anch=1, b=0, a=100, strand="+")
minus_tss <- get_anchored_scores(bed=intersection_table, bw="skew_10_minus.bw", anch=2, b=0, a=100, strand="-")
plus_cpg <- get_anchored_scores(bed=intersection_table, bw="skew_10_plus.bw", anch=4, b=50, a=50, strand="+")
minus_cpg <- get_anchored_scores(bed=intersection_table, bw="skew_10_minus.bw", anch=3, b=50, a=50, strand="-")

plus_tss_index <- get_index(plus_tss)
plus_cpg_index <- get_index(plus_cpg)
plus_index <- plus_tss_index
lt = plus_cpg_index$score
names(lt) = plus_cpg_index$gene_name
plus_index$dist = unname((lt[as.character(plus_index$gene_name)]))
plus_index_noolap <- plus_index %>% filter(!is.na(dist)) %>% filter((cgi_3-gene_5)>=150)

minus_tss_index <- get_index(minus_tss)
minus_cpg_index <- get_index(minus_cpg)
minus_index <- minus_tss_index
lt = minus_cpg_index$score
names(lt) = minus_cpg_index$gene_name
minus_index$dist = unname((lt[as.character(minus_index$gene_name)]))
minus_index_noolap <- minus_index %>% filter(!is.na(dist)) %>% filter((gene_3-cgi_5)>=150)

Dist_signif_plus <- plus_index_noolap %>% filter(dist > score)
TSS_signif_plus <- plus_index_noolap %>% filter(score > dist)
Dist_signif_minus <- minus_index_noolap %>% filter(dist > score)
TSS_signif_minus <- minus_index_noolap %>% filter(score > dist)

write_bed_file(Dist_signif_minus, "dist_skewed_minus.txt")
write_bed_file(Dist_signif_plus, "dist_skewed_plus.txt")
write_bed_file(TSS_signif_minus, "prox_skewed_minus.txt")
write_bed_file(TSS_signif_plus, "prox_skewed_plus.txt")
#For 100 after TSS and 50 before 50 after cpg end
#dist minus is 3489, prox minus is 3870
#dist plus is 3557, prox plus is 3987