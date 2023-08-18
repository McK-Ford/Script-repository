#########################
## 1. Load in packages ##
#########################
library(Rsamtools)
library (tidyverse)
library(rtracklayer)
library(compiler)
library(BSgenome) #make sure Hg38 is locally installed - makes functions much faster for only .7 gB.
library(data.table)
enableJIT(3) #this compiles functions as we go along, which makes them run faster.
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

get_anchored_metaplot_tabs <- function(plus_bw, minus_bw, b, a, class_tab, bs=10){
#INPUT: a bigwig (if data is strand-specific, input strand specific bigwigs, else just put in same bigwig twice), base pairs before and after anchor desired, method to handle multiple bigwig bins per region (sum or mean valid). Could probably add a few more arguments if desired.
#1. Process 5' CpG anchor, 3' CpG anchor, and TSS anchor with get_anchored_scores. Detailed breakdown available in get_anchored_scores function.
  plus_tss = get_anchored_scores(
    bed=intersection_table, bw=plus_bw, anch=5, b=b, a=a, strand="+", bs=bs
  )
   minus_tss = get_anchored_scores(
     bed=intersection_table, bw=minus_bw, anch=6, b=b, a=a, strand="-", bs=bs
   )
  plus_3cpg = get_anchored_scores(
    bed=intersection_table, bw=plus_bw, anch=3, b=b, a=a, strand="+", bs=bs
  )
   minus_3cpg = get_anchored_scores(
     bed=intersection_table, bw=minus_bw, anch=2, b=b, a=a, strand="-", bs=bs
   )
  plus_5cpg = get_anchored_scores(
    bed=intersection_table, bw=plus_bw, anch=2, b=b, a=a, strand="+", bs=bs
  )
   minus_5cpg = get_anchored_scores(
     bed=intersection_table, bw=minus_bw, anch=3, b=b, a=a, strand="-", bs=bs
   )
 # 3. Add a column to all anchor tables telling us what anchor tables they're referencing.. Then bind together.
  plus_tss$anchor = rep_len(c("TSS"), dim(plus_tss)[1])
  minus_tss$anchor = rep_len(c("TSS"), dim(minus_tss)[1])
  plus_5cpg$anchor = rep_len(c("5_cpg"), dim(plus_5cpg)[1])
  minus_5cpg$anchor = rep_len(c("5_cpg"), dim(minus_5cpg)[1])
  plus_3cpg$anchor = rep_len(c("3_cpg"), dim(plus_3cpg)[1])
  minus_3cpg$anchor = rep_len(c("3_cpg"), dim(minus_3cpg)[1])
  plus_anchors = rbind(plus_5cpg, plus_tss, plus_3cpg)
  minus_anchors = rbind(minus_5cpg, minus_tss, minus_3cpg)
  #4. Make lookup tables - retrieve class from class table by TSS, which should be present in both class tables and in anchor tables.
  #5. Get the average score by bin split by anchor and class (ex - bin 1 TSS Dist would be the average score in the first bin after the TSS for all dist classed genes). 
  #6. For graphing - bins start at 1, but if you subtract the amount of base pairs you wanted before the anchor divided by the bin size, then multiply it by bin size, you'll transform the bin numbering into base pairs relative to anchor numbering. This also lets you reorient the minus strand to match the plus strand directionality. Then you can bind the plus and minus strand back together, and summarize their scores too.
  plus_anchors$bin_adjust <- (plus_anchors$bins - (b/10)) * 10
  minus_anchors$bin_adjust <- (minus_anchors$bins - (b/10)) * -10
  merged_anchors <- rbind(plus_anchors, minus_anchors)
  gene_class = class_tab$class
  names(gene_class) = class_tab$gene_name
  merged_anchors$class = unname((gene_class[as.character(merged_anchors$gene_name)]))
   merged_anchors_sum <- merged_anchors %>%
     drop_na() %>% #Na is rare but there's a few genes w/out enough distance between TSS and 3' CpG to be classed.
     group_by(anchor, class, bin_adjust) %>%
     summarise(score = mean(score), .groups = "drop")
  return(merged_anchors_sum)
}

get_anchored_scores_scaled <- function(bed, bw, anch1, anch2, strand, bin_num=40){
  chrs = list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
              "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
              "chr19", "chr20", "chr21", "chr22", "chrX") #for iterating over so you never have too much in memory at once.
  tmp_list = list()
  for (i in seq_along(chrs)) {
    partial_bw = get_bw_by_chrom(bw, chrs[[i]])
    partial_bed = get_bed_by_chrom_n_strand(bed, chrs[[i]], strand)
    partial_bed$reg_width = partial_bed[,anch2] - partial_bed[,anch1]
    binned_partial_bed = partial_bed %>% 
      filter(reg_width>=bin_num)
      binned_partial_bed$bin_width=unlist(floor(binned_partial_bed$reg_width/bin_num))
    #subset the window into bins
    binned_partial_bed$remainder = unlist(binned_partial_bed$reg_width%%bin_num)
    last_bin_start = bin_num-1
    start_bin_string = paste0(seq(from=0, to=last_bin_start, by=1), collapse=",")
    anch_win_w_bins = cbind(binned_partial_bed, start_bin_string)
    anch_win_w_bins_as_cols = separate(
      anch_win_w_bins, col = start_bin_string, into=paste(
        "X", seq(1:bin_num), sep=""
      ), sep=",", remove=TRUE, convert=TRUE
    )
    end_col = paste("X", bin_num, sep="")
    long_anch = pivot_longer(
      anch_win_w_bins_as_cols, cols=X1:end_col, names_to="Bins", names_prefix="X", values_to = "start_adjust"
    ) 
    long_anch$Bins = as.numeric(long_anch$Bins) #apparently if you don't do this it will think bins 5-9 are characters not numbers for some bizarre reason
    #bins are 5' to 3' - take that into account when handling minus down the line
    long_anch$new_start = unlist(long_anch[anch1]+long_anch$bin_width*long_anch$start_adjust)
    #split apply combine
    long_anch_1_thru_39 = long_anch %>% filter(Bins < bin_num)
    long_anch_1_thru_39$new_end = unlist(long_anch_1_thru_39$new_start + long_anch_1_thru_39$bin_width) #technically wouldn't end at 39 if bins weren't 40, but shut up, it works fine as a variable name.
    long_anch_last_bin = long_anch %>% filter(Bins == bin_num)
    long_anch_last_bin$new_end = unlist(long_anch_last_bin$new_start + long_anch_last_bin$bin_width + long_anch_last_bin$remainder)
    long_anch = rbind(long_anch_1_thru_39, long_anch_last_bin)
    long_anch = long_anch %>% arrange (new_end - new_start)
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
    #Actually getting the overlaps + subsetting the granges based on the overlaps table
    olaps = findOverlaps(anchored_granges, partial_bw)
    bw_df = data.frame(partial_bw[subjectHits(olaps)])
    regions_df = data.frame(anchored_granges[queryHits(olaps)])
    regions_w_raw_scores = cbind(bw_df, regions_df)
    names(regions_w_raw_scores) <- c(
      "chrom", "bs", "be", "bw", "star", "score", "chr2",
      "gs", "ge", "gw", "strand", "gn", "g5", "g3", "c5", "c3", "bins"
    )
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
    #olap_df = rbind(olap_df, summarized_regions_w_raw_scores)
    tmp_list[[i]] <- summarized_regions_w_raw_scores
  }
  olap_df <- do.call(rbind, tmp_list)
  olap_df$bins<-as.numeric(olap_df$bins)
  names(olap_df) = c("bins", "gene_5", "chrom", "bin_5", "bin_3", "strand", "gene_name", "gene_3", "cgi_5", "cgi_3", "bins2", "g5", "score")
  return(olap_df)
}

write_bed_file <- function(tab, b_name) {
  write.table(
    tab,
    file=b_name,
    row.names=FALSE,
    col.names=FALSE,
    quote=FALSE,
    sep="\t"
  )
}

metaplotter <- function(tab, method, y_label, color_guide, antisense = FALSE){
  #valid methods- scaled, not_scaled
  #not currently set up to handle different number of bins
  mytheme = theme(
    panel.background=element_rect(fill="white"),
    text=element_text(color="black",face="bold",family="sans"),
    axis.text=element_text(color="black"),
    axis.ticks=element_line(color="black"),
    plot.margin=unit(c(0.25, 0.25, .25, 0.25),"cm"),
    plot.title=element_text(vjust=2),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.title.y=element_text(vjust=0),
    axis.title.x=element_text(vjust=0),
    panel.border = element_blank(),
    axis.line=element_line()
  )
  if (antisense==FALSE) {
    base_geom = ggplot(data=tab,aes(x=bin_adjust, y=score, color=class)) +
      geom_point(size=1.5) +
      geom_line(size=1) +
      mytheme 
  }
  if (antisense==TRUE) {
    base_geom = ggplot(data=tab,aes(x=bin_adjust, y=score, color=class, linetype=direction)) +
      scale_linetype_manual(values= c(sense="solid", antisense="dashed")) +
      geom_line(size=1) +
      mytheme 
  }
     if (method == "scaled") {
       tailored_geom = base_geom +
         geom_vline(xintercept=40,color="grey",size=1, linetype=2) +
         geom_vline(xintercept=80,color="grey",size=1, linetype=2) +
         geom_vline(xintercept=120,color="grey",size=1, linetype=2) +
         xlab("")+
         ylab(y_label)  +
         geom_label(label="-0.8 kb", aes(x=0, y=0), show.legend = FALSE)+
         geom_label(label="5_CpG", aes(x=40, y=0), show.legend = FALSE)+
         geom_label(label="TSS", aes(x=80, y=0), show.legend = FALSE)+
         geom_label(label="3_CpG", aes(x=120, y=0), show.legend = FALSE) +
         geom_label(label="+0.8 kb", aes(x=160, y=0), show.legend = FALSE) +
         scale_x_continuous(breaks=NULL) +
         scale_color_manual(values=color_guide)
     }
   if (method == "not_scaled") {
     tailored_geom = base_geom +
       facet_grid(~fct_relevel(anchor,"5_cpg","TSS","3_cpg",)) + 
       xlab("Distance from Anchor")+
       ylab(y_label)+
       geom_vline(xintercept=0,color="grey",size=1, linetype=2)+
       geom_hline(yintercept = 0, color = "grey", size = 1, linetype=2) +
       theme(panel.grid.major=element_blank())+
       xlim(c(-500,500)) +
       scale_color_manual(values=color_guide) 
   }
   return(tailored_geom)
}

get_classes <- function(file_list, class_names_list, name_col) {
  #you need to return a table that has gene_name and class. All else is just set dressing.
  tmp_list = list()
  for (i in seq_along(file_list)) {
    file = read.delim(file_list[[i]], header=FALSE)
    sub_file = file[,name_col]
    class = rep_len(c(class_names_list[[i]]), length(sub_file))
    sub_file2 = as.data.frame(cbind(sub_file, class))
    tmp_list[[i]] = sub_file2
  }
  file_tab = do.call(rbind, tmp_list)
  names(file_tab) = c("gene_name", "class")
  return(file_tab)
}

get_scaled_cgi_metaplot_tabs <- function(plus_bw, minus_bw, class_tab, around = 800){
  plus_5_tss = get_anchored_scores_scaled(
    bed=intersection_table, bw=plus_bw, anch1=2, anch2=5, strand="+"
  )
  minus_5_tss = get_anchored_scores_scaled(
    bed=intersection_table, bw=minus_bw, anch1=6, anch2=3, strand="-"
  )
  plus_tss_3 = get_anchored_scores_scaled(
    bed=intersection_table, bw=plus_bw, anch1=5, anch2=3, strand="+"
  )
  minus_tss_3 = get_anchored_scores_scaled(
    bed=intersection_table, bw=minus_bw, anch1=2, anch2=6, strand="+"
  )
  #okay i think these are accurate anchors.
  plus_before_cgi = get_anchored_scores(
    bed=intersection_table, bw=plus_bw, anch=2, b=around, a=0, strand="+", bs=20
  )
  minus_before_cgi = get_anchored_scores(
    bed=intersection_table, bw=minus_bw, anch=3, b=0, a=around, strand="-", bs=20
  )
  plus_after_cgi = get_anchored_scores(
    bed=intersection_table, bw=plus_bw, anch=3, b=0, a=around, strand="+", bs=20
  )
  minus_after_cgi = get_anchored_scores(
    bed=intersection_table, bw=minus_bw, anch=2, b=around, a=0, strand="-", bs=20
  )
  # okay. Now I want to get my stuff together into a single DF. How?
  plus_before_cgi$newbin = plus_before_cgi$bins
  plus_5_tss$newbin = plus_5_tss$bins + 40
  plus_tss_3$newbin = plus_tss_3$bins + 80
  plus_after_cgi$newbin = plus_after_cgi$bins + 120
  
  minus_before_cgi$newbin = minus_before_cgi$bins + 120
  minus_5_tss$newbin = minus_5_tss$bins + 80
  minus_tss_3$newbin = minus_tss_3$bins + 40
  minus_after_cgi$newbin = minus_after_cgi$bins
  # 
  plus_bins = rbind(plus_before_cgi, plus_5_tss, plus_tss_3, plus_after_cgi) #so at this stage it's doing what it's supposed to do
  minus_bins = rbind(minus_before_cgi, minus_5_tss, minus_tss_3, minus_after_cgi)
  #
  plus_bins2 = plus_bins %>% filter((cgi_3 - gene_5) >=150)
  minus_bins2 = minus_bins %>% filter((gene_3 - cgi_5) >=150)
  # #Next, summarize by bin.
  plus_bins2$bin_adjust = plus_bins2$newbin
  minus_bins2$bin_adjust = 160 - minus_bins2$newbin +1 #flips to correct direction. Will need to check to make sure it works as desired. Go through all my logic in this document again just in case TBH.
  merged_bins = rbind(plus_bins2, minus_bins2)
  merged_bins = merged_bins %>% group_by(bin_adjust)
  gene_class = class_tab$class
  names(gene_class) = as.character(class_tab$gene_name)
  merged_bins$class = unname((gene_class[as.character(merged_bins$gene_name)]))
  merged_bins_sum <- merged_bins %>%
    drop_na() %>%
    group_by(class, bin_adjust) %>%
    summarise(score = mean(score), .groups = "drop")
}
###################################
## 3. Load in intersection table ##
###################################
intersection_table <- read.delim( #Reading in the CpG/Genes intersection file, from bedtools
  "Refseq_curated_CGI_promoters_filter.hg38.txt",
  header=FALSE,
  quote="")
names(intersection_table) <- c("chrom", "cpg_l", "cpg_h", "cpg_id", "gene_l", "gene_h", "gene_name", "score_filler", "strand")
intersection_table <- intersection_table %>% distinct(gene_name, .keep_all = TRUE)
#######################
## 4. Define classes ##
#######################
proseq_class <- get_classes(file_list = list(
  "dist_propause_plus2.txt",
  "dist_propause_minus2.txt",
  "prox_propause_plus2.txt",
  "prox_propause_minus2.txt",
  "silent_propause_plus2.txt",
  "silent_propause_minus2.txt"), 
  class_names_list = list("dist", "dist", "prox", "prox", "silent", "silent"), name_col = 2)

skew_class = get_classes(file_list = list("~/Lab stuff/R_coding/backedup/prox_skewed_plus.txt",
                                          "~/Lab stuff/R_coding/backedup/prox_skewed_minus.txt",
                                          "~/Lab stuff/R_coding/backedup/dist_skewed_plus.txt",
                                          "~/Lab stuff/R_coding/backedup/dist_skewed_minus.txt"),
                         class_names_list = list(
                           "prox", "prox", "dist", "dist"
                         ),
                         name_col = 2)
cgi_pause_df_3_20_15 <- read.delim(("cgi_pause_df_3_20_15.txt"))
oldsort_newtable <- data.frame(cbind(cgi_pause_df_3_20_15$name2, cgi_pause_df_3_20_15$class))
names(oldsort_newtable) <- c("gene_name", "class")
#####################################
## 5. Get the regions for metaplot ##
#####################################
tab_for_metaplotting <- get_anchored_metaplot_tabs("skew_10_plus.bw", "skew_10_minus.bw", 500, 500, proseq_class)
tbm_graph <- metaplotter(tab = tab_for_metaplotting, method = "not_scaled", y_label="mean skew", color_guide = c(prox="red", dist="blue", silent="green"))
tbm_graph

tab_for_metaplotting_s96 <- get_anchored_metaplot_tabs("Hek_rloop_S96CnT_r2.bw", "Hek_rloop_S96CnT_r2.bw", 500, 500, proseq_class)
rloop_by_proseq_graph <- metaplotter(tab = tab_for_metaplotting_s96, method = "not_scaled", y_label="mean R-loop tags", color_guide = c(prox="red", dist="blue", silent="green"))
rloop_by_proseq_graph

tab_for_metaplotting_hek293tproseq <- get_anchored_metaplot_tabs("GSM4730175_HEK293T-Proseq-rep2.fwd.bw", "GSM4730175_HEK293T-Proseq-rep2.rev.bw", 500, 500, proseq_class)
tab_for_metaplotting_antisense_hek293tproseq <- get_anchored_metaplot_tabs("GSM4730175_HEK293T-Proseq-rep2.rev.bw", "GSM4730175_HEK293T-Proseq-rep2.fwd.bw", 500, 500, proseq_class)
tab_for_metaplotting_hek293tproseq$direction = rep_len(c("sense"), dim(tab_for_metaplotting_hek293tproseq)[1])
tab_for_metaplotting_antisense_hek293tproseq$direction = rep_len(c("antisense"), dim(tab_for_metaplotting_antisense_hek293tproseq)[1])
hekproseq <- rbind(tab_for_metaplotting_hek293tproseq, tab_for_metaplotting__antisense_hek293tproseq)
hekproseq_graph <- metaplotter(tab = hekproseq, method = "not_scaled", y_label="mean proseq tags", color_guide = c(prox="red", dist="blue", silent="green"), antisense = TRUE)
hekproseq_graph


tab_for_metaplotting_groseq <- get_anchored_metaplot_tabs("Kraus.2013.Groseq.MCF7.10_plus.bw", "Kraus.2013.Groseq.MCF7.10_minus.bw", 500, 500, proseq_class)
tab_for_metaplotting_antisense_groseq <- get_anchored_metaplot_tabs("Kraus.2013.Groseq.MCF7.10_minus.bw", "Kraus.2013.Groseq.MCF7.10_plus.bw", 500, 500, proseq_class)
tab_for_metaplotting_groseq$direction = rep_len(c("sense"), dim(tab_for_metaplotting_groseq)[1])
tab_for_metaplotting_antisense_groseq$direction = rep_len(c("antisense"), dim(tab_for_metaplotting_antisense_groseq)[1])
groseq <- rbind(tab_for_metaplotting_groseq, tab_for_metaplotting_antisense_groseq)
groseq_graph <- metaplotter(tab = groseq, method = "not_scaled", y_label="mean groseq tags", color_guide = c(prox="red", dist="blue", silent="green"), antisense = TRUE)
groseq_graph

tab_for_metaplotting_79 <- get_anchored_metaplot_tabs("MCF7.Polyak.G13.H3K79me2.WT.ChIPseq.dedup.rpkm.bw", "MCF7.Polyak.G13.H3K79me2.WT.ChIPseq.dedup.rpkm.bw", 500, 500, proseq_class)
tab_for_metaplotting_79_graph <- metaplotter(tab = tab_for_metaplotting_79, method = "not_scaled", y_label="mean H3K79me2 tags", color_guide = c(prox="red", dist="blue", silent="green"))

tab_for_metaplotting_groseq_graph
tab_for_metaplotting_79_graph

pdf(file = "a_bunch_of_anchored_graphs.pdf")
polII_graph
rloops96_graph
H3K4me3_graph
hek293tproseq_graph
dev.off()

##############################
## 6. Make scaled metaplots ##
##############################
scaled_rloop_by_skew <- get_scaled_cgi_metaplot_tabs("Hek_rloop_S96CnT_r2.bw", "Hek_rloop_S96CnT_r2.bw", skew_class)
scaled_rloop_by_skew_graph <- scaled_graphs(scaled_rloop_by_skew, "mean R-loop tags", y_lim_upper = 10)
scaled_rloop_by_skew_graph

scaled_skew_by_skew <- get_scaled_cgi_metaplot_tabs("skew_10_plus.bw", "skew_10_minus.bw", skew_class)
scaled_skew_by_skew_graph <- scaled_graphs(scaled_skew_by_skew, "mean skew", y_lim_upper = 0.2, y_lim_lower = -0.1)
scaled_skew_by_skew_graph

hek_proseq_1 <- get_scaled_cgi_metaplot_tabs("GSM4730174_HEK293T-Proseq-rep1.fwd.bw", "GSM4730174_HEK293T-Proseq-rep1.rev.bw", skew_class)
scaled_proseq_by_skew_graph <- scaled_graphs(hek_proseq_1, "mean proseq tags", y_lim_upper = 30)
scaled_proseq_by_skew_graph

mcf7_groseq_2013 <- get_scaled_cgi_metaplot_tabs("Kraus.2013.Groseq.MCF7.10_plus.bw", "Kraus.2013.Groseq.MCF7.10_minus.bw", skew_class)
scaled_groseq_by_skew_graph <- scaled_graphs(mcf7_groseq_2013, "mean groseq tags", y_lim_upper = 1.5)
scaled_groseq_by_skew_graph

me3 <- get_scaled_cgi_metaplot_tabs("MCF7.Polyak.G15.H3K4me3.WT.ChIPseq.dedup.rpkm.bw", "MCF7.Polyak.G15.H3K4me3.WT.ChIPseq.dedup.rpkm.bw", proseq_class)
scaled_me3_graph <- scaled_graphs(me3, "mean h3k4me3 tags", y_lim_upper = 75)
scaled_me3_graph
##################################################################
scaled_rloop_by_proseq <- get_scaled_cgi_metaplot_tabs("Hek_rloop_S96CnT_r2.bw", "Hek_rloop_S96CnT_r2.bw", proseq_class)
scaled_rloop_by_proseq_graph <- metaplotter(tab = scaled_rloop_by_proseq, method = "scaled", y_label="mean R-loop tags", color_guide = c(prox="red", dist="blue", silent="green"))
scaled_rloop_by_proseq_graph

scaled_skew_by_proseq <- get_scaled_cgi_metaplot_tabs("skew_10_plus.bw", "skew_10_minus.bw", proseq_class)
scaled_skew_by_proseq_graph <- scaled_graphs(scaled_skew_by_proseq, "mean skew", y_lim_upper = 0.1, y_lim_lower = -0.05)
scaled_skew_by_proseq_graph <- metaplotter(tab = scaled_skew_by_proseq, method = "scaled", y_label="mean skew", color_guide = c(prox="red", dist="blue", silent="green"))
scaled_skew_by_proseq_graph

scaled_prosq_by_proseq <- get_scaled_cgi_metaplot_tabs("GSM4730174_HEK293T-Proseq-rep1.fwd.bw", "GSM4730174_HEK293T-Proseq-rep1.rev.bw", proseq_class)
scaled_prosq_by_proseq_antisense <- get_scaled_cgi_metaplot_tabs("GSM4730174_HEK293T-Proseq-rep1.rev.bw", "GSM4730174_HEK293T-Proseq-rep1.fwd.bw", proseq_class)
scaled_prosq_by_proseq$direction = rep_len(c("sense"), dim(scaled_prosq_by_proseq)[1])
scaled_prosq_by_proseq_antisense$direction = rep_len(c("antisense"), dim(scaled_prosq_by_proseq_antisense)[1])
hekproseq <- rbind(scaled_prosq_by_proseq, scaled_prosq_by_proseq_antisense)
hekproseq_graph <- metaplotter(tab = hekproseq, method = "scaled", y_label="mean proseq tags", color_guide = c(prox="red", dist="blue", silent="green"), antisense = TRUE)
hekproseq_graph

scaled_groseq_by_proseq <- get_scaled_cgi_metaplot_tabs("Kraus.2013.Groseq.MCF7.10_plus.bw", "Kraus.2013.Groseq.MCF7.10_minus.bw", proseq_class)
scaled_groseq_by_proseq_antisense <- get_scaled_cgi_metaplot_tabs("Kraus.2013.Groseq.MCF7.10_minus.bw", "Kraus.2013.Groseq.MCF7.10_plus.bw", proseq_class)
scaled_groseq_by_proseq$direction = rep_len(c("sense"), dim(scaled_groseq_by_proseq)[1])
scaled_groseq_by_proseq_antisense$direction = rep_len(c("antisense"), dim(scaled_groseq_by_proseq_antisense)[1])
groseq <- rbind(scaled_groseq_by_proseq, scaled_groseq_by_proseq_antisense)
groseq_graph <- metaplotter(tab = groseq, method = "scaled", y_label="mean groseq tags", color_guide = c(prox="red", dist="blue", silent="green"), antisense = TRUE)
groseq_graph

k79 <- get_scaled_cgi_metaplot_tabs("MCF7.Polyak.G13.H3K79me2.WT.ChIPseq.dedup.rpkm.bw", "MCF7.Polyak.G13.H3K79me2.WT.ChIPseq.dedup.rpkm.bw", proseq_class)
k79_graph <- metaplotter(tab = k79, method = "scaled", y_label="mean H3K79me2 tags", color_guide = c(prox="red", dist="blue", silent="green"))
k79_graph

##################################
scaled_rloop_by_oldsort <- get_scaled_cgi_metaplot_tabs("Hek_rloop_S96CnT_r2.bw", "Hek_rloop_S96CnT_r2.bw", oldsort_newtable)
scaled_rloop_by_oldsort_graph<- scaled_graphs(scaled_rloop_by_oldsort, "mean R-loop tags", y_lim_upper = 15)
scaled_rloop_by_oldsort_graph

scaled_skew_by_oldsort <- get_scaled_cgi_metaplot_tabs("skew_10_plus.bw", "skew_10_minus.bw", oldsort_newtable)
scaled_skew_by_oldsort_graph <- scaled_graphs(scaled_skew_by_oldsort, "mean skew", y_lim_upper=0.12, y_lim_lower = -0.07)
scaled_skew_by_oldsort_graph

scaled_proseq_by_oldsort <- get_scaled_cgi_metaplot_tabs("GSM4730174_HEK293T-Proseq-rep1.fwd.bw", "GSM4730174_HEK293T-Proseq-rep1.rev.bw", oldsort_newtable)
scaled_proseq_by_oldsort_graph <- scaled_graphs(scaled_proseq_by_oldsort, "mean proseq tags", y_lim_upper = 35)
scaled_proseq_by_oldsort_graph

scaled_groseq_by_oldsort <- get_scaled_cgi_metaplot_tabs("Kraus.2013.Groseq.MCF7.10_plus.bw", "Kraus.2013.Groseq.MCF7.10_minus.bw", oldsort_newtable)
scaled_groseq_by_oldsort_graph <- scaled_graphs(scaled_groseq_by_oldsort, "mean groseq tags", y_lim_upper = 2)
scaled_groseq_by_oldsort_graph
