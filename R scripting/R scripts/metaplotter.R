#this is my metaplotter functions source file now, see application files for how to use
#########################
## 1. Load in packages ##
#########################
library(Rsamtools)
library (tidyverse)
library(rtracklayer)
library(compiler)
library(BSgenome) #make sure Hg38 is locally installed - makes functions much faster for only .7 gB. Script can function without it but is way slower
#BiocManager::install("bsgenome.hsapiens.ucsc.hg38")
library(data.table)
#########################
## 2. Define functions ##
#########################
#Assumptions in these functions: Genome is hg38, input table has chrom column named "chrom" and strand column named "strand",
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
  anch_index=match(anch, names(bed))
  for (i in seq_along(chrs)) { #iterate over chroms
    #get our subsetted starting material
    partial_bw = get_bw_by_chrom(bw, chrs[[i]])
    partial_bed = get_bed_by_chrom_n_strand(bed, chrs[[i]], strand)
    #2. Retrieve desired start coordinates of windows based on our anchors. If plus strand, this is (anchor - the number of base pairs desired before anchor). If minus strand, this is (anchor - the number of base pairs desired after anchor) (since directionality is flipped.) 
    if (strand=="+") {
      partial_bed$win_start = unlist(partial_bed[,anch_index]) - b
    }
    if (strand=="-") {
      partial_bed$win_start = unlist(partial_bed[,anch_index]) - b
    }
    #3. Split window into bins. 1st bin for ea window is win_start to win_start + bs * bin_num (1 for first bin), 2nd is win_start + bs*(bin_num-1) to win_start + bs*(bin_num), etc. So bin coordinates can be determined from bin number. Get last_bin_start from window_width - bs, then generate sequence of all bin starts from 0 to last_bin_start in bs sized steps. Put in string to add to table (table doesn't play nicely w/ it in other formats), then break string apart into bin columns. Pivot into longform table. Calc end coordinates by adding bin width to start coordinates, -1 due to inclusive ranges.
    last_bin_start = b+a-bs
    start_bin_string = paste0(seq(from=0, to=last_bin_start, by=bs), collapse=",") 
    anch_win_w_bins = cbind(partial_bed, start_bin_string)
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
    anchored_granges <- makeGRangesFromDataFrame(long_anch, keep.extra.columns = TRUE, seqnames.field ="chrom", start.field = "new_start", end.field = "new_end") #remember how I said if I found a cleaner way to do my granges function I would fix it? here it is fixed.
    olaps = findOverlaps(anchored_granges, partial_bw)
    bw_df = data.frame(partial_bw[subjectHits(olaps)])
    names(bw_df)=paste("bw_", names(bw_df), sep="")
    regions_df = data.frame(anchored_granges[queryHits(olaps)])
    regions_w_raw_scores = cbind(bw_df, regions_df)

    regions_w_raw_scores$adjustor <- (
      regions_w_raw_scores$bw_width - pmax(
        (
          regions_w_raw_scores$start - regions_w_raw_scores$bw_start
        ),0
      ) - pmax(
        (
          regions_w_raw_scores$bw_end - regions_w_raw_scores$end
        ),0
      )
    )/regions_w_raw_scores$bw_width
 
    regions_w_raw_scores$adjusted_score <- regions_w_raw_scores$bw_score * regions_w_raw_scores$adjustor
    setDT(regions_w_raw_scores) 
    tmp1 = regions_w_raw_scores[, head(.SD, 1), by=.(Bins, start)]
    tmp2 = regions_w_raw_scores[, list((sum(adjusted_score)/sum(adjustor))), by=.(Bins, start)]
    summarized_regions_w_raw_scores = cbind(tmp1, tmp2)
    tmp_list[[i]] <- summarized_regions_w_raw_scores
  }
  olap_df <- do.call(rbind, tmp_list) #bind the list of tables together into a single table.
  olap_df$Bins<-as.numeric(olap_df$Bins)
  return(olap_df)
}

get_anchored_metaplot_tabs <- function(plus_bw, minus_bw, b, a, class_tab, class_col, bs=10, source_tab){
  #function to get anchored table for making metaplot of bigwigs for 5' and 3'  cpg and TSS
  #INPUT: a bigwig (if data is strand-specific, input strand specific bigwigs, else just put in same bigwig twice), base pairs before and after anchor desired, method to handle multiple bigwig bins per region (sum or mean valid).
  #1. Process 5' CpG anchor, 3' CpG anchor, and TSS anchor with get_anchored_scores. Detailed breakdown of this processing available in get_anchored_scores function.
  plus_tss = get_anchored_scores(
    bed=source_tab, bw=plus_bw, anch="gene_l", b=b, a=a, strand="+", bs=bs
  )
  minus_tss = get_anchored_scores(
    bed=source_tab, bw=minus_bw, anch="gene_h", b=b, a=a, strand="-", bs=bs
  )
  plus_3cpg = get_anchored_scores(
    bed=source_tab, bw=plus_bw, anch="cpg_h", b=b, a=a, strand="+", bs=bs
  )
  minus_3cpg = get_anchored_scores(
    bed=source_tab, bw=minus_bw, anch="cpg_l", b=b, a=a, strand="-", bs=bs
  )
  plus_5cpg = get_anchored_scores(
    bed=source_tab, bw=plus_bw, anch="cpg_l", b=b, a=a, strand="+", bs=bs
  )
  minus_5cpg = get_anchored_scores(
    bed=source_tab, bw=minus_bw, anch="cpg_h", b=b, a=a, strand="-", bs=bs
  )
  # 3. Add a column to all anchor tables telling us what anchor tables they're referencing. Then bind together.
  plus_tss$anchor = rep_len(c("TSS"), dim(plus_tss)[1])
  minus_tss$anchor = rep_len(c("TSS"), dim(minus_tss)[1])
  plus_5cpg$anchor = rep_len(c("5_cpg"), dim(plus_5cpg)[1])
  minus_5cpg$anchor = rep_len(c("5_cpg"), dim(minus_5cpg)[1])
  plus_3cpg$anchor = rep_len(c("3_cpg"), dim(plus_3cpg)[1])
  minus_3cpg$anchor = rep_len(c("3_cpg"), dim(minus_3cpg)[1])
  plus_anchors = rbind(plus_5cpg, plus_tss, plus_3cpg)
  minus_anchors = rbind(minus_5cpg, minus_tss, minus_3cpg)
  #4. Convert bin number to base pair location in reference to anchor, including flipping the minus strand directionality. Then, if you're interested in classes add the classes as a column. Then summarize by bin, grouping by class and anchor.
  plus_anchors$bin_adjust <- (plus_anchors$bins - (b/bs)) * bs #here's part of the problem, my adjustment has bins hardcoded in. So (120 - (3000/10))*10 is negative 1800... What about if I BS it? (120-(3000/50))*50 is 3000, then 119 is 2950, 118 is 2900... Cool.
  minus_anchors$bin_adjust <- (minus_anchors$bins - (b/bs)) * -bs + bs
  merged_anchors <- rbind(plus_anchors, minus_anchors)
  if (missing(class_tab)){
    merged_anchors_sum <- merged_anchors %>%
      drop_na() %>% #relic of classification system, can I remove this line?
      group_by(anchor, bin_adjust) %>%
      summarise(score = mean(score), .groups = "drop") 
  } else {
    class_col_num=match(class_col, names(class_tab))
    gene_class = class_tab[,class_col_num]
    names(gene_class) = class_tab$gene_name
    merged_anchors$class = unname((gene_class[as.character(merged_anchors$gene_name)]))
    merged_anchors_sum <- merged_anchors %>%
      drop_na() %>% #Na is rare but there's a few genes w/out enough distance between TSS and 3' CpG to be classed.
      group_by(anchor, class, bin_adjust) %>%
      summarise(score = mean(score), .groups = "drop")
  }
  return(merged_anchors_sum) #This is plottable. Currently, the metaplotter function is only compatible with classed data - see below on plotting unclassed.
}

get_anchored_scores_scaled <- function(bed, bw, anch1, anch2, strand, bin_num=40){
  #The version of get_anchored_scores for regions scaled between two anchors. Much of the logic matches get_anchored_scores logic - see that if not explained here.
  chrs = list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
              "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
              "chr19", "chr20", "chr21", "chr22", "chrX")
  tmp_list = list()
  anch_index1=match(anch1, names(bed))
  anch_index2=match(anch2, names(bed))
  for (i in seq_along(chrs)) {
    partial_bw = get_bw_by_chrom(bw, chrs[[i]])
    partial_bed = get_bed_by_chrom_n_strand(bed, chrs[[i]], strand)
    partial_bed$reg_width = partial_bed[,anch_index2] - partial_bed[,anch_index1]
    binned_partial_bed = partial_bed %>% 
      filter(reg_width>=bin_num) #regions need to be wider than the number of desired bins as you cannot have a bin smaller than 1 bp.
    binned_partial_bed$bin_width=unlist(floor(binned_partial_bed$reg_width/bin_num)) #different width for each window, depending on distance between the two anchors.
    #subset the window into bins
    binned_partial_bed$remainder = unlist(binned_partial_bed$reg_width%%bin_num) #window may not be evenly divisible. chosen handling method is to add remaining bp to final bin. There's got to be a better way tho
    last_bin_start = bin_num-1 #technically these aren't bin starts they're how much to multiply the bin size by to get the bin start
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
    ) #start adjust is 1 less than bin number
    long_anch$Bins = as.numeric(long_anch$Bins) #apparently if you don't do this it will think bins 5-9 are characters not numbers for some bizarre reason
    #bins are 5' to 3' - take that into account when handling minus down the line
    long_anch$new_start = unlist(long_anch[anch1]+long_anch$bin_width*long_anch$start_adjust)
    #split apply combine
    long_anch_1_thru_39 = long_anch %>% filter(Bins < bin_num) #Bins are 1-40, carried over from column names, bin_num is 40, so this gets 1-39 to one table and bin 40 to another.
    long_anch_1_thru_39$new_end = unlist(long_anch_1_thru_39$new_start + long_anch_1_thru_39$bin_width) #technically wouldn't end at 39 if bins weren't 40, but shut up, it works fine as a variable name.
    long_anch_last_bin = long_anch %>% filter(Bins == bin_num)
    long_anch_last_bin$new_end = unlist(long_anch_last_bin$new_start + long_anch_last_bin$bin_width + long_anch_last_bin$remainder)#should = window end, could just use that couldn't I?
    long_anch = rbind(long_anch_1_thru_39, long_anch_last_bin)
    long_anch = long_anch %>% arrange (new_end - new_start)
    anchored_granges <- makeGRangesFromDataFrame(long_anch, keep.extra.columns = TRUE, seqnames.field ="chrom", start.field = anch1, end.field = anch2) #remember how I said if I found a cleaner way to do my granges function I would fix it? here it is fixed.
    #Actually getting the overlaps + subsetting the granges based on the overlaps table
    olaps = findOverlaps(anchored_granges, partial_bw)
    bw_df = data.frame(partial_bw[subjectHits(olaps)])
    names(bw_df)=paste("bw_", names(bw_df), sep="")
    regions_df = data.frame(anchored_granges[queryHits(olaps)])
    regions_w_raw_scores = cbind(bw_df, regions_df)
    #separated out bc the math here is a lot
    regions_w_raw_scores$adjustor <- (
      regions_w_raw_scores$bw_width - pmax(
        (
          regions_w_raw_scores$start - regions_w_raw_scores$bw_start
          ),0
        ) - pmax(
          (
            regions_w_raw_scores$bw_end - regions_w_raw_scores$end
            ),0
          )
      )/regions_w_raw_scores$bw_width
    #
    regions_w_raw_scores$adjusted_score <- regions_w_raw_scores$bw_score * regions_w_raw_scores$adjustor
    setDT(regions_w_raw_scores) 
    tmp1 = regions_w_raw_scores[, head(.SD, 1), by=.(Bins, start)]
    tmp2 = regions_w_raw_scores[, list((sum(adjusted_score)/sum(adjustor))), by=.(Bins, start)]
    summarized_regions_w_raw_scores = cbind(tmp1, tmp2)
    tmp_list[[i]] <- summarized_regions_w_raw_scores
  }
  olap_df <- do.call(rbind, tmp_list)
  olap_df$Bins<-as.numeric(olap_df$Bins)
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
  #You need an appropiately formatted (3 or 6 columns in correct order) table for this but these are good params for making a tab-separated genomic table of any form.
}

metaplotter <- function(tab, method, y_label, color_guide, antisense = FALSE){
  #valid methods- scaled, not_scaled
  #not currently set up to handle different number of bins, not currently set up to handle nonclassed data. Antisense indicates strand specific data.
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

get_scaled_cgi_metaplot_tabs <- function(plus_bw, minus_bw, class_tab, class_col, around = 800){
  #This is not currently compatible with alternate numbers (not 40) of bins - to be worked on.
  plus_5_tss = get_anchored_scores_scaled(
    bed=intersection_table, bw=plus_bw, anch1="cpg_l", anch2="gene_l", strand="+"
  )
  minus_5_tss = get_anchored_scores_scaled(
    bed=intersection_table, bw=minus_bw, anch1="gene_h", anch2="cpg_h", strand="-"
  )
  plus_tss_3 = get_anchored_scores_scaled(
    bed=intersection_table, bw=plus_bw, anch1="gene_l", anch2="cpg_h", strand="+"
  )
  minus_tss_3 = get_anchored_scores_scaled(
    bed=intersection_table, bw=minus_bw, anch1="cpg_l", anch2="gene_h", strand="+"
  )
  plus_before_cgi = get_anchored_scores(
    bed=intersection_table, bw=plus_bw, anch="cpg_l", b=around, a=0, strand="+", bs=20
  )
  minus_before_cgi = get_anchored_scores(
    bed=intersection_table, bw=minus_bw, anch="cpg_h", b=0, a=around, strand="-", bs=20
  )
  plus_after_cgi = get_anchored_scores(
    bed=intersection_table, bw=plus_bw, anch="cpg_h", b=0, a=around, strand="+", bs=20
  )
  minus_after_cgi = get_anchored_scores(
    bed=intersection_table, bw=minus_bw, anch="cpg_l", b=around, a=0, strand="-", bs=20
  )
  # This gets everything into a single dataframe so I can draw pretty scaled graphs.
  plus_before_cgi$newbin = plus_before_cgi$bins
  plus_5_tss$newbin = plus_5_tss$bins + 40
  plus_tss_3$newbin = plus_tss_3$bins + 80
  plus_after_cgi$newbin = plus_after_cgi$bins + 120
  
  minus_before_cgi$newbin = minus_before_cgi$bins + 120
  minus_5_tss$newbin = minus_5_tss$bins + 80
  minus_tss_3$newbin = minus_tss_3$bins + 40
  minus_after_cgi$newbin = minus_after_cgi$bins
  
  plus_bins = rbind(plus_before_cgi, plus_5_tss, plus_tss_3, plus_after_cgi)
  minus_bins = rbind(minus_before_cgi, minus_5_tss, minus_tss_3, minus_after_cgi)
  
  plus_bins2 = plus_bins %>% filter((cgi_3 - gene_5) >=150)
  minus_bins2 = minus_bins %>% filter((gene_3 - cgi_5) >=150) #I should take this out of the general function, I think this was here for proseq pause calling.
  # #Next, summarize by bin.
  plus_bins2$bin_adjust = plus_bins2$newbin
  minus_bins2$bin_adjust = 160 - minus_bins2$newbin +1 #flips to correct direction for graphing purposes
  merged_bins = rbind(plus_bins2, minus_bins2)
  if (missing(class_tab)){
    merged_bins_sum <- merged_bins %>%
      drop_na() %>%
      group_by(bin_adjust) %>%
      summarise(score = mean(score), .groups = "drop") 
  } else{
    class_col_num=match(class_col, names(class_tab))
    gene_class = class_tab[,class_col_num]
    names(gene_class) = class_tab$gene_name
    merged_bins$class = unname((gene_class[as.character(merged_bins$gene_name)]))
    merged_bins_sum <- merged_bins %>%
      drop_na() %>%
      group_by(class, bin_adjust) %>%
      summarise(score = mean(score), .groups = "drop")
  }
  return(merged_bins_sum)
}
