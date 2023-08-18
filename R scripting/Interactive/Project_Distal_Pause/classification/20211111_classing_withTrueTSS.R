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
  
  chrs = list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")
  #chrs = list("chrX")
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
    anchored_granges$gennames <- long_anch$uniq_ID
    anchored_granges$gen_5 <- long_anch$gene_5
    anchored_granges$gen_3 <- long_anch$gene_3
    anchored_granges$cgi_5 <- long_anch$cgi_s
    anchored_granges$cgi_3 <- long_anch$cgi_e
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
    group_by(gene_name) %>%
    summarise(gene_5=dplyr::first(gene_5), strand=dplyr::first(strand), gene_3=dplyr::first(gene_3), cgi_5=dplyr::first(cgi_5), cgi_3=dplyr::first(cgi_3), score=sum(score)) #mean for skew, sum for tags
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

####################
## 3. Import Data ##
####################
load(file="gettingTruTSSes.RData")
col_trimmed_p_silent$tss_s = col_trimmed_p_silent$gene_start
col_trimmed_m_silent$tss_s = col_trimmed_m_silent$gene_start
Plus <- rbind(col_trimmed_p, col_trimmed_p_silent)
Minus <- rbind(col_trimmed_m, col_trimmed_m_silent)

load(file="gettingTruTSSesX.RData")
col_trimmed_p_silent$tss_s = col_trimmed_p_silent$gene_start
col_trimmed_m_silent$tss_s = col_trimmed_m_silent$gene_start
Plus <- rbind(Plus, col_trimmed_p, col_trimmed_p_silent)
Minus <- rbind(Minus, col_trimmed_m, col_trimmed_m_silent)
#########################
## 4. Get skew indices ##
#########################
Plus$gene_5 = Plus$tss_s
Plus$gene_3 = Plus$gene_end
Minus$gene_5 = Minus$gene_start
Minus$gene_3 = Minus$tss_s

plus_tss <- get_anchored_scores(
  bed=Plus, bw="plus_strand_skew_10.bw", anch=18, b=0, a=100, strand="+")
minus_tss <- get_anchored_scores(bed=Minus, bw="minus_strand_skew_10.bw", anch=18, b=0, a=100, strand="-")
plus_cpg <- get_anchored_scores(bed=Plus, bw="plus_strand_skew_10.bw", anch=4, b=50, a=50, strand="+")
minus_cpg <- get_anchored_scores(bed=Minus, bw="minus_strand_skew_10.bw", anch=3, b=50, a=50, strand="-")

plus_tss_index <- get_index(plus_tss)
plus_cpg_index <- get_index(plus_cpg)
plus_index <- plus_tss_index
lt = plus_cpg_index$score
names(lt) = plus_cpg_index$gene_name
plus_index$dist = unname((lt[as.character(plus_index$gene_name)]))
plus_index_skew <- plus_index
#plus_index_noolap <- plus_index %>% filter(!is.na(dist)) %>% filter((cgi_3-gene_5)>=150)

minus_tss_index <- get_index(minus_tss)
minus_cpg_index <- get_index(minus_cpg)
minus_index <- minus_tss_index
lt = minus_cpg_index$score
names(lt) = minus_cpg_index$gene_name
minus_index$dist = unname((lt[as.character(minus_index$gene_name)]))
minus_index_skew <- minus_index

###############################
## 5. Get proseq tag indices ##
###############################
plus_tss <- get_anchored_scores(
  bed=Plus, bw="sorted_plus.bw", anch=18, b=0, a=100, strand="+")
minus_tss <- get_anchored_scores(bed=Minus, bw="sorted_minus.bw", anch=18, b=0, a=100, strand="-")
plus_cpg <- get_anchored_scores(bed=Plus, bw="sorted_plus.bw", anch=4, b=50, a=50, strand="+")
minus_cpg <- get_anchored_scores(bed=Minus, bw="sorted_minus.bw", anch=3, b=50, a=50, strand="-")

plus_tss_index <- get_index(plus_tss)
plus_cpg_index <- get_index(plus_cpg)
plus_index <- plus_tss_index
lt = plus_cpg_index$score
names(lt) = plus_cpg_index$gene_name
plus_index$dist = unname((lt[as.character(plus_index$gene_name)]))
plus_index_tags = plus_index
#plus_index_noolap <- plus_index %>% filter(!is.na(dist)) %>% filter((cgi_3-gene_5)>=150)

minus_tss_index <- get_index(minus_tss)
minus_cpg_index <- get_index(minus_cpg)
minus_index <- minus_tss_index
lt = minus_cpg_index$score
names(lt) = minus_cpg_index$gene_name
minus_index$dist = unname((lt[as.character(minus_index$gene_name)]))
minus_index_tags = minus_index

###################################
## 6. bind this info to the tabs ##
###################################
tag_prox_plus = plus_index_tags$score
names(tag_prox_plus) = plus_index_tags$gene_name
Plus$prox_score_tags = unname((tag_prox_plus[as.character(Plus$uniq_ID)]))

tag_dist_plus = plus_index_tags$dist
names(tag_dist_plus) = plus_index_tags$gene_name
Plus$dist_score_tags = unname((tag_dist_plus[as.character(Plus$uniq_ID)]))

skew_prox_plus = plus_index_skew$score
names(skew_prox_plus) = plus_index_skew$gene_name
Plus$prox_score_skew = unname((skew_prox_plus[as.character(Plus$uniq_ID)]))

skew_dist_plus = plus_index_skew$dist
names(skew_dist_plus) = plus_index_skew$gene_name
Plus$dist_score_skew = unname((skew_dist_plus[as.character(Plus$uniq_ID)]))

Plus = distinct(Plus)

tag_prox_minus = minus_index_tags$score
names(tag_prox_minus) = minus_index_tags$gene_name
Minus$prox_score_tags = unname((tag_prox_minus[as.character(Minus$uniq_ID)]))

tag_dist_minus = minus_index_tags$dist
names(tag_dist_minus) = minus_index_tags$gene_name
Minus$dist_score_tags = unname((tag_dist_minus[as.character(Minus$uniq_ID)]))

skew_prox_minus = minus_index_skew$score
names(skew_prox_minus) = minus_index_skew$gene_name
Minus$prox_score_skew = unname((skew_prox_minus[as.character(Minus$uniq_ID)]))

skew_dist_minus = minus_index_skew$dist
names(skew_dist_minus) = minus_index_skew$gene_name
Minus$dist_score_skew = unname((skew_dist_minus[as.character(Minus$uniq_ID)]))

Minus = distinct(Minus)

########################
## Filtering and labeling (still need to add wendy classes on too)
########################
plus_short = Plus %>% filter(Downdist<200 & chrom!="chrX") #need to fix chrX
plus_goodlen = Plus %>% filter(Downdist>=200 & chrom!="chrX")
minus_short = Minus %>% filter(Downdist<200 & chrom!="chrX") #need to fix chrX
minus_goodlen = Minus %>% filter(Downdist>=200 & chrom!="chrX")

plus_short = Plus %>% filter(Downdist<200) #need to fix chrX
plus_goodlen = Plus %>% filter(Downdist>=200)
minus_short = Minus %>% filter(Downdist<200) #need to fix chrX
minus_goodlen = Minus %>% filter(Downdist>=200)

plus_dist_tab = plus_goodlen %>% filter(prox_score_tags<dist_score_tags)
plus_dist_tab$class_by_tags = rep("dist", dim(plus_dist_tab)[1])
plus_prox_tab = plus_goodlen %>% filter(prox_score_tags>=dist_score_tags)
plus_prox_tab$class_by_tags = rep("prox", dim(plus_prox_tab)[1])
plus_goodlen_c1 = rbind(plus_prox_tab, plus_dist_tab)

plus_dist_skew = plus_goodlen_c1 %>% filter(prox_score_skew<dist_score_skew)
plus_dist_skew$class_by_skew = rep("dist", dim(plus_dist_skew)[1])
plus_prox_skew = plus_goodlen_c1 %>% filter(prox_score_skew>=dist_score_skew)
plus_prox_skew$class_by_skew = rep("prox", dim(plus_prox_skew)[1])
plus_classed = rbind(plus_prox_skew, plus_dist_skew)

minus_dist_tab = minus_goodlen %>% filter(prox_score_tags<dist_score_tags)
minus_dist_tab$class_by_tags = rep("dist", dim(minus_dist_tab)[1])
minus_prox_tab = minus_goodlen %>% filter(prox_score_tags>=dist_score_tags)
minus_prox_tab$class_by_tags = rep("prox", dim(minus_prox_tab)[1])
minus_goodlen_c1 = rbind(minus_prox_tab, minus_dist_tab)

minus_dist_skew = minus_goodlen_c1 %>% filter(prox_score_skew<dist_score_skew)
minus_dist_skew$class_by_skew = rep("dist", dim(minus_dist_skew)[1])
minus_prox_skew = minus_goodlen_c1 %>% filter(prox_score_skew>=dist_score_skew)
minus_prox_skew$class_by_skew = rep("prox", dim(minus_prox_skew)[1])
minus_classed = rbind(minus_prox_skew, minus_dist_skew)

save(minus_classed, plus_classed, minus_short, plus_short, file="TruTSSesClassed.RData")
load(file="TruTSSesClassed.RData")
##############################################
#Import classes
cgi_pause_df_3_20_15 <- read.delim("~/Lab stuff/R_coding/backedup/cgi_pause_df_3_20_15.txt", quote="")
old_w_groups <- cgi_pause_df_3_20_15 %>%
  select(-c(score, cdsStartStat, cdsEndStat)) %>%
  add_count(name2)
old_uniq_name <- old_w_groups %>% filter(n==1)
old_duped_name <- old_w_groups %>% filter(n!=1)
wendy_name_class_tab1 <- old_uniq_name %>% select(c(name2, class)) #this is the first section of this reference being usable.
testing_dupname_dupclass = old_duped_name %>%
  add_count(name2, class)
#if nn == n then we know all the classes are equal by the name, so we can just distinct it out.
names_w_identical_classes = testing_dupname_dupclass %>%
  filter(n==nn) %>%
  select(c(name2, class)) %>%
  distinct()
names_wout_identclasses = testing_dupname_dupclass %>%
  filter(n!=nn)
wendy_name_class_tab2 = rbind(wendy_name_class_tab1, names_w_identical_classes) #usable reference
#then what for the other 700? Pick the class associated with the maximum pausing index since they didn't list dist and prox?
wout_classes_maxpi = names_wout_identclasses %>%
  group_by(name2) %>%
  mutate(maxpi = max(pi)) %>%
  ungroup() %>%
  filter(pi==maxpi)
maxpi_name_class = wout_classes_maxpi %>% select(c(name2, class))
wendy_name_class = rbind(wendy_name_class_tab2, maxpi_name_class)
################################################################
wendylt = wendy_name_class$class
names(wendylt) = wendy_name_class$name2
plus_classed$wendy_class = unname((wendylt[as.character(plus_classed$gene_name)]))
minus_classed$wendy_class = unname((wendylt[as.character(minus_classed$gene_name)]))

write.table(
  plus_classed,
  file="plus_classed_111021.hg38.txt",
  row.names=FALSE,
  col.names=TRUE,
  quote=FALSE,
  sep="\t"
)
write.table(
  minus_classed,
  file="minus_classed_111021.hg38.txt",
  row.names=FALSE,
  col.names=TRUE,
  quote=FALSE,
  sep="\t"
)

##########################################################
#Let's get some information about what these tables look like, then we can heatmap them.
#fudge do I still need to fix chrom X? I may... That's embarrassing but we're probably also alright ignoring it for now and coming back to it later, it is a weird chromosome anyway
plus_classed <- read.delim("~/Lab stuff/R_coding/backedup/plus_classed_111021.hg38.txt", row.names=NULL, quote="")
minus_classed <- read.delim("~/Lab stuff/R_coding/backedup/minus_classed_111021.hg38.txt", row.names=NULL, quote="")
#Okay, information (which will need to be copied to notebook):
#minus classed is 5809 islands. plus classed is 6006 islands. So 11815 total.
#this is filtered for 'intersects with CpG island' and 'downdist is greater than or equal to 200'.
#actually lets pivot this for visualization
plus_classed_long = plus_classed %>% pivot_longer(cols=c(class_by_skew, class_by_tags, wendy_class), names_to = "called_by", values_to = "called_class")
minus_classed_long = minus_classed %>% pivot_longer(cols=c(class_by_skew, class_by_tags, wendy_class), names_to = "called_by", values_to = "called_class")
ggplot(data=plus_classed_long, aes(fill=called_class, x=called_by)) + geom_bar(position="fill")
#oh I need to fix the wendy classes so they match the rest. Also potentially add a silent class to tags. Neutral to skew?
plus_classed$wendy_class = recode(plus_classed$wendy_class, Proximal = "prox", Distal = "dist", Silent = "silent")
minus_classed$wendy_class = recode(minus_classed$wendy_class, Proximal = "prox", Distal = "dist", Silent = "silent")
######################################
summary(plus_classed$prox_score_tags)
summary(plus_classed$dist_score_tags)
quantile(plus_classed$prox_score_tags, .20)
summary(minus_classed$prox_score_tags)
summary(minus_classed$dist_score_tags) #crap I accidentally screwed up minus classed. I need to remake it bc the dist tags are exactly the same as the prox.
#for my information for the plus, the 20th percential is 5 reads. Lets set that as our silent cutoff for now.
plus_classed_sub1 = plus_classed %>% filter(prox_score_tags>5)
plus_classed_sub2 = plus_classed %>% filter(prox_score_tags<=5)
plus_classed_sub2$class_by_tags = rep("silent", dim(plus_classed_sub2)[1])
plus_classed = rbind(plus_classed_sub1, plus_classed_sub2)
##
minus_classed_sub1 = minus_classed %>% filter(prox_score_tags>5)
minus_classed_sub2 = minus_classed %>% filter(prox_score_tags<=5)
minus_classed_sub2$class_by_tags = rep("silent", dim(minus_classed_sub2)[1])
minus_classed = rbind(minus_classed_sub1, minus_classed_sub2)
#####################################
plus_classed_long = plus_classed %>% pivot_longer(cols=c(class_by_skew, class_by_tags, wendy_class), names_to = "called_by", values_to = "called_class")
minus_classed_long = minus_classed %>% pivot_longer(cols=c(class_by_skew, class_by_tags, wendy_class), names_to = "called_by", values_to = "called_class")
temp_for_graphing = rbind(plus_classed_long, minus_classed_long)
ggplot(data=temp_for_graphing, aes(fill=called_class, x=called_by)) + geom_bar(position="fill") + theme_bw()
#########################
# Lets make some beds for trying graphing? Actually, before we do that lets wait to hear back from Paula about the fragmentation step.
