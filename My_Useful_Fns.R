#########################
## 1. Load in packages ##
#########################
library(Rsamtools)
library(tidyverse)
library(compiler)
library(data.table)
library(GenomicAlignments)
library(reshape2)
enableJIT(3)
##################
## 2. Functions ##
##################
#Input files - df in bed format, bam filepath
#method = Regions: peak_centered, single_anch, bi_anch
#Single_anch returns score in bins size bs around anchor point from b bp 
  #upstream to a bp downstream. Plus strand anchor is bed col 2.
  #Col 3 is minus strand anchor. n is not used.
#peak_centered returns score in bins size bs around anchor point from b bp upstream
  #to a bp downstream RELATIVE TO GENOMIC COORDINATES. n is not used.
#bi_anch generates score scaled from start to end of bed region via splitting 
  #said region into n bins. b, a, and bs are not used.
#readsOnly is only relevant for pE mode. If false, insert is included in counts.
#mode = sbp vs fullread. Typically sbp will be used for proseq.
#revcomp ... used for proseq
#debug mode iterates only over 2 small chromosomes so it runs faster.
get_score_matrix <- function(bed, bam, b=NA, a=NA, n=NA, method, bs=10,
                             pairedEnd, mode, revcomp=FALSE,
                             readsOnly=FALSE, debug=FALSE, stranded=FALSE,
                             rnorm=TRUE){
  #getting bam file info
  bam_info = if(pairedEnd) BamFile(bam,asMates=TRUE) else BamFile(bam)
  si=seqinfo(bam_info)
  if (rnorm){
  read_counts= if(pairedEnd) countBam(bam_info, param = ScanBamParam(
        flag=scanBamFlag(isFirstMateRead = TRUE)))$records else
    countBam(bam_info)$records
  } 
  if (debug) print(read_counts)
  #Iterating over chromosomes
  chrs = unlist(list(paste0("chr", rep(1:22)), "chrX")) #got to be cleaner way to do this. No alt chr.
  if (debug){
    chrs = list("chr19", "chr22")
    print("Running in debug mode")
  }
  tmp_list=list()
  for (i in seq_along(chrs)) {
    print(paste("getting regions for chrom ", chrs[[i]], "at", Sys.time()))
    bed_sub = bed %>% filter(bed[[1]]==chrs[[i]])
    if (debug) print(head(bed_sub))
    if (method=="peak_centered") {
      hist = get_center_matrix(bed_sub=bed_sub, a=a, b=b, bs=bs)
      } else if (method=="single_anch") {
        hist = get_single_anch_matrix(bed_sub=bed_sub, a=a, b=b, bs=bs)
      } else {
        #method == bi_anch
        mat_list = get_bi_anch_matrix(bed_sub = bed_sub, n=n)
        hist = mat_list[[1]]
        ends_mat = mat_list[[2]]
        }
    if (debug) print(head(hist))
    hist[hist <= 0] = NA
    long_hist=reshape2::melt(hist, na.rm=TRUE) #longform lets us generate 'all
    #bin starts' and 'all bin ends' vectors for scoring.
    sbp = ScanBamParam(
      which = GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start = 0, end = si@seqlengths[si@seqnames == chrs[[i]]])))
    if (pairedEnd == TRUE & readsOnly == FALSE) {
      bam_aln = granges(readGAlignmentPairs(  #way faster than read GAlignment pairs, like 12x faster... Should I?
        bam, param = sbp)) #using pairs and coercing to GRanges is necessary if you want to keep the insert.
      } else {
      bam_aln = GRanges(readGAlignments(bam_info, param = sbp))
      }
    if (mode=="sbp"){
      bam_aln_tmp = bam_aln
      bam_aln = resize(bam_aln_tmp, width = 1, fix = "end")
    }
    if (revcomp){
        print("getting revcomp")
        temp_plus <- as.character(strand(bam_aln)) == "+"
        strand(bam_aln) <- "+"
        strand(bam_aln)[temp_plus] <- "-"
      }
    if (method=="bi_anch"){
      long_ends=reshape2::melt(ends_mat, na.rm=TRUE)
      test=GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start=as.vector(long_hist[[3]]), end=as.vector(long_ends[[3]]))) #may not need the as.vector, troubleshooting relic
      olap = countOverlaps(test, bam_aln, ignore.strand=stranded) 
      } else {
      olap = countOverlaps(GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start = long_hist[[3]], end = long_hist[[3]]+ bs)), bam_aln, ignore.strand=stranded)
    }
    if (debug) print(head(olap))
    long_hist$value=olap
    hist=acast(long_hist, Var1~Var2) #get the hist back into wideform.
    if (method=="bi_anch" & "-" %in% bed_sub[[6]]) { #flip the - strand scaled stuff
      minus_strand_rows=bed_sub[bed_sub[[6]]=="-",]
      hist[row.names(hist) %in% minus_strand_rows[,4],] =
        hist[row.names(hist) %in% minus_strand_rows[,4],rev(seq_len(ncol(hist)))]
    }
    if (debug) print(head(hist))
    tmp_list[[i]] = hist
  }
  hist = do.call(rbind, tmp_list) #collapse into one matrix
  hist = if(rnorm) hist * 1e6 / read_counts else hist
  hist=as.data.table(hist)
  return(hist)
}
# Common error message if you see the error "aggregation function missing, defaulting to length" this
#means there are duplicates in the name column of the bed, and they need to be unique.

get_center_matrix <- function(bed_sub, a, b, bs){
  print("method: peak_centered")
  ref_point = round((bed_sub[[3]]-bed_sub[[2]])/2)+bed_sub[[2]] #region center
  bins = seq(from=b, to=(a-bs), by=bs)
  hist_rows = bed_sub[[4]]
  hist = matrix(bins, ncol=length(bins), nrow=length(bed_sub[[4]]),
                dimnames=list(hist_rows, bins), byrow = TRUE)
  hist=hist+ref_point #adds center to every bin to get matrix of start pos.
  return(hist)
  }

get_single_anch_matrix <- function(bed_sub, a, b, bs){
  print("method: single_anch")
  ref_point = ifelse(bed_sub[[6]]=="+", bed_sub[[2]], bed_sub[[3]])
  hist_rows=bed_sub[[4]]
  bins = ifelse(
    bed_sub[6]=="+", list(seq(from=b, to=(a-bs), by=bs)),
    list(seq(from=a, to=(b+bs), by=-bs)))
  hist=round(do.call(rbind, bins)) #collapse bins vector into matrix
  bin_names = seq(from=b,to=a-bs,by=bs) #names need to match so can only use +
  hist=hist+ref_point
  dimnames(hist)=list(hist_rows, bin_names)
  return(hist)
  }

get_bi_anch_matrix <- function(bed_sub, n){
  print("method: bi_anch")
  #bed_sub = bed_sub %>% filter(bed_sub[[3]]-bed_sub[[2]]>100)
  ref_point1 = bed_sub[[2]]
  ref_point2 = bed_sub[[3]]
  empty_list = list()
  ends_list = list()
  for (j in seq_along(ref_point1)){
    bin_vec = seq(from=ref_point1[[j]], to=ref_point2[[j]], length.out=n+1)
    #end of region will be included in list, therefore to get bin starts of
    #desired num need to make 1 more seq than desired
    ends_list[[j]] = bin_vec[2]-bin_vec[1]
    empty_list[[j]]=bin_vec[1:n] #then we drop the extra value
  } #is there a -strand problem?
  hist=round(do.call(rbind, empty_list))
  hist_rows = bed_sub[[4]]
  bin_names = seq(from=0,to=n-1,by=1)
  dimnames(hist)=list(hist_rows, bin_names)
  ends_vec = unlist(ends_list)
  hist[hist <= 0] = NA
  ends_mat = round(hist + ends_vec) #this is our bin ends
  mat_list = list(hist, ends_mat)
  return(mat_list)
  }

get_pI <- function(bed, bam, pairedEnd, pause_s, pause_e){
  bed = bed[(bed[[3]] - bed[[2]])>=(pause_e),]
  pause_bed = bed
  body_bed = bed
  pause_bed[pause_bed[[6]] == "+",][2] =  pause_bed[pause_bed[[6]] == "+",][2] + pause_s
  pause_bed[pause_bed[[6]] == "-",][3] =  pause_bed[pause_bed[[6]] == "-",][3] - pause_s
  pause_bed[pause_bed[[6]] == "+",][3] =  pause_bed[pause_bed[[6]] == "+",][2] + pause_e
  pause_bed[pause_bed[[6]] == "-",][2] =  pause_bed[pause_bed[[6]] == "-",][3] - pause_e
  body_bed[body_bed[[6]] == "+",][2] =  body_bed[body_bed[[6]] == "+",][2] + pause_e
  body_bed[body_bed[[6]] == "-",][3] =  body_bed[body_bed[[6]] == "-",][3] - pause_e
  pause_tab = get_score_matrix(bed=pause_bed, bam=bam, n=1,
                               method="bi_anch", mode="sbp",
                               revcomp=TRUE, pairedEnd = pairedEnd, rnorm=FALSE)
  body_tab = get_score_matrix(bed=body_bed, bam=bam, n=1, method="bi_anch", mode= "sbp",
                              revcomp=TRUE, pairedEnd = pairedEnd, rnorm=FALSE)
  joined_tab = cbind(bed, pause_bed[[2]], pause_bed[[3]], pause_tab, body_tab)
  colnames(joined_tab) =
    c(colnames(bed), "pause_start", "pause_end", "pause_counts", "body_counts")
  joined_tab$lengthnorm_pause =
    joined_tab$pause_counts / ( joined_tab$pause_end - joined_tab$pause_start)
  joined_tab$lengthnorm_body =
    joined_tab$body_counts / (joined_tab[[3]]-joined_tab[[2]]-pause_e)
  joined_tab$pI = joined_tab$lengthnorm_pause/joined_tab$lengthnorm_body
  return(joined_tab)
  }
















