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
#possible modes: peak_centered, single_stranded_anchored, bi_stranded_anchored

#single_stranded_anchored returns the score in bins size bs surrounding the
#anchor point from b bp upstream to a bp downstream. Plus strand anchor is
# col 2. Col 3 is minus strand anchor. n is not required for this mode.

#peak_centered returns the score in bins size bs surrounding the center of the
#regions in bed from b bp upstream to a bp downstream relative to genomic
#coordinates. n is not required for this mode.

#bi_stranded_anchored generates score scaled from start to end of bed region
# via splitting said region into n bins. b, a, and bs are not required for this
#mode.

#if you see the error "aggregation function missing, defaulting to length" this
#means there are duplicates in the name column of the bed, and they need to be unique.

get_score_matrix <- function(bed, bam, b=NA, a=NA, n=NA, method, bs=10,
                             pairedEnd=TRUE, readsOnly=TRUE, debug=FALSE){
  if (pairedEnd==TRUE) {
    bam_info=BamFile(bam,asMates=TRUE)
  } else {bam_info=BamFile(bam)}
  si=seqinfo(bam_info)
  read_counts= ifelse(
    pairedEnd==TRUE,
    countBam(
      bam_info, param = ScanBamParam(flag=scanBamFlag(isFirstMateRead = TRUE))
      )$records,
    countBam(
      bam_info)$records
  )
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
    bed_sub = bed %>% filter(bed[[1]]==chrs[[i]])
    if (debug==TRUE){
      print(head(bed_sub))
      }
    if (method=="peak_centered") {
      print("method: peak_centered")
      ref_point = round((bed_sub[[3]]-bed_sub[[2]])/2)+bed_sub[[2]] #region center
      bins = seq(from=b, to=(a-bs), by=bs)
      hist_rows = bed_sub[[4]]
      hist = matrix(bins, ncol=length(bins), nrow=length(bed_sub[[4]]),
                    dimnames=list(hist_rows, bins), byrow = TRUE)
      hist=hist+ref_point #adds center to every bin to get matrix of start pos.
    }
    if (method=="single_stranded_anchored") {
      print("method: single_stranded_anchored")
      ref_point = ifelse(bed_sub[[6]]=="+", bed_sub[[2]], bed_sub[[3]])
      hist_rows=bed_sub[[4]]
      bins = ifelse(
        bed_sub[6]=="+", list(seq(from=b, to=(a-bs), by=bs)),
        list(seq(from=a, to=(b+bs), by=-bs)))
      hist=round(do.call(rbind, bins)) #collapse bins vector into matrix
      bin_names = seq(from=b,to=a-bs,by=bs) #names need to match so can only use +
      hist=hist+ref_point
      dimnames(hist)=list(hist_rows, bin_names)
    }
    if (method=="bi_stranded_anchored") {
      print("method: bi_stranded_anchored")
      bed_sub = bed_sub %>% filter(bed_sub[[3]]-bed_sub[[2]]>100)
      ref_point1 = bed_sub[[2]] #Have to fix minus strand directionality later
      ref_point2 = bed_sub[[3]] #granges doesn't like backwards regions
      empty_list = list()
      ends_list = list()
      for (j in seq_along(ref_point1)){
        bin_vec = seq(from=ref_point1[[j]], to=ref_point2[[j]], length.out=n+1)
        #end of region will be included in list, therefore to get bin starts of
        #desired num need to make 1 more seq than desired
        ends_list[[j]] = bin_vec[2]-bin_vec[1]
        empty_list[[j]]=bin_vec[1:n] #then we drop the extra value
      }
      hist=round(do.call(rbind, empty_list))
      hist_rows = bed_sub[[4]]
      bin_names = seq(from=0,to=n-1,by=1)
      dimnames(hist)=list(hist_rows, bin_names)
      ends_vec = unlist(ends_list)
      hist[hist <= 0] = NA
      ends_mat = round(hist + ends_vec) #this is our bin ends
    }
    if (debug==TRUE){
      print(head(hist))
    }
    hist[hist <= 0] = NA
    sbp = ScanBamParam(
      which = GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start = 0, end = si@seqlengths[si@seqnames == chrs[[i]]])))
    if (pairedEnd==TRUE & readsOnly==FALSE) {
        bam_aln = granges(readGAlignmentPairs(  #way faster than read GAlignment pairs, like 12x faster... Should I?
          bam, param = sbp)) #using pairs and coercing to GRanges is necessary if you want to keep the insert.
        ### Something I want to test 1:
      #bam_aln = bam_aln[width(bam_aln)>=200]
      #bam_aln=resize(bam_aln, width=50, fix="center")
        } #12 secs if you want readGAlignmentPairs, 1-2 if you read GAlignments
    else {
      bam_aln = readGAlignments(bam_info, param = sbp)
      #bam_aln=resize(granges(bam_aln), width=1, fix="start")
      #bam_aln = readGAlignmentPairs(bam_info, param = sbp)
    }
    long_hist=reshape2::melt(hist, na.rm=TRUE) #longform lets us generate 'all
    #bin starts' and 'all bin ends' vectors for scoring.
    if (method=="bi_stranded_anchored"){
      long_ends=reshape2::melt(ends_mat, na.rm=TRUE)
      test=GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start=as.vector(long_hist[[3]]), end=as.vector(long_ends[[3]]))) #may not need the as.vector, troubleshooting relic
      olap = countOverlaps(test, bam_aln, ignore.strand=TRUE) 
    }
    else {
      olap = countOverlaps(GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start = long_hist[[3]], end = long_hist[[3]]+ bs)), bam_aln, ignore.strand=TRUE)
    }
    if (debug==TRUE){
      print(head(olap))
    }
    long_hist$value=olap
    hist=acast(long_hist, Var1~Var2) #get the hist back into wideform.
    if (method=="bi_stranded_anchored" & "-" %in% bed_sub[[6]]) { #flip the - strand scaled stuff
      minus_strand_rows=bed_sub[bed_sub[[6]]=="-",]
      hist[row.names(hist) %in% minus_strand_rows[,4],] = hist[row.names(hist) %in% minus_strand_rows[,4],rev(seq_len(ncol(hist)))]
    }
    if (debug==TRUE){
      print(head(hist))
    }
    tmp_list[[i]] = hist
  }
  hist = do.call(rbind, tmp_list) #collapse into one matrix
  hist = hist * 1e6 / read_counts #normalize
  hist=as.data.table(hist, keep.rownames=TRUE)
  return(hist)
}
