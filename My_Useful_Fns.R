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
#Args: bed = table with first 6 col in bed format (regions), bam= string indicating bam filepath (reads)
#method = "peak_centered", "single_anch", "bi_anch"
#Single_anch returns score in bins size bs (int) around anchor point from b (int) bp 
  #upstream to a (int) bp downstream. Plus strand anchor is bed col 2.
  #Col 3 is minus strand anchor. n is not used.
#peak_centered returns score in bins size bs (int) around anchor point from b (int) bp upstream
  #to a (int) bp downstream RELATIVE TO GENOMIC COORDINATES. n is not used.
#bi_anch generates score scaled from start to end of bed region via splitting 
  #said region into n (int) bins. b, a, and bs are not used.
#readsOnly is only relevant for pE mode. If false, insert is included in counts.
#pairedEnd and rnorm are bools, indicating whether data is paired end and whether to rpkm normalize.
#mode = "sbp" vs "fullread". Typically sbp will be used for proseq.
#revcomp used for proseq. bool
#debug mode iterates only over 2 small chromosomes so it runs faster, and prints out extra troubleshooting info.
score_matrix <- function(bed, bam, b=NA, a=NA, n=NA, method, bs=10,
                             pairedEnd, mode, revcomp=FALSE,
                             readsOnly=FALSE, debug=FALSE, ignorestrand=FALSE,
                             rnorm=TRUE){
  #getting bam file info
  bam_info = if(pairedEnd) BamFile(bam,asMates=TRUE) else BamFile(bam)
  si=seqinfo(bam_info)
  readcounts = if (rnorm) readcounts(bam_info, pairedEnd, debug) else NA
  #Iterating over chromosomes
  chrs = list_chrs()
  if (debug){
    chrs = list("chr19", "chr22")
    print("Running in debug mode")
  }
 # chrs = list("chr14", "chr19")
  print(paste0("method is ", method, " and mode is ", mode))
  tmp_list=list=chrs
  tmp_ref_list=chrs
  for (i in seq_along(chrs)) {
    print(paste0("getting regions for chrom ", chrs[[i]], " at ", Sys.time()))
    bed_sub = bed %>% filter(bed[[1]]==chrs[[i]])
    if (method=="peak_centered") {
      hist = center_mat(bed_sub=bed_sub, a=a, b=b, bs=bs)
      } else if (method=="single_anch") {
        hist = single_anch_mat(bed_sub=bed_sub, a=a, b=b, bs=bs)
        strandvec = rep(bed_sub[[6]], each=(a-b)/bs)
      } else {
        mat_list = bi_anch_mat(bed_sub = bed_sub, n=n)
        hist = mat_list[[1]]
        ends_mat = mat_list[[2]]
        strandvec = rep(bed_sub[[6]], each=n)
        }
    if (debug) print(head(hist))
    long_hist=reshape2::melt(hist, na.rm=TRUE) #longform lets us generate 'all
    #bin starts' and 'all bin ends' vectors for scoring.
    sbp = ScanBamParam(
      which = GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start = 0, end = si@seqlengths[si@seqnames == chrs[[i]]])))
    if (pairedEnd == TRUE & readsOnly == FALSE) {
      bam_aln = granges(readGAlignmentPairs(
        bam, param = sbp)) #using pairs and coercing to GRanges is necessary if you want to keep the insert.
      if (debug) print("getting bam pe")
      } else {
        if (debug) print("getting bam se")
        bam_aln = GRanges(readGAlignments(bam_info, param = sbp))
        }
    if (debug) print(head(bam_aln))
    if (revcomp){
      if (debug) print("getting revcomp")
        temp_plus <- as.character(strand(bam_aln)) == "+"
        strand(bam_aln) <- "+"
        strand(bam_aln)[temp_plus] <- "-"
      }
    if (mode=="sbp"){
      if (debug) print("getting sbp")
      bam_aln_tmp = bam_aln
      bam_aln = resize(bam_aln_tmp, width = 1, fix = "end")
      }
    if (method=="bi_anch"){
      long_ends=reshape2::melt(ends_mat, na.rm=TRUE)
      test=GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start=long_hist[[3]], end=long_ends[[3]]), strand=strandvec) 
      olap = countOverlaps(test, bam_aln, ignore.strand=ignorestrand) 
      } else {
        olap = countOverlaps(GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start = long_hist[[3]], end = long_hist[[3]]+ bs), strand=strandvec), bam_aln, ignore.strand=ignorestrand)
        }
    if (debug) print(head(olap))
    long_hist$value=olap
    hist=acast(long_hist, Var1~Var2) #get the hist back into wideform.
    if (method=="bi_anch" & "-" %in% bed_sub[[6]]) { #flip the - strand scaled stuff
      minus_strand_rows=bed_sub[bed_sub[[6]]=="-",]
      hist[row.names(hist) %in% minus_strand_rows[,4],] =
        hist[row.names(hist) %in% minus_strand_rows[,4],rev(seq_len(ncol(hist)))]
      if (debug) print(head(hist))
      if (debug) print("flipped - hist")
    }
    if (debug) print(head(hist))
    tmp_list[[i]] = hist
    tmp_ref_list[[i]] = bed_sub
  }
  hist = do.call(rbind, tmp_list) #collapse into one matrix
  hist = if(rnorm) hist * 1e6 / readcounts else hist
  bed=do.call(rbind,tmp_ref_list)
  hist=cbind(as.data.table(bed), as.data.table(hist))
  return(hist)
}
# Common error message if you see the error "aggregation function missing, defaulting to length" this
#means there are duplicates in the name column of the bed, and they need to be unique.

center_mat <- function(bed_sub, a, b, bs){
  ref_point = round((bed_sub[[3]]-bed_sub[[2]])/2)+bed_sub[[2]] #region center
  bins = seq(from=b, to=(a-bs), by=bs)
  hist_rows = bed_sub[[4]]
  hist = matrix(bins, ncol=length(bins), nrow=length(bed_sub[[4]]),
                dimnames=list(hist_rows, bins), byrow = TRUE)
  hist=hist+ref_point #adds center to every bin to get matrix of start pos.
  hist[hist <= 0] = NA
  return(hist)
  }

single_anch_mat <- function(bed_sub, a, b, bs){
  ref_point = ifelse(bed_sub[[6]]=="+", bed_sub[[2]], bed_sub[[3]])
  hist_rows=bed_sub[[4]]
  bins = ifelse(
    bed_sub[6]=="+", list(seq(from=b, to=(a-bs), by=bs)),
    list(seq(from=a, to=(b+bs), by=-bs)))
  hist=round(do.call(rbind, bins)) #collapse bins vector into matrix
  bin_names = seq(from=b,to=a-bs,by=bs)
  hist=hist+ref_point
  dimnames(hist)=list(hist_rows, bin_names)
  hist[hist <= 0] = NA
  return(hist)
  }

bi_anch_mat <- function(bed_sub, n){
  ref_point1 = bed_sub[[2]]
  ref_point2 = bed_sub[[3]]
  empty_list = list()
  ends_list = list()
  for (j in seq_along(ref_point1)){
    bin_vec = seq(from=ref_point1[[j]], to=ref_point2[[j]], length.out=n+1)
    ends_list[[j]] = bin_vec[2]-bin_vec[1]
    empty_list[[j]]=bin_vec[1:n] #drop the extra value
  }
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

pI <- function(bed, bam, pairedEnd, pause_s, pause_e, body_s=pause_e, mode="sbp"){
  #Args: bed is table with 1st 6 cols in bed format, bam is string to bam filepath,
  #pairedEnd is bool, pause_s and pause_e are integer positions indicating bp relative
  #to TSS. gene_body implicitly defined as pause_e through gene_end. Switch mode to fullread for groseq.
  bed = bed[(bed[[3]] - bed[[2]])>=(pause_e),]
  bed = order_bed_by_chrom(bed)
  pause_bed = bed
  body_bed = bed
  pause_bed[pause_bed[[6]] == "+",][2] =  pause_bed[pause_bed[[6]] == "+",][2] + pause_s
  pause_bed[pause_bed[[6]] == "-",][3] =  pause_bed[pause_bed[[6]] == "-",][3] - pause_s
  pause_bed[pause_bed[[6]] == "+",][3] =  pause_bed[pause_bed[[6]] == "+",][2] + pause_e
  pause_bed[pause_bed[[6]] == "-",][2] =  pause_bed[pause_bed[[6]] == "-",][3] - pause_e
  body_bed[body_bed[[6]] == "+",][2] =  body_bed[body_bed[[6]] == "+",][2] + body_s
  body_bed[body_bed[[6]] == "-",][3] =  body_bed[body_bed[[6]] == "-",][3] - body_s
  pause_tab = score_matrix(bed=pause_bed, bam=bam, n=1,
                               method="bi_anch", mode=mode,
                               revcomp=TRUE, pairedEnd = pairedEnd, rnorm=FALSE, ignorestrand=FALSE)
  body_tab = score_matrix(bed=body_bed, bam=bam, n=1, method="bi_anch", mode=mode,
                              revcomp=TRUE, pairedEnd = pairedEnd, rnorm=FALSE, ignorestrand=FALSE)
  joined_tab = cbind(bed, pause_bed[,2:3], pause_tab[[ncol(pause_tab)]], body_tab[[ncol(body_tab)]])
  colnames(joined_tab) = c(colnames(bed), "pause_start", "pause_end", "pause_counts", "body_counts")
  bam_info = if(pairedEnd) BamFile(bam,asMates=TRUE) else BamFile(bam)
  #normalizing factor
  rcm=readcounts(bam_info, pairedEnd, debug=FALSE)/1e6
  tot_counts = joined_tab$pause_counts + joined_tab$body_counts
  gene_lens = (joined_tab[[3]]-joined_tab[[2]])/1e3
  joined_tab$total_norm = (tot_counts/gene_lens)/rcm
  #pI_calc
  joined_tab$lengthnorm_pause = joined_tab$pause_counts / (joined_tab$pause_end - joined_tab$pause_start)
  joined_tab$lengthnorm_body = joined_tab$body_counts / (joined_tab[[3]]-joined_tab[[2]]-pause_e)
  joined_tab$pI = joined_tab$lengthnorm_pause/joined_tab$lengthnorm_body
  return(joined_tab)
  }


readcounts <- function(bam_info, pairedEnd, debug){
  read_counts= if(pairedEnd) countBam(bam_info, param = ScanBamParam(
    flag=scanBamFlag(isFirstMateRead = TRUE)))$records else
      countBam(bam_info)$records
  if (debug) print(paste0("rnorm ", read_counts))
  return(read_counts)
  }

#order bed by chrom?
order_bed_by_chrom <- function(bed){
  chrs=list_chrs()
  el=list()
  for (k in seq_along(chrs)){
    bed_sub = bed %>% filter(bed[[1]]==chrs[[k]])
    el[[k]] = bed_sub
    }
  bed = do.call(rbind, el)
  return(bed)
}

list_chrs <- function(){
  chrs = as.list(paste0("chr", rep(1:22)))
  chrs = append(chrs, "chrX")
  return(chrs)
}









