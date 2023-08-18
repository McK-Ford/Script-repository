source("/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/NewMetaplotterJoshBased.R")

###############
## Get files ##
###############
refseq_curated_longest_cgi_minus <- read.delim("/Users/kayle/Box/Vertinolab/McKayla Ford/Data/refseq_curated_longest_cgi_minus.hg38.txt")
refseq_curated_longest_cgi_plus <- read.delim("/Users/kayle/Box/Vertinolab/McKayla Ford/Data/refseq_curated_longest_cgi_plus.hg38.txt")
nascentDirectory="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/"
plus_strand_list=list("MCF7_H9_plus.bam",
                      "MCF7_B7_plus.bam",
                      "MCF7_G11_plus.bam",
                      "MCF7_C11_plus.bam")
plus_strand_list=paste0(nascentDirectory, plus_strand_list)
minus_strand_list=list("MCF7_H9_minus.bam",
                      "MCF7_B7_minus.bam",
                      "MCF7_G11_minus.bam",
                      "MCF7_C11_minus.bam")
minus_strand_list=paste0(nascentDirectory, minus_strand_list)
refseq_curated_longest_cgi_plus$name=paste0(refseq_curated_longest_cgi_plus$name,
                                            "_",
                                            refseq_curated_longest_cgi_plus$gene_s)
TSS_CpG_end_plus= refseq_curated_longest_cgi_plus %>%
  select(1,5,3,8,2,7) #needs to be a bed
refseq_curated_longest_cgi_minus$name=paste0(refseq_curated_longest_cgi_minus$name,
                                            "_",
                                            refseq_curated_longest_cgi_minus$gene_e)
TSS_CpG_end_minus= refseq_curated_longest_cgi_minus %>%
  select(1,2,6,8,3,7) #cgi end placeholder
#########################
## Get promoter scores ##
#########################
#Goal is to classify expression
plus_mats = lapply(plus_strand_list, get_score_matrix, bed=TSS_CpG_end_plus, b=0, a=500,
                  method="single_stranded_anchored", bs=500, pairedEnd=FALSE)
minus_mats = lapply(minus_strand_list, get_score_matrix, bed=TSS_CpG_end_minus, b=0, a=500,
                   method="single_stranded_anchored", bs=500, pairedEnd=FALSE)
merged_mat_plus=cbind(plus_mats[[1]], plus_mats[[2]][[2]], plus_mats[[3]][[2]],
                 plus_mats[[4]][[2]])
merged_mat_minus=cbind(minus_mats[[1]], minus_mats[[2]][[2]], minus_mats[[3]][[2]],
                      minus_mats[[4]][[2]])
merged_mat=rbind(merged_mat_plus, merged_mat_minus)
#names
colnames(merged_mat)=c("Gene_ID", "H9_TR", "B7_TS", "G11_TR", "C11_TS")
summary(merged_mat$H9_TR)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.0000   0.0000   0.7899   7.5820  10.0571 546.8181 
summary(merged_mat$G11_TR)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.0000   0.1513   2.0805  10.5770  14.4838 804.1535 
summary(merged_mat$B7_TS)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.0000   0.1937   2.1303   8.2630  11.1122 596.6663 
summary(merged_mat$C11_TS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   0.123   2.214   9.121  12.243 655.496

#lets first say the 'basically silent' has a score lower than 2. Prob fairly consistent but pay attention anyway
Silent=merged_mat %>% filter(C11_TS<2)
#7285, about half our genes.
not_silent=merged_mat %>% filter(C11_TS>2)
#7671
summary(not_silent$C11_TS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.030   5.043  11.883  17.371  22.086 655.496 
mat_2_2_10 <- not_silent %>% filter(C11_TS<10                               )
#3389
mat_10_2_20 <- not_silent %>% filter(C11_TS>10 & C11_TS<20)
#2075
mat_20_2_50 <- not_silent %>% filter(C11_TS>20 & C11_TS<50)
#1846
mat_50_n_up <- not_silent %>% filter(C11_TS>50)
#361
#########################################
#those will be my categories, then.
refseq_c = rbind(refseq_curated_longest_cgi_plus, refseq_curated_longest_cgi_minus)
merged_mat=as.data.frame(merged_mat)
refseq_c2 = merge(x=refseq_c, y=merged_mat, by.y="Gene_ID", by.x="name")

merged_mat <- read.csv("/Users/kayle/Box/Vertinolab/McKayla Ford/Data/merged_mat.txt.gz", sep="")

#ugh of course I used a slightly dif gene ID for that... ugh.
refseq_c2$gene_ID = sub("_[^.]*$", "", refseq_c2$name)
refseq_c2$gene_ID = ifelse(refseq_c2$strand=="+",
                           paste0(refseq_c2$gene_ID, "_",
                                  refseq_c2$cpg_e),
                           paste0(refseq_c2$gene_ID,
                                  "_", refseq_c2$cpg_s))
refseq_c3 = merge(x=refseq_c2, y=merged_mat, by="gene_ID", all=TRUE)
##########################################################
refseq_c3$Ts_class = "silent"
refseq_c3$Ts_class[refseq_c3$C11_TS>2] = "sig_02-10"
refseq_c3$Ts_class[refseq_c3$C11_TS>10] = "sig_10-20"
refseq_c3$Ts_class[refseq_c3$C11_TS>20] = "sig_20-50"
refseq_c3$Ts_class[refseq_c3$C11_TS>50] = "sig_50+"
mm = refseq_c3 %>%
  group_by(gene_ID) %>%
  mutate(tot_H3K4me3 = sum(H3K4me3_EV)) %>%
  ungroup()


h1 <- ggplot(data=mm, mapping=aes(
  x=Bins, y=reorder(gene_ID, tot_H3K4me3),fill=H3K4me3_EV)) + geom_raster() +
  scale_fill_gradient(name="H3K4me3_EV", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,15)) +
  theme(axis.text.y = element_blank()) +
  facet_wrap(~Ts_class, scales="free")
h1

mm = refseq_c3 %>%
  group_by(gene_ID) %>%
  mutate(tot_H2AZ = sum(H2AZ)) %>%
  ungroup()

h1 <- ggplot(data=mm, mapping=aes(
  x=Bins, y=reorder(gene_ID, tot_H2AZ),fill=H2AZ)) + geom_raster() +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,2)) +
  theme(axis.text.y = element_blank()) +
  facet_wrap(~Ts_class, scales="free")
h1

mm = refseq_c3 %>%
  group_by(gene_ID) %>%
  mutate(tot_KDM5B = sum(KDM5B_EV)) %>%
  ungroup()

h1 <- ggplot(data=mm, mapping=aes(
  x=Bins, y=reorder(gene_ID, tot_KDM5B),fill=KDM5B_EV)) + geom_raster() +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,2)) +
  theme(axis.text.y = element_blank()) +
  facet_wrap(~Ts_class, scales="free")
h1

write.table(refseq_c3, "merged_mat_w_labs.txt")

merged_mat_w_labs <- read.csv("~/merged_mat_w_labs.txt", sep="")
################################################################################
#skew has to be a bw
library(rtracklayer)
skew_matrix <- function(bed, bw, b=NA, a=NA, n=NA, method, bs=10){
 chrs = list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
             "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
             "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
             "chr19", "chr20", "chr21", "chr22", "chrX")
  tmp_list=list()
  for (i in seq_along(chrs)) {
    bed_sub = bed %>% filter(bed[[1]]==chrs[[i]])
    if (method=="peak_centered") {
      ref_point = round((bed_sub[[3]]-bed_sub[[2]])/2)+bed_sub[[2]] #region center
      bins = seq(from=b, to=(a-bs), by=bs)
      hist_rows = bed_sub[[4]]
      hist = matrix(bins, ncol=length(bins), nrow=length(bed_sub[[4]]),
                    dimnames=list(hist_rows, bins), byrow = TRUE)
      hist=hist+ref_point #adds center to every bin to get matrix of start pos.
    }
    if (method=="single_stranded_anchored") {
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
    hist[hist <= 0] = NA
    print(paste(chrs[[i]], Sys.time()))
    partial_bw = import(bw, selection=GenomicSelection("hg38",
                                                       chrom=chrs[[i]],
                                                       colnames="score"))
    long_hist=reshape2::melt(hist, na.rm=TRUE) #longform lets us generate 'all
    #bin starts' and 'all bin ends' vectors for scoring.
    if (method=="bi_stranded_anchored"){
      long_ends=reshape2::melt(ends_mat, na.rm=TRUE)
      test=GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start=as.vector(long_hist[[3]]), end=as.vector(long_ends[[3]]))) #may not need the as.vector, troubleshooting relic
      olap = findOverlaps(test, partial_bw, ignore.strand=TRUE) 
    } #does count overlaps still work? no...
    else {
      olap = findOverlaps(GRanges( seqnames = chrs[[i]], ranges = IRanges(
        start = long_hist[[3]], end = long_hist[[3]]+ bs)), partial_bw, ignore.strand=TRUE)
    }
     bw_df = data.frame(partial_bw[subjectHits(olap)])
     names(bw_df)=paste("bw_", names(bw_df), sep="")
     regions_df = data.frame(long_hist[queryHits(olap),])
     t=list(bw_df, long_hist, olap, regions_df)
    regions_w_raw_scores = cbind(bw_df, regions_df)
    regions_w_raw_scores$adjustor = abs(regions_w_raw_scores$value-regions_w_raw_scores$bw_start)/10
    regions_w_raw_scores$score=regions_w_raw_scores$bw_score*regions_w_raw_scores$adjustor
      setDT(regions_w_raw_scores) 
      tmp1 = regions_w_raw_scores[c(7,8,11)]
      tmp2 = regions_w_raw_scores[, mean(score), by=.(Var1, Var2)]
      tmp3 = acast(tmp2, Var1~Var2)
      tmp_list[[i]] <- tmp3
  }
  hist = do.call(rbind, tmp_list) #collapse into one matrix
  hist=as.data.table(hist, keep.rownames=TRUE)
}

t = skew_matrix(TSS_CpG_end_plus, "/Users/kayle/Box/Vertinolab/McKayla Ford/Data/skew/track_files/plus_strand_skew_10.bw",
                b=-1000, a=1000, method="single_stranded_anchored")
mat_melt = data.table::melt(t, measure.vars=2:201,
                            variable.name="Bins", value.name = "Score")
mm2 = mat_melt %>%
  group_by(rn) %>%
  mutate(tot_skew = sum(Score)) %>%
  ungroup()
mm2$Bins=as.numeric(mm2$Bins)


h1 <- ggplot(data=mm2, mapping=aes(
  x=Bins, y=reorder(rn, tot_skew),fill=Score)) + geom_raster() +
  scale_fill_gradient2(name="H3K4me3_EV", high = "#A90E02",
                      low = "#012345", midpoint = 0, mid = "#FFFFFF", 
                      oob=scales::squish, limits=c(-0.25,0.25)) +
  theme(axis.text.y = element_blank())
h1
