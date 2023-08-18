source("/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/NewMetaplotterJoshBased.R")
enableJIT(3)
##############
## plotting ##
##############
bamDirectory="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/"
refseq_curated_longest_cgi_minus <- read.delim(paste0(
  bamDirectory,
  "refseq_curated_longest_cgi_minus.hg38.txt"))
refseq_curated_longest_cgi_plus <- read.delim(paste0(
  bamDirectory,
  "refseq_curated_longest_cgi_plus.hg38.txt"))

TSS_CpG_end_plus= refseq_curated_longest_cgi_plus %>%
  select(1,5,3,8,2,7) #using CpG start as score placeholder
TSS_CpG_end_plus$name=paste0(TSS_CpG_end_plus$name, "_", TSS_CpG_end_plus$cpg_e)
TSS_CpG_end_minus= refseq_curated_longest_cgi_minus %>% select(1,2,6,8,3,7) #cgi end placeholder
TSS_CpG_end_minus$name=paste0(TSS_CpG_end_minus$name, "_", TSS_CpG_end_minus$cpg_s)

name_vec = c("chrom", "start", "end", "name", "score", "strand")
colnames(TSS_CpG_end_minus) = name_vec
colnames(TSS_CpG_end_plus) = name_vec
TSS_CpG_end = rbind(TSS_CpG_end_minus, TSS_CpG_end_plus)
bamDirectory1="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/bams/hg38_dedup/"

test1 = get_score_matrix(bam=paste0(
  bamDirectory1, "ZS1_19_MCF7_EV_H3K4me3.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE)

test2 = get_score_matrix(bam=paste0(
  bamDirectory1, "ZS1_19_MCF7_EV_H3K4me3.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE,
readsOnly = FALSE)

test3 = get_score_matrix(bam=paste0(    ###added in line "bam_aln = bam_aln[width(bam_aln)>=200]"
  bamDirectory1, "ZS1_19_MCF7_EV_H3K4me3.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE)

test4 = get_score_matrix(bam=paste0(    ###added in line "bam_aln = bam_aln[width(bam_aln)>=200]"
  bamDirectory1, "ZS1_19_MCF7_EV_H3K4me3.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE,
readsOnly = FALSE)

test5 = get_score_matrix(bam=paste0(    ###added in line "bam_aln=resize(bam_aln, width=50, fix="center")"
  bamDirectory1, "ZS1_19_MCF7_EV_H3K4me3.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE,
readsOnly = FALSE)

test6 = get_score_matrix(bam=paste0(    ###added in lines "bam_aln = bam_aln[width(bam_aln)>=250]" and "bam_aln=resize(bam_aln, width=50, fix="center")"
  bamDirectory1, "ZS1_19_MCF7_EV_H3K4me3.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE,
readsOnly = FALSE)

testlist=list(test1,test2,test3,test4,test5, test6)
testlist2 <- lapply(testlist, data.table::melt, measure.vars=c(2:201),
                    variable.name="Bins", value.name="Score")
mergedtests=cbind(testlist2[[1]], testlist2[[2]][[3]], testlist2[[3]][[3]],
                  testlist2[[4]][[3]], testlist2[[5]][[3]], testlist2[[6]][[3]])
colnames(mergedtests)=c("gene_ID", "Bins", "noinsert", "insert", "noinsert>=200",
                        "insert>=200", "center_50", "center_50_insert>=200"
)

mergedtests[,gene_ID:=NULL]
mean_mergedtests = mergedtests[,lapply(.SD, mean), by=Bins]
mean_mergedtests$Bins = as.numeric(mean_mergedtests$Bins)*10-1000
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

tests_melt = data.table::melt(mean_mergedtests, measure.vars=2:7,
                              variable.name="Test", value.name = "Score")

geom = ggplot(data=tests_melt, aes(x=Bins, y=Score, color=Test)) +
  mytheme +
  geom_line(size=1.5) +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed")
geom

summary(mean_mergedtests$noinsert)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.112   2.373   3.188   3.445   4.756   5.771 
mean_mergedtests$noinsert=mean_mergedtests$noinsert/5.771

summary(mean_mergedtests$insert)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.840   3.870   4.933   5.428   7.633   8.714  
mean_mergedtests$insert=mean_mergedtests$insert/8.714

summary(mean_mergedtests$`insert>=200`)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.462   3.113   3.846   4.255   5.893   6.853 
mean_mergedtests$`insert>=200`=mean_mergedtests$`insert>=200`/6.853

summary(mean_mergedtests$center_50)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.4115  0.8684  1.2119  1.2601  1.7202  2.3769 
mean_mergedtests$center_50=mean_mergedtests$center_50/2.3769

summary(mean_mergedtests$`center_50_insert>=200`)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2583  0.5521  0.7082  0.7825  1.0924  1.4899 
mean_mergedtests$`center_50_insert>=200`=mean_mergedtests$`center_50_insert>=200`/1.4899

mean_mergedtests=mean_mergedtests[,c(1:3, 5:7)]

tests_melt = data.table::melt(mean_mergedtests, measure.vars=2:6,
                              variable.name="Test", value.name = "Score")

geom = ggplot(data=tests_melt, aes(x=Bins, y=Score, color=Test)) +
  mytheme +
  geom_line(size=1.5) +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed")
geom

svg(file="test_with_Y_Scaling_H3K4me3.svg")
geom
dev.off()

#################################################

test1 = get_score_matrix(bam=paste0(
  bamDirectory1, "ZS1_19_MCF7_EV_H2AZ.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE)

test2 = get_score_matrix(bam=paste0(
  bamDirectory1, "ZS1_19_MCF7_EV_H2AZ.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE,
readsOnly = FALSE)

test3 = get_score_matrix(bam=paste0(    ###added in line "bam_aln = bam_aln[width(bam_aln)>=200]"
  bamDirectory1, "ZS1_19_MCF7_EV_H2AZ.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE)

test4 = get_score_matrix(bam=paste0(    ###added in line "bam_aln = bam_aln[width(bam_aln)>=200]"
  bamDirectory1, "ZS1_19_MCF7_EV_H2AZ.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE,
readsOnly = FALSE)

test5 = get_score_matrix(bam=paste0(    ###added in line "bam_aln=resize(bam_aln, width=50, fix="center")"
  bamDirectory1, "ZS1_19_MCF7_EV_H2AZ.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE,
readsOnly = FALSE)

test6 = get_score_matrix(bam=paste0(    ###added in lines "bam_aln = bam_aln[width(bam_aln)>=250]" and "bam_aln=resize(bam_aln, width=50, fix="center")"
  bamDirectory1, "ZS1_19_MCF7_EV_H2AZ.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE,
readsOnly = FALSE)

testlist=list(test1,test2,test3,test4,test5, test6)
testlist2 <- lapply(testlist, data.table::melt, measure.vars=c(2:201),
                    variable.name="Bins", value.name="Score")
mergedtests=cbind(testlist2[[1]], testlist2[[2]][[3]], testlist2[[3]][[3]],
                  testlist2[[4]][[3]], testlist2[[5]][[3]], testlist2[[6]][[3]])
colnames(mergedtests)=c("gene_ID", "Bins", "noinsert", "insert", "noinsert>=200",
                        "insert>=200", "center_50", "center_50_insert>=200"
)

mergedtests[,gene_ID:=NULL]
mean_mergedtests = mergedtests[,lapply(.SD, mean), by=Bins]
mean_mergedtests$Bins = as.numeric(mean_mergedtests$Bins)*10-1000

tests_melt = data.table::melt(mean_mergedtests, measure.vars=2:7,
                              variable.name="Test", value.name = "Score")

geom = ggplot(data=tests_melt, aes(x=Bins, y=Score, color=Test)) +
  mytheme +
  geom_line(size=1.5) +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed")
geom

summary(mean_mergedtests$`insert>=200`)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.3337  0.4929  0.6293  0.6229  0.7728  0.8703 
mean_mergedtests$`insert>=200`=mean_mergedtests$`insert>=200`/0.8703

summary(mean_mergedtests$noinsert)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.2953  0.4719  0.6256  0.6373  0.7709  1.0931 
mean_mergedtests$noinsert=mean_mergedtests$noinsert/1.0931

summary(mean_mergedtests$insert)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.4526  0.6783  0.8925  0.8914  1.1011  1.2968  
mean_mergedtests$insert=mean_mergedtests$insert/1.2968

summary(mean_mergedtests$center_50)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.1112  0.1609  0.2297  0.2346  0.3032  0.4108 
mean_mergedtests$center_50=mean_mergedtests$center_50/0.4108

summary(mean_mergedtests$`center_50_insert>=200`)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.06241 0.08878 0.11779 0.11951 0.14784 0.17916 
mean_mergedtests$`center_50_insert>=200`=mean_mergedtests$`center_50_insert>=200`/0.17916

mean_mergedtests=mean_mergedtests[,c(1:3, 5:7)]

tests_melt = data.table::melt(mean_mergedtests, measure.vars=2:6,
                              variable.name="Test", value.name = "Score")

geom = ggplot(data=tests_melt, aes(x=Bins, y=Score, color=Test)) +
  mytheme +
  geom_line(size=1.5) +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed")
geom

svg(file="test_with_Y_Scaling_H2AZ.svg")
geom
dev.off()

########################################################

test1 = get_score_matrix(bam=paste0(
  bamDirectory1, "ZS1_19_MCF7_EV_KDM5B.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE)

test2 = get_score_matrix(bam=paste0(
  bamDirectory1, "ZS1_19_MCF7_EV_KDM5B.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE,
readsOnly = FALSE)

test3 = get_score_matrix(bam=paste0(    ###added in line "bam_aln = bam_aln[width(bam_aln)>=200]"
  bamDirectory1, "ZS1_19_MCF7_EV_KDM5B.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE)

test4 = get_score_matrix(bam=paste0(    ###added in line "bam_aln = bam_aln[width(bam_aln)>=200]"
  bamDirectory1, "ZS1_19_MCF7_EV_KDM5B.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE,
readsOnly = FALSE)

test5 = get_score_matrix(bam=paste0(    ###added in line "bam_aln=resize(bam_aln, width=50, fix="center")"
  bamDirectory1, "ZS1_19_MCF7_EV_KDM5B.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE,
readsOnly = FALSE)

test6 = get_score_matrix(bam=paste0(    ###added in lines "bam_aln = bam_aln[width(bam_aln)>=250]" and "bam_aln=resize(bam_aln, width=50, fix="center")"
  bamDirectory1, "ZS1_19_MCF7_EV_KDM5B.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE,
readsOnly = FALSE)

testlist=list(test1,test2,test3,test4,test5, test6)
testlist2 <- lapply(testlist, data.table::melt, measure.vars=c(2:201),
                    variable.name="Bins", value.name="Score")
mergedtests=cbind(testlist2[[1]], testlist2[[2]][[3]], testlist2[[3]][[3]],
                  testlist2[[4]][[3]], testlist2[[5]][[3]], testlist2[[6]][[3]])
colnames(mergedtests)=c("gene_ID", "Bins", "noinsert", "insert", "noinsert>=200",
                        "insert>=200", "center_50", "center_50_insert>=200"
)

mergedtests[,gene_ID:=NULL]
mean_mergedtests = mergedtests[,lapply(.SD, mean), by=Bins]
mean_mergedtests$Bins = as.numeric(mean_mergedtests$Bins)*10-1000

tests_melt = data.table::melt(mean_mergedtests, measure.vars=2:7,
                              variable.name="Test", value.name = "Score")

geom = ggplot(data=tests_melt, aes(x=Bins, y=Score, color=Test)) +
  mytheme +
  geom_line(size=1.5) +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed")
geom

summary(mean_mergedtests$`insert>=200`)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.1921  0.2994  0.4018  0.4422  0.6011  0.7019  
mean_mergedtests$`insert>=200`=mean_mergedtests$`insert>=200`/0.7019

summary(mean_mergedtests$noinsert)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.1612  0.2524  0.3428  0.4435  0.5346  1.0975
mean_mergedtests$noinsert=mean_mergedtests$noinsert/1.0975

summary(mean_mergedtests$insert)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.2478  0.3839  0.5170  0.6026  0.8460  1.0643  
mean_mergedtests$insert=mean_mergedtests$insert/1.0643

summary(mean_mergedtests$center_50)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.05924 0.09465 0.13027 0.16637 0.22790 0.38696 
mean_mergedtests$center_50=mean_mergedtests$center_50/0.38696

summary(mean_mergedtests$`center_50_insert>=200`)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.03567 0.05450 0.06946 0.08173 0.10044 0.19525 
mean_mergedtests$`center_50_insert>=200`=mean_mergedtests$`center_50_insert>=200`/0.19525

mean_mergedtests=mean_mergedtests[,c(1:3, 5:7)]

tests_melt = data.table::melt(mean_mergedtests, measure.vars=2:6,
                              variable.name="Test", value.name = "Score")

geom = ggplot(data=tests_melt, aes(x=Bins, y=Score, color=Test)) +
  mytheme +
  geom_line(size=1.5) +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed")
geom

svg(file="test_with_Y_Scaling_KDM5B.svg")
geom
dev.off()
