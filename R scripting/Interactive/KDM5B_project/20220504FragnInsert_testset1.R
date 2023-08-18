#so what was the goal of this script? 
# I took the CGI TSSes, and metaplotted the 1 kb around the TSS
#for H3K4me3, H2AZ, and KDM5B with and without insert, with a filter applied
#to insert size, and with it resized to just the center 50 reads.

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

test3 = get_score_matrix(bam=paste0(    ###added in line "bam_aln = bam_aln[width(bam_aln)>=100]"
  bamDirectory1, "ZS1_19_MCF7_EV_H3K4me3.dedup.bam"
),
bed=TSS_CpG_end,
b=-1000,
a=1000,
method="single_stranded_anchored",
debug=TRUE,
readsOnly = FALSE)

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
colnames(mergedtests)=c("gene_ID", "Bins", "noinsert", "insert", "insert>=100",
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

svg(file="test_no_Y_Scaling2.svg")
geom
dev.off()


summary(mean_mergedtests$noinsert)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.112   2.373   3.188   3.445   4.756   5.771 
mean_mergedtests$noinsert=mean_mergedtests$noinsert/5.771

summary(mean_mergedtests$insert)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.840   3.870   4.933   5.428   7.633   8.714  
mean_mergedtests$insert=mean_mergedtests$insert/8.714

summary(mean_mergedtests$`insert>=100`)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.793   3.773   4.703   5.257   7.400   8.480
mean_mergedtests$`insert>=100`=mean_mergedtests$`insert>=100`/8.480

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


tests_melt = data.table::melt(mean_mergedtests, measure.vars=2:7,
                              variable.name="Test", value.name = "Score")

geom = ggplot(data=tests_melt, aes(x=Bins, y=Score, color=Test)) +
  mytheme +
  geom_line(size=1.5) +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed")
geom

svg(file="test_with_Y_Scaling2.svg")
geom
dev.off()

tests_melt2 = mean_mergedtests[,c(1:3,6:7)]
tests_melt3 = data.table::melt(tests_melt2, measure.vars=2:5,
                               variable.name="Test", value.name = "Score")
geom = ggplot(data=tests_melt3, aes(x=Bins, y=Score, color=Test)) +
  mytheme +
  geom_line(size=1.5) +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(n.breaks = 11)
geom


svg(file="test_with_Y_Scaling3.svg")
geom
dev.off()
