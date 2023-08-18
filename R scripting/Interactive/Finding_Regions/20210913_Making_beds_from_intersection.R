###################
## Initial setup ##
###################
library(dplyr)
intersection_table <- read.delim( #Reading in the CpG/Genes intersection file, from bedtools
  "Refseq_curated_CGI_promoters_filter.hg38.txt",
  header=FALSE,
  quote="")
#16108 genes.
#############################################3
plus_strand_only <- intersection_table %>%
  filter(V9=="+")
#8114 genes

remove_outliers_p <- plus_strand_only %>% filter (
  (V3-V5)<3000 & (V3-V5) > 20
) #removes longer than 3k and shorter than 20 distance btwn TSS and end of CpG island.
#7969 genes

TSS_cpg_end_p <- select(remove_outliers_p, V1, V5, V3, V7, V8, V9)
#selects only the columns used in a bed6 file. Here we want our start to be TSS and end to be end
#of CpG island, but we could also make gene files by selecting gene end or CpG files w/ CpG start.

TSS_cpg_end_p <- TSS_cpg_end_p %>%
  rename(
    start=`V5`,
    end=`V3`
  )
# 
# write.table(
#   TSS_cpg_end_p,
#   file="TSS_CPG_end_p.hg38.bed",
#   row.names=FALSE,
#   col.names=FALSE,
#   quote=FALSE,
#   sep="\t"
# )
###############################################

minus_strand_only <- intersection_table %>%
  filter(V9=="-")
#7994 genes
remove_outliers_m <- minus_strand_only %>% filter (
  (V6-V2)<3000 & (V6-V2) > 20
)
#7860 genes

TSS_cpg_end_m <- select(remove_outliers_m, V1, V2, V6, V7, V8, V9)
TSS_cpg_end_m <- TSS_cpg_end_m %>%
  rename(
    start=`V2`,
    end=`V6`
  )
# 
# write.table(
#   TSS_cpg_end_m,
#   file="TSS_CPG_end_m.hg38.bed",
#   row.names=FALSE,
#   col.names=FALSE,
#   quote=FALSE,
#   sep="\t"
# )
###################################3
merged_bed6 <- rbind.data.frame(
  TSS_cpg_end_p,
  TSS_cpg_end_m
)

write.table(
  merged_bed6,
  file="TSS_CPG_end_curated.hg38.bed",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)

######################################
#windows
# ##########################################
plus_strand_only <- intersection_table %>%
filter(V9=="+")
#8114 genes
#summary stats for CpG start to TES are
summary(plus_strand_only$V6-plus_strand_only$V2)
#basically median is 24k, 3rd quartile is 59k, and max is over 2 million. But that's okay
#we're just going 1000 before and 3k after, no more.
t1 <- plus_strand_only %>% filter(V6 > V3)
t2 <- plus_strand_only %>% filter(V6 < V3)
t1$win_end <- t1$V6 + 3000
t2$win_end <- t2$V3 + 3000
t3 <- rbind(t1, t2)
t3$win_start <- t3$V5-3000
winp_1 <- select(t3, V1, win_start, win_end, V7, V8, V9)

win_end <- plus_strand_only$V6+3000
win_end2 <- plus_strand_only$V3 + 3000 #test
winp_1 <- select(plus_strand_only, V1, V7, V8, V9)
winp_2 <- cbind.data.frame(winp_1, win_start, win_end, win_end2) #test variant
tst1 <- winp_2 %>% filter(win_end > win_end2)
tst2 <- winp_2 %>% filter(win_end2 > win_end)
winp_3 <- winp_2 %>% relocate(win_start, win_end, .after = V1)
winp_4 <- winp_1 %>%
  filter(!grepl('chrY|([\\w_]+)alt|random|fix', V1))
#filters around 600 out.
write.table(
  winp_4,
  file="plus_win_for_skew.bed",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)

minus_strand_only <- intersection_table %>%
  filter(V9=="-")
summary(minus_strand_only$V3-minus_strand_only$V5)
minus_strand_only$win_end <- minus_strand_only$V3+3000
t1 <- minus_strand_only %>% filter(V5 > V2)
t2 <- minus_strand_only %>% filter(V5 < V2)
t1$win_start <- t1$V2 - 3000
t2$win_start <- t2$V5 - 3000
t3 <- rbind(t1, t2)
winm_1 <- select(t3, V1, win_start, win_end, V7, V8, V9)

win_start <- minus_strand_only$V5-3000
winm_1 <- select(minus_strand_only, V1, V7, V8, V9)
winm_2 <- cbind.data.frame(winm_1, win_start, win_end)
winm_3 <- winm_2 %>% relocate(win_start, win_end, .after = V1)
winm_4 <- winm_1 %>%
  filter(!grepl('chrY|([\\w_]+)alt|random|fix', V1))
#filters around 600 out.
write.table(
  winm_4,
  file="minus_win_for_skew.bed",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)

