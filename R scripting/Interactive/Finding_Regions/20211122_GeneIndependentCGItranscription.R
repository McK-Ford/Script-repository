library(tidyverse)
library(data.table)

ptabX <- read.delim(
  "intersection_files/11162021_cpgs_proseq_noanno_intersect_plus_chrX.txt",
  header=FALSE,
  quote=""
)
p_chrX_dedup1 <- ptabX %>%
  group_by(V2) %>%
  summarise(
    chrom=first(V1),
    end=first(V3),
    uniq_ID=first(V4),
    score=sum(V8)
  ) #this lets us get the total score in each window.
#Okay now that I've tested with plus X here's the loop for remaining chromosomes
chromlist = paste0("chr", seq(from=1, to=22))
out_list = list()
for (i in seq_along(chromlist)) {
  filename = paste0(
    "intersection_files/11162021_cpgs_proseq_noanno_intersect_plus_",
    chromlist[[i]],
    ".txt"
  )
  file = read.delim(
    filename,
    header=FALSE,
    quote=""
  )
  p_chr_dedup1 <- file %>%
    group_by(V2) %>%
    summarise(
      chrom=first(V1),
      end=first(V3),
      uniq_ID=first(V4),
      score=sum(V8)
    )
  out_list[[i]] = p_chr_dedup1
}
plus_vals_0 = do.call(rbind, out_list)
plus_vals = rbind(plus_vals_0, p_chrX_dedup1)
plus_vals=unique(plus_vals)
summary(plus_vals$score)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 0
#1.000    1.000    3.000    7.283    6.000 5855.000 

mtabX <- read.delim(
  "intersection_files/11162021_cpgs_proseq_noanno_intersect_minus_chrX.txt",
  header=FALSE,
  quote=""
)

m_chrX_dedup1 <- mtabX %>%
  group_by(V2) %>%
  summarise(
    chrom=first(V1),
    end=first(V3),
    uniq_ID=first(V4),
    score=sum(V8)
  )

out_list = list()
for (i in seq_along(chromlist)) {
  filename = paste0(
    "intersection_files/11162021_cpgs_proseq_noanno_intersect_minus_",
    chromlist[[i]],
    ".txt"
  )
  file = read.delim(
    filename,
    header=FALSE,
    quote=""
  )
  m_chr_dedup1 <- file %>%
    group_by(V2) %>%
    summarise(
      chrom=first(V1),
      end=first(V3),
      uniq_ID=first(V4),
      score=sum(V8)
    )
  out_list[[i]] = m_chr_dedup1
}
minus_vals_0 = do.call(rbind, out_list)
minus_vals = rbind(minus_vals_0, m_chrX_dedup1)

p_vals_dedup_maxscore <- plus_vals %>%
  ungroup() %>%
  group_by(uniq_ID) %>%
  mutate(score_max=max(score)) %>%
  filter(score==score_max) %>%
  summarize(
    chrom=first(chrom),
    start=first(V2),
    end=max(end),
    score=first(score)
  )
summary(p_vals_dedup_maxscore$score)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    3.00   11.00   48.63   45.00 5855.00 
p_vals_dedup_5prime <- plus_vals %>%
  ungroup() %>%
  filter(score>3) %>%
  group_by(uniq_ID) %>%
  mutate(fiveprime =min(V2)) %>%
  filter(V2==fiveprime)
summary(p_vals_dedup_5prime$score)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4.000   4.000   5.000   5.861   6.000 360.000

m_vals_dedup_maxscore <- minus_vals %>%
  ungroup() %>%
  group_by(uniq_ID) %>%
  mutate(score_max=max(score)) %>%
  filter(score==score_max) %>%
  summarize(
    chrom=first(chrom),
    start=first(V2),
    end=max(end),
    score=first(score)
  )
summary(m_vals_dedup_maxscore$score)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    3.00   11.00   51.69   46.00 4740.00 
m_vals_dedup_5prime <- minus_vals %>%
  ungroup() %>%
  filter(score>3) %>%
  group_by(uniq_ID) %>%
  mutate(fiveprime =max(V2)) %>%
  filter(V2==fiveprime)
summary(m_vals_dedup_5prime$score)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4.000   4.000   4.000   5.294   6.000 82.000
save(p_vals_dedup_5prime, p_vals_dedup_maxscore, m_vals_dedup_5prime, m_vals_dedup_maxscore, file="cgi_tss_noanno_intersect_deduped.RData")
load(file="cgi_tss_noanno_intersect_deduped.RData")

#get beds ... oh these aren't matched up with islands of course :(
cpgIslands <- read.delim("~/Lab stuff/R_coding/backedup/cpgIsland.hg38.txt", header=FALSE, quote="")
no_alt_chroms <- cpgIslands %>%
  filter(!grepl('chrY|([\\w_]+)alt|random|fix|v1', V1))
no_alt_chroms$uniq_ID <- paste0(no_alt_chroms$V4, no_alt_chroms$V2)

ltTruStart = no_alt_chroms$V2
ltTruEnd = no_alt_chroms$V3
ltTruID = no_alt_chroms$V4

names(ltTruStart) = no_alt_chroms$uniq_ID
names(ltTruEnd) = no_alt_chroms$uniq_ID
names(ltTruID) = no_alt_chroms$uniq_ID

#wait do I even need the islands? keep this code ofc

mvals5_bed = m_vals_dedup_5prime %>% select(c(chrom, V2, end, uniq_ID, score))
pvals5_bed = p_vals_dedup_5prime %>% select(c(chrom, V2, end, uniq_ID, score))
mvalsmax_bed = m_vals_dedup_maxscore %>% select(c(chrom, start, end, uniq_ID, score))
pvalsmax_bed = p_vals_dedup_maxscore %>% select(c(chrom, start, end, uniq_ID, score))

mvals5_bed$strand = rep("-", dim(mvals5_bed)[1])
pvals5_bed$strand = rep("+", dim(pvals5_bed)[1])
mvalsmax_bed$strand = rep("-", dim(mvalsmax_bed)[1])
pvalsmax_bed$strand = rep("+", dim(pvalsmax_bed)[1])

write.table(
  mvals5_bed,
  file="mvals5.bed",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)
write.table(
  pvals5_bed,
  file="pvals5.bed",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)
write.table(
  mvalsmax_bed,
  file="mvalsmax.bed",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)
write.table(
  pvalsmax_bed,
  file="pvalsmax.bed",
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t"
)


testing_max = merge(pvalsmax_bed, mvalsmax_bed, by="uniq_ID", all=TRUE)
testing_5 = merge(pvals5_bed, mvals5_bed, by="uniq_ID", all=TRUE)
#total cgi are 31k
summary(testing_max$start.x - testing_max$start.y)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
#-6944.000  -214.000    86.000    -6.249   191.000 11489.000      6631 
#length is 20178
summary(testing_max$score.x - testing_max$score.y)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
#-4615.000   -23.000     0.000    -3.858    20.000  5376.000      6631 
summary(testing_5$V2.x - testing_5$V2.y) #V2 is start
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#-12977.0   -734.0   -412.0   -507.3    -97.0   2005.0     5365 
#length is 15147
summary(testing_5$score.x - testing_5$score.y)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-75.000  -1.000   0.000   0.518   1.000 216.000    5365

unidir_5 = testing_5 %>% filter(is.na(chrom.x)|is.na(chrom.y))
#is 5365 observations
unidir_max = testing_max %>% filter(is.na(chrom.x)|is.na(chrom.y))
#is 6631 obs
bidir_5 = testing_5 %>% filter(!is.na(chrom.x) & !is.na(chrom.y))
#is 9782 obs (bc less than 3 score is filtered out)
bidir_max = testing_max %>% filter(!is.na(chrom.x) & !is.na(chrom.y))
#is 13547 obs
smoothScatter(x=bidir_5$score.x, y=bidir_5$score.y, xlim=c(0,100), ylim=c(0,100))
smoothScatter(x=bidir_max$score.x, y=bidir_max$score.y, xlim=c(0,1000), ylim=c(0,1000))
OverLapping = bidir_max %>% filter((bidir_max$start.x - bidir_max$start.y)<=0)
# is 5428 obs
summary(OverLapping$score.x - OverLapping$score.y)
NotOverLapping = bidir_max %>% filter((bidir_max$start.x - bidir_max$start.y)>0)
# is 8119 obs
summary(NotOverLapping$score.x - NotOverLapping$score.y)
smoothScatter(x=OverLapping$score.x, y=OverLapping$score.y, xlim=c(0,1000), ylim=c(0,1000))
smoothScatter(x=OverLapping$score.x, y=OverLapping$score.y, xlim=c(0,1000), ylim=c(0,1000))
