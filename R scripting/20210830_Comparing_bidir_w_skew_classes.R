### Load in modules ###
library(tidyverse)

### read in files ###
Dist_skew <- read.delim(
  "qnd_Dist_skew.bed",
  header=FALSE,
  quote="")
Prox_skew <- read.delim(
  "qnd_TSS_skew.bed",
  header=FALSE,
  quote="")
#5813, 5850 (11663)

neg_neu <- read.delim(
  "neg_neu.bed",
  header=FALSE,
  quote="")
Dist_signif <- read.delim(
  "dist_signif.bed",
  header=FALSE,
  quote="")
Prox_signif <- read.delim(
  "TSS_signif.bed",
  header=FALSE,
  quote="")
#3468, 4225, 4365 (12058) #I need to figure out why something makes it into this bed that doesn't make it into the
#other bed... 


unidir <- read.delim(
  "bidirectional islands/unidir.bed",
  header=FALSE,
  quote="")
cbi <- read.delim(
  "bidirectional islands/close_bidir.bed",
  header=FALSE,
  quote="")
fbi <- read.delim(
  "bidirectional islands/far_bidir.bed",
  header=FALSE,
  quote="")
#14851 total genes.

#1. add labels of what each thing is to it.
unidir$isle_class <- rep_len(c("unidir"), dim(unidir)[1])
fbi$isle_class <- rep_len(c("far_bi"), dim(fbi)[1])
cbi$isle_class <- rep_len(c("close_bi"), dim(cbi)[1])
Dist_skew$skew_class <- rep_len(c("prox"), dim(Dist_skew)[1])
Prox_skew$skew_class <- rep_len(c("dist"), dim(Prox_skew)[1])
neg_neu$skew_class <- rep_len(c("neg_neu"), dim(neg_neu)[1])
Dist_signif$skew_class <- rep_len(c("prox"), dim(Dist_signif)[1])
Prox_signif$skew_class <- rep_len(c("dist"), dim(Prox_signif)[1])

#2. get what the labels of the other thing would be.
bidir <- rbind(unidir, fbi, cbi)
skew_3way <- rbind(Dist_signif, Prox_signif, neg_neu)
skew_2way <- rbind(Dist_skew, Prox_skew)

tmp_lookup = bidir$isle_class
names(tmp_lookup)=bidir$V4
skew_3way$isle_class=unname(tmp_lookup[skew_3way$V4])
skew_2way$isle_class=unname(tmp_lookup[skew_2way$V4])

ggplot(data=skew_3way) + geom_bar(aes(x=skew_class, fill=isle_class), position="dodge") +theme_bw()
ggplot(data=skew_3way) + geom_bar(aes(x=isle_class, fill=skew_class), position="dodge") +theme_bw()
ggplot(data=skew_2way) + geom_bar(aes(x=skew_class, fill=isle_class), position="dodge") +theme_bw()
ggplot(data=skew_2way) + geom_bar(aes(x=isle_class, fill=skew_class), position="dodge") +theme_bw()
for_testing <- skew_3way %>% count(skew_class, isle_class)

dis_row <- c(420, 556, 3241, 148)
prox_row <- c(473, 480, 3173, 161)
neg_neu_row <-c(429, 341, 2555, 143)
mymat <- as.matrix(rbind(dis_row, prox_row, neg_neu_row))
chisq <- chisq.test(mymat)
chisq
chisq$expected
chisq$observed
chisq$residuals
#Neg_neutral is more likely to be associated with close-bidirectional islands,
#while distal is most likely to be associated with far bidirectional islands. There is also a small association
#between neg_neutral and NA, but I'm not sure what that means.

tst <- skew_3way %>% filter(is.na(isle_class))
