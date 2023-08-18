library(tidyverse)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
H2AZ_driven_intersection <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/H2AZ_driven_intersection.txt",
  header=FALSE, na.strings='.')
KDM5B_driven_intersection <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/KDM5B_driven_intersection.txt",
  header=FALSE, na.strings='.')
H3K27me3_driven_intersection <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/H3K27me3_driven_intersection.txt",
  header=FALSE, na.strings='.')
H3K4me3_driven_intersection <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/H3K4me3_driven_intersection.txt",
  header=FALSE, na.strings='.')
##### get the data in ####
H2AZ_driven_intersection$V7[is.na(H2AZ_driven_intersection$V7)] = "placeholder"
H2AZ_single_peak_per_peak = H2AZ_driven_intersection %>%
  group_by(V6, V7) %>% mutate(max_col = max(V11)) %>%
  ungroup() %>% filter(V11==max_col)

collapsed_H2AZ <- data.frame(
  do.call(paste, c(H2AZ_single_peak_per_peak[1:6], sep=",")),
  H2AZ_single_peak_per_peak$V7,
  do.call(paste, c(H2AZ_single_peak_per_peak[8:13], sep=",")))
names(collapsed_H2AZ) <- c("col_H2AZ", "bed_id", "col_info")
collapsed_H2AZ$bed_id <- str_remove_all(
  collapsed_H2AZ$bed_id,
  "ZS1_19_MCF7_EV_|ZS1_19_MCF7_0025EV_|\\.stringent\\.bed")
#collapsed_H2AZ$bed_id[is.na(bed_id)] = "empty_placeholder"
wide_collapsed_H2AZ <- collapsed_H2AZ %>%
  pivot_wider(
    names_from = bed_id, values_from = col_info,
    names_prefix = "samp_")
######
KDM5B_driven_intersection$V7[is.na(KDM5B_driven_intersection$V7)] = "placeholder"
KDM5B_single_peak_per_peak = KDM5B_driven_intersection %>%
  group_by(V6, V7) %>% mutate(max_col = max(V11)) %>%
  ungroup() %>% filter(V11==max_col)
collapsed_KDM5B <- data.frame(
  do.call(paste, c(KDM5B_single_peak_per_peak[1:6], sep=",")),
  KDM5B_single_peak_per_peak$V7,
  do.call(paste, c(KDM5B_single_peak_per_peak[8:13], sep=",")))
names(collapsed_KDM5B) <- c("col_KDM5B", "bed_id", "col_info")
collapsed_KDM5B$bed_id <- str_remove_all(
  collapsed_KDM5B$bed_id,
  "ZS1_19_MCF7_EV_|ZS1_19_MCF7_0025EV_|\\.stringent\\.bed")
wide_collapsed_KDM5B <- collapsed_KDM5B %>%
  pivot_wider(
    names_from = bed_id, values_from = col_info,
    names_prefix = "samp_")
######
H3K27me3_driven_intersection$V7[is.na(H3K27me3_driven_intersection$V7)] = "placeholder"
H3K27me3_single_peak_per_peak = H3K27me3_driven_intersection %>%
  group_by(V6, V7) %>% mutate(max_col = max(V11)) %>%
  ungroup() %>% filter(V11==max_col)

collapsed_H3K27me3 <- data.frame(
  do.call(paste, c(H3K27me3_single_peak_per_peak[1:6], sep=",")),
  H3K27me3_single_peak_per_peak$V7,
  do.call(paste, c(H3K27me3_single_peak_per_peak[8:13], sep=",")))
names(collapsed_H3K27me3) <- c("col_H3K27me3", "bed_id", "col_info")
collapsed_H3K27me3$bed_id <- str_remove_all(
  collapsed_H3K27me3$bed_id,
  "ZS1_19_MCF7_EV_|ZS1_19_MCF7_0025EV_|\\.stringent\\.bed")
wide_collapsed_H3K27me3 <- collapsed_H3K27me3 %>%
  pivot_wider(
    names_from = bed_id, values_from = col_info,
    names_prefix = "samp_")
###
H3K4me3_driven_intersection$V7[is.na(H3K4me3_driven_intersection$V7)] = "placeholder"
H3K4me3_single_peak_per_peak = H3K4me3_driven_intersection %>%
  group_by(V6, V7) %>% mutate(max_col = max(V11)) %>%
  ungroup() %>% filter(V11==max_col)

collapsed_H3K4me3 <- data.frame(
  do.call(paste, c(H3K4me3_single_peak_per_peak[1:6], sep=",")),
  H3K4me3_single_peak_per_peak$V7,
  do.call(paste, c(H3K4me3_single_peak_per_peak[8:13], sep=",")))
names(collapsed_H3K4me3) <- c("col_H3K4me3", "bed_id", "col_info")
collapsed_H3K4me3$bed_id <- str_remove_all(
  collapsed_H3K4me3$bed_id,
  "ZS1_19_MCF7_EV_|ZS1_19_MCF7_0025EV_|\\.stringent\\.bed")
wide_collapsed_H3K4me3 <- collapsed_H3K4me3 %>%
  pivot_wider(
    names_from = bed_id, values_from = col_info,
    names_prefix = "samp_")

######
#now I need to arrange these all identically
wide_collapsed_KDM5B = wide_collapsed_KDM5B[c(1,2,3,4)]
wide_collapsed_H2AZ = wide_collapsed_H2AZ[c(3,1,2,4)]
wide_collapsed_H3K4me3 = wide_collapsed_H3K4me3[c(3,2,1,4)]
wide_collapsed_H3K27me3 = wide_collapsed_H3K27me3[c(3,2,4,1)]

names(wide_collapsed_KDM5B) = c("KDM5B", "H2AZ", "H3K4me3", "H3K27me3")
names(wide_collapsed_H2AZ) = c("KDM5B", "H2AZ", "H3K4me3", "H3K27me3")
names(wide_collapsed_H3K4me3) = c("KDM5B", "H2AZ", "H3K4me3", "H3K27me3")
names(wide_collapsed_H3K27me3) = c("KDM5B", "H2AZ", "H3K4me3", "H3K27me3")

intersection_stats = unique(rbind(
  wide_collapsed_KDM5B, wide_collapsed_H2AZ,
  wide_collapsed_H3K4me3, wide_collapsed_H3K27me3))
#51902 preunique, 11128 after once I fixed the scientific notation error ...
# ... that seems unreasonable.. No wait I wasn't looking at the wide-collapsed
#version before, that was silly

intersection_stats$KDM5B_stats = intersection_stats$KDM5B
intersection_stats$KDM5B_stats[is.na(intersection_stats$KDM5B)]=FALSE
intersection_stats$KDM5B_stats[!is.na(intersection_stats$KDM5B)]=TRUE

intersection_stats$H3K27me3_stats = intersection_stats$H3K27me3
intersection_stats$H3K27me3_stats[is.na(intersection_stats$H3K27me3)] = FALSE
intersection_stats$H3K27me3_stats[!is.na(intersection_stats$H3K27me3)] = TRUE

intersection_stats$H3K4me3_stats = intersection_stats$H3K4me3
intersection_stats$H3K4me3_stats[is.na(intersection_stats$H3K4me3)] = FALSE
intersection_stats$H3K4me3_stats[!is.na(intersection_stats$H3K4me3)] = TRUE

intersection_stats$H2AZ_stats = intersection_stats$H2AZ
intersection_stats$H2AZ_stats[is.na(intersection_stats$H2AZ)] = FALSE
intersection_stats$H2AZ_stats[!is.na(intersection_stats$H2AZ)] = TRUE

#####
intersection_stats = separate(
  intersection_stats, col=KDM5B,
  into=paste0("KDM5B_", c("chr", "start", "end", "tot_score", "max_score",
                          "max_reg")), sep=",", remove=TRUE, convert=TRUE)
intersection_stats = separate(
  intersection_stats, col=H3K4me3,
  into=paste0("H3K4me3_", c("chr", "start", "end","tot_score", "max_score",
                          "max_reg")), sep=",", remove=TRUE, convert=TRUE)
intersection_stats = separate(
  intersection_stats, col=H3K27me3,
  into=paste0("H3K27me3_", c("chr", "start", "end","tot_score", "max_score",
                          "max_reg")), sep=",", remove=TRUE, convert=TRUE)
intersection_stats = separate(
  intersection_stats, col=H2AZ,
  into=paste0("H2AZ_", c("chr", "start", "end","tot_score", "max_score",
                          "max_reg")), sep=",", remove=TRUE, convert=TRUE)
######
# Now I have my table (which I need to write to file).
#Steps to repeat from KDM5B peaks stats include:
#Annotate. Can do that with Ben's pipeline instead this time, and just define
#regions based on that info.
#get peakdifs
#summary of peakdifs
#additionally summary of categories
#don't forget CpG islands
#also I need to graph heatmaps at peaks not just at genes
#im flagging though, I need caffiene
intersection_stats$peak_id <- paste0("peak_", seq(1:dim(intersection_stats)[1]))
save_state = intersection_stats #now if I mess up I can load it from this
intersection_stats$flattened_peak_start <- apply(
  intersection_stats[,c(2,8,14,20)], 1, min, na.rm=TRUE)
intersection_stats$flattened_peak_end <- apply(
  intersection_stats[,c(3,9,15,21)], 1, max, na.rm=TRUE)
intersection_stats$flattened_chrom <- apply(
  intersection_stats[,c(1,7,13,19)], 1, max, na.rm=TRUE)
#by all logic max shouldn't work here, but it does, look I just need a chrom

intersection_grange = makeGRangesFromDataFrame(
  intersection_stats,
  start.field = "flattened_peak_start",
  end.field = "flattened_peak_end",
  seqnames.field = "flattened_chrom",
  keep.extra.columns = TRUE)

tx = TxDb.Hsapiens.UCSC.hg38.knownGene
hs = org.Hs.eg.db
source("seqTools.annotBedtoTxDbGene.R")

annot = annotBedtoTxDbGene(intersection_grange, tx, prefix = "hg38", org = hs)
intersection_anno = as.data.frame(annot)

intersection_anno <- intersection_anno %>% dplyr::select(-c(7,12,18,24))

write.table(
  intersection_stats,
  file="intersection_stats.txt",
  row.names=FALSE,
  quote=FALSE,
  sep="\t"
)
write.table(
  intersection_anno,
  file="intersection_anno.txt",
  row.names=FALSE,
  quote=FALSE,
  sep="\t"
)

intersection_anno <- read.delim("~/intersection_anno.txt")
intersection_anno_promoter <- intersection_anno %>% filter(abs(hg38.tssDist)<=1000)
