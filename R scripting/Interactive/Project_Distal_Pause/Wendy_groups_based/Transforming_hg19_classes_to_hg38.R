library (tidyverse)
library(data.table)
write_bed_file <- function(tab, b_name) {
  write.table(
    tab,
    file=b_name,
    row.names=FALSE,
    col.names=FALSE,
    quote=FALSE,
    sep="\t"
  )
}

cgi_pause_df_3_20_15 <- read.delim("cgi_pause_df_3_20_15.txt")
cgi_pause_df_3_20_15 = unique(cgi_pause_df_3_20_15)
no_alt <- cgi_pause_df_3_20_15 %>%
  filter(!grepl('chrY|([\\w_]+)alt|random|fix|v1', chrom)) #1196
#I'm starting to think part of the problem is names aren't unique so technically I'm not joining by a unique column... Will UCSC liftover care if I'm not using a real gene name? if not lets make one up based on whatever the most unique factor is.
#cgi start is 9856 uniques, start is 11404, end is 11409, oooh TSS is 11922, given we have 1196 vars I probably won't get better than that, yeah name is only 11787
no_alt_uniq_id <- no_alt %>% distinct(TSS, .keep_all = TRUE)
uniq_id <- paste0("TSS", no_alt_uniq_id$TSS)
no_alt_uniq_id$uniq_id <- uniq_id

bed1_cgi = no_alt_uniq_id %>% select(1:3, 39, 24, 16)
bed2_gene = no_alt_uniq_id %>% select(4:6, 39, 24, 16)
bed3_tx = no_alt_uniq_id %>% select(15, 17, 18, 39, 24, 16)
bed4_cds = no_alt_uniq_id %>% select(15, 19, 20, 39, 24, 16)
minus = no_alt_uniq_id %>% filter(strand == "-")
plus = no_alt_uniq_id %>% filter(strand == "+")
bed5_exon1p = plus %>% select(15, 29, 30, 39, 24, 16) #may need to split exon file by strand, yep that's the issue there.
bed5_exon1m = minus %>% select(15, 30, 29, 39, 24, 16)

write_bed_file(bed1_cgi, "bed1_cgi.bed")
write_bed_file(bed2_gene, "bed2_gene.bed")
write_bed_file(bed3_tx, "bed3_tx.bed")
write_bed_file(bed4_cds, "bed4_cds.bed")
write_bed_file(bed5_exon1p, "bed5_exon1p.bed")
write_bed_file(bed5_exon1m, "bed5_exon1m.bed")

#Submit to UCSC liftover tool, https://genome.ucsc.edu/cgi-bin/hgLiftOver
bed1_cgi <- unique(read.delim("bed1_edit.bed", header = FALSE))
bed2_gene <- unique(read.delim("bed2_edit.bed", header = FALSE))
bed3_tx <- unique(read.delim("bed3_edit.bed", header = FALSE))
bed4_cds <- unique(read.delim("bed4_edit.bed", header = FALSE))
bed5_exon1p <- unique(read.delim("bed5p_edit.bed", header = FALSE))
bed5_exon1m <- unique(read.delim("bed5m_edit.bed", header = FALSE))

bed5_exon1 <- rbind(bed5_exon1p, bed5_exon1m)
names(bed1_cgi) <- c("chrom1", "start1", "end1", "uniq_id", "score1", "strand1")
names(bed2_gene) <- c("chrom2", "start2", "end2", "uniq_id", "score2", "strand2")
names(bed3_tx) <- c("chrom3", "start3", "end3", "uniq_id", "score3", "strand3")
names(bed4_cds) <- c("chrom4", "start4", "end4", "uniq_id", "score4", "strand4")
names(bed5_exon1) <- c("chrom5", "start5", "end5", "uniq_id", "score5", "strand5")

merge1 <- merge(no_alt_uniq_id, bed1_cgi, by = "uniq_id")
merge2 <- merge(merge1, bed2_gene, by = "uniq_id")
merge3 <- merge(merge2, bed3_tx, by = "uniq_id")
merge4 <- merge(merge3, bed4_cds, by = "uniq_id")
merge5 <- merge(merge4, bed5_exon1, by = "uniq_id")
#Beautiful, with the new uniq id the merging actually does what it is supposed to!!!! We do lose about 100 genes in the process that did not properly liftover in one of the files, though.
setDT(merge5)
merge5[, c("score1", "score2", "score3", "score4", "score5") := NULL]
merge5_chromsubset = merge5[chrom==chr & chrom == chrom.1 & chrom==chrom1 & chrom == chrom2 & chrom ==chrom3 & chrom == chrom4 & chrom == chrom5] #lose 2 observations. Not especially concerned about that - probably were observations of genes moved to alt chromosomes. Yep, just ran a unique on chrom5 - there's a couple genes in an alt chrom22.
merge5_chromsubset[, c("chr", "chrom.1", "chrom1", "chrom2", "chrom3", "chrom4", "chrom5") := NULL]
merge5_strandsubset = merge5_chromsubset[strand1 == strand2 & strand1 ==strand3 & strand1 == strand4 & strand1 == strand5] #this loses 1 row.
merge5_strandsubset[, c("strand", "strand2", "strand3", "strand4", "strand5") := NULL]
#okay, other columns to remove? Then start renaming and make columns that are of the same thing as in original dataset, unless do not know what they are. But rn it's 52 variables, could use cleaning up.
merge5_strandsubset[, c("sstrand", "Updist", "Downdist", "TSS", "X3.CGI", "exonStarts", "exonEnds", "exonFrames", "TSS_1", "first_ex_len", "downdist", "diff", "cgi_length", "gene_length") := NULL] #TSS and TSS.1 are the same thing. TSS should equal transcription start.
#remember if any of these columns are desired it would be very quick to remake this table with this code and the source files.
cleaned_table <- merge5_strandsubset %>% select(c(2, 23:33, 7:9, 15, 16))
names(cleaned_table) <- c("chrom", "cgi_s", "cgi_e", "strand", "gene_s", "gene_e", "tx_s", "tx_e", "cds_s", "cds_e", "exon1_s", "exon1_e", "pi", "class", "name", "score", "name2")
cleaned_table$gene_len <- cleaned_table$gene_e - cleaned_table$gene_s
cleaned_table$cgi_len <- cleaned_table$cgi_e - cleaned_table$cgi_s #dataset already filtered to greater than 200 long
cleaned_table$exon1_len <- cleaned_table$exon1_e - cleaned_table$exon1_s
cleaned_table$p_gene_in_cgi <- (cleaned_table$gene_len - pmax((cleaned_table$cgi_s - cleaned_table$gene_s),0) - pmax((cleaned_table$gene_e-cleaned_table$cgi_e),0))/cleaned_table$gene_len
cleaned_table$pex_in_cgi <- (cleaned_table$exon1_len - pmax((cleaned_table$cgi_s - cleaned_table$exon1_s),0) - pmax((cleaned_table$exon1_e-cleaned_table$cgi_e),0))/cleaned_table$exon1_len
write.table(
  cleaned_table,
  file="cgi_pause_df_hg38_10_8_21.txt",
  row.names=FALSE,
  quote=FALSE,
  sep="\t"
)
