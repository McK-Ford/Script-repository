merged_mat <- read.csv("~/merged_mat.txt", sep="")
library(tidyverse)
merged_mat$Bins=as.numeric(merged_mat$Bins)
h1 <- ggplot(data=merged_mat, mapping=aes(
  x=Bins, y=gene_ID,fill=H3K4me3_EV)) + geom_raster()
h1

merged_mat$sqrtH3K4me3=sqrt(merged_mat$H3K4me3_EV)
h1 <- ggplot(data=merged_mat, mapping=aes(
  x=Bins, y=gene_ID,fill=H3K4me3_EV)) + geom_raster(interpolate=TRUE) +
  scale_fill_gradient(name="H3K4me3_EV", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(1,10)) +theme(axis.text.y = element_blank())
h1

h1 <- ggplot(data=merged_mat, mapping=aes(
  x=Bins, y=gene_ID,fill=H2AZ)) + geom_raster(interpolate=TRUE) +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1)) +theme(axis.text.y = element_blank())
h1

# Version 1 - reorder by max scores
mm = merged_mat %>%
  group_by(gene_ID) %>%
  mutate(tot_H3K4me3 = sum(H3K4me3_EV)) %>%
  ungroup() %>%
  arrange(desc(tot_H3K4me3), Bins)

h1 <- ggplot(data=mm, mapping=aes(
  x=Bins, y=reorder(gene_ID, tot_H3K4me3),fill=H3K4me3_EV)) + geom_raster() +
  scale_fill_gradient(name="H3K4me3_EV", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,10)) +theme(axis.text.y = element_blank())
h1

h1 <- ggplot(data=mm, mapping=aes(
  x=Bins, y=reorder(gene_ID, tot_H3K4me3),fill=H3K4me3_SH2)) + geom_raster(interpolate=TRUE) +
  scale_fill_gradient(name="H3K4me3_SH2", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,10)) +theme(axis.text.y = element_blank())
h1

#version 2 - by CpG island size
TSS_CpG_end$len = TSS_CpG_end$end - TSS_CpG_end$start
lt=TSS_CpG_end$len
names(lt)=TSS_CpG_end$name
mm$len = unname((lt[as.character(mm$gene_ID)]))

h1 <- ggplot(data=mm, mapping=aes(
  x=Bins, y=reorder(gene_ID, len),fill=H3K4me3_EV)) + geom_raster() +
  scale_fill_gradient(name="H3K4me3_EV", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,10)) +theme(axis.text.y = element_blank())
h1

h1 <- ggplot(data=mm, mapping=aes(
  x=Bins, y=reorder(gene_ID, len),fill=H2AZ)) + geom_raster() +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.25)) +theme(axis.text.y = element_blank())
h1

#####################
h1 <- ggplot(data=mm, mapping=aes(
  x=Bins, y=reorder(gene_ID, len),fill=H3K4me3_pWT)) + geom_raster() +
  scale_fill_gradient(name="H3K4me3_pWT", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,5)) +theme(axis.text.y = element_blank())
h1