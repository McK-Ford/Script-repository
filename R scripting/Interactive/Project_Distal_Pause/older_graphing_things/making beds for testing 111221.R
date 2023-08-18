load(file="TruTSSesClassed.RData")
load(file="gettingTruTSSes.RData")
library(tidyverse)
verified_genes_long_plus = plus_classed %>% select(c(chrom, tss_s, gene_end))
verified_genes_long_minus = minus_classed %>% select(c(chrom, gene_start, tss_s))
verified_genes_short_plus = plus_short %>% select(c(chrom, tss_s, gene_end))
verified_genes_short_minus = minus_short %>% select(c(chrom, gene_start, tss_s))
not_verified_plus = col_trimmed_p_silent %>% select(c(chrom, gene_start, gene_end))
not_verified_minus = col_trimmed_m_silent %>% select(c(chrom, gene_start, gene_end))
#add labels then merge? bed is name score strand?
naming=c("chr", "gs", "ge")
names(verified_genes_long_plus)=naming
names(verified_genes_long_minus)=naming
names(verified_genes_short_plus)=naming
names(verified_genes_short_minus)=naming
names(not_verified_plus)=naming
names(not_verified_minus)=naming
vgl=rbind(verified_genes_long_minus, verified_genes_long_plus)
vgl$gn = rep("ver_long", dim(vgl)[1])
vgs=rbind(verified_genes_short_minus, verified_genes_short_plus)
vgs$gn = rep("ver_short", dim(vgs)[1])
nvg=rbind(not_verified_plus, not_verified_plus)
nvg$gn = rep("not_ver", dim(nvg)[1])
gene_bed1 = rbind(vgl, vgs, nvg)
gene_bed1$score = rep(0, dim(gene_bed1)[1])
gene_bed1$strand = rep("+", dim(gene_bed1)[1])

write.table(
  gene_bed1, "ver_tss_genes.bed", quote=FALSE,
  row.names=FALSE, col.names=FALSE, sep="\t"
)

verified_cgi_long_plus = plus_classed %>% select(c(chrom, cgi_s, cgi_e))
verified_cgi_long_minus = minus_classed %>% select(c(chrom, cgi_s, cgi_e))
verified_cgi_short_plus = plus_short %>% select(c(chrom, cgi_s, cgi_e))
verified_cgi_short_minus = minus_short %>% select(c(chrom, cgi_s, cgi_e))
not_ver_cgi_plus = col_trimmed_p_silent %>% select(c(chrom, cgi_s, cgi_e))
not_ver_cgi_minus = col_trimmed_m_silent %>% select(c(chrom, cgi_s, cgi_e))
#add labels then merge? bed is name score strand?
naming=c("chr", "cs", "ce")
names(verified_cgi_long_plus)=naming
names(verified_cgi_long_minus)=naming
names(verified_cgi_short_plus)=naming
names(verified_cgi_short_minus)=naming
names(not_ver_cgi_plus)=naming
names(not_ver_cgi_minus)=naming
vcl=rbind(verified_cgi_long_minus, verified_cgi_long_plus)
vcl$gn = rep("ver_long", dim(vcl)[1])
vcs=rbind(verified_cgi_short_minus, verified_cgi_short_plus)
vcs$gn = rep("ver_short", dim(vcs)[1])
nvc=rbind(not_ver_cgi_plus, not_ver_cgi_minus)
nvc$gn = rep("not_ver", dim(nvc)[1])
cgi_bed1 = rbind(vcl, vcs, nvc)
cgi_bed1$score = rep(0, dim(cgi_bed1)[1])
cgi_bed1$strand = rep("+", dim(cgi_bed1)[1])

write.table(
  cgi_bed1, "ver_tss_cgis.bed", quote=FALSE,
  row.names=FALSE, col.names=FALSE, sep="\t"
)
