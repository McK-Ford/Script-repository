load(file="gettingTruTSSes.RData")
library(tidyverse)

plus_tab = rbind (col_trimmed_p, col_trimmed_p_silent)
plus_tab$start_skew = plus_tab$cgi_s - 1000
long_isle = plus_tab %>% filter(cgi_e>gene_end)
std_isle = plus_tab %>% filter(cgi_e<=gene_end)
long_isle$end_skew = long_isle$cgi_e + 1000
std_isle$end_skew = std_isle$gene_end + 1000
plus_w_skew_regs = rbind(long_isle, std_isle)

minus_tab = rbind (col_trimmed_m, col_trimmed_m_silent)
minus_tab$end_skew = minus_tab$cgi_e + 1000
long_isle_m = minus_tab %>% filter(cgi_s<gene_start)
std_isle_m = minus_tab %>% filter (cgi_s>=gene_start)
long_isle_m$start_skew = long_isle_m$cgi_s - 1000
std_isle_m$start_skew = std_isle_m$gene_start - 1000
minus_w_skew_regs = rbind(long_isle_m, std_isle_m)

tstbed_p = cbind(
  plus_w_skew_regs$chrom, plus_w_skew_regs$start_skew, plus_w_skew_regs$end_skew,
  plus_w_skew_regs$gene_name, plus_w_skew_regs$score, plus_w_skew_regs$strand
)
tstbed_m = cbind(
  minus_w_skew_regs$chrom, minus_w_skew_regs$start_skew, minus_w_skew_regs$end_skew,
  minus_w_skew_regs$gene_name, minus_w_skew_regs$score, minus_w_skew_regs$strand
)
write.table(
  tstbed_m, "minus_win.bed", quote=FALSE,
  row.names = FALSE, col.names = FALSE
)
write.table(
  tstbed_p, "plus_win.bed", quote=FALSE,
  row.names = FALSE, col.names = FALSE
)
