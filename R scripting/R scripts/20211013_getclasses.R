library (tidyverse)
library(compiler)

get_classes <- function(file_list, class_names_list) {
  #you need to return a table that has gene_name and class. All else is just set dressing.
  tmp_list = list()
  for (i in seq_along(file_list)) {
    file = read.delim(file_list[[i]])
    class = rep_len(c(class_names_list[[i]]), dim(file)[1])
    sub_file = as.data.frame(cbind(file, class))
    tmp_list[[i]] = sub_file
  }
  file_tab = do.call(rbind, tmp_list)
  return(file_tab)
}


proseq_class <- get_classes(file_list = list(
  "dist_propause_plus2.txt",
  "dist_propause_minus2.txt",
  "prox_propause_plus2.txt",
  "prox_propause_minus2.txt",
  "silent_propause_plus2.txt",
  "silent_propause_minus2.txt"), 
  class_names_list = list("dist", "dist", "prox", "prox", "silent", "silent"))

skew_class = get_classes(file_list = list("~/Lab stuff/R_coding/backedup/prox_skewed_plus.txt",
                                          "~/Lab stuff/R_coding/backedup/prox_skewed_minus.txt",
                                          "~/Lab stuff/R_coding/backedup/dist_skewed_plus.txt",
                                          "~/Lab stuff/R_coding/backedup/dist_skewed_minus.txt"),
                         class_names_list = list(
                           "prox", "prox", "dist", "dist"
                         ))

cgi_pause <- read.delim( #Reading in the CpG/Genes intersection file, from bedtools
  "cgi_pause_df_hg38_10_8_21.txt",
  quote="")
split_strand_p <- cgi_pause %>% filter(strand == "+")
split_strand_m <- cgi_pause %>% filter(strand == "-")
p_dedup_cgi <- split_strand_p %>% group_by(cgi_s) %>% mutate(tss=min(gene_s)) %>% filter(tss==gene_s)
m_dedup_cgi <- split_strand_m %>% group_by(cgi_s) %>% mutate(tss=max(gene_e)) %>% filter(tss==gene_e)
dedup_cgi <- rbind(p_dedup_cgi, m_dedup_cgi)

gencode_table <- read.delim( #Reading in the CpG/Genes intersection file, from bedtools
  "Gencode_CGI_genes_intersection_table.hg38.txt",
  quote="")

uniq_id <- paste0(gencode_table$tss, gencode_table$ensembl_ID)
gencode_table$gene_name <- uniq_id

heklt = proseq_class$class
names(heklt) = proseq_class$gene_name
gencode_table$hek_proseq_class = unname((heklt[as.character(gencode_table$gene_name)]))

skewlt = skew_class$class
names(skewlt) = skew_class$gene_name
gencode_table$skew_class = unname((skewlt[as.character(gencode_table$gene_name)]))

gt_w_count <- gencode_table %>% group_by(geneName) %>% add_count()
dedup_cgi_w_counts <- dedup_cgi %>% group_by(name2) %>% add_count()
gt_w_count_dedup <- gt_w_count %>% filter(n == 1)
cgi_pause_1name <- dedup_cgi_w_counts %>% filter(n == 1)

duped_me <- gt_w_count %>% filter(n != 1) #these are genes with multiple TSSes that intersect unique CGIs.
duped_josh <- dedup_cgi_w_counts %>% filter(n != 1) 

duped_me_p <- duped_me %>% filter(strand == "+") %>% group_by(gene_s) %>% mutate(min_end = min(gene_e)) %>% filter(min_end == gene_e) %>% ungroup() %>% group_by(geneName) %>% add_count()

duped_me_m <- duped_me %>% filter(strand == "-") %>% group_by(gene_e) %>% mutate(min_end = max(gene_s)) %>% filter(min_end == gene_s) %>% ungroup() %>% group_by(geneName) %>% add_count()

now_single_p <- duped_me_p %>% filter(nn == 1) %>% select(!c(nn, min_end))
now_single_m <- duped_me_m %>% filter(nn == 1) %>% select(!c(nn, min_end))
now_single <- rbind(now_single_p, now_single_m)

gt_w_count_dedup <- rbind (gt_w_count_dedup, now_single) #gives us 11028 total.

multi_tss_n_islands = duped_me_p %>% filter(nn != 1)
multi_tss_n_islands_m = duped_me_m %>% filter(nn != 1)
multi_tss_n_islands = rbind(multi_tss_n_islands, multi_tss_n_islands_m)
#unique genes= 819 lost to this
length(unique(multi_tss_n_islands$geneName))
write.table(multi_tss_n_islands, "multi_tss_n_islands.txt")

#lookup by name
wendylt = cgi_pause_1name$class
names(wendylt) = cgi_pause_1name$name2
gt_w_count_dedup$wendy_class = unname((wendylt[as.character(gt_w_count_dedup$geneName)]))

write.table(gt_w_count_dedup, "classed_cgi_gene_table.txt")

super_prox <- gt_w_count_dedup %>% filter(hek_proseq_class == "prox" & skew_class == "prox" & wendy_class == "Proximal")
dim(super_prox) #1962 genes of the 11028
super_dist <- gt_w_count_dedup %>% filter(hek_proseq_class == "dist" & skew_class == "dist" & wendy_class == "Distal")
dim(super_dist) #259 genes of the 11028
all_na <- gt_w_count_dedup %>% filter(is.na(hek_proseq_class) & is.na(skew_class) & is.na(wendy_class))
dim(all_na) #377 genes of the 11028
write.table(super_prox, "super_prox.txt")
write.table(super_dist, "super_dist.txt")
write.table(all_na, "all_na.txt")
