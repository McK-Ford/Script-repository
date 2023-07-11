library(tidyverse)
project_cell_lines_pI <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/project_cell_lines_pI_for_laura/project_cell_lines_pI.txt")
pI_tab_longform = project_cell_lines_pI %>%
  select(-c(13:18, 20:25, 27:32,34:39, 41:46)) %>%
  pivot_longer(cols=c(13:17), names_to = c(".value", "cell_line"), names_pattern = "(.+)\\.(.+)")
library(umap)
#version 1: remove NA. Version 2 will be NA as 0.
pI_tab_longform_noNa = pI_tab_longform %>%
  filter(!is.na(pI)) %>%
  mutate(row_id = row_number())
pi.data = pI_tab_longform_noNa[,14]
pi_umap = umap(pi.data)

umap_df <- pi_umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(row_id = row_number()) %>%
  inner_join(pI_tab_longform_noNa, by="row_id")
head(umap_df)         

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

ggplot(data=umap_df, aes(
   x = UMAP1, y = UMAP2, color=cell_line
   )) +
   geom_point() +
   labs(x = "UMAP1", y = "UMAP2",
        title = "Project cell lines pI UMAP") +
  mytheme

#that cant be right
pI_tab_NaAs0 = project_cell_lines_pI %>%
  select(-c(13:18, 20:25, 27:32,34:39, 41:46)) %>%
  pivot_longer(cols=c(13:17),
               names_to = c(".value", "cell_line"),
               names_pattern = "(.+)\\.(.+)") %>%
  filter(!is.na(pI)) %>%
  pivot_wider(values_from=c(14),
              names_from = cell_line, names_sep=".",
              names_vary="slowest", values_fill = 0) %>%
  mutate(row_id = row_number())

#breaks if not force NA to 0.
pi.data = pI_tab_NaAs0[,13:17]
pi_umap = umap(pi.data, n_neighbors=5) 

umap_df <- pi_umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(row_id = row_number()) %>%
  inner_join(pI_tab_longform_noNa, by="row_id")
head(umap_df)         

UMAP = ggplot(data=umap_df, aes(
  x = UMAP1, y = UMAP2, color=cell_line
)) +
  geom_point() +
  labs(x = "UMAP1", y = "UMAP2",
       title = "Project cell lines pI UMAP") +
  mytheme
pdf(file="all_genes_umap.pdf")
UMAP
dev.off()
#looks better when cuts rownames out and nearest neighbor to 5 but still weird

pca_t1 = prcomp(pi.data, scale. = TRUE, center = TRUE)
#if use option 1, get pc1 = 1 bc there's only one column. So probably not that.
summary(pca_t1)
pca_t1_df = as.data.frame(pca_t1$x) 
test = cbind(pca_t1_df, pI_tab_NaAs0)
ggplot(data=pca_t1_df, aes(x=PC1, y=PC2)) +
    geom_point() + mytheme
pca_t1$x

pI_tab_NaAs0_2 = pI_tab_NaAs0 %>% filter((MCF7+MB231+MCF10A+SUM159+T47D)>0) #I didn't realize there were around 400 still all 0s...
transposed_pca = t(pI_tab_NaAs0_2)
pi.data = as.data.frame(transposed_pca[13:17,]) %>%
  type_convert()
#why are these character?
pi_umap = umap(pi.data, n_neighbors=5) 

umap_df<- pi_umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(row_id = row_number()) %>%
  inner_join(pI_tab_longform_noNa, by="row_id")
head(umap_df)         

ggplot(data=umap_df, aes(
  x = UMAP1, y = UMAP2, color=cell_line
)) +
  geom_point() +
  labs(x = "UMAP1", y = "UMAP2",
       title = "Project cell lines pI UMAP") +
  mytheme

pca_t1 = prcomp(pi.data, scale. = TRUE, center = TRUE)
pca_t1
#if use option 1, get pc1 = 1 bc there's only one column. So probably not that.
summary(pca_t1)
pca_t1_df = as.data.frame(pca_t1$x) 
test = cbind(pca_t1_df, pI_tab_NaAs0)
ggplot(data=pca_t1_df, aes(x=PC1, y=PC2)) +
  geom_point()
mytheme
pca_t1$x

#Okay based on these results, I don't think we're having the rotation problem I was worried about.
quantile(pI_tab_NaAs0$MCF7, 0.9)
#90th quantile by pI is 54
quantile(pI_tab_NaAs0$MB231, 0.9)
#90th quant is 13
quantile(pI_tab_NaAs0$MCF10A, 0.9)
#90th quant is 29
quantile(pI_tab_NaAs0$SUM159, 0.9)
#90th quant is 20
quantile(pI_tab_NaAs0$T47D, 0.9)
#90th quant 27
pI_tab_high_sig = pI_tab_NaAs0 %>% filter(MCF7>quantile(MCF7, 0.9) |
                                          MB231>quantile(MB231, 0.9) |
                                          MCF10A>quantile(MCF10A, 0.9) |
                                          SUM159>quantile(SUM159, 0.9) |
                                          T47D>quantile(T47D, 0.9)
                                          )
pi.data = pI_tab_high_sig[,13:17]
pi_umap = umap(pi.data, n_neighbors=5) 

umap_df <- pi_umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(row_id = row_number()) %>%
  inner_join(pI_tab_longform_noNa, by="row_id")
head(umap_df)         

UMAP2 = ggplot(data=umap_df, aes(
  x = UMAP1, y = UMAP2, color=cell_line
)) +
  geom_point() +
  labs(x = "UMAP1", y = "UMAP2",
       title = "Project cell lines pI UMAP") +
  mytheme
pdf(file="top10percentbyline_umap.pdf")
UMAP2
dev.off()

#top 10percent total?
quantile((pI_tab_NaAs0$T47D+pI_tab_NaAs0$SUM159 +
          pI_tab_NaAs0$MCF10A + pI_tab_NaAs0$MB231 +
            pI_tab_NaAs0$MCF7), 0.9)
#90th quant 129
pI_tab_high_sig2 = pI_tab_NaAs0 %>% filter((T47D+SUM159 +
                                             MCF10A + MB231 +
                                             MCF7) > 
                                            quantile((T47D+SUM159 +
                                                     MCF10A + MB231 +
                                                     MCF7), 0.9))
pi.data = pI_tab_high_sig2[,13:17]
pi_umap = umap(pi.data, n_neighbors=5) 

umap_df <- pi_umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(row_id = row_number()) %>%
  inner_join(pI_tab_longform_noNa, by="row_id")
head(umap_df)         

UMAP3 = ggplot(data=umap_df, aes(
  x = UMAP1, y = UMAP2, color=cell_line
)) +
  geom_point() +
  labs(x = "UMAP1", y = "UMAP2",
       title = "Project cell lines pI UMAP") +
  mytheme
pdf(file="top10percentoverall_umap.pdf")
UMAP3
dev.off()
##
pca_t1 = prcomp(pi.data, scale. = TRUE, center = TRUE)
summary(pca_t1)
pca_t1_df = as.data.frame(pca_t1$x) 
test = cbind(pI_tab_high_sig2, pca_t1_df)
ggplot(data=pca_t1_df, aes(x=PC1, y=PC2)) +
  geom_point() + mytheme
#maybe try an outgroup cell line or two for compare?


source("~/Script repository/My_Useful_Fns.R")

gc <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/Regions/Untouched source data/Gencode39_knowngene.hg38.txt.gz",
  header=TRUE
)
#remove irrelevant columns
gc2 <- gc %>% dplyr::select(-c(7:17, 20, 22, 24:26)) #unless Ching-Hua used detect transcripts?
gc3 <- gc2[,c(1:3, 7, 5:6, 4, 8:10)]
colnames(gc3) <- c("chrom", "start", "end", "symbol",
                   "score", "strand", "ENST_ID", "ID3", "class", "type")
gc_genes_only <- gc3 %>%
  filter(type=="protein_coding" & !grepl('chrY|([\\w_]+)alt|random|fix|v1', chrom)) ###86K out of 266K
gc4p <- gc_genes_only %>%
  group_by(symbol) %>%
  filter(strand == "+") %>%
  mutate(start_site=min(start)) %>%
  filter(start==start_site) %>%
  dplyr::select(-c(11)) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  ungroup()
gc4m <- gc_genes_only %>%
  group_by(symbol) %>%
  filter(strand == "-") %>%
  mutate(start_site=max(end)) %>%
  filter(end==start_site) %>%
  dplyr::select(-c(11)) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  ungroup()

gc4=rbind(gc4p, gc4m)
Danko_MCF7 <- pI(bed=gc4,
          bam="C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/PublicData/GSE93229_Danko_2018_MCF7_PROseq/bamdeduped/MCF7_B7.deduped.bam",
          pairedEnd=FALSE, pause_s=0, pause_e=150)
Fang_293T <- pI(bed=gc4,
           bam="C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/PublicData/GSE156400_Fang_2021_HEK293T_PROseq_Rloops/Liang2021HEK_Rloops/Pro-seq/20211207 realigned/bam/Liang_hek_proseq_2021_1.bam",
           pairedEnd = FALSE, pause_s=0, pause_e = 150)
Sawarker_293T <- pI(bed=gc4,
             bam="C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/PublicData/GSE112379_Sawarkar_2019_HEK293T_PROseq/20211207 Realigned/bam/Garcia_hek_proseq_2021_1.bam",
             pairedEnd = FALSE, pause_s=0, pause_e = 150)



##
pI_tab2 <- cbind(Danko_MCF7, Fang_293T[,13:18], Sawarker_293T[,13:18])
colnames(pI_tab2) = c(colnames(pI_tab2[,1:12]),
                     paste0("Danko_MCF7",           ".", colnames(pI_tab2[,13:18])),
                     paste0("Fang_293T",           ".", colnames(pI_tab2[,13:18])),
                     paste0("Sawarker_293T",         ".", colnames(pI_tab2[,13:18])))
pI_tab_longform2 = pI_tab2 %>%
  pivot_longer(cols=c(13:30), names_to = c("cell_line", ".value"), names_pattern = "(.+)\\.(.+)")
pI_tab_longform3 = pI_tab_longform2 %>%
  filter(!is.na(pI) & pI!="Inf" & total_norm>=1)
#get shared.
outgroups_longform_trimmed = pI_tab_longform3 %>%
  select(-c(14:18)) %>%
  filter(symbol %in% pI_tab_NaAs0$symbol) %>%
  filter(!is.na(pI))
pI_tab_NaAs0_2 = pI_tab_NaAs0 %>%
  pivot_longer(cols=c(13:17),
               names_to = "cell_line",
               values_to = "pI")
project_lines_w_outgroups = rbind(pI_tab_NaAs0_2[,c(1:12, 14:15)],
                                  outgroups_longform_trimmed)  %>%
  mutate(row_id = row_number())
project_lines_w_outgroups_wide = project_lines_w_outgroups %>%
  pivot_wider(values_from=c(14),
              names_from = cell_line, names_sep=".",
              names_vary="slowest", values_fill = 0)

#breaks if not force NA to 0.
pi.data = project_lines_w_outgroups_wide[,13:20]
pi_umap = umap(pi.data, n_neighbors=8) 

umap_df <- pi_umap$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(row_id = row_number()) %>%
  inner_join(project_lines_w_outgroups, by="row_id")
head(umap_df)         

UMAP = ggplot(data=umap_df, aes(
  x = UMAP1, y = UMAP2, color=cell_line
)) +
  geom_point() +
  labs(x = "UMAP1", y = "UMAP2",
       title = "Project cell lines pI UMAP") +
  mytheme
pdf(file="all_genes_umap_w_outgroups.pdf")
UMAP
dev.off()
