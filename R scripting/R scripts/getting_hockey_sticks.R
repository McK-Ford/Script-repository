library(tidyverse)
library(data.table)
project_cell_lines_pI <- read.delim("C:/Users/kayle/Box/Vertinolab/GenomicSandbox/project_cell_lines/ProSeq/pI/V2/project_cell_lines_pI.txt")
setDT(project_cell_lines_pI)
project_cell_lines_pI[, ':='(pI_ro.T47D=frank(pI.T47D, na.last="keep"),  #this is data.table rank order function, bc much faster than base R
                             pI_ro.MB231=frank(pI.MB231, na.last="keep"),
                             pI_ro.SUM159=frank(pI.SUM159, na.last="keep"),
                             pI_ro.MCF7=frank(pI.MCF7, na.last="keep"),
                             pI_ro.MCF10A=frank(pI.MCF10A, na.last="keep"))]
mytheme = theme(
  panel.background=element_blank(),
  text=element_text(color="black",face="bold",family="sans"),
  axis.text=element_text(color="black"),
  axis.ticks=element_line(color="black"),
  plot.margin=unit(c(0.25, 0.25, .25, 0.25),"cm"),
  plot.title=element_text(vjust=2, hjust=0.5),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  axis.title.y=element_text(vjust=0),
  axis.title.x=element_text(vjust=0),
  panel.border = element_blank(),
  axis.line=element_line()
)

quantile(project_cell_lines_pI$pI_ro.T47D, 0.9, na.rm=TRUE)
#90% 
#5005 
T47_hs = ggplot(data=project_cell_lines_pI) +
  geom_point(aes(y=pI.T47D, x=pI_ro.T47D)) +
  mytheme + 
  xlab("pI_rank") +
  ylab("pI") +
  ggtitle("T47D")+
  geom_vline(aes(xintercept=5005), linetype="dotted")

quantile(project_cell_lines_pI$pI_ro.MB231, 0.9, na.rm=TRUE)
#90% 
#4891.6 
MB231_hs = ggplot(data=project_cell_lines_pI) +
  geom_point(aes(y=pI.MB231, x=pI_ro.MB231)) +
  mytheme + 
  xlab("pI_rank") +
  ylab("pI") +
  ggtitle("MB231")+
  geom_vline(aes(xintercept=4891), linetype="dotted")

quantile(project_cell_lines_pI$pI_ro.SUM159, 0.9, na.rm=TRUE)
#90% 
#4607.2  
SUM159_hs = ggplot(data=project_cell_lines_pI) +
  geom_point(aes(y=pI.SUM159, x=pI_ro.SUM159)) +
  mytheme + 
  xlab("pI_rank") +
  ylab("pI") +
  ggtitle("SUM159")+
  geom_vline(aes(xintercept=4607), linetype="dotted")

quantile(project_cell_lines_pI$pI_ro.MCF7, 0.9, na.rm=TRUE)
#90% 
#4645 
MCF7_hs = ggplot(data=project_cell_lines_pI) +
  geom_point(aes(y=pI.MCF7, x=pI_ro.MCF7)) +
  mytheme + 
  xlab("pI_rank") +
  ylab("pI") +
  ggtitle("MCF7")+
  geom_vline(aes(xintercept=4645), linetype="dotted")

quantile(project_cell_lines_pI$pI_ro.MCF10A, 0.9, na.rm=TRUE)
#90% 
#4380.4 
MCF10A_hs = ggplot(data=project_cell_lines_pI) +
  geom_point(aes(y=pI.MCF10A, x=pI_ro.MCF10A)) +
  mytheme + 
  xlab("pI_rank") +
  ylab("pI") +
  ggtitle("MCF10A")+
  geom_vline(aes(xintercept=4380), linetype="dotted")

pdf(file="project_cell_lines_hockey_sticks.pdf", width = 7)
T47_hs
MB231_hs
SUM159_hs
MCF7_hs
MCF10A_hs
dev.off()

top_10perc_SUM159 = project_cell_lines_pI %>% filter(pI_ro.SUM159 >= 4607)
top_10perc_T47D = project_cell_lines_pI %>% filter(pI_ro.T47D >= 5005)
top_10perc_MB231 = project_cell_lines_pI %>% filter(pI_ro.MB231 >= 4891)
top_10perc_MCF7 = project_cell_lines_pI %>% filter(pI_ro.MCF7 >= 4645)
top_10perc_MCF10A = project_cell_lines_pI %>% filter(pI_ro.MCF10A >= 4380)

#NAs being dif for each change number slightly but all are around 500.

MCF7_list = top_10perc_MCF7[[4]]
library(AnnotationDbi)
library(org.Hs.eg.db)
library(enrichR)
dbs <-listEnrichrDbs() #I think I want GO and KEGG terms right?
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", "KEGG_2021_Human")
options(enrichR.base.address="https://amp.pharm.mssm.edu/Enrichr/")
enriched <- enrichr(MCF7_list, dbs)
#head(enriched[["GO_Molecular_Function_2021"]])
#plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
SUM159_list = top_10perc_SUM159[[4]]
SUM159_enriched <- enrichr(SUM159_list, dbs)
T47D_list = top_10perc_T47D[[4]]
T47D_enriched <- enrichr(T47D_list, dbs)
MB231_list = top_10perc_MB231[[4]]
MB231_enriched <- enrichr(MB231_list, dbs)
MCF10A_list = top_10perc_MCF10A[[4]]
MCF10A_enriched <- enrichr(MCF10A_list, dbs)
#alright, what next?
plotEnrich(enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enriched[[4]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
#based on those P values cancer cell line encyclopedia will be useless. Can I also filter final result tab by P value?
#recreating this graphic on my own - y axis terms, x axis count, gradient P value is pretty easy. Geom bar.
#LoL = list of lists, stop judging me
LoL = list(enriched, SUM159_enriched, T47D_enriched, MB231_enriched, MCF10A_enriched)
dbs = as.list(dbs)
line_ids = list("MCF7", "SUM159", "T47D", "MB231", "MCF10A")
#for loop is going to be clearest, I don't want to nest lapply or do more copy paste
#probably should have looped some of the earlier stuff too tbh
#but sometimes copy paste is genuinely quicker for interactive work
for (i in seq_along(LoL)){
  tmp = Map(cbind, LoL[[i]], db_id=dbs)
  tmp2 = do.call(rbind, tmp)
  LoL[[i]] = tmp2
  }
LoL2 = Map(cbind, LoL, line_id=line_ids)
enrichment_tab = do.call(rbind, LoL2)

#safe to remove columns old ajusted p value and old p value
#probably only want examples with p value less than 0.05, we'll start with that then look at documentation for adjusted

enrichment_tab = enrichment_tab %>%
  dplyr::select(-c(5,6))

#reran without P value filtering bc information suggests true point of P value in enrichment analysis is to get a rank
#order, not to filter. Exploratory data analysis.

write.table(enrichment_tab, file="ontology.txt", quote=FALSE, sep="\t", row.names = FALSE)
####
ontology <- read.delim("~/ontology.txt")
library(tidyverse)

ont = ontology %>%
  separate(col=2, into=c("num_genes_match", "num_genes_set"),
           sep = "/", remove=TRUE, convert=TRUE)

#okay lets split this out by set type otherwise I'm making things more difficult for myself.
ont_GO_molec = ont[ont$db_id == "GO_Molecular_Function_2021",]
ont_KEGG = ont[ont$db_id == "KEGG_2021_Human",]
ont_GO_biol = ont[ont$db_id == "GO_Biological_Process_2021",]
ont_GO_cell_comp = ont[ont$db_id == "GO_Cellular_Component_2021",]

subprocesses = list(ont_KEGG, ont_GO_molec, ont_GO_cell_comp, ont_GO_biol)
for (i in seq_along(subprocesses)) {
  subprocesses[[i]] = subprocesses[[i]] %>% 
    group_by(line_id) %>%
    slice_min(P.value, n=10) %>%
    ungroup()
  }

mytheme = theme(
  panel.background=element_blank(),
  text=element_text(color="black",face="bold",family="sans"),
  axis.text=element_text(color="black", size=6),
  axis.ticks=element_line(color="black"),
  plot.margin=unit(c(0.25, 0.25, .25, 0.25),"cm"),
  plot.title=element_text(vjust=2, hjust=0.5),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  axis.title.y=element_text(vjust=0),
  axis.title.x=element_text(vjust=0),
  panel.border = element_blank(),
  axis.line=element_line()
)
ont_KEGG_plot = ggplot(data=subprocesses[[1]]) +
  geom_col(aes(y=reorder(Term, P.value, decreasing=TRUE), x=num_genes_match, fill=P.value)) +
  mytheme + 
  xlab("pI_rank") +
  ylab("pI") +
  ggtitle("KEGG") +
  facet_wrap(vars(line_id), scales="free")
ont_GO_molec_plot = ggplot(data=subprocesses[[2]]) +
  geom_col(aes(y=reorder(Term, P.value, decreasing=TRUE), x=num_genes_match, fill=P.value)) +
  mytheme + 
  xlab("pI_rank") +
  ylab("pI") +
  ggtitle("GO_molec") +
  facet_wrap(vars(line_id), scales="free")
ont_GO_cell_comp_plot = ggplot(data=subprocesses[[3]]) +
  geom_col(aes(y=reorder(Term, P.value, decreasing=TRUE), x=num_genes_match, fill=P.value)) +
  mytheme + 
  xlab("pI_rank") +
  ylab("pI") +
  ggtitle("GO_cell_comp") +
  facet_wrap(vars(line_id), scales="free")
ont_GO_biol_plot = ggplot(data=subprocesses[[4]]) +
  geom_col(aes(y=reorder(Term, P.value, decreasing=TRUE), x=num_genes_match, fill=P.value)) +
  mytheme + 
  xlab("pI_rank") +
  ylab("pI") +
  ggtitle("GO_biol") +
  facet_wrap(vars(line_id), scales="free")


pdf(file="project_cell_lines_ontology.pdf", width = 15)
ont_KEGG_plot
ont_GO_molec_plot
ont_GO_cell_comp_plot
ont_GO_biol_plot
dev.off()
