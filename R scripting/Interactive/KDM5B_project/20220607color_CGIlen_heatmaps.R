### The goal of this script was to:
#Import my metaplotter
#make heatmaps of annotated TSSes at CGIs for H3K4me3, H2AZ, and KDM5B (VER0046)
#I wanted them to be sorted by distance to the end of the CpG island, in both directions.
#but my first attempt went wrong bc I didn't carry over the distance in the metadata.
#so I ran over it several times.


# need to make: with and without insert, H2AZ, H3K4me3, KDM5B, sorted by CGI size, TSS downstream edge and upstream edge to TSS.
source("/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/NewMetaplotterJoshBased.R")
enableJIT(3)

bamDirectory="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/"
refseq_curated_longest_cgi_minus <- read.delim(paste0(
  bamDirectory,
  "refseq_curated_longest_cgi_minus.hg38.txt"))
refseq_curated_longest_cgi_plus <- read.delim(paste0(
  bamDirectory,
  "refseq_curated_longest_cgi_plus.hg38.txt"))

TSS_CpG_end_plus= refseq_curated_longest_cgi_plus %>%
  dplyr::select(1,5,3,8,2,7) #using CpG start as score placeholder
TSS_CpG_end_plus$name=paste0(TSS_CpG_end_plus$name, "_", TSS_CpG_end_plus$cpg_e)
TSS_CpG_end_minus= refseq_curated_longest_cgi_minus %>% dplyr::select(1,2,6,8,3,7) #cgi end placeholder
TSS_CpG_end_minus$name=paste0(TSS_CpG_end_minus$name, "_", TSS_CpG_end_minus$cpg_s)

name_vec = c("chrom", "start", "end", "name", "score", "strand")
colnames(TSS_CpG_end_minus) = name_vec
colnames(TSS_CpG_end_plus) = name_vec
TSS_CpG_end = rbind(TSS_CpG_end_minus, TSS_CpG_end_plus)

CpG_start_TSS_plus= refseq_curated_longest_cgi_plus %>%
  dplyr::select(1,2,5,8,3,7)
CpG_start_TSS_plus$name=paste0(CpG_start_TSS_plus$name, "_", CpG_start_TSS_plus$cpg_e)
CpG_start_TSS_minus= refseq_curated_longest_cgi_minus %>% dplyr::select(1,6,3,8,2,7) #cgi end placeholder
CpG_start_TSS_minus$name=paste0(CpG_start_TSS_minus$name, "_", CpG_start_TSS_minus$cpg_s)

colnames(CpG_start_TSS_plus) = name_vec
colnames(CpG_start_TSS_minus) = name_vec
CpG_start_TSS = rbind(CpG_start_TSS_plus, CpG_start_TSS_minus)

bamDirectory1="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/bams/hg38_dedup/"


bam_list = list(paste0(
  bamDirectory1, "ZS1_19_MCF7_EV_H3K4me3.dedup.bam"
),
paste0(
  bamDirectory1, "ZS1_19_MCF7_EV_H2AZ.dedup.bam"
),
paste0(
  bamDirectory1, "ZS1_19_MCF7_EV_KDM5B.dedup.bam"
)
)

m_TCE_NI = lapply(bam_list, get_score_matrix, bed=TSS_CpG_end, b=-3000, a=3000,
                  method="single_stranded_anchored")
#m_CST_NI = lapply(bam_list, get_score_matrix, bed=CpG_start_TSS, b=-1000, a=1000,
     #             method="single_stranded_anchored")
# m_TCE_I = lapply(bam_list, get_score_matrix, bed=TSS_CpG_end, b=-1000, a=1000,
#                   method="single_stranded_anchored", readsOnly=FALSE)
# m_CST_I = lapply(bam_list, get_score_matrix, bed=CpG_start_TSS, b=-1000, a=1000,
#                   method="single_stranded_anchored", readsOnly=FALSE)

m_TCE_NI2 <- lapply(m_TCE_NI, data.table::melt, measure.vars=c(2:601),
                    variable.name="Bins", value.name="Score")
#m_CST_NI2 <- lapply(m_CST_NI, data.table::melt, measure.vars=c(2:201),
        #            variable.name="Bins", value.name="Score")
# m_TCE_I2 <- lapply(m_TCE_I, data.table::melt, measure.vars=c(2:201),
#                     variable.name="Bins", value.name="Score")
# m_CST_I2 <- lapply(m_CST_I, data.table::melt, measure.vars=c(2:201),
#                     variable.name="Bins", value.name="Score")

#I also want to see if I can strand the CGIs at enhancers.


merged_mat=cbind(m_TCE_NI2[[1]], m_TCE_NI2[[2]][[3]], m_TCE_NI2[[3]][[3]])
colnames(merged_mat)=c("gene_ID", "Bins", "TCE_NI_H3K4me3", "TCE_NI_H2AZ",
                       "TCE_NI_KDM5B"
)
merged_mat$Bins = as.numeric(merged_mat$Bins)*10-3000
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


ggplot(data=merged_mat, mapping=aes(
  x=Bins, y=gene_ID,fill=TCE_NI_H3K4me3)) + geom_raster() +
  scale_fill_gradient(name="H3K4me3_EV_insert", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,11)) + theme(axis.text.y = element_blank())

# 
# pdf(file="test_TSS_CpG_end.pdf", width = 3.5)
# ggplot(data=merged_mat, mapping=aes(
#   x=Bins, y=gene_ID,fill=TCE_NI_H3K4me3)) + geom_raster() +
#   scale_fill_gradient(name="H3K4me3_EV", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,14)) + theme(axis.text.y = element_blank())
# ggplot(data=merged_mat, mapping=aes(
#   x=Bins, y=gene_ID,fill=TCE_NI_H2AZ)) + geom_raster() +
#   scale_fill_gradient(name="H2AZ_EV", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,2)) + theme(axis.text.y = element_blank())
# ggplot(data=merged_mat, mapping=aes(
#   x=Bins, y=gene_ID,fill=TCE_NI_KDM5B)) + geom_raster() +
#   scale_fill_gradient(name="KDM5B_EV", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,2)) + theme(axis.text.y = element_blank())
# ggplot(data=merged_mat, mapping=aes(
#   x=Bins, y=gene_ID,fill=TCE_I_H3K4me3)) + geom_raster() +
#   scale_fill_gradient(name="H3K4me3_EV_insert", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,22)) + theme(axis.text.y = element_blank())
# ggplot(data=merged_mat, mapping=aes(
#   x=Bins, y=gene_ID,fill=TCE_I_H2AZ)) + geom_raster() +
#   scale_fill_gradient(name="H2AZ_EV_insert", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,3)) + theme(axis.text.y = element_blank())
# ggplot(data=merged_mat, mapping=aes(
#   x=Bins, y=gene_ID,fill=TCE_I_KDM5B)) + geom_raster() +
#   scale_fill_gradient(name="KDM5B_EV_insert", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,2)) + theme(axis.text.y = element_blank())
# dev.off()

save(m_TCE_NI, file="mats.RData")
load("mats.RData")


#oh that's dumb, I ordered the 5' end in the order of the TSS_CpG end bc of what I merged them to...
#I also technically annotated them to the wrong anchor... So take the TSS_CpG end and order according to 5' end?
m_TCE_NI2 <- lapply(m_TCE_NI, data.table::melt, measure.vars=c(2:601),
                    variable.name="Bins", value.name="Score")
mm1 = cbind(m_TCE_NI2[[1]], m_TCE_NI2[[2]][[3]], m_TCE_NI2[[3]][[3]])
mm1$Bins=as.numeric(mm1$Bins)

colnames(mm1)=c("gene_ID", "Bins", "NI_H3K4me3", "NI_H2AZ",
                       "NI_KDM5B")

bamDirectory="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/"

refseq_curated_longest_cgi_minus <- read.delim(paste0(
  bamDirectory,
  "refseq_curated_longest_cgi_minus.hg38.txt"))
refseq_curated_longest_cgi_plus <- read.delim(paste0(
  bamDirectory,
  "refseq_curated_longest_cgi_plus.hg38.txt"))

refseq_curated_longest_cgi_plus$uniq_id=paste0(
  refseq_curated_longest_cgi_plus$name, "_",
  refseq_curated_longest_cgi_plus$cpg_e)
refseq_curated_longest_cgi_minus$uniq_id=paste0(
  refseq_curated_longest_cgi_minus$name, "_",
  refseq_curated_longest_cgi_minus$cpg_s)
refseq_curated_longest_cgi_plus$dist_tss_3prime =
  refseq_curated_longest_cgi_plus$cpg_e - refseq_curated_longest_cgi_plus$gene_s
refseq_curated_longest_cgi_plus$dist_5prime_tss=
  refseq_curated_longest_cgi_plus$gene_s - refseq_curated_longest_cgi_plus$cpg_s
refseq_curated_longest_cgi_minus$dist_tss_3prime=
  refseq_curated_longest_cgi_minus$gene_e - refseq_curated_longest_cgi_minus$cpg_s
refseq_curated_longest_cgi_minus$dist_5prime_tss=
  refseq_curated_longest_cgi_minus$cpg_e - refseq_curated_longest_cgi_minus$gene_e
refseq_curated_longest_cgi = rbind(refseq_curated_longest_cgi_plus, refseq_curated_longest_cgi_minus)

mm2 = merge(x=refseq_curated_longest_cgi, y=mm1, by.x="uniq_id", by.y="gene_ID", all=TRUE)
mm2$Bins = as.numeric(as.character(mm2$Bins))*10-3000

pdf(file="test.pdf", height = 8, width = 4)
ggplot(data=mm2, mapping=aes(
   x=Bins, y=reorder(uniq_id, dist_tss_3prime),fill=NI_H3K4me3)) + geom_raster() +
   scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#800000", oob=scales::squish, limits=c(0,16)) +
   theme(axis.text.y = element_blank())
dev.off()

pdf(file="test2.pdf", height = 8, width = 4)
ggplot(data=mm2, mapping=aes(
  x=Bins, y=reorder(uniq_id, dist_tss_3prime),fill=NI_H2AZ)) + geom_raster() +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#114945", oob=scales::squish, limits=c(0,2)) +
  theme(axis.text.y = element_blank())
dev.off()

pdf(file="test3.pdf", height = 8, width = 4)
ggplot(data=mm2, mapping=aes(
  x=Bins, y=reorder(uniq_id, dist_tss_3prime),fill=NI_KDM5B)) + geom_raster() +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#711971", oob=scales::squish, limits=c(0,2.4)) +
  theme(axis.text.y = element_blank())
dev.off()


pdf(file="test4.pdf", height = 8, width = 4)
ggplot(data=mm2, mapping=aes(
  x=Bins, y=reorder(uniq_id, dist_5prime_tss),fill=NI_H3K4me3)) + geom_raster() +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#800000", oob=scales::squish, limits=c(0,16)) +
  theme(axis.text.y = element_blank())
dev.off()

pdf(file="test5.pdf", height = 8, width = 4)
ggplot(data=mm2, mapping=aes(
  x=Bins, y=reorder(uniq_id, dist_5prime_tss),fill=NI_H2AZ)) + geom_raster() +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#114945", oob=scales::squish, limits=c(0,2)) +
  theme(axis.text.y = element_blank())
dev.off()

pdf(file="test6.pdf", height = 8, width = 4)
ggplot(data=mm2, mapping=aes(
  x=Bins, y=reorder(uniq_id, dist_5prime_tss),fill=NI_KDM5B)) + geom_raster() +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#711971", oob=scales::squish, limits=c(0,2.4)) +
  theme(axis.text.y = element_blank())
dev.off()

#still looks kind of bad, why???

mm3 = mm2 %>% filter(Bins == -990)

pdf(file="test7.pdf", height = 8, width = 2)
ggplot(data=mm3, mapping=aes(
  x=dist_tss_3prime, y=reorder(uniq_id, dist_tss_3prime))) + geom_point(size=0.25) +
  theme(axis.text = element_blank()) + xlim(0,3000)
dev.off()

pdf(file="test8.pdf", height = 8, width = 4)
ggplot(data=mm3, mapping=aes(
  x=dist_5prime_tss, y=reorder(uniq_id, dist_5prime_tss))) + geom_point(size=0.25) +
  theme(axis.text = element_blank()) + xlim(0,3000)
dev.off()

 # mm4 = mm3[c(1:10),]
# ggplot(data=mm4, mapping=aes(
#   x=Bins, y=gene_ID,fill=I_H3K4me3)) + geom_raster() +
#   scale_fill_gradient(name="H3K4me3_EV_insert", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,22))
# #wtf it's not preserving y axis order, fine we can force it to (and fix the other one too!)
# bamDirectory="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/"
# refseq_curated_longest_cgi_minus <- read.delim(paste0(
#   bamDirectory,
#   "refseq_curated_longest_cgi_minus.hg38.txt"))
# refseq_curated_longest_cgi_plus <- read.delim(paste0(
#   bamDirectory,
#   "refseq_curated_longest_cgi_plus.hg38.txt"))
# 
# TSS_CpG_end_plus= refseq_curated_longest_cgi_plus %>%
#   select(1,5,3,8,2,7) #using CpG start as score placeholder
# TSS_CpG_end_plus$name=paste0(TSS_CpG_end_plus$name, "_", TSS_CpG_end_plus$cpg_e)
# TSS_CpG_end_minus= refseq_curated_longest_cgi_minus %>% select(1,2,6,8,3,7) #cgi end placeholder
# TSS_CpG_end_minus$name=paste0(TSS_CpG_end_minus$name, "_", TSS_CpG_end_minus$cpg_s)
# 
# name_vec = c("chrom", "start", "end", "name", "score", "strand")
# colnames(TSS_CpG_end_minus) = name_vec
# colnames(TSS_CpG_end_plus) = name_vec
# TSS_CpG_end = rbind(TSS_CpG_end_minus, TSS_CpG_end_plus)
# TSS_CpG_end = arrange(TSS_CpG_end, desc(end-start))
# 
# CpG_start_TSS_plus= refseq_curated_longest_cgi_plus %>%
#   select(1,2,5,8,3,7)
# CpG_start_TSS_plus$name=paste0(CpG_start_TSS_plus$name, "_", CpG_start_TSS_plus$cpg_e)
# CpG_start_TSS_minus= refseq_curated_longest_cgi_minus %>% select(1,6,3,8,2,7) #cgi end placeholder
# CpG_start_TSS_minus$name=paste0(CpG_start_TSS_minus$name, "_", CpG_start_TSS_minus$cpg_s)
# 
# colnames(CpG_start_TSS_plus) = name_vec
# colnames(CpG_start_TSS_minus) = name_vec
# CpG_start_TSS = rbind(CpG_start_TSS_plus, CpG_start_TSS_minus)
# 
# #okay at this rate it might be easier to start from scratch... Ugh. still, try it for now
# mm0 = m_TCE_NI[[1]]
# mm0$len_3prime = TSS_CpG_end$end - TSS_CpG_end$start
# #how to get 5'
# 
# CpG_start_TSS_plus$len_3prime = CpG_start_TSS_plus$score - CpG_start_TSS_plus$end #bc score was after all CpG end...
# CpG_start_TSS_minus$len_3prime = CpG_start_TSS_minus$start - CpG_start_TSS_minus$score
# CpG_start_TSS = rbind(CpG_start_TSS_plus, CpG_start_TSS_minus)
# CpG_start_TSS = arrange(CpG_start_TSS, desc(len_3prime))
# 
# mm0$len_5prime = CpG_start_TSS$end - CpG_start_TSS$start
# mm0 = data.table::melt(mm0, measure.vars=c(2:201), variable.name="Bins", value.name="Score")
# mm1 = cbind(mm0, m_TCE_NI2[[2]][[3]], m_TCE_NI2[[3]][[3]],
#             m_TCE_I2[[1]][[3]], m_TCE_I2[[2]][[3]], m_TCE_I2[[3]][[3]])
# 
# colnames(mm1)=c("gene_ID", "len_3prime", "len_5prime", "Bins", "NI_H3K4me3", "NI_H2AZ",
#                 "NI_KDM5B", "I_H3K4me3", "I_H2AZ",
#                 "I_KDM5B")
# 
# mm1$Bins = as.numeric(as.character(mm1$Bins))
# mm1$end_3prime = mm1$len_3prime/10
# 
# geom = ggplot() +
#   geom_raster(data=mm1, aes(x=Bins, y=reorder(gene_ID, len_3prime), fill=I_H3K4me3)) +
#   scale_fill_gradient(name="H3K4me3_EV_insert", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,22)) +
#   theme(axis.text.y = element_blank())
# 
# geom 
# mm2 = mm1 %>% filter(Bins==0)
# ggplot(data=mm2) + geom_point(aes(y=reorder(gene_ID, len_3prime), x=end_3prime)) +
#   xlim(-1000, 1000)
# 
# 
# mm4 = mm1[c(1:10),]
# ggplot(data=mm4, mapping=aes(
#   x=Bins, y=reorder(gene_ID, len_5prime), fill=I_H3K4me3)) + geom_raster() +
#   scale_fill_gradient(name="H3K4me3_EV_insert", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,22)) #so it is reordering now...
# 
# pdf(file="CpG_start_TSS.pdf", width = 3.5)
# ggplot(data=mm1, mapping=aes(
#   x=Bins, y=reorder(gene_ID, len_5prime), fill=NI_H3K4me3)) + geom_raster() +
#   scale_fill_gradient(name="H3K4me3_EV", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,14)) + theme(axis.text.y = element_blank())
# ggplot(data=mm1, mapping=aes(
#   x=Bins, y=reorder(gene_ID, len_5prime), fill=NI_H2AZ)) + geom_raster() +
#   scale_fill_gradient(name="H2AZ_EV", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,2)) + theme(axis.text.y = element_blank())
# ggplot(data=mm1, mapping=aes(
#   x=Bins, y=reorder(gene_ID, len_5prime), fill=NI_KDM5B)) + geom_raster() +
#   scale_fill_gradient(name="KDM5B_EV", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,2)) + theme(axis.text.y = element_blank())
# ggplot(data=mm1, mapping=aes(
#   x=Bins, y=reorder(gene_ID, len_5prime), fill=I_H3K4me3)) + geom_raster() +
#   scale_fill_gradient(name="H3K4me3_EV_insert", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,22)) + theme(axis.text.y = element_blank())
# ggplot(data=mm1, mapping=aes(
#   x=Bins, y=reorder(gene_ID, len_5prime), fill=I_H2AZ)) + geom_raster() +
#   scale_fill_gradient(name="H2AZ_EV_insert", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,3)) + theme(axis.text.y = element_blank())
# ggplot(data=mm1, mapping=aes(
#   x=Bins, y=reorder(gene_ID, len_5prime), fill=I_KDM5B)) + geom_raster() +
#   scale_fill_gradient(name="KDM5B_EV_insert", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,2)) + theme(axis.text.y = element_blank())
# dev.off()
# 
# pdf(file="test_TSS_CpG_end.pdf", width = 3.5)
# ggplot(data=mm1, mapping=aes(
#   x=Bins, y=reorder(gene_ID, len_3prime), fill=NI_H3K4me3)) + geom_raster() +
#   scale_fill_gradient(name="H3K4me3_EV", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,14)) + theme(axis.text.y = element_blank())
# ggplot(data=mm1, mapping=aes(
#   x=Bins, y=reorder(gene_ID, len_3prime), fill=NI_H2AZ)) + geom_raster() +
#   scale_fill_gradient(name="H2AZ_EV", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,2)) + theme(axis.text.y = element_blank())
# ggplot(data=mm1, mapping=aes(
#   x=Bins, y=reorder(gene_ID, len_3prime), fill=NI_KDM5B)) + geom_raster() +
#   scale_fill_gradient(name="KDM5B_EV", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,2)) + theme(axis.text.y = element_blank())
# ggplot(data=mm1, mapping=aes(
#   x=Bins, y=reorder(gene_ID, len_3prime), fill=I_H3K4me3)) + geom_raster() +
#   scale_fill_gradient(name="H3K4me3_EV_insert", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,22)) + theme(axis.text.y = element_blank())
# ggplot(data=mm1, mapping=aes(
#   x=Bins, y=reorder(gene_ID, len_3prime), fill=I_H2AZ)) + geom_raster() +
#   scale_fill_gradient(name="H2AZ_EV_insert", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,3)) + theme(axis.text.y = element_blank())
# ggplot(data=mm1, mapping=aes(
#   x=Bins, y=reorder(gene_ID, len_3prime), fill=I_KDM5B)) + geom_raster() +
#   scale_fill_gradient(name="KDM5B_EV_insert", low = "#FFFFFF",
#                       high = "#012345", oob=scales::squish, limits=c(0,2)) + theme(axis.text.y = element_blank())
# dev.off()


##########################################
NonCGI_genes <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NonCGI_genes.txt", header=FALSE)
NonCGI_genes1 = NonCGI_genes %>% filter(V8==-1)
#that's 18375. Don't forget to deduplicate these...
test_inter_p = NonCGI_genes %>% filter(V8!=-1 & V5=="+")
#that's 22.5k
test_inter_m = NonCGI_genes %>% filter(V8!=-1 & V5=="-")
#that's 21.5k
test_inter_p = test_inter_p %>%
  group_by(V6) %>%
  mutate(testing_val = min(V8)) %>%
  filter(V8==testing_val) %>%
  mutate(testing_val2 = min(V2)) %>%
  filter(V2==testing_val2) %>%
  ungroup() %>%
  distinct(V6, .keep_all = TRUE) #was this just for picking 1 cgi per? and 1 tss
#now its 8943
tru_nocgi_p = test_inter_p %>% filter(V8>V2) %>% select(c(1:3, 6, 4:5))
#okay 2022 true ones then, that happen to have gene body cgis
test_inter_m = test_inter_m %>%
  group_by(V6) %>%
  mutate(testing_val = max(V9)) %>%
  filter(V9==testing_val) %>%
  mutate(testing_val2 = max(V3)) %>%
  filter(V3==testing_val2) %>%
  ungroup() %>%
  distinct(V6, .keep_all = TRUE)
#now its 8624
tru_nocgi_m = test_inter_m %>% filter(V3>V9) %>% select(c(1:3, 6, 4:5))
#1893
#no inter cleanup
no_inter_p = NonCGI_genes1 %>%
  filter(V5=="+") %>%
  group_by(V6) %>%
  mutate(testing_val = min(V2)) %>%
  filter(V2==testing_val) %>%
  ungroup() %>%
  distinct(V6, .keep_all = TRUE) %>% select(c(1:3, 6, 4:5))
#6566 more
no_inter_m = NonCGI_genes1 %>%
  filter(V5=="-") %>%
  group_by(V6) %>%
  mutate(testing_val = max(V3)) %>%
  filter(V3==testing_val) %>%
  ungroup() %>%
  distinct(V6, .keep_all = TRUE) %>% select(c(1:3, 6, 4:5))
#6377 more
nointer_tab = rbind(tru_nocgi_p, tru_nocgi_m, no_inter_p, no_inter_m) #16858
nointer_tab <-nointer_tab %>% filter(!grepl("chrY|([\\w_]+)alt|random|fix|v1", V1)) %>% distinct(V6, .keep_all = TRUE) #15413
#########################
nointer_mat = lapply(bam_list, get_score_matrix, bed=nointer_tab, b=-1000, a=1000,
                     method="single_stranded_anchored")
nointer_mat3 <- lapply(nointer_mat, data.table::melt, measure.vars=c(2:201),
                       variable.name="Bins", value.name="Score")
merged_nointer_mat=cbind(nointer_mat3[[1]], nointer_mat3[[2]][[3]], nointer_mat3[[3]][[3]]
)
colnames(merged_nointer_mat)=c("gene_ID", "Bins", "H2AZ", "H3K4me3", "KDM5B")

merged_nointer_mat[,gene_ID:=NULL]
mean_merged_nointer_mat = merged_nointer_mat[,lapply(.SD, mean), by=Bins]
mean_merged_nointer_mat$Bins = as.numeric(mean_merged_nointer_mat$Bins)*10-1000

mat_melt2 = data.table::melt(mean_merged_nointer_mat, measure.vars=2:4,
                             variable.name="Antibody", value.name = "Score")



geom = ggplot(data=mat_melt2, aes(
  x=Bins, y=Score, color=Antibody)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_color_manual(values=c("H2AZ"="goldenrod3", "H3K4me3_EV"="darkseagreen4",
                              "KDM5B_EV"="firebrick4"))
geom
