#### get stranded CGIs ####
##### Setup #####
cpgIslands <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/cpgIsland.hg38.txt",
  header=FALSE
) #standard no modification islands text
source(
  "/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/NewMetaplotterJoshBased.R"
)
enableJIT(3)

nascentDirectory="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/"
plus_strand_list=list("MCF7_H9_plus.bam",
                      "MCF7_B7_plus.bam",
                      "MCF7_G11_plus.bam",
                      "MCF7_C11_plus.bam")
plus_strand_list=paste0(nascentDirectory, plus_strand_list)
minus_strand_list=list("MCF7_H9_minus.bam",
                       "MCF7_B7_minus.bam",
                       "MCF7_G11_minus.bam",
                       "MCF7_C11_minus.bam")
minus_strand_list=paste0(nascentDirectory, minus_strand_list)

#make the CGI file into a bed
cpgIslands$score = 0
cpgIslands$strand = "+"
cpgIslands$V4 = paste0("CGI_", cpgIslands$V2, "_", cpgIslands$V3) #unique ID

##### Get plus & minus transcription across island #####
plus_mat = lapply(plus_strand_list, get_score_matrix, bed=cpgIslands, n=1,
                  method="bi_stranded_anchored", pairedEnd=FALSE)
minus_mat = lapply(minus_strand_list, get_score_matrix, bed=cpgIslands, n=1,
                   method="bi_stranded_anchored", pairedEnd=FALSE)

cpgIslands_stranded = cpgIslands[1:4] #just cut false strand and score off bed
cpgIslands_stranded = cpgIslands_stranded %>%
  filter(!grepl("chrY|([\\w_]+)alt|random|fix|v1|v2", V1)) # cut out alt chroms
pm = cbind(plus_mat[[1]],
           plus_mat[[2]][[2]],
           plus_mat[[3]][[2]],
           plus_mat[[4]][[2]])
cpgIslands_stranded2 = merge(x=cpgIslands_stranded,
                             y=pm, by.x="V4",
                             by.y="rn", all=TRUE)
mm = cbind(minus_mat[[1]],
           minus_mat[[2]][[2]],
           minus_mat[[3]][[2]],
           minus_mat[[4]][[2]])
cpgIslands_stranded3 = merge(x=cpgIslands_stranded2,
                             y=mm, by.x="V4", by.y="rn",
                             all=TRUE)
colnames(cpgIslands_stranded3) = c("ID", "chrom", "start",
                                   "end", "p1","p2","p3",
                                   "p4","m1","m2","m3","m4")
cpgIslands_stranded3$plus =
  cpgIslands_stranded3$p1 + cpgIslands_stranded3$p2 +
  cpgIslands_stranded3$p3 + cpgIslands_stranded3$p4
cpgIslands_stranded3$minus =
  cpgIslands_stranded3$m1 + cpgIslands_stranded3$m2 +
  cpgIslands_stranded3$m3 + cpgIslands_stranded3$m4
#now we have a table of cgis and associated transcriptional scores across entire

##### assign these CGI a strand based on this transcription #####
summary(cpgIslands_stranded3$plus)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.     
#   0.00    0.00    2.10   31.80   33.87 2949.83  
summary(cpgIslands_stranded3$minus)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.       
#   0.0000    0.0615    3.4093   32.5872   41.9621 2226.0639   
summary(cpgIslands_stranded3$minus+cpgIslands_stranded3$plus)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.  
#   0.000    0.149    6.783   64.385   81.469 5175.509   
#lets say greater than 6 reads in it total, because we can't strand things that
#aren't expressed?
CpGIslands_silent = cpgIslands_stranded3 %>% filter(minus+plus<6)
CpGIslands_plus = cpgIslands_stranded3 %>% filter(minus+plus>6 & plus>minus)
# a good portion of these are likely genic of course
CpGIslands_minus = cpgIslands_stranded3 %>% filter(minus+plus>6 & plus<minus)
#okay next I would need to see if they're associated with genes.
#Strength of dif is also an ordering option... Let's calculate that.
CpGIslands_plus$strand = "+"
CpGIslands_plus$fc = (CpGIslands_plus$plus+1)/(CpGIslands_plus$minus+1)
CpGIslands_minus$strand = "-"
CpGIslands_minus$fc = (CpGIslands_minus$minus+1)/(CpGIslands_minus$plus+1)
cpgIslands4 = rbind(CpGIslands_minus, CpGIslands_plus)
cgi_bed = cpgIslands4 %>% select(c(2:4, 1, 16, 15)) #bed format for future use.
write.table(cgi_bed, file="strandedCGIs.bed",
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


#### get loc of maximal transcription ####
##### get location of maximal transcription to be anchor point ####
source("/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/NewMetaplotterJoshBased.R")
enableJIT(3)
nascentDirectory="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/"
plus_strand_list=list("MCF7_H9_plus.bam",
                      "MCF7_B7_plus.bam",
                      "MCF7_G11_plus.bam",
                      "MCF7_C11_plus.bam")
plus_strand_list=paste0(nascentDirectory, plus_strand_list)
minus_strand_list=list("MCF7_H9_minus.bam",
                       "MCF7_B7_minus.bam",
                       "MCF7_G11_minus.bam",
                       "MCF7_C11_minus.bam")
minus_strand_list=paste0(nascentDirectory, minus_strand_list)
strandedCGIs <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/strandedCGIs.bed", header=FALSE)

summary(strandedCGIs$V3-strandedCGIs$V2)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#201.0   479.0   753.0   938.9  1147.5 32637.0 
quantile(strandedCGIs$V3-strandedCGIs$V2, .95)
#95% 
#2165.6
#maybe we could stick with the ones that are less than 25 hundred?
#Then cut them to the actual cgi len.
strandedCGIs = strandedCGIs %>% filter(strandedCGIs$V3-strandedCGIs$V2<2500)
#loses 469 CGI.

plusCGI = strandedCGIs %>% filter(V6=="+")
minusCGI = strandedCGIs %>% filter(V6=="-")
#most of these are probably genes, annotation will be necessary
#goal is to get what bin for each of these CGI has maximal transcription
#and set that as 'location of transcription initiation.' Though it's probably
#the pause really.
plus_mat = lapply(
  plus_strand_list, get_score_matrix, bed=plusCGI,
  b=0, a=2500, bs = 25, method="single_stranded_anchored", pairedEnd=FALSE)
minus_mat = lapply(
  minus_strand_list, get_score_matrix, bed=minusCGI,
  b=0, a=2500, bs = 25, method="single_stranded_anchored", pairedEnd=FALSE)

plus_mat_long <- lapply(plus_mat, data.table::melt, measure.vars=c(2:101),
                        variable.name="Bins", value.name="Score")
minus_mat_long <- lapply(minus_mat, data.table::melt, measure.vars=c(2:101),
                         variable.name="Bins", value.name="Score")


pm = cbind(plus_mat_long[[1]], plus_mat_long[[2]][[3]],
           plus_mat_long[[3]][[3]], plus_mat_long[[4]][[3]])
cpgIslands_stranded2 = merge(x=plusCGI, y=pm, by.x="V4",
                             by.y="rn", all=TRUE)

mm = cbind(minus_mat_long[[1]], minus_mat_long[[2]][[3]],
           minus_mat_long[[3]][[3]], minus_mat_long[[4]][[3]])
cpgIslands_stranded3 = merge(x=minusCGI, y=mm, by.x="V4",
                             by.y="rn", all=TRUE)
colnames(cpgIslands_stranded2) = c("ID", "chrom", "start", "end",
                                   "fc", "strand", "bin", "p1","p2","p3","p4")
colnames(cpgIslands_stranded3) = c("ID", "chrom", "start", "end",
                                   "fc", "strand", "bin", "m1","m2","m3","m4")
cpgIslands_stranded2$bin = as.numeric(cpgIslands_stranded2$bin)
cpgIslands_stranded3$bin = as.numeric(cpgIslands_stranded3$bin)
##### get maximum score, + in case of ties, minimum bin #####
cpgIslands_stranded2$max_end_bin =
  (cpgIslands_stranded2$end - cpgIslands_stranded2$start)/25
cpgIslands_stranded3$max_end_bin =
  (cpgIslands_stranded3$end - cpgIslands_stranded3$start)/25
plus_maxreg = cpgIslands_stranded2 %>%
  pivot_longer(cols=8:11, names_to="source", values_to = "score") %>%
  group_by(ID) %>%
  filter(bin<max_end_bin+1) %>%
  ungroup()#rounding up
plus_maxreg = plus_maxreg %>%
  group_by(ID, bin) %>%
  mutate(sum_score = sum(score)) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(max_score=max(sum_score)) %>%
  filter(sum_score==max_score & source=="p1") %>%
  mutate(minbin=min(bin)) %>%
  filter(bin==minbin)

minus_maxreg = cpgIslands_stranded3 %>%
  pivot_longer(cols=8:11, names_to="source", values_to = "score") %>%
  group_by(ID) %>%
  filter(bin<max_end_bin+1) %>%
  ungroup()
minus_maxreg = minus_maxreg %>%
  group_by(ID, bin) %>%
  mutate(sum_score = sum(score)) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(max_score=max(sum_score)) %>%
  filter(sum_score==max_score & source=="m1") %>%
  mutate(minbin=min(bin)) %>%
  filter(bin==minbin)

plus_maxreg$TS = plus_maxreg$start + (plus_maxreg$bin*25)
minus_maxreg$TS = minus_maxreg$end - (minus_maxreg$bin*25)

maxreg = rbind(plus_maxreg, minus_maxreg)
maxreg_trim = maxreg %>% select(c(2:4,1,5:6, 12, 14))
write.table(maxreg_trim, file="strandedCGIs.bed", row.names=FALSE,
            col.names=FALSE, quote=FALSE, sep="\t")
#### annotating ####
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
tx = TxDb.Hsapiens.UCSC.hg38.knownGene
hs = org.Hs.eg.db
source("C:/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/seqTools.annotBedtoTxDbGene_ALTERED.R")
CGI_grange = makeGRangesFromDataFrame(strandedCGIs, keep.extra.columns = TRUE, seqnames.field = c("V1"), start.field = c("V2"), end.field = c("V3"))

annot_gr = annotgrtoTxDbGene(org=hs, gr=CGI_grange, tx=tx)

annot_df = data.frame(annot_gr)
write.table(annot_df, "annotatedcgis.txt", quote=FALSE, sep="\t", row.names = FALSE)


#### getting symbols ####
library(tidyverse)
library(biomaRt)
annotatedcgis <- read.delim("~/annotatedcgis.txt")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
conversion = getBM(attributes = c('ensembl_transcript_id', 'hgnc_symbol'),
                   mart = ensembl)
#apparently after dot is just version number so I can just trim it?
annotatedcgis$fixed_txnames = sapply(                #okay not really my style
  strsplit(annotatedcgis$TXNAME, ".", fixed=T),   #got it online
  function(x) x[1])  #but it works and that's what's important
annotatedcgis$fixed_tssnames = sapply(
  strsplit(annotatedcgis$tss, ".", fixed=T),
  function(x) x[1])
lt = conversion[,2]
names(lt) = conversion[,1]
annotatedcgis$txsymbol = unname((
  lt[as.character(annotatedcgis$fixed_txnames)]))
annotatedcgis$tsssymbol = unname((
  lt[as.character(annotatedcgis$fixed_tssnames)]))

annotatedcgis = annotatedcgis %>% dplyr::select(-c(5,10,23,24))
colnames(annotatedcgis) = c("chrom", "CGI_start", "CGI_end", "CGI_width",
                            "CGI_ID", "proseq_fc", "CGI_imputed_strand", 
                            "CGI_maxscore", "maxscore_loc",
                            colnames(annotatedcgis)[10:22])

annotatedcgis$maxscore_loc2 = annotatedcgis$maxscore_loc + 1

intergenic_cgis = annotatedcgis %>% filter(CGI_end<=txStart | CGI_start>=txEnd) 
#only 628. How many did josh find, well were all his transcribed this much?
#reminder these are only transcribed CGIs.

intergenic_cgis = annotatedcgis %>% filter(txDist!=0) #lets go with this one,
#it's easier. 606 here.
genebody_cgis = annotatedcgis %>% filter(txDist==0 & tssDist!=0) # then 2239. 
#Some of these may be intronic enhancers.
tss_cgis = annotatedcgis %>% filter(tssDist==0) #10841. makes sense.

#so what will happen if we order them by the strength of the stranded 
#transcriptional difference? bc pulling
#up the top intergenic in igv it's stranded... Or by strength of maxloc for 
#that matter. Make beds from them, first of all.
write.table(annotatedcgis, file="annotatedcgis.txt",
            row.names=FALSE, quote=FALSE, sep="\t")
##### Getting H3K4me3/H2AZ/KDM5B for these locs ####
iCGIs = intergenic_cgis %>% dplyr::select(c(1, 9, 23, 5:7)) #need these as beds
gbCGIs = genebody_cgis %>% dplyr::select(c(1, 9, 23, 5:7))
tssCGIs = tss_cgis %>% dplyr::select(c(1, 9, 23, 5:7))

source(
  "/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/NewMetaplotterJoshBased.R"
  )
enableJIT(3)

bamDirectory=
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/bams/hg38_dedup/"
bam_list=list("ZS1_19_MCF7_EV_H3K4me3.dedup.bam",
              "ZS1_19_MCF7_EV_H2AZ.dedup.bam",
              "ZS1_19_MCF7_EV_KDM5B.dedup.bam")
bam_list=paste0(bamDirectory, bam_list)

iCGIs_mat = lapply(
  bam_list, get_score_matrix, bed=iCGIs, b=-3000, a=3000,
  bs = 10, method="single_stranded_anchored"
  )
gbCGIs_mat = lapply(
  bam_list, get_score_matrix, bed=gbCGIs, b=-3000, a=3000,
  bs = 10, method="single_stranded_anchored"
  )
tssCGIs_mat = lapply(
  bam_list, get_score_matrix, bed=tssCGIs, b=-3000, a=3000,
  bs = 10, method="single_stranded_anchored"
  )

iCGIs_mat <- lapply(iCGIs_mat, data.table::melt, measure.vars=c(2:601),
                    variable.name="Bins", value.name="Score")
gbCGIs_mat <- lapply(gbCGIs_mat, data.table::melt, measure.vars=c(2:601),
                     variable.name="Bins", value.name="Score")
tssCGIs_mat <- lapply(tssCGIs_mat, data.table::melt, measure.vars=c(2:601),
                      variable.name="Bins", value.name="Score")

iCGIs_merged_mat=cbind(
  iCGIs_mat[[1]], iCGIs_mat[[2]][[3]], iCGIs_mat[[3]][[3]]
  ) #should I tag on the metadata? I could, use key-value pairs
gbCGIs_merged_mat=cbind(
  gbCGIs_mat[[1]], gbCGIs_mat[[2]][[3]], gbCGIs_mat[[3]][[3]]
  )
tssCGIs_merged_mat=cbind(
  tssCGIs_mat[[1]], tssCGIs_mat[[2]][[3]], tssCGIs_mat[[3]][[3]]
  )

name_vec =c("gene_ID", "Bins", "H3K4me3", "H2AZ",
            "KDM5B")

colnames(iCGIs_merged_mat) = name_vec
colnames(gbCGIs_merged_mat) = name_vec
colnames(tssCGIs_merged_mat) = name_vec

iCGIs_merged_mat$Bins = as.numeric(as.character(iCGIs_merged_mat$Bins))
gbCGIs_merged_mat$Bins = as.numeric(as.character(gbCGIs_merged_mat$Bins))
tssCGIs_merged_mat$Bins = as.numeric(as.character(tssCGIs_merged_mat$Bins))

mytheme = theme(
  panel.background=element_blank(),
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


quantile(iCGIs_merged_mat$H3K4me3,0.95)
#95% 
#3.142425 

quantile(gbCGIs_merged_mat$H3K4me3,0.95)
#95% 
#4.109324

quantile(tssCGIs_merged_mat$H3K4me3,0.95)
#95% 
#11.6028

##

quantile(iCGIs_merged_mat$H2AZ,0.95)
#95% 
#1.661545 

quantile(gbCGIs_merged_mat$H2AZ,0.95)
#95% 
#1.083616 

quantile(tssCGIs_merged_mat$H2AZ,0.95)
#95% 
#1.733786 

##

quantile(iCGIs_merged_mat$KDM5B,0.95)
#95% 
#1.20822 

quantile(gbCGIs_merged_mat$KDM5B,0.95)
#95% 
#0.9061648 

quantile(tssCGIs_merged_mat$KDM5B,0.95)
#95% 
#1.510275 

ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=gene_ID,fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.75)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

#I think they might all be backwards? idk. Unless the point of maximal
#transcription is downstream of the TSS, which is entirely possible, in which 
#case I'd expect a band of H2K4me3 followed by a band of H2AZ and KDM5B.
#ugh this computer really doesn't want to do large heatmaps... sigh.

lt_fc = intergenic_cgis$proseq_fc
names(lt_fc) = intergenic_cgis$CGI_ID
iCGIs_merged_mat$fc = unname((lt_fc[as.character(iCGIs_merged_mat$gene_ID)]))

lt_maxscore = intergenic_cgis$CGI_maxscore
names(lt_maxscore) = intergenic_cgis$CGI_ID
iCGIs_merged_mat$maxscore = unname((
  lt_maxscore[as.character(iCGIs_merged_mat$gene_ID)]))

ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,3.25)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,3.25)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.25)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

#########################################################################testing
test = iCGIs_merged_mat %>%
  filter(
    gene_ID == "CGI_100192759_100193274" |
      gene_ID =="CGI_100562075_100562638"|
      gene_ID == "CGI_39824286_39824504" |
      gene_ID == "CGI_61694582_61694802" |
      gene_ID == 'CGI_181169863_181170124'
    )
ggplot(data=test) +
  geom_raster(aes(
    x=Bins, y=gene_ID,fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.25)) +
  mytheme +
  geom_vline(aes(xintercept=0), linetype="dashed")
ggplot(data=test) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.25)) +
  mytheme +
  geom_vline(aes(xintercept=0), linetype="dashed")
#returned to this. Fold change seems to be sorting correctly when I don't say descending = TRUE.
ggplot(data=test) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.25)) +
  mytheme +
  geom_vline(aes(xintercept=0), linetype="dashed")
ggplot(data=test) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.25)) +
  mytheme +
  geom_vline(aes(xintercept=0), linetype="dashed")

#why is the signal always way upstream? Oh. 
#It's something to do with how I handled the minus strand genes... 
#They're the only ones trending horribly upstream. 
#Probably their anchor location is in a bad place??? 
#Did I do the math wrong when I was doing the whole maxbin thing?

#what does it look like if I only look at plus for now and fix later? 
#Oh I don't have strands in this dataset, ofc...

######################Better this time, don't need to strand the data.
ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,3.5)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,3.5)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.75)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")
summary(iCGIs_merged_mat$H2AZ)
ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.75)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.25)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.75)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")
###
lt_strand = genebody_cgis$CGI_imputed_strand
names(lt_strand) = genebody_cgis$CGI_ID
gbCGIs_merged_mat$strand = unname((
  lt_strand[as.character(gbCGIs_merged_mat$gene_ID)]))
gbplus = gbCGIs_merged_mat %>% filter(strand == "+")

lt_fc = genebody_cgis$proseq_fc
names(lt_fc) = genebody_cgis$CGI_ID
gbCGIs_merged_mat$fc = unname((lt_fc[as.character(gbCGIs_merged_mat$gene_ID)]))

lt_fc = genebody_cgis$CGI_maxscore
names(lt_fc) = genebody_cgis$CGI_ID
gbCGIs_merged_mat$maxscore = unname((lt_fc[as.character(gbCGIs_merged_mat$gene_ID)]))

quantile(gbCGIs_merged_mat$KDM5B, 0.95)
ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,4)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,4)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

lt_fc = tss_cgis$proseq_fc
names(lt_fc) = tss_cgis$CGI_ID
tssCGIs_merged_mat$fc = unname((lt_fc[as.character(tssCGIs_merged_mat$gene_ID)]))

lt_maxscore = tss_cgis$CGI_maxscore
names(lt_maxscore) = tss_cgis$CGI_ID
tssCGIs_merged_mat$maxscore = unname((lt_maxscore[as.character(tssCGIs_merged_mat$gene_ID)]))

quantile(tssCGIs_merged_mat$H3K4me3, 0.95)

ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,12)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,12)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.8)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.8)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.5)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.5)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed")

#I'm guessing my foldchange is biased towards the small ones. Let's keep with maxscore then.

quantile(iCGIs_merged_mat$H3K4m,0.9)





pdf(file="Stranded_CGIS2.pdf", width=4, height = 8)
#
ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,3.5)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (intergenic)")
#
ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.75)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (intergenic)")
#
ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.75)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (intergenic)")
#
ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,4)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (gene body)")
#
ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (gene body)")
#
ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (gene body)")
#
ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,12)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (TSS)")
#
ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.8)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (TSS)")
#
ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.5)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (TSS)")
dev.off()


pdf(file="Stranded_CGIS2.pdf", width=4, height = 8)
#
ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.5)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (intergenic)")
#
ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.1)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (intergenic)")
#
ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,0.6)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (intergenic)")
#
ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.7)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (gene body)")
#
ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,0.7)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (gene body)")
#
ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,0.6)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (gene body)")
#
ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,8)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (TSS)")
#
ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.25)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (TSS)")
#
ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, maxscore),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,0.9)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (TSS)")
dev.off()

pdf(file="Stranded_CGIS3.pdf", width=4, height = 8)
#
ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.5)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (intergenic)")
#
ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.1)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (intergenic)")
#
ggplot(data=iCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,0.6)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (intergenic)")
#
ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.7)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (gene body)")
#
ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,0.7)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (gene body)")
#
ggplot(data=gbCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,0.6)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (gene body)")
#
ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H3K4me3))  +
  scale_fill_gradient(name="H3K4me3", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,8)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (TSS)")
#
ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=H2AZ))  +
  scale_fill_gradient(name="H2AZ", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,1.25)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (TSS)")
#
ggplot(data=tssCGIs_merged_mat) +
  geom_raster(aes(
    x=Bins, y=reorder(gene_ID, fc),fill=KDM5B))  +
  scale_fill_gradient(name="KDM5B", low = "#FFFFFF",
                      high = "#012345", oob=scales::squish, limits=c(0,0.9)) +
  mytheme +
  theme(axis.text.y = element_blank()) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  xlab("distance from max transcription (bp)") +
  ylab("CgIs (TSS)")
dev.off()
