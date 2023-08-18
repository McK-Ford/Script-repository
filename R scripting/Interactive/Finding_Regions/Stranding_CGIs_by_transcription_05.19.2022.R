#### Setup ####
cpgIslands <- read.delim(
  "C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/regions/cpgIsland.hg38.txt",
  header=FALSE
  ) #standard no modification islands text
source(
  "/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/NewMetaplotterJoshBased.R"
  )
enableJIT(3)

nasDir="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/PublicData/GSE93229_Danko_2018_MCF7_PROseq/"
plus_strand_list=list("MCF7_H9_plus.bam",
                      "MCF7_B7_plus.bam",
                      "MCF7_G11_plus.bam",
                      "MCF7_C11_plus.bam")
plus_strand_list=paste0(nasDir, plus_strand_list)
minus_strand_list=list("MCF7_H9_minus.bam",
                       "MCF7_B7_minus.bam",
                       "MCF7_G11_minus.bam",
                       "MCF7_C11_minus.bam")
minus_strand_list=paste0(nasDir, minus_strand_list)

#make the CGI file into a bed
cpgIslands$score = 0
cpgIslands$strand = "+"
cpgIslands$V4 = paste0("CGI_", cpgIslands$V2, "_", cpgIslands$V3) #unique ID

#### Get plus & minus transcription across island ####
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

#### assign these CGI a strand based on this transcription ####
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


##### get location of maximal transcription to be anchor point ####
source("/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/NewMetaplotterJoshBased.R")
enableJIT(3)
nasDir="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/PublicData/GSE93229_Danko_2018_MCF7_PROseq/"
plus_strand_list=list("MCF7_H9_plus.bam",
                      "MCF7_B7_plus.bam",
                      "MCF7_G11_plus.bam",
                      "MCF7_C11_plus.bam")
plus_strand_list=paste0(nasDir, plus_strand_list)
minus_strand_list=list("MCF7_H9_minus.bam",
                       "MCF7_B7_minus.bam",
                       "MCF7_G11_minus.bam",
                       "MCF7_C11_minus.bam")
minus_strand_list=paste0(nasDir, minus_strand_list)
strandedCGIs <- read.delim("strandedCGIs.bed", header=FALSE)

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
#### get maximum score, + in case of ties, minimum bin ####
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

##################################
geneDir="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/Regions/"
refseq_curated_longest_cgi_minus <- read.delim(paste0(
  geneDir,
  "refseq_curated_longest_cgi_minus.hg38.txt"))
refseq_curated_longest_cgi_plus <- read.delim(paste0(
  geneDir,
  "refseq_curated_longest_cgi_plus.hg38.txt"))

TSS_CpG_end_plus= refseq_curated_longest_cgi_plus %>%
  select(1,5,3,8,2,7) #using CpG start as score placeholder
TSS_CpG_end_plus$name=paste0(TSS_CpG_end_plus$name, "_",
                             TSS_CpG_end_plus$cpg_e)
TSS_CpG_end_minus= refseq_curated_longest_cgi_minus %>% select(1,2,6,8,3,7)
#cgi end placeholder
TSS_CpG_end_minus$name=paste0(TSS_CpG_end_minus$name, "_",
                              TSS_CpG_end_minus$cpg_s)

name_vec = c("chrom", "start", "end", "name", "score", "strand")
colnames(TSS_CpG_end_minus) = name_vec
colnames(TSS_CpG_end_plus) = name_vec
TSS_CpG_end_minus$start=TSS_CpG_end_minus$end-1
TSS_CpG_end_plus$end=TSS_CpG_end_plus$end+1
TSS_CpG_end = rbind(TSS_CpG_end_minus, TSS_CpG_end_plus)
write.table(TSS_CpG_end, file="TSS_only.bed", row.names=FALSE,
            col.names=FALSE, quote=FALSE, sep="\t")
