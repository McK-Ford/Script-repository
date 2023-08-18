cpgIslands <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/Data/cpgIsland.hg38.txt", header=FALSE)

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

cpgIslands$score = 0
cpgIslands$strand = "+"
cpgIslands$V4 = paste0("CGI_", cpgIslands$V2, "_", cpgIslands$V3)
plus_mat = lapply(plus_strand_list, get_score_matrix, bed=cpgIslands, n=1,
                     method="bi_stranded_anchored", pairedEnd=FALSE)
minus_mat = lapply(minus_strand_list, get_score_matrix, bed=cpgIslands, n=1,
                  method="bi_stranded_anchored", pairedEnd=FALSE)

cpgIslands_stranded = cpgIslands[1:4]
cpgIslands_stranded = cpgIslands_stranded %>% filter(!grepl("chrY|([\\w_]+)alt|random|fix|v1|v2", V1))
pm = cbind(plus_mat[[1]], plus_mat[[2]][[2]], plus_mat[[3]][[2]], plus_mat[[4]][[2]])
cpgIslands_stranded2 = merge(x=cpgIslands_stranded, y=pm, by.x="V4", by.y="rn", all=TRUE)

mm = cbind(minus_mat[[1]], minus_mat[[2]][[2]], minus_mat[[3]][[2]], minus_mat[[4]][[2]])
cpgIslands_stranded3 = merge(x=cpgIslands_stranded2, y=mm, by.x="V4", by.y="rn", all=TRUE)
colnames(cpgIslands_stranded3) = c("ID", "chrom", "start", "end", "p1","p2","p3","p4","m1","m2","m3","m4")
cpgIslands_stranded3$plus = cpgIslands_stranded3$p1 + cpgIslands_stranded3$p2 +
  cpgIslands_stranded3$p3 + cpgIslands_stranded3$p4
cpgIslands_stranded3$minus = cpgIslands_stranded3$m1 + cpgIslands_stranded3$m2 +
  cpgIslands_stranded3$m3 + cpgIslands_stranded3$m4
#now how am I going to filter these?
summary(cpgIslands_stranded3$plus)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.     
#   0.00    0.00    2.10   31.80   33.87 2949.83  
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.       
#   0.0000    0.0615    3.4093   32.5872   41.9621 2226.0639   
#what are the NAs?
summary(cpgIslands_stranded3$minus+cpgIslands_stranded3$plus)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.  
#   0.000    0.149    6.783   64.385   81.469 5175.509   
#lets say greater than 6 reads in it?
CpGIslands_silent = cpgIslands_stranded3 %>% filter(minus+plus<6)
CpGIslands_plus = cpgIslands_stranded3 %>% filter(minus+plus>6 & plus>minus)
# a good portion of these are likely genic of course
CpGIslands_minus = cpgIslands_stranded3 %>% filter(minus+plus>6 & plus<minus)
#okay next I would need to see if they're associated with genes. Strength of dif is also an ordering option...
CpGIslands_plus$strand = "+"
CpGIslands_plus$fc = (CpGIslands_plus$plus+1)/(CpGIslands_plus$minus+1)
CpGIslands_minus$strand = "-"
CpGIslands_minus$fc = (CpGIslands_minus$minus+1)/(CpGIslands_minus$plus+1)
cpgIslands4 = rbind(CpGIslands_minus, CpGIslands_plus)
cgi_bed = cpgIslands4 %>% select(c(2:4, 1, 16, 15))
write.table(cgi_bed, file="strandedCGIs.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
