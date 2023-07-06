###################
## Initial setup ##
###################
library(tidyverse)

intersection <- read.delim("curated_refseq_cpg_intersect_raw.txt", header=FALSE)
no_alt_chroms <- intersection %>% filter(!grepl("chrY|([\\w_]+)alt|random|fix|v1", V1)) %>% unique()
plus_strand_only <- no_alt_chroms %>%
  filter(V6=="+")
plus_TSS_intersect_only <- plus_strand_only %>%
  filter(V2>V8)
unique_tsses_plus <- plus_TSS_intersect_only %>% distinct(V2, .keep_all = TRUE)
one_tss_per_isle_plus <- unique_tsses_plus %>% group_by(V8) %>% add_tally() %>% filter(n==1)
more_than_one_tss_per_island <- unique_tsses_plus %>% group_by(V8) %>% add_tally() %>% filter(n!=1)
#info on num islands affected:
dim(more_than_one_tss_per_island %>% distinct(V8))
#this problem affects 1111 CpG islands on the plus strand. (if using refseq curated.)


minus = intersection %>% filter(V6=="-")
minus_TSS_intersect_only <- minus %>%
  filter(V9>V3) %>% unique()
unique_tsses_minus <- minus_TSS_intersect_only %>% distinct(V3, .keep_all = TRUE)
one_tss_per_isle_minus <- unique_tsses_minus %>% group_by(V8) %>% add_tally() %>% filter(n==1)
more_than_one_tss_per_island_minus <- unique_tsses_minus %>% group_by(V8) %>% add_tally() %>% filter(n!=1)
#info on num islands affected:
dim(more_than_one_tss_per_island_minus %>% distinct(V8))
#this problem affects 1201 CpG islands on the minus strand. (if using refseq curated.)

t1 = more_than_one_tss_per_island_minus %>% group_by(V8) %>% mutate(sd_tss = sd(V3))
t2 = t1 %>% distinct(V8, .keep_all=true)
summary(t2$sd_tss)
 #Min.  1st Qu.   Median     Mean  3rd Qu.     Max. #0.707   71.418  206.713  298.560  389.580 6030.207

###################################33
#Ensembl
ensembl_trim1 <- read.delim("C:/Users/kayle/Box/Vertinolab/McKayla Ford/ensembl_trim1.txt", header=FALSE)
#remove unnecessary columns:
#thickStart, thickEnd, reserved, blockCounts, blockSizes, chromStarts, oChromStart, oChromEnd, oStrand, oChromSize, oChromStarts, oSequence, oCDS, chromSize, match, misMatch, repMatch, nCount, seqType, srcChromStart, srcChromEnd, srcIdent, srcAligned
library(tidyverse)
ensembl_trim2 <- ensembl_trim1 %>% select(-c(7:25, 29:31))
no_alt_chroms <- ensembl_trim2 %>% filter(!grepl("chrY|([\\w_]+)alt|random|fix|v1", V1)) %>% unique()


plus_strand_only <- no_alt_chroms %>%
  filter(V6=="+")
unique_tsses_plus <- plus_strand_only %>% distinct(V2, .keep_all = TRUE)

minus_strand_only = no_alt_chroms %>% filter(V6=="-")
unique_tsses_minus <- minus_strand_only %>% distinct(V3, .keep_all = TRUE)

one_tss_per_isle_plus <- unique_tsses_plus %>% group_by(V42) %>% add_tally() %>% filter(n==1)
more_than_one_tss_per_island <- unique_tsses_plus %>% group_by(V42) %>% add_tally() %>% filter(n!=1)
#info on num islands affected:
dim(more_than_one_tss_per_island %>% distinct(V42))
#this problem affects 1111 CpG islands on the plus strand. (if using refseq curated.)


one_tss_per_isle_minus <- unique_tsses_minus %>% group_by(V42) %>% add_tally() %>% filter(n==1)
more_than_one_tss_per_island_minus <- unique_tsses_minus %>% group_by(V42) %>% add_tally() %>% filter(n!=1)
#info on num islands affected:
dim(more_than_one_tss_per_island_minus %>% distinct(V42))
#this problem affects 1201 CpG islands on the minus strand. (if using refseq curated.)

t1 = more_than_one_tss_per_island_minus %>% group_by(V42) %>% mutate(sd_tss = sd(V3))
t2 = t1 %>% distinct(V42, .keep_all=TRUE)
summary(t2$sd_tss)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.71     34.65     98.47   1054.48    244.04 168783.02