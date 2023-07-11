danko_mcf7_subclones_proseq_short <- read.delim("~/danko_mcf7_subclones_proseq_short.txt")
#first filter
library(tidyverse)
danko_trimmed_1 <- danko_mcf7_subclones_proseq_short %>% filter(B7_dist>=0.3 | H9_dist>=0.3 | C11_dist>=0.3 | G11_dist>=0.3) #7507
danko_trimmed_2 <- danko_mcf7_subclones_proseq_short %>% filter(B7_dist>=0.1 | H9_dist>=0.1 | C11_dist>=0.1 | G11_dist>=0.1) #9522
danko_trimmed_2_B7 <- danko_mcf7_subclones_proseq_short %>% filter(B7_score>=0.1 | B7_dist>=0.1)
danko_trimmed_2_H9 <- danko_mcf7_subclones_proseq_short %>% filter(H9_score>=0.1 | H9_dist>=0.1)
danko_trimmed_2_C11 <- danko_mcf7_subclones_proseq_short %>% filter(C11_score>=0.1 | C11_dist>=0.1)
danko_trimmed_2_G11 <- danko_mcf7_subclones_proseq_short %>% filter(G11_score>=0.1 | G11_dist>=0.1)

danko_trimmed_3 <- danko_mcf7_subclones_proseq_short %>% filter(B7_dist>=0.1 | H9_dist>=0.1 | C11_dist>=0.1 | G11_dist>=0.1 | B7_score>= 0.1 | H9_score>=0.1 | C11_score>=0.1 | G11_score>=0.1)
#okay, 10640 of the 13382 genes have some signal.

B7_distal_classed <- danko_trimmed_3 %>% filter(B7_dist>=0.1 & B7_dist>=B7_score)
B7_proximal_classed <- danko_trimmed_3 %>% filter(B7_score>=0.1 & B7_dist<=B7_score)
B7_silent_classed <- danko_trimmed_3 %>% filter(B7_score<0.1 & B7_dist<0.1)
B7_distal_classed$B7_class <- rep("dist", dim(B7_distal_classed)[1])
B7_proximal_classed$B7_class <- rep("prox", dim(B7_proximal_classed)[1])
B7_silent_classed$B7_class <- rep("silent", dim(B7_silent_classed)[1])
B7_classed <- rbind(B7_proximal_classed, B7_distal_classed, B7_silent_classed)

H9_distal_classed <- danko_trimmed_3 %>% filter(H9_dist>=0.1 & H9_dist>=H9_score)
H9_proximal_classed <- danko_trimmed_3 %>% filter(H9_score>=0.1 & H9_dist<=H9_score)
H9_silent_classed <- danko_trimmed_3 %>% filter(H9_score<0.1 & H9_dist<0.1)
H9_distal_classed$H9_class <- rep("dist", dim(H9_distal_classed)[1])
H9_proximal_classed$H9_class <- rep("prox", dim(H9_proximal_classed)[1])
H9_silent_classed$H9_class <- rep("silent", dim(H9_silent_classed)[1])
H9_classed <- rbind(H9_proximal_classed, H9_distal_classed, H9_silent_classed)

C11_distal_classed <- danko_trimmed_3 %>% filter(C11_dist>=0.1 & C11_dist>=C11_score)
C11_proximal_classed <- danko_trimmed_3 %>% filter(C11_score>=0.1 & C11_dist<=C11_score)
C11_silent_classed <- danko_trimmed_3 %>% filter(C11_score<0.1 & C11_dist<0.1)
C11_distal_classed$C11_class <- rep("dist", dim(C11_distal_classed)[1])
C11_proximal_classed$C11_class <- rep("prox", dim(C11_proximal_classed)[1])
C11_silent_classed$C11_class <- rep("silent", dim(C11_silent_classed)[1])
C11_classed <- rbind(C11_proximal_classed, C11_distal_classed, C11_silent_classed)

G11_distal_classed <- danko_trimmed_3 %>% filter(G11_dist>=0.1 & G11_dist>=G11_score)
G11_proximal_classed <- danko_trimmed_3 %>% filter(G11_score>=0.1 & G11_dist<=G11_score)
G11_silent_classed <- danko_trimmed_3 %>% filter(G11_score<0.1 & G11_dist<0.1)
G11_distal_classed$G11_class <- rep("dist", dim(G11_distal_classed)[1])
G11_proximal_classed$G11_class <- rep("prox", dim(G11_proximal_classed)[1])
G11_silent_classed$G11_class <- rep("silent", dim(G11_silent_classed)[1])
G11_classed <- rbind(G11_proximal_classed, G11_distal_classed, G11_silent_classed)

B7_lt = B7_classed$B7_class
names(B7_lt)=B7_classed$gene_name
danko_trimmed_3$B7_class =unname((B7_lt[as.character(danko_trimmed_3$gene_name)]))

H9_lt = H9_classed$H9_class
names(H9_lt)=H9_classed$gene_name
danko_trimmed_3$H9_class =unname((H9_lt[as.character(danko_trimmed_3$gene_name)]))

C11_lt = C11_classed$C11_class
names(C11_lt)=C11_classed$gene_name
danko_trimmed_3$C11_class =unname((C11_lt[as.character(danko_trimmed_3$gene_name)]))

G11_lt = G11_classed$G11_class
names(G11_lt)=G11_classed$gene_name
danko_trimmed_3$G11_class =unname((G11_lt[as.character(danko_trimmed_3$gene_name)]))

danko_trimmed_silent <- danko_mcf7_subclones_proseq_short %>% filter(B7_dist<=0.1 & H9_dist<=0.1 & C11_dist<=0.1 & G11_dist<=0.1 & B7_score<= 0.1 & H9_score<=0.1 & C11_score<=0.1 & G11_score<=0.1)
danko_trimmed_silent$B7_class <- rep("silent", dim(danko_trimmed_silent)[1])
danko_trimmed_silent$H9_class <- rep("silent", dim(danko_trimmed_silent)[1])
danko_trimmed_silent$C11_class <- rep("silent", dim(danko_trimmed_silent)[1])
danko_trimmed_silent$G11_class <- rep("silent", dim(danko_trimmed_silent)[1])

classed = rbind(danko_trimmed_3, danko_trimmed_silent)

write.table(
  classed,
  file="danko_mcf7_subclone_proseq_classes.txt",
  row.names=FALSE,
  quote=FALSE,
  sep="\t"
)
