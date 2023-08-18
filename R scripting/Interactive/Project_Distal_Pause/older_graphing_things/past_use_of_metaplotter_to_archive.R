source("C:/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/metaplotter.R"); #note this requires "data.table", "BSgenome", "compiler", "rtracklayer", "tidyverse", and "Rsamtools" packages

###executing metaplotter

enableJIT(3) #this compiles functions as we go along, which makes them run faster.
###################################
## 3. Load in intersection table ##
###################################
groseq_tab_classes_new_short <- read.delim("~/danko_mcf7_subclone_proseq_classes.txt")
#####################################
## 5. Get the regions for metaplot ##
#####################################
#Vertino_classed
intersection_table = groseq_tab_classes_new_short
B7 = get_anchored_metaplot_tabs("MCF7_B7_fwd_10.bw", "MCF7_B7_rev_10.bw", 500, 500, class_tab=intersection_table, class_col=intersection_table$B7_class, bs=25)
H9_b7 = get_anchored_metaplot_tabs("MCF7_H9_fwd_10.bw", "MCF7_H9_rev_10.bw", 500, 500, class_tab=intersection_table, class_col=intersection_table$B7_class, bs=25)
C11_b7 = get_anchored_metaplot_tabs("MCF7_C11_fwd_10.bw", "MCF7_C11_rev_10.bw", 500, 500, class_tab=intersection_table, class_col=intersection_table$B7_class, bs=25)
G11_b7 = get_anchored_metaplot_tabs("MCF7_G11_fwd_10.bw", "MCF7_G11_rev_10.bw", 500, 500, class_tab=intersection_table, class_col=intersection_table$B7_class, bs=25)

H9 = get_anchored_metaplot_tabs("MCF7_H9_fwd_10.bw", "MCF7_H9_rev_10.bw", 500, 500, class_tab=intersection_table, class_col=intersection_table$H9_class, bs=25)
C11 = get_anchored_metaplot_tabs("MCF7_C11_fwd_10.bw", "MCF7_C11_rev_10.bw", 500, 500, class_tab=intersection_table, class_col=intersection_table$C11_class, bs=25)
G11 = get_anchored_metaplot_tabs("MCF7_G11_fwd_10.bw", "MCF7_G11_rev_10.bw", 500, 500, class_tab=intersection_table, class_col=intersection_table$G11_class, bs=25)

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

B7$subclone = rep("B7", dim(B7)[1])
H9_b7$subclone = rep("H9", dim(H9_b7)[1])
C11_b7$subclone = rep("C11", dim(C11_b7)[1])
G11_b7$subclone = rep("G11", dim(G11_b7)[1])

B7_classed = rbind(B7, H9_b7, C11_b7, G11_b7)
base_geom = ggplot(data=B7_classed,aes(x=bin_adjust, y=score, color=class, linetype=subclone)) +
  #geom_point(size=1) +
  geom_line(size=1) +
  mytheme 
tailored_geom = base_geom +
  facet_grid(~fct_relevel(anchor,"5_cpg","TSS","3_cpg",)) + 
  xlab("Distance from Anchor")+
  ylab("Mean proseq tags, B7 classes")+
  geom_vline(xintercept=0,color="grey",size=1, linetype=2)+
  geom_hline(yintercept = 0, color = "grey", size = 1, linetype=2) +
  theme(panel.grid.major=element_blank())+
  xlim(c(-500,500)) +
  scale_color_manual(values=c("dist"="plum4", "prox"="darkorange3", "silent"="darkgreen")) +
  scale_linetype_manual(values=c("B7"="solid", "C11"="dotted", "G11"="longdash", "H9"="dotdash"))
tailored_geom


#okay now a scaled version
B7_scaled = get_scaled_cgi_metaplot_tabs("MCF7_B7_fwd_10.bw", "MCF7_B7_rev_10.bw", class_tab=intersection_table, class_col=intersection_table$B7_class)
H9_b7_scaled = get_scaled_cgi_metaplot_tabs("MCF7_H9_fwd_10.bw", "MCF7_H9_rev_10.bw", class_tab=intersection_table, class_col=intersection_table$B7_class)
C11_b7_scaled = get_scaled_cgi_metaplot_tabs("MCF7_C11_fwd_10.bw", "MCF7_C11_rev_10.bw", class_tab=intersection_table, class_col=intersection_table$B7_class)
G11_scaled = get_scaled_cgi_metaplot_tabs("MCF7_G11_fwd_10.bw", "MCF7_G11_rev_10.bw", class_tab=intersection_table, class_col=intersection_table$B7_class)

B7_scaled$subclone = rep("B7", dim(B7_scaled)[1])
H9_b7_scaled$subclone = rep("H9", dim(H9_b7_scaled)[1])
C11_b7_scaled$subclone = rep("C11", dim(C11_b7_scaled)[1])
G11_scaled$subclone = rep("G11", dim(G11_scaled)[1])
B7_classed_scaled = rbind(B7_scaled, H9_b7_scaled, C11_b7_scaled, G11_scaled)

base_geom = ggplot(data=B7_classed_scaled,aes(x=bin_adjust, y=score, color=class, linetype=subclone)) +
  #geom_point(size=1.5) +
  geom_line(size=1) +
  mytheme 
tailored_geom_scaled = base_geom +
  geom_vline(xintercept=40,color="grey",size=1, linetype=2) +
  geom_vline(xintercept=80,color="grey",size=1, linetype=2) +
  geom_vline(xintercept=120,color="grey",size=1, linetype=2) +
  xlab("Scaled Islands")+
  ylab("Mean Tags")  +
  geom_label(label="-0.8 kb", aes(x=0, y=0), show.legend = FALSE)+
  geom_label(label="5_CpG", aes(x=40, y=0), show.legend = FALSE)+
  geom_label(label="TSS", aes(x=80, y=0), show.legend = FALSE)+
  geom_label(label="3_CpG", aes(x=120, y=0), show.legend = FALSE) +
  geom_label(label="+0.8 kb", aes(x=160, y=0), show.legend = FALSE) +
  scale_x_continuous(breaks=NULL)+
  scale_linetype_manual(values=c("B7"="solid", "C11"="dotted", "G11"="longdash", "H9"="dotdash"))+
  scale_color_manual(values=c("dist"="plum4", "prox"="darkorange3", "silent"="darkgreen"))
tailored_geom_scaled


pdf(file = "test.pdf")
tailored_geom_scaled
dev.off()


B7_bimodal_classed <- intersection_table %>% filter(B7_dist>=0.3 & B7_score>=0.3)
B7_bimodal_prox <- intersection_table %>% filter(B7_dist>=0.3 & B7_score>=0.3 & B7_score >= 2*B7_dist)
B7_bimodal_dist <- intersection_table %>% filter(B7_dist>=0.3 & B7_score>=0.3 & B7_dist >= 2*B7_score)
B7_tru_bimodal <- intersection_table %>% filter(B7_dist>=0.3 & B7_score>=0.3 & B7_dist < 2*B7_score & B7_score < 2 *B7_dist)
B7_prox <- intersection_table %>% filter(B7_score>0.3 & B7_dist<=B7_score & B7_dist<0.3)
B7_dist <- intersection_table %>% filter(B7_dist>0.3 & B7_dist>B7_score & B7_score<0.3)
B7_silent <- intersection_table %>% filter(B7_score<0.3 & B7_dist<0.3)

B7_dist$B7_class2 <- rep("dist", dim(B7_dist)[1])
B7_prox$B7_class2 <- rep("prox", dim(B7_prox)[1])
B7_silent$B7_class2 <- rep("silent", dim(B7_silent)[1])
B7_bimodal_prox$B7_class2 <- rep("bimod_prox", dim(B7_bimodal_prox)[1])
B7_bimodal_dist$B7_class2 <- rep("bimod_dist", dim(B7_bimodal_dist)[1])
B7_tru_bimodal$B7_class2 <- rep("tru_bimod", dim(B7_tru_bimodal)[1])
B7_classed <- rbind(B7_dist, B7_prox, B7_silent, B7_bimodal_prox, B7_bimodal_dist, B7_tru_bimodal)

######################################
B7 = get_anchored_metaplot_tabs("MCF7_B7_fwd_10.bw", "MCF7_B7_rev_10.bw", 250, 250, class_tab=B7_classed, class_col="B7_class2", bs=25)
B7_antisense = get_anchored_metaplot_tabs("MCF7_B7_rev_10.bw", "MCF7_B7_fwd_10.bw", 250, 250, class_tab=B7_classed, class_col="B7_class2", bs=25)
B7_antisense$score = B7_antisense$score*-1
B7$dir <- rep("sense", dim(B7)[1])
B7_antisense$dir <- rep("antisense", dim(B7_antisense)[1])
B7_graph = rbind(B7, B7_antisense)
base_geom = ggplot(data=B7_graph,aes(x=bin_adjust, y=score, color=class, linetype=dir)) +
  #geom_point(size=1) +
  geom_line(size=1) +
  mytheme 
tailored_geom = base_geom +
  facet_grid(~fct_relevel(anchor,"5_cpg","TSS","3_cpg",)) + 
  xlab("Distance from Anchor")+
  ylab("Mean proseq tags, B7 classes")+
  geom_vline(xintercept=0,color="grey",size=1, linetype=2)+
  geom_hline(yintercept = 0, color = "grey", size = 1, linetype=2) +
  theme(panel.grid.major=element_blank())+
  xlim(c(-250,250)) +
  scale_linetype_manual(values=c("sense"="solid", "antisense"="dotted"))
tailored_geom

B7_scaled = get_scaled_cgi_metaplot_tabs("MCF7_B7_fwd_10.bw", "MCF7_B7_rev_10.bw", class_tab=B7_classed, class_col="B7_class2")
B7_scaled_antisense = get_scaled_cgi_metaplot_tabs("MCF7_B7_rev_10.bw", "MCF7_B7_fwd_10.bw", class_tab=B7_classed, class_col="B7_class2")

B7_scaled$dir <- rep("sense", dim(B7_scaled)[1])
B7_scaled_antisense$dir <- rep("antisense", dim(B7_scaled_antisense)[1])
B7_scaled_antisense$score = B7_scaled_antisense$score*-1
B7_scaled_graph = rbind(B7_scaled, B7_scaled_antisense)

base_geom = ggplot(data=B7_scaled_graph,aes(x=bin_adjust, y=score, color=class, linetype=dir)) +
  #geom_point(size=0.5) +
  geom_line(size=1) +
  mytheme 
tailored_geom_scaled = base_geom +
  geom_vline(xintercept=40,color="grey",size=1, linetype=2) +
  geom_vline(xintercept=80,color="grey",size=1, linetype=2) +
  geom_vline(xintercept=120,color="grey",size=1, linetype=2) +
  xlab("Scaled Islands")+
  ylab("Mean Tags")  +
  geom_label(label="-0.8 kb", aes(x=0, y=0), show.legend = FALSE)+
  geom_label(label="5_CpG", aes(x=40, y=0), show.legend = FALSE)+
  geom_label(label="TSS", aes(x=80, y=0), show.legend = FALSE)+
  geom_label(label="3_CpG", aes(x=120, y=0), show.legend = FALSE) +
  geom_label(label="+0.8 kb", aes(x=160, y=0), show.legend = FALSE) +
  scale_x_continuous(breaks=NULL)+
  scale_linetype_manual(values=c("sense"="solid", "antisense"="dotted"))
tailored_geom_scaled
