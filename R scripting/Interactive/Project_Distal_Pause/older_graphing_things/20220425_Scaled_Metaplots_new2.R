#Running NewMetaplotterJoshBased.R
source("/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/NewMetaplotterJoshBased.R")
enableJIT(3)
##############
## plotting ##
##############
bamDirectory="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/CutnTag_n_chip/Vertino lab data/VER00046 2021.11.17 Cut and Tag/bams/hg38_dedup/"

refseq_curated_longest_cgi_minus <- read.delim("/Users/kayle/Box/Vertinolab/McKayla Ford/Data/refseq_curated_longest_cgi_minus.hg38.txt")
refseq_curated_longest_cgi_plus <- read.delim("/Users/kayle/Box/Vertinolab/McKayla Ford/Data/refseq_curated_longest_cgi_plus.hg38.txt")

TSS_CpG_end_plus= refseq_curated_longest_cgi_plus %>%
  select(1,5,3,8,2,7) #using CpG start as score placeholder
CpG_start_TSS_plus= refseq_curated_longest_cgi_plus %>% select(1,2,5,8,3,7) #cgi end placeholder
TSS_CpG_end_plus$name=paste0(TSS_CpG_end_plus$name, "_", TSS_CpG_end_plus$cpg_e)
CpG_start_TSS_plus$name=paste0(CpG_start_TSS_plus$name, "_", CpG_start_TSS_plus$cpg_e)

TSS_CpG_end_minus= refseq_curated_longest_cgi_minus %>% select(1,2,6,8,3,7) #cgi end placeholder
CpG_start_TSS_minus= refseq_curated_longest_cgi_minus %>% select(1,6,3,8,2,7) #cgi start placeholder
TSS_CpG_end_minus$name=paste0(TSS_CpG_end_minus$name, "_", TSS_CpG_end_minus$cpg_s)
CpG_start_TSS_minus$name=paste0(CpG_start_TSS_minus$name, "_", CpG_start_TSS_minus$cpg_s)

name_vec = c("chrom", "start", "end", "name", "score", "strand")
colnames(TSS_CpG_end_minus) = name_vec
colnames(TSS_CpG_end_plus) = name_vec
TSS_CpG_end = rbind(TSS_CpG_end_minus, TSS_CpG_end_plus)

colnames(CpG_start_TSS_minus) = name_vec
colnames(CpG_start_TSS_plus) = name_vec
CpG_start_TSS = rbind(CpG_start_TSS_minus, CpG_start_TSS_plus)

#I wonder if it's a length issue - what if I picked a number of bins that would lead to a bin size closer to what's in the flanking regions?
# Ex - median TSS-CpG end dist is about 500, if I did 50 bins that would make them average 10 bp each instead of the 5 they're currently doing.
# Meanwhile for CpG_start - TSS, the median is 300. So if we want to hit something that'll work for both of them, what about 40 bins?

notstranded_bam_list = list(paste0(
  bamDirectory, "ZS1_19_MCF7_EV_H2AZ.dedup.bam"),
  paste0(
    bamDirectory, "ZS1_19_MCF7_EV_H3K4me3.dedup.bam"),
  paste0(
    bamDirectory, "ZS1_19_MCF7_EV_KDM5B.dedup.bam")
)

bed_list=list(CpG_start_TSS, TSS_CpG_end)

cgi_plot_nostranding = function(bam_list, bed,
                                n=40, pairedEnd=TRUE){ #still a terrible function
  mat_list = lapply(bam_list, get_score_matrix, bed=bed, n=n,
                    method="bi_stranded_anchored", pairedEnd=pairedEnd)
  metaplot_list <- lapply(mat_list, data.table::melt, measure.vars=c(2:41),
                          variable.name="Bins", value.name="Score")
  }

region_lists=lapply(bed_list, cgi_plot_nostranding, bam_list=notstranded_bam_list)

cpg_start_anchor_plus= refseq_curated_longest_cgi_plus %>%
  select(1,2,3,8,5,7)
cpg_end_anchor_plus= refseq_curated_longest_cgi_plus %>% select(1,3,6,8,5,7)
cpg_start_anchor_plus$name=paste0(cpg_start_anchor_plus$name, "_", cpg_start_anchor_plus$cpg_e)
cpg_end_anchor_plus$name=paste0(cpg_end_anchor_plus$name, "_", cpg_end_anchor_plus$cpg_e)

cpg_start_anchor_minus= refseq_curated_longest_cgi_minus %>%
  select(1,2,3,8,5,7)
cpg_end_anchor_minus= refseq_curated_longest_cgi_minus %>% select(1,5,2,8,6,7)
cpg_start_anchor_minus$name=paste0(cpg_start_anchor_minus$name, "_", cpg_start_anchor_minus$cpg_s)
cpg_end_anchor_minus$name=paste0(cpg_end_anchor_minus$name, "_", cpg_end_anchor_minus$cpg_s)


bed_list_b=list(cpg_start_anchor_plus, cpg_start_anchor_minus)
bed_list_a=list(cpg_end_anchor_plus, cpg_end_anchor_minus)

cgi_plot_nostranding_noscale = function(bam_list, bed, dist=800, b_or_a, pairedEnd=TRUE){
  if (b_or_a=="b"){
    mat_list = lapply(bam_list, get_score_matrix, bed=bed, b=-dist, a=0,
                      method="single_stranded_anchored", pairedEnd=pairedEnd)
  }
  else {
    mat_list = lapply(bam_list, get_score_matrix, bed=bed, b=0, a=dist,
                      method="single_stranded_anchored", pairedEnd=pairedEnd)
  }
  metaplot_list <- lapply(mat_list, data.table::melt, measure.vars=c(2:81),
                          variable.name="Bins", value.name="Score")
}

###############################################################################################
region_lists_noscale_b=lapply(bed_list_b, cgi_plot_nostranding_noscale, b_or_a="b", bam_list=notstranded_bam_list)
region_lists_noscale_a=lapply(bed_list_a, cgi_plot_nostranding_noscale, b_or_a="a", bam_list=notstranded_bam_list)
# Alright, what's the landscape here? We've got startplus (bams), startminus (bams). Same for end. The other
# section is instead cpg_start_tss (bams), tss cpgend (bams).

upstream_p = cbind(region_lists_noscale_b[[1]][[1]],
                   region_lists_noscale_b[[1]][[2]][[3]],
                   region_lists_noscale_b[[1]][[3]][[3]])
upstream_m = cbind(region_lists_noscale_b[[2]][[1]],
                   region_lists_noscale_b[[2]][[2]][[3]],
                   region_lists_noscale_b[[2]][[3]][[3]])
upstream = rbind(upstream_m, upstream_p)

downstream_p = cbind(region_lists_noscale_a[[1]][[1]],
                   region_lists_noscale_a[[1]][[2]][[3]],
                   region_lists_noscale_a[[1]][[3]][[3]])
downstream_m = cbind(region_lists_noscale_a[[2]][[1]],
                   region_lists_noscale_a[[2]][[2]][[3]],
                   region_lists_noscale_a[[2]][[3]][[3]])
downstream = rbind(downstream_m, downstream_p)

CpG_start_TSS_scores = cbind(region_lists[[1]][[1]],
                             region_lists[[1]][[2]][[3]],
                             region_lists[[1]][[3]][[3]])
TSS_CpG_end_scores = cbind(region_lists[[2]][[1]],
                           region_lists[[2]][[2]][[3]],
                           region_lists[[2]][[3]][[3]])
upstream_flat = upstream[,.(H2AZ=mean(Score), H3K4me3=mean(V2),
                            KDM5B=mean(V3)),by=Bins]
downstream_flat = downstream[,.(H2AZ=mean(Score), H3K4me3=mean(V2),
                            KDM5B=mean(V3)),by=Bins]
CpG_start_TSS_scores_flat = CpG_start_TSS_scores[,.(H2AZ=mean(Score), H3K4me3=mean(V2),
                            KDM5B=mean(V3)),by=Bins]
TSS_CpG_end_scores_flat = TSS_CpG_end_scores[,.(H2AZ=mean(Score), H3K4me3=mean(V2),
                            KDM5B=mean(V3)),by=Bins]

colnames(TSS_CpG_end_scores_flat)=c("Bins", "H2AZ", "H3K4me3",
                               "KDM5B")
colnames(CpG_start_TSS_scores_flat)=c("Bins", "H2AZ", "H3K4me3",
                                 "KDM5B")
colnames(upstream_flat)=c("Bins", "H2AZ", "H3K4me3",
                               "KDM5B")
colnames(downstream_flat)=c("Bins", "H2AZ", "H3K4me3",
                                 "KDM5B")
#total bins are 180 on each side.
#so bin_n.
### ****
upstream_flat$binN = seq(from=-120, to=-41, by=1)
CpG_start_TSS_scores_flat$binN = seq(from=-40, to=-1, by=1)
TSS_CpG_end_scores_flat$binN = seq(from=0, to=39, by=1)
downstream_flat$binN = seq(from=40, to=119, by=1)
CGI_scores=rbind(upstream_flat, CpG_start_TSS_scores_flat,
                 TSS_CpG_end_scores_flat, downstream_flat)
#CGI_scores=CGI_scores %>%   ##divide by max to normalize y-axis
#  mutate_at(vars(-Bins), ~(./max(.)))

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

tailored_geom2 = ggplot(data=CGI_scores,aes(
  x=binN, y=H2AZ)) +
  mytheme +
  geom_point(size=1) +
  xlab("Distance from Anchor")+
  ylab("H2AZ Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept = -41, linetype="dashed") +
  geom_vline(xintercept = 39, linetype="dashed")
tailored_geom2

tailored_geom2 = ggplot(data=CGI_scores,aes(
  x=binN, y=H3K4me3)) +
  mytheme +
  geom_point(size=1) +
  xlab("Distance from Anchor")+
  ylab("H3K4me3 Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept = -41, linetype="dashed") +
  geom_vline(xintercept = 40, linetype="dashed")
tailored_geom2

tailored_geom2 = ggplot(data=CGI_scores,aes(
  x=binN, y=KDM5B)) +
  mytheme +
  geom_point(size=1) +
  xlab("Distance from Anchor")+
  ylab("KDM5B Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept = -41, linetype="dashed") +
  geom_vline(xintercept = 40, linetype="dashed")
tailored_geom2

CGI_scores_melt = CGI_scores %>% 
  select(-c(1)) %>% pivot_longer(cols=c(H2AZ, H3K4me3, KDM5B))

tailored_geom2 = ggplot(data=CGI_scores_melt, aes(
  x=binN, y=value, color=name)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept = -41, linetype="dashed") +
  geom_vline(xintercept = 40, linetype="dashed")
tailored_geom2

tailored_geom2 = ggplot(data=CGI_scores_melt, aes(
  x=binN, y=value, color=name)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  ylim(0.3,1.5)+
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept = -41, linetype="dashed") +
  geom_vline(xintercept = 40, linetype="dashed")
tailored_geom2
#############################
notstranded_bam_list = list(paste0(
  bamDirectory, "ZS1_19_MCF7_EV_H2AZ.dedup.bam"),
  paste0(
    bamDirectory, "ZS1_19_MCF7_EV_H3K4me3.dedup.bam"),
  paste0(
    bamDirectory, "ZS1_19_MCF7_EV_KDM5B.dedup.bam"),
  paste0(
    bamDirectory, "ZS1_19_MCF7_SH1_H3K4me3.dedup.bam"
  ),
  paste0(
    bamDirectory, "ZS1_19_MCF7_SH2_H3K4me3.dedup.bam"
  )
) #***

bed_list=list(CpG_start_TSS, TSS_CpG_end)

cgi_plot_nostranding = function(bam_list, bed,
                                n=100, pairedEnd=TRUE){ #still a terrible function
  mat_list = lapply(bam_list, get_score_matrix, bed=bed, n=n,
                    method="bi_stranded_anchored", pairedEnd=pairedEnd)
  metaplot_list <- lapply(mat_list, data.table::melt, measure.vars=c(2:101),
                          variable.name="Bins", value.name="Score")
}

region_lists=lapply(bed_list, cgi_plot_nostranding, bam_list=notstranded_bam_list)


cpg_start_anchor_plus= refseq_curated_longest_cgi_plus %>%
  select(1,2,3,8,5,7)
cpg_end_anchor_plus= refseq_curated_longest_cgi_plus %>% select(1,3,6,8,5,7)
cpg_start_anchor_plus$name=paste0(cpg_start_anchor_plus$name, "_", cpg_start_anchor_plus$cpg_e)
cpg_end_anchor_plus$name=paste0(cpg_end_anchor_plus$name, "_", cpg_end_anchor_plus$cpg_e)

cpg_start_anchor_minus= refseq_curated_longest_cgi_minus %>%
  select(1,2,3,8,5,7)
cpg_end_anchor_minus= refseq_curated_longest_cgi_minus %>% select(1,5,2,8,6,7)
cpg_start_anchor_minus$name=paste0(cpg_start_anchor_minus$name, "_", cpg_start_anchor_minus$cpg_s)
cpg_end_anchor_minus$name=paste0(cpg_end_anchor_minus$name, "_", cpg_end_anchor_minus$cpg_s)


bed_list_b=list(cpg_start_anchor_plus, cpg_start_anchor_minus)
bed_list_a=list(cpg_end_anchor_plus, cpg_end_anchor_minus)

cgi_plot_nostranding_noscale = function(bam_list, bed, dist=800, b_or_a, pairedEnd=TRUE){
  if (b_or_a=="b"){
    mat_list = lapply(bam_list, get_score_matrix, bed=bed, b=-dist, a=0,
                      method="single_stranded_anchored", pairedEnd=pairedEnd)
  }
  else {
    mat_list = lapply(bam_list, get_score_matrix, bed=bed, b=0, a=dist,
                      method="single_stranded_anchored", pairedEnd=pairedEnd)
  }
  metaplot_list <- lapply(mat_list, data.table::melt, measure.vars=c(2:81),
                          variable.name="Bins", value.name="Score")
}

###############################################################################################
region_lists_noscale_b=lapply(bed_list_b, cgi_plot_nostranding_noscale, b_or_a="b", bam_list=notstranded_bam_list)
region_lists_noscale_a=lapply(bed_list_a, cgi_plot_nostranding_noscale, b_or_a="a", bam_list=notstranded_bam_list)
# Alright, what's the landscape here? We've got startplus (bams), startminus (bams). Same for end. The other
# section is instead cpg_start_tss (bams), tss cpgend (bams).

upstream_p = cbind(region_lists_noscale_b[[1]][[1]],
                   region_lists_noscale_b[[1]][[2]][[3]],
                   region_lists_noscale_b[[1]][[3]][[3]],
                   region_lists_noscale_b[[1]][[4]][[3]],
                   region_lists_noscale_b[[1]][[5]][[3]])
upstream_m = cbind(region_lists_noscale_b[[2]][[1]],
                   region_lists_noscale_b[[2]][[2]][[3]],
                   region_lists_noscale_b[[2]][[3]][[3]],
                   region_lists_noscale_b[[2]][[4]][[3]],
                   region_lists_noscale_b[[2]][[5]][[3]])
upstream = rbind(upstream_m, upstream_p)

downstream_p = cbind(region_lists_noscale_a[[1]][[1]],
                     region_lists_noscale_a[[1]][[2]][[3]],
                     region_lists_noscale_a[[1]][[3]][[3]],
                     region_lists_noscale_a[[1]][[4]][[3]],
                     region_lists_noscale_a[[1]][[5]][[3]])
downstream_m = cbind(region_lists_noscale_a[[2]][[1]],
                     region_lists_noscale_a[[2]][[2]][[3]],
                     region_lists_noscale_a[[2]][[3]][[3]],
                     region_lists_noscale_a[[2]][[4]][[3]],
                     region_lists_noscale_a[[2]][[5]][[3]])
downstream = rbind(downstream_m, downstream_p)

CpG_start_TSS_scores = cbind(region_lists[[1]][[1]],
                             region_lists[[1]][[2]][[3]],
                             region_lists[[1]][[3]][[3]],
                             region_lists[[1]][[4]][[3]],
                             region_lists[[1]][[5]][[3]])
TSS_CpG_end_scores = cbind(region_lists[[2]][[1]],
                           region_lists[[2]][[2]][[3]],
                           region_lists[[2]][[3]][[3]],
                           region_lists[[2]][[4]][[3]],
                           region_lists[[2]][[5]][[3]])
upstream_flat = upstream[,.(H2AZ=mean(Score), H3K4me3=mean(V2),
                            KDM5B=mean(V3), H3K4me3_SH1=mean(V4), H3K4me3_SH2=mean(V5)),by=Bins]
downstream_flat = downstream[,.(H2AZ=mean(Score), H3K4me3=mean(V2),
                                KDM5B=mean(V3), H3K4me3_SH1=mean(V4), H3K4me3_SH2=mean(V5)),by=Bins]
CpG_start_TSS_scores_flat = CpG_start_TSS_scores[,.(H2AZ=mean(Score), H3K4me3=mean(V2),
                                                    KDM5B=mean(V3), H3K4me3_SH1=mean(V4), H3K4me3_SH2=mean(V5)),by=Bins]
TSS_CpG_end_scores_flat = TSS_CpG_end_scores[,.(H2AZ=mean(Score), H3K4me3=mean(V2),
                                                KDM5B=mean(V3), H3K4me3_SH1=mean(V4), H3K4me3_SH2=mean(V5)),by=Bins]

colnames(TSS_CpG_end_scores_flat)=c("Bins", "H2AZ", "H3K4me3",
                                    "KDM5B", "H3K4me3_SH1", "H3K4me3_SH2")
colnames(CpG_start_TSS_scores_flat)=c("Bins", "H2AZ", "H3K4me3",
                                      "KDM5B", "H3K4me3_SH1", "H3K4me3_SH2")
colnames(upstream_flat)=c("Bins", "H2AZ", "H3K4me3",
                          "KDM5B", "H3K4me3_SH1", "H3K4me3_SH2")
colnames(downstream_flat)=c("Bins", "H2AZ", "H3K4me3",
                            "KDM5B", "H3K4me3_SH1", "H3K4me3_SH2")
#total bins are 180 on each side.
#so bin_n.
### ****
upstream_flat$binN = seq(from=-180, to=-101, by=1)
CpG_start_TSS_scores_flat$binN = seq(from=-100, to=-1, by=1)
TSS_CpG_end_scores_flat$binN = seq(from=0, to=99, by=1)
downstream_flat$binN = seq(from=100, to=179, by=1)
CGI_scores=rbind(upstream_flat, CpG_start_TSS_scores_flat,
                 TSS_CpG_end_scores_flat, downstream_flat)
#CGI_scores=CGI_scores %>%   ##divide by max to normalize y-axis
#  mutate_at(vars(-Bins), ~(./max(.)))

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


CGI_scores_melt = CGI_scores %>% 
  select(-c(1)) %>% pivot_longer(cols=c(H2AZ, H3K4me3, KDM5B, H3K4me3_SH1, H3K4me3_SH2))

tailored_geom2 = ggplot(data=CGI_scores_melt, aes(
  x=binN, y=value, color=name)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept = -100, linetype="dashed") +
  geom_vline(xintercept = 100, linetype="dashed")
tailored_geom2

tailored_geom2 = ggplot(data=CGI_scores_melt, aes(
  x=binN, y=value, color=name)) +
  mytheme +
  geom_line(size=1) +
  xlab("Distance from Anchor")+
  ylab("Tags") +
  ylim(0.3,1.5)+
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept = -41, linetype="dashed") +
  geom_vline(xintercept = 40, linetype="dashed")
tailored_geom2
