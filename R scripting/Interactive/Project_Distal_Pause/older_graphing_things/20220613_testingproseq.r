#am I doing single basepair or full read pileup? Bc I think
#I might be doing full read which wouldn't necessarily be
#correct for proseq, since this is based off bams and bams are the full read.

source("/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/NewMetaplotterJoshBased.R")
enableJIT(3)

geneDir="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/Regions/"
refseq_curated_longest_cgi_minus <- read.delim(paste0(
  geneDir,
  "refseq_curated_longest_cgi_minus.hg38.txt"))
refseq_curated_longest_cgi_plus <- read.delim(paste0(
  geneDir,
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

nasDir="/Users/kayle/Box/Vertinolab/McKayla Ford/Data/NascentSeq/PublicData/GSE93229_Danko_2018_MCF7_PROseq/"
plus_strand=paste0(nasDir, "MCF7_B7_plus.bam")
minus_strand=paste0(nasDir, "MCF7_B7_minus.bam")

V1_fulllen_plus = get_score_matrix(bam=plus_strand, bed=TSS_CpG_end_plus,
                                   b=-1000, a=1000, method="single_stranded_anchored",
                                   pairedEnd = FALSE)
V1_fulllen_minus = get_score_matrix(bam=minus_strand, bed=TSS_CpG_end_minus,
                                   b=-1000, a=1000, method="single_stranded_anchored",
                                   pairedEnd = FALSE)
"C:/Users/kayle/Box/"
#what is the alteration I need to make?
#"bam_aln=resize(granges(bam_aln), width=1, fix="end")"

V2_endanchored_plus = get_score_matrix(bam=plus_strand, bed=TSS_CpG_end_plus,
                                   b=-1000, a=1000, method="single_stranded_anchored",
                                   pairedEnd = FALSE)
V2_endanchored_minus = get_score_matrix(bam=minus_strand, bed=TSS_CpG_end_minus,
                                    b=-1000, a=1000, method="single_stranded_anchored",
                                    pairedEnd = FALSE)

#just in case we'll do start as well bc I'm easily confused when it comes to proseq
V3_startanchored_plus = get_score_matrix(bam=plus_strand, bed=TSS_CpG_end_plus,
                                       b=-1000, a=1000, method="single_stranded_anchored",
                                       pairedEnd = FALSE)
V3_startanchored_minus = get_score_matrix(bam=minus_strand, bed=TSS_CpG_end_minus,
                                        b=-1000, a=1000, method="single_stranded_anchored",
                                        pairedEnd = FALSE)
tests = list(V1_fulllen_plus, V1_fulllen_minus, V2_endanchored_plus,
             V2_endanchored_minus, V3_startanchored_plus, V3_startanchored_minus)
tests2 <- lapply(tests, data.table::melt, measure.vars=c(2:201),
                    variable.name="Bins", value.name="Score")

merged_tests_plus=cbind(tests2[[1]], tests2[[3]][[3]], tests2[[5]][[3]])
colnames(merged_tests_plus)=c("gene_ID", "Bins", "fulllen", "end",
                       "start"
)
merged_tests_plus$Bins = as.numeric(merged_tests_plus$Bins)*10-1000
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

merged_tests_plus[,gene_ID:=NULL]
mean_merged_tests_plus = merged_tests_plus[,lapply(.SD, mean), by=Bins]


tests_melt = data.table::melt(mean_merged_tests_plus, measure.vars=2:4,
                              variable.name="Test", value.name = "Score")

geom = ggplot(data=tests_melt, aes(x=Bins, y=Score, color=Test)) +
  mytheme +
  geom_line(size=1.5) +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed")
geom

svg(file="test_proseq_metaplotting_plus.svg")
geom
dev.off()

merged_tests_minus=cbind(tests2[[2]], tests2[[4]][[3]], tests2[[6]][[3]])
colnames(merged_tests_minus)=c("gene_ID", "Bins", "fulllen", "end",
                              "start"
)
merged_tests_minus$Bins = as.numeric(merged_tests_minus$Bins)*10-1000

merged_tests_minus[,gene_ID:=NULL]
mean_merged_tests_minus = merged_tests_minus[,lapply(.SD, mean), by=Bins]


tests_melt = data.table::melt(mean_merged_tests_minus, measure.vars=2:4,
                              variable.name="Test", value.name = "Score")

geom = ggplot(data=tests_melt, aes(x=Bins, y=Score, color=Test)) +
  mytheme +
  geom_line(size=1.5) +
  xlab("Distance from Anchor") +
  ylab("Tags") +
  geom_vline(xintercept = 0, linetype="dashed")
geom

svg(file="test_proseq_metaplotting_minus.svg")
geom
dev.off()
