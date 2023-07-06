source("C:/Users/kayle/Box/Vertinolab/McKayla Ford/Scripts/R scripts/metaplotter.R") #note this requires "data.table", "BSgenome", "compiler", "rtracklayer", "tidyverse", and "Rsamtools" packages

###executing metaplotter
enableJIT(3) #this compiles functions as we go along, which makes them run faster.
newtab <- read.delim("~/newtab.txt") #chromosome column must be named chrom, strand column must be named strand.
example_notscaled_classed = get_anchored_metaplot_tabs("fwd.bw", "rev.bw", 500, 500, class_tab=newtab, class_col=newtab$class)
H2AZall_genenotscaled_plus = get_anchored_scores(bed=Hg38_refseq_genes, bw="ZS1_19_MCF7_EV_H2AZ.dedup.rpkm.bw", anch="txStart", b=500, a=1000, strand="+", bs=50)
example_notscaled_notclassed = get_anchored_metaplot_tabs("fwd.bw", "rev.bw", 500, 500)
B7_scaled = get_scaled_cgi_metaplot_tabs("MCF7_B7_fwd_10.bw", "MCF7_B7_rev_10.bw")