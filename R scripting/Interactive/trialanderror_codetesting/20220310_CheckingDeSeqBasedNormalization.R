Hg38_refseq_genes.txt <- read.delim("~/Hg38_refseq_genes.txt.gz", header=FALSE, comment.char="#")
ZS1_19_MCF7_EV_H2AZ.countmatrix <- read.delim("~/ZS1_19_MCF7_EV_H2AZ.countmatrix.bedGraph", header=FALSE)
ZS1_19_MCF7_SH1_H2AZ.countmatrix <- read.delim("~/ZS1_19_MCF7_SH1_H2AZ.countmatrix.bedGraph", header=FALSE)
ZS1_19_MCF7_SH2_H2AZ.countmatrix <- read.delim("~/ZS1_19_MCF7_SH2_H2AZ.countmatrix.bedGraph", header=FALSE)

reg <- paste0(ZS1_19_MCF7_EV_H2AZ.countmatrix$V1,
              ZS1_19_MCF7_EV_H2AZ.countmatrix$V2,
              ZS1_19_MCF7_EV_H2AZ.countmatrix$V3)
raw.counts <- cbind(ZS1_19_MCF7_EV_H2AZ.countmatrix$V4,
                    ZS1_19_MCF7_SH1_H2AZ.countmatrix$V4,
                    ZS1_19_MCF7_SH2_H2AZ.countmatrix$V4)

#wtf why do they have dif row numbers... that makes no sense.
library(DESeq2)
library(edgeR)
BiocManager::install("edgeR")
NormFactor <- calcNormFactors(object=raw.counts, method="RLE")
LibSize <- colSums(raw.counts)
sizeFactors <- NormFactor*LibSize/1000000
sizeFactors.reciprocal <- 1/sizeFactors
sizeFactors
#[1] 28.78717 26.97790 20.30271
#that would be the number it would divide by.
#that's really close to the numbers it actually used most likely...

#by edgeR
NormFactor <- calcNormFactors(object=raw.counts, method="TMM")
LibSize <- colSums(raw.counts)
sizeFactors <- NormFactor*LibSize/1000000
sizeFactors
#[1] 28.48768 27.03486 20.47292

RPKM_factors <- LibSize/1000000*0.01
RPKM_factors
#[1] 0.2796612 0.2730083 0.2065157

#Yeah. Once you correct for bin size they're very similar.
#EV is a little higher, SH1 and SH2 are a little lower.