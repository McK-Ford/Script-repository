library(tidyverse)

KDM5B_peaks_stats <- read.delim("~/KDM5B_peaks_stats.txt")
KDM5B_peaks_stats_groupcheck <- KDM5B_peaks_stats %>% group_by(start) %>% add_count() %>% mutate(n>1)
KDM5B_bed_for_annotation <- distinct(KDM5B_peaks_stats[1:3])

source("seqTools.annotBedtoTxDbGene.R")
bedgr = GRanges(seqnames = KDM5B_bed_for_annotation[[1]], ranges = IRanges(start = KDM5B_bed_for_annotation[[2]], end = KDM5B_bed_for_annotation[[3]]))

library("TxDb.Hsapiens.UCSC.hg38.knownGene")
tx = TxDb.Hsapiens.UCSC.hg38.knownGene

library("org.Hs.eg.db")
hs = org.Hs.eg.db

annot = annotBedtoTxDbGene(bedgr, tx, prefix = "hg38", org = hs)
#'select()' returned many:1 mapping between keys and columns
#'select()' returned many:1 mapping between keys and columns
#1655 genes were dropped because they have exons located on both strands of the same reference sequence or
#on more than one reference sequence, so cannot be represented by a single genomic range.
#Use 'single.strand.genes.only=FALSE' to get all the genes in a GRangesList object, or use
#suppressMessages() to suppress this message.
#'select()' returned many:1 mapping between keys and columns
#'select()' returned many:1 mapping between keys and columns
#Warning messages:
#  1: In .set_group_names(ans, use.names, x, "tx") :
#  some group names are NAs or duplicated
#2: In .set_group_names(ans, use.names, txdb, "tx") :
#  some group names are NAs or duplicated
#3: In .set_group_names(ans, use.names, txdb, "tx") :
#  some group names are NAs or duplicated

#his fwrite isn't working for me :(
tst = as.data.frame(annot)

gtfFile = "Homo_sapiens.GRCh38.gtf.gz"
tx = makeTxDbFromGFF(gtfFile, format = "gtf")
#Import genomic features from the file as a GRanges object ... OK
#Prepare the 'metadata' data frame ... OK
#Make the TxDb object ... OK
#Warning message:
#  In .get_cds_IDX(mcols0$type, mcols0$phase) :
#  The "phase" metadata column contains non-NA values for features of type stop_codon. This information was
#ignored.
gtfAnnot = annotateGtf(gtfFile, outFile = gsub(".gtf.gz", ".gtf.flatAnnot.txt", gtfFile))
annot = annotBedtoTxDbGene(bedgr, tx, prefix = "GRCh38")
#'select()' returned many:1 mapping between keys and columns
#'select()' returned many:1 mapping between keys and columns
annot = as.data.table(annot)
annot$GRCh38.tsSymbol  = gtfAnnot$symbol[match(annot$GRCh38.tsKg, gtfAnnot$gene_id)]
annot$GRCh38.tssSymbol = gtfAnnot$symbol[match(annot$GRCh38.tssKg, gtfAnnot$gene_id)]


#so I've got tst and annot, and should join them to my main file, but do they have gene vs nongene in their annotation? Localization? looks like not, so I'm probably better off leaving it off for now.