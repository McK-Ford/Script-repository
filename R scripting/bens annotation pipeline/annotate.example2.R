########### Benjamin G. Barwick
## Example to annotate to custom TxDbs (e.g. Ensembl / Gencode)

#set working directory
wd = "/home/bbarwick/Documents/Vertino/annotation/"
setwd(wd)

#load R code for annotation
source("seqTools.annotBedtoTxDbGene.R")
# note this requires "data.table", "BSgenome", "GenomicRanges", "rtracklayer",
# and "GenomicFeatures" packages

#fast read in bed file and convert to GRanges object 
#Note replace 'chr1' with '1'

bedFile <- "bedFiles.diffBind/regionRef.fdr_0001.dedup.rpkm.9.bed"
bed     <- fread(bedFile)
bedgr   = GRanges(seqnames = gsub("chr", "", bed[[1]]), ranges = IRanges(start = bed[[2]], end = bed[[3]]))

#Build txDb from a GTF file 
gtfFile = "Homo_sapiens.GRCh38.104.gtf.gz"
tx = makeTxDbFromGFF(gtfFile, format = "gtf")

#compile flat annotation file from gtf - write out just in case
gtfAnnot = annotateGtf(gtfFile, outFile = gsub(".gtf.gz", ".gtf.flatAnnot.txt", gtfFile)) 

#annotate. note the TxDb and bed chromosomes have to have the same notation as the (e.g. 'chr1' or '1'). also note output is a GRanges object
#Note that TxDb often do not contain gene symbols. 
annot = annotBedtoTxDbGene(bedgr, tx, prefix = "GRCh38")

### Add gene symbol after the fact
annot = as.data.table(annot)
annot$GRCh38.tsSymbol  = gtfAnnot$symbol[match(annot$GRCh38.tsKg, gtfAnnot$gene_id)]
annot$GRCh38.tssSymbol = gtfAnnot$symbol[match(annot$GRCh38.tssKg, gtfAnnot$gene_id)]

#fast write out
fwrite(annot, gsub(".bed", ".annot2.bed", bedFile), sep = "\t", row.names = F)


### Session info
sessionInfo()
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] GenomicFeatures_1.38.2 AnnotationDbi_1.48.0   Biobase_2.46.0        
 [4] data.table_1.12.8      BSgenome_1.54.0        rtracklayer_1.46.0    
 [7] Biostrings_2.54.0      XVector_0.26.0         GenomicRanges_1.38.0  
[10] GenomeInfoDb_1.22.1    IRanges_2.20.2         S4Vectors_0.24.4      
[13] BiocGenerics_0.32.0   

loaded via a namespace (and not attached):
 [1] SummarizedExperiment_1.16.1 progress_1.2.2             
 [3] tidyselect_1.1.0            purrr_0.3.4                
 [5] lattice_0.20-41             vctrs_0.3.1                
 [7] generics_0.0.2              BiocFileCache_1.10.2       
 [9] blob_1.2.1                  XML_3.99-0.3               
[11] rlang_0.4.6                 pillar_1.4.4               
[13] glue_1.4.1                  DBI_1.1.0                  
[15] rappdirs_0.3.1              BiocParallel_1.20.1        
[17] bit64_0.9-7                 dbplyr_1.4.4               
[19] matrixStats_0.56.0          GenomeInfoDbData_1.2.2     
[21] lifecycle_0.2.0             stringr_1.4.0              
[23] zlibbioc_1.32.0             memoise_1.1.0              
[25] biomaRt_2.42.1              curl_4.3                   
[27] Rcpp_1.0.5                  openssl_1.4.2              
[29] DelayedArray_0.12.3         bit_1.1-15.2               
[31] Rsamtools_2.2.3             hms_0.5.3                  
[33] askpass_1.1                 digest_0.6.25              
[35] stringi_1.4.6               dplyr_1.0.0                
[37] grid_3.6.3                  tools_3.6.3                
[39] bitops_1.0-6                magrittr_1.5               
[41] RCurl_1.98-1.2              RSQLite_2.2.0              
[43] tibble_3.0.2                crayon_1.3.4               
[45] pkgconfig_2.0.3             ellipsis_0.3.1             
[47] Matrix_1.2-18               prettyunits_1.1.1          
[49] assertthat_0.2.1            httr_1.4.1                 
[51] R6_2.4.1                    GenomicAlignments_1.22.1   
[53] compiler_3.6.3             
>
