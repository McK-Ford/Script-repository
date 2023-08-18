########### Benjamin G. Barwick
## Example to annotate to UCSC TxDb w/ ENTREZ defined genes 
## See Example 2 for annotating to other TxDbs (e.g. Ensembl / Gencode)

#set working directory
wd = "/home/bbarwick/Documents/Vertino/annotation/"
setwd(wd)

#load R code for annotation
source("seqTools.annotBedtoTxDbGene.R"); #note this requires "data.table", "BSgenome", "GenomicRanges", "rtracklayer", and "GenomicFeatures" packages

#fast read in bed file and convert to GRanges object 
bedFile = "bedFiles.diffBind/regionRef.fdr_0001.dedup.rpkm.9.bed"
bed = fread(bedFile)
bedgr = GRanges(seqnames = bed[[1]], ranges = IRanges(start = bed[[2]], end = bed[[3]]))

#Use a prebuilt txDb object for annotation
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
tx = TxDb.Hsapiens.UCSC.hg38.knownGene

#load organism to capture gene symbols (most TxDbs only contain transcript IDs, but not gene symbols). 
#Note this functionality is only supported for UCSC TxDbs that use ENTREZID to reference genes. For Ensembl / Gencode annotations please see example 2 
library("org.Hs.eg.db")
hs = org.Hs.eg.db;

#annotate. note the TxDb and bed chromosomes have to have the same notation as the (e.g. 'chr1' or '1').
annot = annotBedtoTxDbGene(bedgr, tx, prefix = "hg38", org = hs)

#fast write out
fwrite(as.data.table(annot), gsub(".bed", ".annot1.bed", bedFile), sep = "\t", row.names = F)

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
 [1] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
 [2] GenomicFeatures_1.38.2                 
 [3] AnnotationDbi_1.48.0                   
 [4] Biobase_2.46.0                         
 [5] data.table_1.12.8                      
 [6] BSgenome_1.54.0                        
 [7] rtracklayer_1.46.0                     
 [8] Biostrings_2.54.0                      
 [9] XVector_0.26.0                         
[10] GenomicRanges_1.38.0                   
[11] GenomeInfoDb_1.22.1                    
[12] IRanges_2.20.2                         
[13] S4Vectors_0.24.4                       
[14] BiocGenerics_0.32.0                    

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

