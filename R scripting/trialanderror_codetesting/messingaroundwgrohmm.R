BiocManager::install("groHMM")
BiocManager::install("GenomicFeatures")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("edgeR")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(groHMM)
options(mc.cores=getCores(4))
# S0mR1 <- as(readGAlignments(system.file("extdata", "S0mR1.bam", package="groHMM")), "GRanges")
# S0mR2 <- as(readGAlignments(system.file("extdata", "S0mR1.bam", package="groHMM")), "GRanges") # Use R1 as R2
# S40mR1 <- as(readGAlignments(system.file("extdata", "S40mR1.bam", package="groHMM")), "GRanges")
# S40mR2 <- as(readGAlignments(system.file("extdata", "S40mR1.bam", package="groHMM")), "GRanges") # Use R1 as R2
# 
# # Combine replicates
# S0m <- c(S0mR1, S0mR2)
# S40m <- c(S40mR1, S40mR2)
# 
# writeWiggle(reads=S0m, file="S0m_Plus.wig", fileType="wig", strand="+", reverse=FALSE)
# writeWiggle(reads=S0m, file="S0m_Minus.wig", fileType="wig", strand="-", reverse=TRUE)
# 
# # For BigWig file:
# #library(BSgenome.Hsapiens.UCSC.hg19)
# #si <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
# #writeWiggle(reads=S0m, file="S0m_Plus.wig", fileType="BigWig", strand="+", reverse=FALSE, seqinfo=si)
# 
# # Normalized wiggle files
# expCounts <- mean(c(NROW(S0m), NROW(S40m)))
# writeWiggle(reads=S0m, file="S0m_Plus_Norm.wig", fileType="wig", strand="+", normCounts=expCounts/NROW(S0m), reverse=FALSE)
# 
# Sall <- sort(c(S0m, S40m))
# # hmmResult <- detectTranscripts(Sall, LtProbB=-200, UTS=5, threshold=1)
#  # Load hmmResult from the saved previous run
# load(system.file("extdata", "Vignette.RData", package="groHMM"))
# txHMM <- hmmResult$transcripts
# head(txHMM)
# 
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# kgdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# library(GenomicFeatures)
# # For refseq annotations:
# # rgdb <- makeTxDbFromUCSC(genome="hg19", tablename="refGene")
# # saveDb(hg19RGdb, file="hg19RefGene.sqlite")
# # rgdb <- loadDb("hg19RefGene.sqlite")
# kgChr7 <- transcripts(kgdb, filter=list(tx_chrom = "chr7"), columns=c("gene_id", "tx_id", "tx_name"))
# seqlevels(kgChr7) <- seqlevelsInUse(kgChr7)
# 
# # Collapse overlapping annotations
# kgConsensus <- makeConsensusAnnotations(kgChr7, keytype="gene_id", mc.cores=getOption("mc.cores"))
# library(org.Hs.eg.db)
# map <- select(org.Hs.eg.db, keys=unlist(mcols(kgConsensus)$gene_id), columns=c("SYMBOL"), keytype=c("ENTREZID"))
# mcols(kgConsensus)$symbol <- map$SYMBOL
# mcols(kgConsensus)$type <- "gene"
# 
# e <- evaluateHMMInAnnotations(txHMM, kgConsensus)
# e$eval
# 
# getExpressedAnnotations <- function(features, reads) {
#   fLimit <- limitToXkb(features)
#   count <- countOverlaps(fLimit, reads)
#   features <- features[count!=0,]
#   return(features[(quantile(width(features), .05) < width(features)) & (width(features) < quantile(width(features), .95)),])}
# conExpressed <- getExpressedAnnotations(features=kgConsensus,reads=Sall)
# 
# td <- getTxDensity(txHMM, conExpressed, mc.cores=getOption("mc.cores"))
# #Merged annotations: 51
# #Dissociated a single annotation: 35
# #Overlaps between transcript and annotation:
#  # Total = 494 Used for density = 389
# u <- par("usr")
# lines(c(u[1], 0, 0, 1000, 1000, u[2]), c(0,0,u[4]-.04,u[4]-.04,0,0), col="red")
# legend("topright", lty=c(1,1), col=c(2,1), c("ideal", "groHMM"))
# text(c(-500,500), c(.05,.5), c("FivePrimeFP", "TP"))
# 
# bPlus <- breakTranscriptsOnGenes(txHMM, kgConsensus, strand="+")
# #19 transcripts are broken into 46
# bMinus <- breakTranscriptsOnGenes(txHMM, kgConsensus, strand="-")
# #14 transcripts are broken into 32
# txBroken <- c(bPlus, bMinus)
# txFinal <- combineTranscripts(txBroken, kgConsensus)
# #87 transcripts are combined to 34
# tdFinal <- getTxDensity(txFinal, conExpressed, mc.cores=getOption("mc.cores"))
# 
# #I'm not going to worry about the expression section at the moment
# 
# # okay, pausing index
# features2 <- GRanges("chr7", IRanges(2394074,2420377), strand="+")
# reads <- as(readGAlignments(system.file("extdata", "S0mR1.bam", package="groHMM")), "GRanges")
# pitst2 <- pausingIndex(features2, reads, up=400)
# features2 <- GRanges("chr7", IRanges(2393274,2420477), strand="+")
# ptst3 <- pausingIndex(features2, reads)


#What if we go over stuff thing by thing... Let's do it using a subset of my own stuff to better understand where everything comes from.
intersection_table <- read.delim( #Reading in the CpG/Genes intersection file, from bedtools
  "Refseq_curated_CGI_promoters_filter.hg38.txt",
  header=FALSE,
  quote="")
names(intersection_table) <- c("chrom", "cpg_l", "cpg_h", "cpg_id", "gene_l", "gene_h", "gene_name", "score_filler", "strand")
library(tidyverse)
intersection_table <- intersection_table %>% distinct(gene_name, .keep_all = TRUE)
test_set <- intersection_table[sample(nrow(intersection_table), 20), ]
my_reads <- as(readGAlignments("Kraus.2013.Groseq.MCF7.10_QC.sort.bam"), "GRanges")
my_reads #okay why is my IRanges a real range and theirs isn't? I suppose I only got the single BP after the bam stage?
my_features <- GRanges(seqnames = test_set$chrom, ranges = IRanges(start=test_set$gene_l, end=test_set$gene_h), strand = test_set$strand)
mytst <- pausingIndex(my_features2, my_reads)
#so it does run though we get a NaN and we don't know what matches up with what genes. I ran until I got a no NaN version to make life easier.
#what about when you go line by line through the document? Might have problems at the C part but we should try and see

reads <- my_reads[order(as.character(seqnames(my_reads)), start(my_reads)), ]
f <- data.frame(chrom=as.character(seqnames(my_features)), 
                start=as.integer(start(my_features)), end=as.integer(end(my_features)), 
                strand=as.character(strand(my_features)))
#lol I like that, you make the granges back into a data frame anyway, why require features to be a DF?
if ("symbol" %in% names(mcols(my_features))){
  f <- cbind(f, symbol=my_features$symbol) 
} else {
  f <- cbind(f, symbol=GeneID <- as.character(seq_len(NROW(f))))
}
#this is what adds the 1-whatever list. Couldn't I do 'gene_name' instead of symbol if I didn't want to get numbers? lets try... or just as.data.frame the GRanges, yes I didn't type specify but it should still be fine dammit
my_features2 <- GRanges(seqnames = test_set$chrom, ranges = IRanges(start=test_set$gene_l, end=test_set$gene_h), strand = test_set$strand)
my_features2$symbol <- test_set$gene_name
f2 <- data.frame(chrom=as.character(seqnames(my_features2)), 
                 start=as.integer(start(my_features2)), end=as.integer(end(my_features2)), 
                 strand=as.character(strand(my_features2)), gene_name=my_features2$gene_name)
f=f2
p <- as.data.frame(reads) #changed this too to just an as.data.frame nm not explicitly calling types causes issues downstream
C <- sort(as.character(unique(f[[1]]))) #yep that's the chrom list
#make a bunch of placeholder vectors
Pause <- rep(0,NROW(f))
Body  <- rep(0,NROW(f))
Fish  <- rep(0,NROW(f))
GeneID <- rep("", NROW(f))
CIl  <- rep(0,NROW(f))
CIu  <- rep(0,NROW(f))

## Pass back information for the fisher test...
PauseCounts <- rep(0, NROW(f))
BodyCounts  <- rep(0, NROW(f))
UpCounts    <- rep(0, NROW(f))
UgCounts    <- rep(0, NROW(f))

size = 50
up = 1000
down = 1000

size        <- as.integer(size)
up          <- as.integer(up)
down        <- as.integer(down)

#we'll do it their way even though mine is just as quick -  might have problems with type of table my way after all
PLUS_INDX <- which(f[[4]] == "+")
MINU_INDX <- which(f[[4]] == "-")

debug = TRUE

if(debug) {
  message("Calculating TSS and gene ends for each gene based 
            on strand information.")
}
c_tss_indx <- rep(0,NROW(f))
c_tss_indx[PLUS_INDX] <- 2
c_tss_indx[MINU_INDX] <- 3
c_tss <- unlist(lapply(c(1:NROW(f)), function(x) { f[x, c_tss_indx[x]] })) #okay when going line by line I can understand this I think. It just grabs the tss, respective to strand, so it just grabs the relevant column of f.

###### Now calculate left and right position for gene body, based 
### on '+' or '-'.
### Calculate gene end.  Gene start is contiguous with the coordinates 
###for the promoter.
c_gene_end_indx <- rep(0,NROW(f))
c_gene_end_indx[PLUS_INDX] <- 3
c_gene_end_indx[MINU_INDX] <- 2
c_gene_end <- unlist(lapply(c(1:NROW(f)), function(x) { 
  f[x,c_gene_end_indx[x]] }))
#that was the same thing just flipped


### Assign left and right.
gLEFT   <- rep(0,NROW(c_tss))
gRIGHT  <- rep(0,NROW(c_tss))
gLEFT[PLUS_INDX]    <- c_tss[PLUS_INDX] + down
gRIGHT[PLUS_INDX]   <- c_gene_end[PLUS_INDX]
gLEFT[MINU_INDX]    <- c_gene_end[MINU_INDX]
gRIGHT[MINU_INDX]   <- c_tss[MINU_INDX] - down

## Run parallel version.
#### mcp <- mclapply(c(1:NROW(C)), pausingIndex_foreachChrom, C=C, f=f, p=p, 
     #######           gLEFT=gLEFT, gRIGHT=gRIGHT, c_tss=c_tss, 
      ######          size=size, up=up, down=down, UnMAQ=UnMAQ, debug=debug, ...)
#obv that's iterative - lets try chr12 instead cause it has both directions, to look at closer
indxF   <- which(as.character(f[[1]]) == C[3])
indxPrb <- which(as.character(p[[1]]) == C[3])
#that was selecting for rows that correspond to the chroms of interest
#if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
  # Order -- Make sure, b/c this is one of our main assumptions.  
  # Otherwise violated for DBTSS.
  Ford <- order(f[indxF,2])
  Pord <- order(p[indxPrb,2])

  # Type coersions.
  FeatureStart    <- as.integer(gLEFT[indxF][Ford])
  FeatureEnd  <- as.integer(gRIGHT[indxF][Ford])
  FeatureStr  <- as.character(f[indxF,4][Ford])
  FeatureTSS  <- as.integer(c_tss[indxF][Ford])
  
  PROBEStart  <- as.integer(p[indxPrb,2][Pord])
  PROBEEnd    <- as.integer(p[indxPrb,3][Pord])
  PROBEStr    <- as.character(p[indxPrb,4][Pord])  

  # Set dimensions.
  dim(FeatureStart)   <- c(NROW(FeatureStart), NCOL(FeatureStart))
  dim(FeatureEnd)     <- c(NROW(FeatureEnd),   NCOL(FeatureEnd))
  dim(FeatureTSS)     <- c(NROW(FeatureTSS),   NCOL(FeatureTSS))
  dim(FeatureStr)     <- c(NROW(FeatureStr),   NCOL(FeatureStr))
  dim(PROBEStart)     <- c(NROW(PROBEStart),   NCOL(PROBEStart))
  dim(PROBEEnd)       <- c(NROW(PROBEEnd),     NCOL(PROBEEnd))
  dim(PROBEStr)       <- c(NROW(PROBEStr),     NCOL(PROBEStr)) 

  ## Run the calculations on the gene.
  ## Calculate the maximal 50 bp window.
  if(debug) {
    message(C[3],": Counting reads in pause peak.")
  }
  HPause <- .Call("NumberOfReadsInMaximalSlidingWindow", 
                  FeatureTSS, FeatureStr, PROBEStart, PROBEEnd, 
                  PROBEStr, size, up, down, PACKAGE = "groHMM") #all these values exist, why are they garbage?
  #just returns an int string from the aforementioned C function ... Which is screwing up for some damn reason... Dammit I don't want to have to eff with the C function too...
  ## Run the calculate on the gene body...
  if(debug) {
    message(C[3],": Counting reads in gene.")
  }
  HGeneBody <- .Call("CountReadsInFeatures", FeatureStart, 
                     FeatureEnd, FeatureStr, PROBEStart, PROBEEnd, PROBEStr, 
                     PACKAGE = "groHMM")
  
  ## Get size of gene body
  Difference <- FeatureEnd-FeatureStart
  Difference[Difference < 0] <- 0 
  ## Genes < 1kb, there is no suitable area in the body of the gene.
  
  #Not concerned at the moment about mapQ, gonna skip that
  ## Now use Fisher's Exact.
  if(debug) {
    message(C[3],": Using Fisher's exact.")
  }
  # Make uniform reads.
  Up <- round(HPause + HGeneBody)*(size)/(size+Difference) 
  ## Expted in pause == 
  ##((total reads)/ (total size) [reads/ base]) * 
  ## size [reads/ pause window]
  Ug <- round(HPause + HGeneBody)*(Difference)/(size+Difference) 
  ## Expted reads in body == 
  ## ((total reads)/ (total size) [reads/ base]) * 
  ## (gene size) [reads/ gene body]
  HFish <- unlist(lapply(c(1:NROW(Ford)), function(x) {
    fisher.test(
      data.frame(
        c(HPause[x], round(Up[x])), 
        c(HGeneBody[x], round(Ug[x]))
      )
    )$p.value
  } ))
  ## Make return values.
  Pause_c     <- as.double(HPause/size)
  Body_c  <- as.double((HGeneBody+1)/Difference) 
  ## 6-5-2012 -- Add a pseudocount here, 
  ## forcing at least 1 read in the gene body.   
  Fish_c  <- as.double(HFish)
  GeneID_c    <- as.character(f[indxF,5][Ford])
  
  PauseCounts_c <- HPause
  BodyCounts_c  <- HGeneBody
  UpCounts_c    <- round(Up)
  UgCounts_c    <- round(Ug)
  
  approx.ratio.CI <- function(x1, x2, alpha=0.05) {
    t = qnorm(1 - alpha/2)
    n = x1 + x2
    zp = (t^2/2 + x1 + 1/2)^2 - ((n + t^2)/n) * (x1 + 1/2)^2
    zn = (t^2/2 + x1 - 1/2)^2 - ((n + t^2)/n) * (x1 - 1/2)^2
    
    a = (t^2/2 + x1 + 1/2 + sqrt(zp)) / (n + t^2/2 - x1 - 1/2 - sqrt(zp))
    b = (t^2/2 + x1 - 1/2 - sqrt(zn)) / (n + t^2/2 - x1 + 1/2 + sqrt(zn))
    
    return(c(b, a))
  }
  
  approx.ratios.CI <- function(num.counts, denom.counts, alpha=0.05) {
    stopifnot(length(num.counts) == length(denom.counts))
    N = length(num.counts)
    
    result = matrix(data=0, nrow=N, ncol=2)
    
    for (i in 1:N)
      result[i, ] = approx.ratio.CI(num.counts[i], denom.counts[i], alpha)
    
    return(result)
  }
  
  aCI <- approx.ratios.CI(HPause, HGeneBody)
  scaleFactor <- Difference/size ## Body/ pause, must be 1/ PI units.
  CIl_c <- as.double(aCI[,1]*scaleFactor)
  CIu_c <- as.double(aCI[,2]*scaleFactor)
  
  foo = list(Pause= Pause_c, Body= Body_c, Fish= Fish_c, 
              GeneID= GeneID_c, PauseCounts= PauseCounts_c, 
              BodyCounts= BodyCounts_c, UpCounts= UpCounts_c, 
              UgCounts= UgCounts_c, CIl= CIl_c, CIu= CIu_c, ord= Ford)
#then there's something weird going on with the list stuff
  #for the sake of my sanity I'm just gonna as.data.frame it
 bar = as.data.frame(foo)
 #unfortunately something a bit odd happens. The pause function didn't retrieve any values dammit while 'mytest' did so it should have. It also didn't return body counts?
 
 mytst2 <- pausingIndex2(my_features2, my_reads)
 pausingIndex2 <- function(features, reads, size=50, up=1000, down=1000, 
                          UnMAQ=NULL, debug=FALSE, ...) {
   # make sure reads are sorted
   reads <- reads[order(as.character(seqnames(reads)), start(reads)), ]
   f <- data.frame(chrom=as.character(seqnames(features)), 
                   start=as.integer(start(features)), end=as.integer(end(features)), 
                   strand=as.character(strand(features)))
   if ("symbol" %in% names(mcols(features))){
     f <- cbind(f, symbol=features$symbol) 
   } else {
     f <- cbind(f, symbol=GeneID <- as.character(seq_len(NROW(f))))
   }
   p <- data.frame(chrom=as.character(seqnames(reads)), 
                   start=as.integer(start(reads)), end=as.integer(end(reads)), 
                   strand=as.character(strand(reads)))
   
   C <- sort(as.character(unique(f[[1]])))
   Pause <- rep(0,NROW(f))
   Body  <- rep(0,NROW(f))
   Fish  <- rep(0,NROW(f))
   GeneID <- rep("", NROW(f))
   CIl  <- rep(0,NROW(f))
   CIu  <- rep(0,NROW(f))
   
   ## Pass back information for the fisher test...
   PauseCounts <- rep(0, NROW(f))
   BodyCounts  <- rep(0, NROW(f))
   UpCounts    <- rep(0, NROW(f))
   UgCounts    <- rep(0, NROW(f))
   
   size        <- as.integer(size)
   up          <- as.integer(up)
   down        <- as.integer(down)
   
   ###### Calculate PLUS and MINUS index, for DRY compliance.
   PLUS_INDX <- which(f[[4]] == "+")
   MINU_INDX <- which(f[[4]] == "-")
   
   ###### Identify TSS -- Start for '+' strand, End for '-' strand.
   if(debug) {
     message("Calculating TSS and gene ends for each gene based 
            on strand information.")
   }
   c_tss_indx <- rep(0,NROW(f))
   c_tss_indx[PLUS_INDX] <- 2
   c_tss_indx[MINU_INDX] <- 3
   c_tss <- unlist(lapply(c(1:NROW(f)), function(x) { f[x, c_tss_indx[x]] }))
   
   ###### Now calculate left and right position for gene body, based 
   ### on '+' or '-'.
   ### Calculate gene end.  Gene start is contiguous with the coordinates 
   ###for the promoter.
   c_gene_end_indx <- rep(0,NROW(f))
   c_gene_end_indx[PLUS_INDX] <- 3
   c_gene_end_indx[MINU_INDX] <- 2
   c_gene_end <- unlist(lapply(c(1:NROW(f)), function(x) { 
     f[x,c_gene_end_indx[x]] }))
   
   ### Assign left and right.
   gLEFT   <- rep(0,NROW(c_tss))
   gRIGHT  <- rep(0,NROW(c_tss))
   gLEFT[PLUS_INDX]    <- c_tss[PLUS_INDX] + down
   gRIGHT[PLUS_INDX]   <- c_gene_end[PLUS_INDX]
   gLEFT[MINU_INDX]    <- c_gene_end[MINU_INDX]
   gRIGHT[MINU_INDX]   <- c_tss[MINU_INDX] - down
   
   ## Run parallel version.
   mcp <- mclapply(c(1:NROW(C)), pausingIndex_foreachChrom, C=C, f=f, p=p, 
                   gLEFT=gLEFT, gRIGHT=gRIGHT, c_tss=c_tss, 
                   size=size, up=up, down=down, UnMAQ=UnMAQ, debug=debug, ...)
   
   ## Unlist and re-order values for printing in a nice data.frame.
   for(i in 1:NROW(C)) {
     # Which KG?  prb?
     indxF   <- which(as.character(f[[1]]) == C[i])
     indxPrb <- which(as.character(p[[1]]) == C[i])
     
     if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
       Pause[indxF][mcp[[i]][["ord"]]]     <- mcp[[i]][["Pause"]]
       Body[indxF][mcp[[i]][["ord"]]]  <- mcp[[i]][["Body"]]
       Fish[indxF][mcp[[i]][["ord"]]]  <- mcp[[i]][["Fish"]]
       GeneID[indxF][mcp[[i]][["ord"]]]    <- mcp[[i]][["GeneID"]]
       
       PauseCounts[indxF][mcp[[i]][["ord"]]] <- mcp[[i]][["PauseCounts"]]
       BodyCounts[indxF][mcp[[i]][["ord"]]]  <- mcp[[i]][["BodyCounts"]]
       UpCounts[indxF][mcp[[i]][["ord"]]]    <- mcp[[i]][["UpCounts"]]
       UgCounts[indxF][mcp[[i]][["ord"]]]    <- mcp[[i]][["UgCounts"]]
       
       CIl[indxF][mcp[[i]][["ord"]]] <- mcp[[i]][["CIl"]]
       CIu[indxF][mcp[[i]][["ord"]]] <- mcp[[i]][["CIu"]]
     }
   }
   
   return(data.frame(Pause= Pause, Body= Body, Fisher= Fish, GeneID= GeneID, 
                     CIlower=CIl, CIupper=CIu, PauseCounts= PauseCounts, 
                     BodyCounts= BodyCounts, uPCounts= UpCounts, uGCounts= UgCounts))
 }
 
 pausingIndex_foreachChrom <- function(i, C, f, p, gLEFT, gRIGHT, c_tss, size, 
                                       up, down, UnMAQ, debug) {
   if(debug) {
     message(C[i])
   }
   
   # Which KG?  prb?
   indxF   <- which(as.character(f[[1]]) == C[i])
   indxPrb <- which(as.character(p[[1]]) == C[i])
   
   if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
     # Order -- Make sure, b/c this is one of our main assumptions.  
     # Otherwise violated for DBTSS.
     Ford <- order(f[indxF,2])
     Pord <- order(p[indxPrb,2])
     
     # Type coersions.
     FeatureStart    <- as.integer(gLEFT[indxF][Ford])
     FeatureEnd  <- as.integer(gRIGHT[indxF][Ford])
     FeatureStr  <- as.character(f[indxF,4][Ford])
     FeatureTSS  <- as.integer(c_tss[indxF][Ford])
     
     PROBEStart  <- as.integer(p[indxPrb,2][Pord])
     PROBEEnd    <- as.integer(p[indxPrb,3][Pord])
     PROBEStr    <- as.character(p[indxPrb,4][Pord])
     
     # Set dimensions.
     dim(FeatureStart)   <- c(NROW(FeatureStart), NCOL(FeatureStart))
     dim(FeatureEnd)     <- c(NROW(FeatureEnd),   NCOL(FeatureEnd))
     dim(FeatureTSS)     <- c(NROW(FeatureTSS),   NCOL(FeatureTSS))
     dim(FeatureStr)     <- c(NROW(FeatureStr),   NCOL(FeatureStr))
     dim(PROBEStart)     <- c(NROW(PROBEStart),   NCOL(PROBEStart))
     dim(PROBEEnd)       <- c(NROW(PROBEEnd),     NCOL(PROBEEnd))
     dim(PROBEStr)       <- c(NROW(PROBEStr),     NCOL(PROBEStr))
     
     ## Run the calculations on the gene.
     ## Calculate the maximal 50 bp window.
     if(debug) {
       message(C[i],": Counting reads in pause peak.")
     }
     HPause <- .Call("NumberOfReadsInMaximalSlidingWindow", 
                     FeatureTSS, FeatureStr, PROBEStart, PROBEEnd, 
                     PROBEStr, size, up, down, PACKAGE = "groHMM")
     
     print(HPause) #It's not FeatureTSS or FeatureStr that's the problem...
     ## Run the calculate on the gene body...
     if(debug) {
       message(C[i],": Counting reads in gene.")
     }
     HGeneBody <- .Call("CountReadsInFeatures", FeatureStart, 
                        FeatureEnd, FeatureStr, PROBEStart, PROBEEnd, PROBEStr, 
                        PACKAGE = "groHMM")
     ## Get size of gene body
     Difference <- FeatureEnd-FeatureStart
     Difference[Difference < 0] <- 0 
     ## Genes < 1kb, there is no suitable area in the body of the gene.
     
     ## Calculate UN-MAQable regions...
     if(!is.null(UnMAQ)) {
       
       ## Count start index.
       chr_indx <- which(UnMAQ[[1]][[1]] == C[i])
       CHRSIZE <- as.integer(UnMAQ[[1]][[2]][chr_indx])
       CHRSTART <- as.integer(0)
       if(chr_indx > 1) {  ## Running on 1:0 gives c(1, 0)
         CHRSTART <- as.integer( 
           sum(UnMAQ[[1]][[2]][
             c(1:(chr_indx-1))
           ]) +1)
       }
       
       if(debug) {
         message(C[i],": Counting unMAQable regions.")
         message("CHRSIZE:", CHRSIZE, "CHRSTART:", CHRSTART)
       }
       
       ## Count unMAQable regions, and size of everything ... 
       nonmappable <- .Call("CountUnMAQableReads", 
                            FeatureStart, FeatureEnd, as.integer(UnMAQ[[2]]), 
                            CHRSTART, CHRSIZE, PACKAGE = "groHMM")
       
       ## Adjust size of gene body.
       Difference <- Difference - nonmappable + 1 
       ## Otherwise, get -1 for some.
       
       if(debug) {
         print(head(nonmappable))
         print(as.integer(head(Difference)))
       }
     }
     
     ## Now use Fisher's Exact.
     if(debug) {
       message(C[i],": Using Fisher's exact.")
     }
     # Make uniform reads.
     Up <- round(HPause + HGeneBody)*(size)/(size+Difference) 
     ## Expted in pause == 
     ##((total reads)/ (total size) [reads/ base]) * 
     ## size [reads/ pause window]
     Ug <- round(HPause + HGeneBody)*(Difference)/(size+Difference) 
     ## Expted reads in body == 
     ## ((total reads)/ (total size) [reads/ base]) * 
     ## (gene size) [reads/ gene body]
     HFish <- unlist(lapply(c(1:NROW(Ford)), function(x) {
       fisher.test(
         data.frame(
           c(HPause[x], round(Up[x])), 
           c(HGeneBody[x], round(Ug[x]))
         )
       )$p.value
     } ))
     
     ## Make return values.
     Pause_c     <- as.double(HPause/size)
     Body_c  <- as.double((HGeneBody+1)/Difference) 
     ## 6-5-2012 -- Add a pseudocount here, 
     ## forcing at least 1 read in the gene body.   
     Fish_c  <- as.double(HFish)
     GeneID_c    <- as.character(f[indxF,5][Ford])
     
     PauseCounts_c <- HPause
     BodyCounts_c  <- HGeneBody
     UpCounts_c    <- round(Up)
     UgCounts_c    <- round(Ug)
     
     aCI <- approx.ratios.CI(HPause, HGeneBody)
     scaleFactor <- Difference/size ## Body/ pause, must be 1/ PI units.
     CIl_c <- as.double(aCI[,1]*scaleFactor)
     CIu_c <- as.double(aCI[,2]*scaleFactor)
     
     return(list(Pause= Pause_c, Body= Body_c, Fish= Fish_c, 
                 GeneID= GeneID_c, PauseCounts= PauseCounts_c, 
                 BodyCounts= BodyCounts_c, UpCounts= UpCounts_c, 
                 UgCounts= UgCounts_c, CIl= CIl_c, CIu= CIu_c, ord= Ford))
   }
   return(integer(0))
 }
 