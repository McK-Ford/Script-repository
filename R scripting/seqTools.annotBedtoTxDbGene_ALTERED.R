#### setup ####
library( "BSgenome" )
library( "data.table" )
library( "GenomicFeatures" )
library( "rtracklayer" )

annotgrtoTxDbGene <- function( gr, tx, org = c( NULL ) ) {
  ts <- transcripts( tx, columns = c( "TXID", "TXNAME", "GENEID" ) ) #get transcripts
  #Determine the the closest transcript
  closeTs <- distanceToNearest( gr, ts, ignore.strand = TRUE ) #returns a hits object
  txqh <- queryHits( closeTs ) #bc it's cleaner if we only call it once
  txsh <- subjectHits( closeTs )
  ##assign metadata
  gr$TXNAME[ txqh ] <- ts$TXNAME[ txsh ]
  gr$txStart[ txqh ] <- start( ts )[ txsh ]
  gr$txEnd[ txqh ] <- end( ts )[ txsh ]
  gr$txStrand[ txqh ] <- as.character( strand( ts )[ txsh ] )
#skip distance stranding, I have ways I trust more for that.
  #Assign ts dist	
  gr$txDist[ txqh ] <- closeTs@elementMetadata$distance
 
#### TSS portion
  tss <- transcripts( tx, columns = c( "TXID", "TXNAME", "GENEID" ) ); #, ...);	
  #set txs to only tsses
  end( tss )[ as.character( strand( tss ) ) == "+" ] <- start( tss )[ as.character( strand( tss ) ) == "+" ]
  start( tss )[ as.character( strand( tss ) ) == "-" ] <- end( tss )[ as.character( strand( tss ) ) == "-" ]
  #Determine the relative location
  closeTss <- distanceToNearest( gr, tss, ignore.strand = T )
  tssqh <- queryHits( closeTss ) #bc it's cleaner if we only call it once
  tsssh <- subjectHits( closeTss )
  
  #Assign tss name to annotation column (AnnotCol) "ts"
  gr$tss[ tssqh ] <- tss$TXNAME[ tsssh ]
  gr$tssStart[ tssqh ] <- start(ts)[ tsssh ]
  gr$tssEnd[ tssqh ] <- end( ts )[ tsssh ]
  gr$tssStrand[ tssqh ] <- as.character( strand( ts )[ tsssh ] )
  #Assign tss dist	
  gr$tssDist[ tssqh ] <- closeTss@elementMetadata$distance
  
#actual genes ... how does txdb for ucsc knowngene differentiate between gene and transcript?
  #Assign gene overlap
  genes <- genes( tx, columns = c( "GENEID" ) ) #this is where genes may be dropped
  overlaps <- findOverlaps( gr, genes )
  overlaps <- cbind( as.data.frame( overlaps ), id = as.character( genes$GENEID[ subjectHits( overlaps ) ]))
  overlaps <- data.table( unique( overlaps[ , c( 1, 3 ) ] ) )
  overlaps <- overlaps[ !is.na( overlaps$id ), ]
  ag <- overlaps[ , list( id = paste( id, collapse = "," ) ), by = c( 'queryHits' ) ]
  gr$gene[ ag$queryHits ] <- ag$id 
  return( gr );
}
