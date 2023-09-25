obj_name <- function( obj ) { deparse( substitute( obj ) ) }

list_chrs <- function( inc_X = FALSE, inc_Y = FALSE, ...) {
  #does not include non-standard chroms
  chrs <- as.list( paste0( "chr", rep( 1:22 ) ) )
  if ( inc_X ) { chrs <- append( chrs, "chrX" ) }
  if ( inc_Y ) { chrs <- append( chrs, "chrY" ) }
  return( chrs )
}

order_bed_by_chrom <- function( bed, ... ){
  #This function will order a bed or anything where FIRST COLUMN is chromosomes.
  chrs <- list_chrs( ... )
  bedsub <- lapply( chrs, function( chr ) bed[ bed[[1]] == chr, ] )
  bed <- ls_to_mat( bedsub )
  return( bed )
}

readcounts <- function( bam, pairedEnd, debug = FALSE, ... ) {
  #pairedEnd is a bool, bam a BamFile obj
  require( Rsamtools )
  read_counts <- if( pairedEnd ) countBam(
    bam,
    param = ScanBamParam( flag = scanBamFlag( isProperPair = TRUE, isFirstMateRead = TRUE, isSecondaryAlignment = FALSE, isDuplicate = FALSE ) ) 
  )$records else
    countBam( bam,
              param = ScanBamParam( flag = scanBamFlag( isSecondaryAlignment = FALSE, isDuplicate = FALSE ) ) )$records
  if ( debug ) print( paste0( "rnorm ", read_counts ) )
  return( read_counts )
}

ls_to_mat <- function( ls ) { do.call( rbind, ls ) }

randsamp <- function( tab, s = 10 ) { tab[ sample( nrow( tab ), size = s ), ] } #bc I never remember how to do this...

invert_list <- function( LoL ) {
  tmp_ls <- lapply( LoL, '[', names( LoL[[1]] ) )
  apply( do.call( rbind, tmp_ls ), 2, as.list )
}

write_genetsv <- function(tab, name, coln = FALSE){
  write.table(tab, name, quote=FALSE,
  row.names=FALSE, col.names=coln, sep = "\t", eol="\n")
}
