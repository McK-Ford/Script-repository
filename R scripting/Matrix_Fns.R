#########################
## 1. Load in packages ##
#########################
library( Rsamtools ) #def need
library( data.table )
library( GenomicAlignments ) # need this for count / find overlaps
library( rtracklayer ) #need this for bws
source( "~/Script repository/R scripting/Utility_Fns.R" )
#####
# Functions
####

non_scaled_mat <- function( bed_sub, a, b, method, bs = 10, ... ) {
  if ( missing( a ) | missing( b ) ){
    stop( "a and/or b are missing and must be specified for this method" )
  }
  if ( method == "peak_centered" ) {
    ref_point <- round( ( bed_sub[[3]] - bed_sub[[2]] ) / 2 ) + bed_sub[[2]] #region center
  } else {
    ref_point <- ifelse( bed_sub[[6]] == "+", bed_sub[[2]], bed_sub[[3]] )
  }
  hist_rows <- bed_sub[[4]]
  bins <- ifelse(
    bed_sub[[6]] == "+",
    list( seq( from = -b, to = a - bs, by = bs ) ),
    list( seq( from = b - bs, to = -a, by = -bs ) )
  )
  hist <- round( ls_to_mat( bins ) )
  bin_names <- seq( from = -b, to = a - bs, by = bs )
  hist <- hist + ref_point
  dimnames( hist ) <- list( hist_rows, bin_names )
  hist[ hist <= 0 ] <- NA
  strandvec <- rep( bed_sub[[6]], times = ( a + b ) / bs )
  return( list( hist, strandvec ) )
}

bi_anch_mat <- function( bed_sub, n, ... ){
  if ( missing( n ) ){
    stop( "n is missing and must be specified for this method" )
  }
  ref_point1 <- bed_sub[[2]]
  ref_point2 <- bed_sub[[3]]
  empty_list <- list()
  ends_list <- list()
  for ( j in seq_along( ref_point1 ) ) {
    bin_vec <- seq(
      from = ref_point1[[j]], to = ref_point2[[j]], length.out = n + 1 
      )
    ends_list[[j]] <- bin_vec[2] - bin_vec[1]
    empty_list[[j]] <- bin_vec[ 1:n ] #drop the extra value
  }
  hist <- round( ls_to_mat( empty_list ) )
  hist_rows <- bed_sub[[4]]
  bin_names <- seq( from = 0,to = n - 1, by = 1 )
  dimnames( hist ) <- list( hist_rows, bin_names )
  ends_vec <- unlist( ends_list )
  hist[ hist <= 0 ] <- NA
  ends_mat <- round( hist + ends_vec ) #this is our bin ends
  strandvec <- rep( bed_sub[[6]], times = n )
  mat_list <- list( hist, ends_mat, strandvec, ends_vec )
  return( mat_list )
}

sbp_from_bam <- function( bam, debug = FALSE, sbp = FALSE, ... ) {
  if (sbp) {
  if ( debug ) print( "getting sbp" )
  bam_tmp <- resize( bam, width = 1, fix = "end" )
  return( bam_tmp )
  } else {
    return( bam )
  }
}

bam_revcomp <- function( bam, debug = FALSE, revcomp = FALSE, ... ) {
  if ( revcomp ) {
  if ( debug ) print( "getting revcomp" )
  temp_plus <- as.character( strand( bam ) ) == "+"
  strand( bam ) <- "+"
  strand( bam )[ temp_plus ] <- "-"
  return( bam ) 
  } else {
    return( bam )
  }
}

bam_vars <- function( bam = bam, pairedEnd = FALSE, rnorm = TRUE, debug = FALSE, ... ) {
  bam_info <- if( pairedEnd ) BamFile( bam, asMates = TRUE ) else BamFile( bam )
  si <- seqinfo( bam_info )
  readcounts <- if ( rnorm ) readcounts( bam = bam_info, pairedEnd, debug ) else 1e6 #as normalization is RPM, this is a scaling factor of 1.
  return( list("bam_info" = bam_info, "si" = si, "readcounts" = readcounts) )
}

bw_olaps <- function( gr, bw_sub, ignorestrand = FALSE, debug = FALSE, len_adj = TRUE, ... ){
    olap <-  findOverlaps( gr, bw_sub, ignore.strand = ignorestrand )
    bw_df <-  data.frame( bw_sub[ subjectHits( olap ) ] )
    regions_df <-  data.frame( gr[ queryHits( olap ) ] )
    scored_regions <-  cbind( bw_df, regions_df )
    names( scored_regions ) <- c(
      "chrom", "bs", "be", "bin_w", "star", "score",
      names( scored_regions[ , 7:13 ] )
    )
    scored_regions$bs <- scored_regions$bs - 1 #0 index
    gs_before_bs <-  pmax( ( scored_regions$start - scored_regions$bs ), 0 )
    ge_after_be <-  pmax( ( scored_regions$be - scored_regions$end ), 0 )
    scored_regions$adjustor <- ( 
      scored_regions$bin_w - gs_before_bs - ge_after_be ) / scored_regions$bin_w
    scored_regions$adjusted_score <- scored_regions$score * scored_regions$adjustor
    setDT( scored_regions ) 
    tmp2 <- scored_regions[ ,
                            list( ( sum( adjusted_score ) / sum( adjustor ) ) ),
                            by = .( ID, BI )
                            ]
    hist <- dcast( tmp2, ID~BI, value.var = "V1" ) 
    if (len_adj) {hist <- length_adjuster(hist, ...)}
}

bam_olaps <- function( gr, bam_sub, long_hist, ignorestrand = FALSE, debug = FALSE, len_adj = TRUE, ... ){
  olap <- as.numeric(countOverlaps( query = gr, subject = bam_sub, ignore.strand = ignorestrand))
  if ( debug ) print( head( olap ) )
  if ( debug ) print( tail( olap ) )
  long_hist$value <- olap
  hist <- dcast( long_hist, rn~variable ) #get the hist back into wideform.
  colnames( hist ) <- c( "ID", colnames( hist[, 2:ncol( hist ) ] ) )
  if ( debug ) print( head( hist ) )
  if ( debug ) print( tail( hist ) )
  if (len_adj) {hist <- length_adjuster(hist, ...)}
  return( hist )
}

length_adjuster <- function( tab, bs = 10, ... ){ 
  tab[ , 2:ncol( tab ) ] <- tab[ , 2:ncol( tab ) ] / bs 
  return(tab)
  }

flipped_hist <- function( scoremat, method, ... ) {
  if (method == "bi_anch" & "-" %in% scoremat$strand ) { #flip the - strand scaled stuff
    #I can't do this quite the same way I used to - used a matrix previously and this is a data.table.
    # So I needed to find a new way. This is pretty similar tho.
    plus_tab <- scoremat[ scoremat$strand == "+" ]
    minus_tab <- scoremat[ scoremat$strand == "-" ]
    binvar <- rev( seq( from = 7, to = ncol( minus_tab ), by = 1 ) )
    minus_tab[ , 7:ncol( minus_tab ) ] <- minus_tab[ , binvar, with = FALSE ]
    scoremat <- rbind( plus_tab, minus_tab)
    return( scoremat )
  } else {
    return( scoremat )
  }
}

#Args:
#bed = table with first 6 col in bed format (regions) 
#readsource = string indicating bam/bw filepath (reads)
#method = "peak_centered", "single_anch", "bi_anch"
#Single_anch returns score in bins size bs around anchor point from b bp upstream to a bp downstream. 
#Plus strand anchor is bed col 2. Col 3 is minus strand anchor.
#peak_centered returns score in bins size bs around center of region from b bp upstream
#to a bp downstream.
#bi_anch generates score scaled from start to end of bed region via splitting 
#said region into n bins.
#
#If method is peak_centered or single_anch, need these arguments:
##b - number of basepairs before anchor point
##a - number of basepairs after anchor point
##bs - Optional, desired bin size, defaults to 10 bp.
#
#If method is bi_anch, need this argument:
##n - number of bins to split region evenly into.
#
#bam_or_bw: options are 'bam' and 'bw'. Pretty obvious function.
#
# Optional arguments:
## revcomp - bool TRUE/FALSE. Determines if you should flip bam. Enter TRUE for pro/gro seq. Default FALSE.
## debug - bool TRUE/FALSE. Only looks at Chr 19 and 22, to run faster and prints out extra troubleshooting info. Default FALSE.
## rnorm - bool TRUE/FALSE. RPM normalizes signal. Only works for bam files. Default FALSE.
## sbp - bool TRUE/FALSE. Defaults to FALSE, use TRUE for proseq to only get the bp at the end. Only works for bam files.
## readsOnly - bool TRUE/FALSE. Only relevant for pairedEnd = TRUE with bam files. Determines whether or not insert is included. Defaults to FALSE.
##  ##pairedEnd indicating whether data is paired end.
##ignorestrand
score_matrix <- function( bed, readsource, bam_or_bw, method, debug = FALSE, ... ){
  if ( debug ) {
    print( paste0( "Reads file is located at " , readsource ) )
    print( paste0( "Filetype is", bam_or_bw ))
    print( paste0( "Method is ", method ) )
  }
  if ( bam_or_bw == "bam" )  { #if bam, get readcounts for normalization
    bv <- bam_vars( bam = readsource, ... )
    readcounts <- 1e6 / bv[[3]]
  } else {
      bv = NA
      readcounts <- 1  #this way I don't need an if statement below, will just be unchanged if rnorm false. May play with, not 100% happy with this solution.
    }
  #Iterating over chromosomes
  chrs <- list_chrs( ... )
  if ( debug ){
    chrs <- list( "chr19", "chr22" ) #two small chrs - some errors only show up w/ more than one
    print("Running in debug mode" )
  }
  tmp_LoL <- lapply(
    chrs,
    single_chrom_mat,
    debug = debug,
    bed = bed,
    readsource = readsource,
    method = method,
    bam_or_bw = bam_or_bw,
    bv = bv,
    ...
    )
  names( tmp_LoL ) <- chrs
  tmp_LoL_flip <- invert_list( tmp_LoL )
  hist <- as.data.table( ls_to_mat( tmp_LoL_flip[[2]] ) ) #collapse into one matrix
  hist[ , 2:ncol( hist ) ] <- hist[ , 2:ncol( hist ) ] * readcounts
  bed <- as.data.table( ls_to_mat( tmp_LoL_flip[[1]] ) )
  scoremat <- merge( bed, hist, by.x = "uniq_id", by.y = "ID", all = TRUE ) #sorted slightly differently so can't just do cbind
  scoremat <- flipped_hist( scoremat, method, ... )
  return( scoremat )
}
# Common error message if you see the error "aggregation function missing, defaulting to length" this
#means there are duplicates in the name column of the bed, and they need to be unique.

single_chrom_mat <- function( 
    chr, bed, method, readsource, bam_or_bw, pairedEnd = FALSE, readsOnly = FALSE,
    ignorestrand = FALSE, debug = FALSE, bs, bv, ... ) {
    if ( debug ) { 
      print( paste0( "Analyzing chrom ", chr, " at ", Sys.time() ) )
    }
    bed_sub <- bed[ bed[[1]] == chr, ]
    if ( method == "peak_centered" | method == "single_anch") {
      hist_ls <- non_scaled_mat( bed_sub = bed_sub, method = method, bs = bs,... )
      hist <- hist_ls[[1]]
      strandvec <- hist_ls[[2]]
    } else {
      mat_list <- bi_anch_mat( bed_sub = bed_sub, ... )
      hist <- mat_list[[1]]
      ends_mat <- mat_list[[2]]
      strandvec <- mat_list[[3]]
      bs <- mat_list[[4]]
    }
    if ( debug ) print( head( hist ) )
    long_hist <- melt(
      as.data.table( hist, keep.rownames = TRUE ),
      na.rm = TRUE,
      id.vars = "rn"
    )
    if ( method == "peak_centered" | method == "single_anch" ) {
      long_ends <- long_hist[[3]] + bs
    } else {
      long_ends <- melt(
        as.data.table( ends_mat, keep.rownames = TRUE),
        na.rm=TRUE,
        id.vars = "rn")
      long_ends <- long_ends[[3]]
    }
    if (bam_or_bw == "bw") { 
      bw_sub <- import(
      readsource,
      selection = GenomicSelection( "hg38", chrom = chr, colnames = "score" )
      )
      } else {
        sbp <- ScanBamParam( 
          which = GRanges(
            seqnames = chr,
            ranges = IRanges( start = 0, end = bv$si@seqlengths[ bv$si@seqnames == chr ] )
          )
        )
        if ( pairedEnd == TRUE & readsOnly == FALSE ) {
          bam_sub <- GRanges( readGAlignmentPairs( bv$bam_info, param = sbp ) ) 
          #using pairs and coercing to GRanges is necessary if you want to keep the insert.
          if ( debug ) print( "getting bam pe" )
        } else {
          if ( debug ) print( "getting bam se" )
          bam_sub <- GRanges( readGAlignments( bv$bam_info, param = sbp ) )
        }
        if ( debug ) print( head( bam_sub ) )
        bam_sub <- bam_revcomp( bam = bam_sub, debug = debug, ... )
        bam_sub <- sbp_from_bam( bam = bam_sub, debug = debug, ... )
      }
      gr <- GRanges( 
        seqnames = chr,
        ranges = IRanges( start = long_hist[[3]] + 1, end = long_ends ),
        strand = strandvec
      )
    gr$ID <-  long_hist[[1]]
    gr$BI <- long_hist[[2]]
    if (bam_or_bw == "bam") { tmp_mat <- bam_olaps( gr = gr, bam_sub = bam_sub, ignorestrand = ignorestrand, long_hist = long_hist, debug = debug, bs = bs, ... ) }
    if (bam_or_bw == "bw") { tmp_mat <- bw_olaps( gr, bw_sub = bw_sub, ignorestrand = ignorestrand, debug = debug, bs = bs, ... ) }
    return( list( "region" = bed_sub, "scores" = tmp_mat ) )
  }
