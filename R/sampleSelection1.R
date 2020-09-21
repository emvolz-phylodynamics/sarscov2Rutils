#' Sample sequence IDs from within a region 
#'
#' @param md data frame with gisaid metadata
#' @param n sample size within region 
#' @param inclusion_rules  Required list of rules of the form c( <meta data column> , <regular expression> ). Where this pattern matches, a sequence will be included in the sample
#' @param dedup If TRUE identical sequences will not be included 
#' @param time_stratify If TRUE, will collect a time stratified sample from within the region instead of simple random sample 
#' @examples
#'\dontrun{
#' #This will get a sample from King County in Washington and exclude sequences labelled Washington from the exog sample, since we aren't sure if they are from King County or not
#'	ipatt = '^KingCounty$'
#'	epatt = '.*Washington.*'
#'	regiontips = region_sampler1( md, n = 10  , inclusion_rules = list( c('CityOrCounty', ipatt) ))
#'	exogtips = exog_sampler1( md, 20, D, s, exclusion_rules = list( c('CityOrCounty', epatt) )  )
#'}
#' @export
region_sampler1 <- function( md
 , n = 50 
 , inclusion_rules = list() #list ( c( 'CityOrCounty' , '^Manhattan$') )
 , dedup = TRUE 
 , time_stratify = TRUE 
){
	stopifnot( length( inclusion_rules ) > 0 )
	include = do.call( c, lapply( inclusion_rules, function(r){
		f = as.character( md[[r[1]]] )
		i = grepl( x = f , pattern = r[2] , perl=TRUE)
		if ( dedup )
			i <- i & (md$inNoDups==1)
		md$seqName[ i ]
	}))
	include = unique( include )
	if ( n > length( include ))
		return(include)
	
	mdr = md[ match( include, md$seqName ), ]
	sts <- lubridate::ymd( as.character( mdr$sampleDate))
	names( sts ) <- include 
	
	if ( time_stratify )
	{
	    ssts = sort(sts)
		N <- length(sts)
		ibins = unique(floor(seq(2, N, length = n-1)))
      		keep <- c(names(ssts)[1])
     		for (k in 2:(n-1)) {
        		keep <- c(keep, sample(names(ssts)[ibins[k - 1]:(ibins[k]-1)], 
                               size = 1))
      		}
		
		
		keep <- c( keep, tail(names(ssts),1) )
	} else {
		keep <- sample( include, replace=FALSE, size = n)
	}
	return (keep)
}


#' Sample exogenous sequence IDs using a time stratified design and including close genetic distance to a regional sample 
#'
#' @param md data frame with gisaid metadata 
#' @param n sample size 
#' @param D distance matrix between sequences
#' @param region_sample vector of sequence IDs from region
#' @param exclusion_rules Optional list of rules of the form c( <meta data column> , <regular expression> ). Where this pattern matches, a sequence will not be included in the sample
#' @param dedup If TRUE will not include identical sequences
#' @examples
#'\dontrun{
#' #This will get a sample from King County in Washington and exclude sequences labelled Washington from the exog sample, since we aren't sure if they are from King County or not
#'	ipatt = '^KingCounty$'
#'	epatt = '.*Washington.*'
#'	regiontips = region_sampler1( md, n = 10  , inclusion_rules = list( c('CityOrCounty', ipatt) ))
#'	exogtips = exog_sampler1( md, 20, D, s, exclusion_rules = list( c('CityOrCounty', epatt) )  )
#'}
#' @export 
exog_sampler1 <- function( 
	md
	, n
	, D 
	, region_sample 
	, exclusion_rules = list() # list ( c( 'CityOrCounty' , '.*NewYork.*') )
	, dedup = TRUE
){
	exclude <- region_sample
	if ( length( exclusion_rules ) > 0){
		exclude0 = do.call( c, lapply( exclusion_rules, function(r){
			f = as.character( md[[r[1]]] )
			i = grepl( x = f , pattern = r[2] , perl=TRUE)
			md$seqName[ i ]
		}))
		exclude <- unique( c( exclude, exclude0 ))
	}
	
	mde = md[ !(md$seqName %in% exclude) , ]
	if ( dedup )
		mde <- mde[ mde$inNoDups == 1, ]
	
	# find close matches 
	D <- as.matrix(D)
	Drr <- D[region_sample, region_sample ]
	Drnr <- D[region_sample, mde$seqName]
	Drnr[is.na( Drnr ) ] <- Inf 
	keep <- c()
	nregion <- length( region_sample )
	for ( i in 1:nregion){
		keep <- c( keep, colnames(Drnr)[ which.min( Drnr[i,] ) ] )
	}
	keep <- unique( keep )
	
	# time strat sample in exog 
	sts <- lubridate::ymd( as.character( mde$sampleDate))
	names(sts) <- mde$seqName 
	ssts = sort( sts )
	N <- length(ssts)
    ibins = floor(seq(N/n, N, length = n))
    i0 = 1
    keepstrat <- c()
    for (k in 1:n) {
        keepstrat <- c(keepstrat, sample(names(ssts)[i0:ibins[k]], size = 1))
        i0 <- ibins[k] + 1
    }
    keepstrat <- unique(keepstrat)

	c( keepstrat, keep )
}


#' Sample exogenous sequence IDs using a time stratified design and including close genetic distance to a regional sample 
#'
#' @param md data frame with gisaid metadata 
#' @param n sample size 
#' @param smallGDpairs A data frame or path to csv containing <ID1 ID2 Distance> such as produced by the 'tn93' tool 
#' @param region_sample vector of sequence IDs from region
#' @param exclusion_rules Optional list of rules of the form c( <meta data column> , <regular expression> ). Where this pattern matches, a sequence will not be included in the sample
#' @param dedup If TRUE will not include identical sequences
#' @examples
#'\dontrun{
#' #This will get a sample from King County in Washington and exclude sequences labelled Washington from the exog sample, since we aren't sure if they are from King County or not
#'	ipatt = '^KingCounty$'
#'	epatt = '.*Washington.*'
#'	regiontips = region_sampler1( md, n = 10  , inclusion_rules = list( c('CityOrCounty', ipatt) ))
#'	exogtips = exog_sampler1( md, 20, D, s, exclusion_rules = list( c('CityOrCounty', epatt) )  )
#'}
#' @export 
exog_sampler2 <- function( 
	md
	, n
	, smallGDpairs
	, region_sample 
	, exclusion_rules = list() # list ( c( 'CityOrCounty' , '.*NewYork.*') )
	, dedup = FALSE
){
	if ( is.character( smallGDpairs ))
		smallGDpairs = read.csv( smallGDpairs, stringsAs=FALSE )
	
	exclude <- region_sample
	if ( length( exclusion_rules ) > 0){
		exclude0 = do.call( c, lapply( exclusion_rules, function(r){
			f = as.character( md[[r[1]]] )
			i = grepl( x = f , pattern = r[2] , perl=TRUE)
			md$seqName[ i ]
		}))
		exclude <- unique( c( exclude, exclude0 ))
	}
	
	mde = md[ !(md$seqName %in% exclude) , ]
	if ( dedup )
		mde <- mde[ mde$inNoDups == 1, ]
	
	# find close matches 
	nregion <- length( region_sample )
	rsgdp = smallGDpairs[ smallGDpairs$ID1 %in% region_sample, ]
	rsgdp = rsgdp[ rsgdp$ID2 %in% mde$seqName, ] # do not match to exclude , filter out dup's if applicable 
	keep <- c() 
	us = intersect( rsgdp$ID1, region_sample )
	if ( length( us ) > 0 ){
		#		keep0 <- sapply( us , function(u) {
		#			rsgdp1 = rsgdp[ rsgdp$ID1==u, ]
		#			rsgdp1$ID2[ which.min( rsgdp1$Distance )[1] ]
		#		})
		closeexog <- lapply( us , function(u) {
			rsgdp1 = rsgdp[ rsgdp$ID1==u, ]
			rsgdp1$ID2[ which.min( rsgdp1$Distance ) ]
		})
		names( closeexog ) <- us 
		# prioritize exog that are close to more than one regional 
		ce2 <- do.call( c, closeexog )
		tce2 = table( ce2 )
		keep = names( tce2 )[ tce2 > 1 ]
		
		for ( u in us ){
			if (!any( closeexog[[u]] %in% keep )){
				keep <- c( keep, closeexog[[u]][1] )
			}
		}
	}
	
	keep <- unique( keep )
	
	# time strat sample in exog 
	sts <- lubridate::ymd( as.character( mde$sampleDate))
	names(sts) <- mde$seqName 
	ssts = sort( sts )
	N <- length(ssts)
    ibins = floor(seq(N/n, N, length = n))
    i0 = 1
    keepstrat <- c()
    for (k in 1:n) {
        keepstrat <- c(keepstrat, sample(names(ssts)[i0:ibins[k]], size = 1))
        i0 <- ibins[k] + 1
    }
    keepstrat <- unique(keepstrat)

	unique( c( keepstrat, keep ) )
}




#' Prepare an alignment of regional and exogenous sequences with tip labels for SEIJR in PhyDyn
#'
#' @param algnfn path to gisaid aligment 
#' @param outfn path to save analysis alignment 
#' @param regiontips vector of sequence id
#' @param exogtips vector of sequence id
#' @param metadata data frame with gisaid metadata 
#' @param useBiostrings set TRUE to use Biostrings::readDNAStringSet - more efficient for large alignments 
#' @export
prep_tip_labels_seijr <- function( algnfn , outfn , regiontips, exogtips, metadata, useBiostrings=FALSE ){
	md = metadata 
	if(useBiostrings==TRUE){
	library(Biostrings)
	d <- as.DNAbin(readDNAStringSet(algnfn))
	s = intersect(c(regiontips, exogtips), names(d))
	dd=d[s]
	} else {
	d = read.dna( algnfn, 'fasta' )
	s= intersect( c( regiontips, exogtips ) , rownames(d))
	dd=d[s, ]
	}
	nms = names(dd)
	.md <- md [ match( s, md$seqName ) , ]
	sts <- setNames( lubridate::decimal_date( lubridate::ymd( as.character( .md$sampleDate))), .md$seqName ) 
	demes <- setNames( rep( 'exog', length( nms ) ), nms )
	demes [ nms %in% regiontips ] <- 'Il'
	nms =  paste(sep = "|", nms, sts[nms], paste0("_",  demes[nms]))
	names(dd)  <- nms
	write.dna( dd, file=outfn, format='fasta' )
	invisible( dd )
}
