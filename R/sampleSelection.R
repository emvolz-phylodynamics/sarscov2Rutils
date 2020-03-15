library( ape ) 
library( treedater )
library( lubridate )


#' Selects a background reference set of sequences with good quality  
#'
#" Given an aligment with dates in tip labels, will remove 
#' - bat/pangolin/canine
#' - sequences with many gaps 
#' - sequences within incomplete dates
#' - sequences with poor relationship with a molecular clock
#'
#" @param path_to_align A DNAbin alignment *or* a system path (type character) where original alignment can be found, such as /gisaid/gisaid_cov2020_sequences_March14_aligned.fas
#' @param path_to_save Where to store (as fasta) the filtered alignment
#' @param q_threshold Clock outlier threshold 
#' @param minEdge minimum branch length (substitutions per site) to stabilize clock inference 
#'
#' @return Writes a new fasta aligment. Return value is a treedater tree and dnabin
#' @export 
filter_quality0 <- function( path_to_align, path_to_save = NULL , q_threshold=.05, minEdge=1/29e3/10, ... )
{
	if ( inherits( path_to_align, 'DNAbin' ) )
		d = path_to_align
	else
		d = read.dna( path_to_align, format = 'fasta')
	
	n0 <- nrow(d)
	
	# remove bat and pangolin and canine if present 
	i = grepl( pattern = 'bat' , rownames(d)) | grepl( pattern = 'pangolin' , rownames(d)) | grepl( pattern = 'canine' , rownames(d))
	d = d[!i, ]

	# remove any with incomplete dates if present
	dts = lapply( strsplit( rownames(d), '\\|' ) , function(x){
		suppressWarnings( ymd( tail(x,1)  ) )
	})
	keep <- sapply( dts, function(x) !is.na(x) )
	d <- d[ keep , ]

	# remove seq with more than 10% missing 
	ngaps = sapply( 1:nrow(d) , function(i) sum( as.character(d[i,])=='-' ) )
	keep <- ngaps < (.1*ncol(d))
	d = d[keep, ]
	
	# remove outliers according to a strict molecular clock 
	D <- dist.dna( d, 'F84', pairwise=TRUE)
	sts <- sapply( strsplit( rownames(d), '\\|' ) , function(x){
		decimal_date( ymd( tail(x,1)))
	})
	names(sts) <- rownames(d)
	tr <- nj( D )
	tr$edge.length <- pmax( tr$edge.length , minEdge)
	td <- suppressWarnings( dater( unroot( tr ), sts, s = 29e3, omega0=.001, ... )  )
	ot0 = outlierTips(td) 
	outs0 = as.character( ot0$taxon )[ ot0$q < q_threshold ]
	tr1 <- unroot( drop.tip( tr, outs0 ) )
	td1 <- suppressWarnings( dater( unroot( tr1 ), sts, s = 29e3, omega0 = .001, ...) )
	
	d1 = d[ tr1$tip.label, ]
	
	if ( !is.null( path_to_save ))
		write.dna( d1, file = path_to_save, format = 'fasta' )
	
	n1 <- Ntip( tr1 )
	cat( paste('Removed', n0-n1, 'sequences.', n1, 'remaining. Aligment saved to', path_to_save, '\n') )
	
	list( alignment = d1, time_tree = td1 )
}

#' Select a random sample from aligment stratified through time 
#'
#' Select samples based on quantile of sample time distribution. Requires date to be at and of sequence label
#'
#" @param path_to_align A DNAbin alignment *or* a system path (type character) where original alignment can be found, such as /gisaid/gisaid_cov2020_sequences_March14_aligned.fas
#' @param path_to_save Where to store (as fasta) the filtered alignment
#' @param q_threshold Clock outlier threshold 
#' @param minEdge minimum branch length (substitutions per site) to stabilize clock inference 
#' 
#' @return A DNAbin alignment. Will also save to path_to_save 
#' @export 
time_stratified_sample <- function(n, path_to_align, path_to_save = NULL ) {
	if ( inherits( path_to_align, 'DNAbin' ) )
		d = path_to_align
	else
		d = read.dna( path_to_align, format = 'fasta')
	
	sts <- sapply( strsplit( rownames(d), '\\|' ) , function(x){
		decimal_date( ymd( tail(x,1)))
	})
	names( sts ) <- rownames( d )
	
	ssts = sort( sts )
	N <- length( sts )
	ibins = floor( seq( N/n , N , length = n ) )
	i0 = 1 
	keep <- c() 
	for ( k in 1:n){
		keep <- c( keep, sample( names(ssts)[i0:ibins[k]], size = 1 )  )
		i0 <- ibins[ k ] + 1
	}
	keep <- unique( keep )
	
	d1 = d[keep, ]
	
	if ( !is.null( path_to_save ))
		write.dna( d1, file = path_to_save, format = 'fasta' )
	
	list( 
		alignment = d1 
		, sids = keep 
		, sts= sts[keep ]
	) 
}

#' Select a random sample from aligment stratified through time 
#' AND including close distance matches to a subset of sequences from a particular location
#' All sequences from given location will be included 
#'
#' Select samples based on quantile of sample time distribution. Requires date to be at and of sequence label
#'
#' @param region_regex Sample names matching this regular expression will be retained and closest matches also retained
#" @param path_to_align A DNAbin alignment *or* a system path (type character) where original alignment can be found, such as /gisaid/gisaid_cov2020_sequences_March14_aligned.fas
#' @param path_to_save Where to store (as fasta) the filtered alignment
#' @param q_threshold Clock outlier threshold 
#' @param minEdge minimum branch length (substitutions per site) to stabilize clock inference 
#' 
#' @return A DNAbin alignment. Will also save to path_to_save 
#' @export 
region_time_stratified_sample <- function(region_regex, n, path_to_align, path_to_save = NULL ) {
	if ( inherits( path_to_align, 'DNAbin' ) )
		d = path_to_align
	else
		d = read.dna( path_to_align, format = 'fasta')
	
	dr = d[ grepl(pattern= region_regex, x = rownames(d) ) , ]
	dnr = d[ setdiff( rownames(d), rownames(dr)), ]
	D <- dist.dna( d, 'F84' , pairwise.deletion=TRUE )
	Drnr <- as.matrix( D )[ rownames(dr), rownames(dnr ) ]
	keep <- c()
	for ( i in 1:nrow(dr )){
		keep <- c( keep, colnames(Drnr)[ which.min( Drnr[i,] ) ] )
	}
	keep <- unique( keep )
	keep <- c( rownames( dr ), keep )
	
	dr2 = d[ keep, ]
	dnr2 = d[ setdiff( rownames(d), keep ), ]
	o2 = time_stratified_sample(n, dnr2, path_to_save = NULL )
	d3 = rbind( dr2, o2$alignment )
	if ( !is.null( path_to_save ))
		write.dna( d3, file = path_to_save, format = 'fasta' )
	
	d3
}

# debug 
if (FALSE){
	p0 = '/home/erikvolz/git/sarscov2-phylodynamics/gisaid/gisaid_cov2020_sequences_aligned_March14_noGaps.fasta'
	o = filter_quality0( p0, NULL, ncpu = 6 )
#~ 	o1 = time_stratified_sample( 50, o$alignment,  '/home/erikvolz/git/sarscov2-phylodynamics/gisaid/gisaid_cov2020_sequences_March14_aligned_filter0_tempsamp50.fas')
	o1 = region_time_stratified_sample( '.*/Netherlands/.*' , 25, o$alignment,  '/home/erikvolz/git/sarscov2-phylodynamics/gisaid/gisaid_cov2020_sequences_aligned_March14_noGaps_filter0_netherlands.fas')
	o1 = region_time_stratified_sample( '.*/USA/WA.*' , 25, o$alignment,  '/home/erikvolz/git/sarscov2-phylodynamics/gisaid/gisaid_cov2020_sequences_aligned_March14_noGaps_filter0_WA.fas')
}
