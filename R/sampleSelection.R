
# prototytpe: does not work 
.cutree_hclust <- function( D ){
	D <- as.matrix(D)
	g = -1
	for ( k in 1:nrow(D)){
		i <- which( D[k, ] == 0 )
		#if ( length(i) > 0 ) browser()
		{
			D[i,i] <- g
			g <- g - 1
		}
	}
	D [ D > 0 ] <- 0
	D <- abs(D) 
	g <- abs(g) - 1
	ct = apply( D, MAR=1, function( x ) max(unique(x)) )
	ct
}
#~ readRDS('tmp.rds') -> D
#~ ct = .cutree_hclust( D )


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
#' @param deduplicate_identical if TRUE, will only include one sequence (most recent) from among sets of identical sequences
#'
#' @return Writes a new fasta aligment. Return value is a treedater tree and dnabin
#' @export 
filter_quality0 <- function( path_to_align, path_to_save = NULL , q_threshold=.05, minEdge=1/29e3/10, deduplicate_identical=TRUE,  ... )
{
	library( ape ) 
	library( treedater )
	library( lubridate )

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
	D <- as.matrix( dist.dna( d, 'F84', pairwise=TRUE) )
	sts <- sapply( strsplit( rownames(d), '\\|' ) , function(x){
		decimal_date( ymd( tail(x,1)))
	})
	names(sts) <- rownames(d)
	# deduplicate identical sequences if indicated 
#~ browser()
	if (deduplicate_identical){
		drop <- c() 
		ct = cutree( hclust( as.dist(D) ) , h = 1e-6  )
		for ( k in unique( ct )){
			s = names( ct )[ ct == k ]
			if ( length( s ) > 1 ){
				u = s[ which.max(sts[s])[1]  ]
				drop <- c( drop , setdiff( s, u ) )
			}
		}
		if ( length( drop  ) > 1 ){
			keep <- setdiff( rownames(D), drop )
			D <- D[keep, keep]
			sts <- sts[ rownames(D)]
		}
	}
	tr <- nj( D )
	tr$edge.length <- pmax( tr$edge.length , minEdge)
	# treedater designed to be fast rather than accurate for outlier detection: 
#~ browser()
	td <- suppressWarnings( dater( unroot( tr ), sts, s = 29e3, omega0=.001, numStartConditions=1, maxit=1,meanRateLimits=c(.00099,.00101), ... )  )
	ot0 = outlierTips(td) 
	outs0 = as.character( ot0$taxon )[ ot0$q < q_threshold ]
	tr1 <- unroot( drop.tip( tr, outs0 ) )
	#td1 <- suppressWarnings( dater( unroot( tr1 ), sts, s = 29e3, omega0 = .001,numStartConditions=1, maxit=1,meanRateLimits=c(.00099,.00101),  ...) )
	
	d1 = d[ tr1$tip.label, ]
	
	if ( !is.null( path_to_save ))
		write.dna( d1, file = path_to_save, format = 'fasta' )
	
	n1 <- Ntip( tr1 )
	cat( paste('Removed', n0-n1, 'sequences.', n1, 'remaining. Aligment saved to', path_to_save, '\n') )
	
	list( alignment = d1, time_tree = td )
}
#~ path_to_align = '../../../gisaid/gisaid_cov2020_sequences_aligned_March14_noGaps.fasta'
#~ a = filter_quality0( path_to_align, path_to_save = NULL , q_threshold=.05, minEdge=1/29e3/10, deduplicate_identical=TRUE )






#' Remove samples from non-human hosts and sequences with incomplete date information and sequences with more than 10 per cent missing 
#'
#" @param path_to_align A DNAbin alignment OR a path to a fasta on disk 
#' @param path_to_save Path to save fasta 
#' @return DNAbin alignment
#' @export 
filter_hostAndDates <- function( path_to_align, path_to_save = NULL )
{
	library( ape ) 
	library( treedater )
	library( lubridate )

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
	
	rownames(d) <- gsub( rownames(d), pattern=' ', replacement='_' )
	
	f0fasta = path_to_save 
	if ( is.null( f0fasta ))
		f0fasta = tempfile() 
	
	write.dna(d, file=f0fasta, format='fasta')
	d
}

#' Selects a background reference set of sequences with good quality  
#'
#" Given an aligment with dates in tip labels, will remove 
#' - sequences with poor relationship with a molecular clock
#' - optionally deduplicate identical sequences 
#'
#" @param path_to_align A path (type character) where original alignment can be found, such as /gisaid/gisaid_cov2020_sequences_March14_aligned.fas
#' @param fn_tree Optional path to a maximum likelihood tree. If not provided, will estimate using fasttree. 
#' @param path_to_save Where to store (as RDS) the filtered alignment etc
#' @param q_threshold Clock outlier threshold 
#' @param minEdge minimum branch length (substitutions per site) to stabilize clock inference 
#' @param deduplicate_identical if TRUE, will only include one sequence (most recent) from among sets of identical sequences
#'
#' @return Writes a RDS file containing filtered aligment, time tree, fasttree, and sample times. Return value is a list with the same 
#' @export 
filter_quality1 <- function(path_to_align, fn_tree=NULL, path_to_save = NULL , q_threshold=.05, minEdge=1/29e3/10, deduplicate_identical=TRUE, ncpu = 6,  ... )
{
	library( ape ) 
	library( treedater )
	library( lubridate )
	
	d = read.dna( path_to_align, format = 'fasta')
	n0 <- nrow( d)
	
	# run fasttree 
	if ( is.null( fn_tree )){
		fn_tree = tempfile()
		system( paste('fasttreeMP -nt', path_to_align, ' > ', fn_tree ) )
		cat ( paste( 'Fasttree saved to ', fn_tree, '\n'))
		tr0 = read.tree( fn_tree)
		D <- as.matrix( cophenetic.phylo( tr0 ) ) 
	}
	
	# remove outliers according to a strict molecular clock 
	sts <- sapply( strsplit( rownames(d), '\\|' ) , function(x){
		decimal_date( ymd( tail(x,1)))
	})
	names(sts) <- rownames(d)
	# deduplicate identical sequences if indicated 
	if (deduplicate_identical){
		drop <- c() 
		keep <- rownames(D)
		ct = cutree( hclust( as.dist(D) ) , h = 1e-6  )
		for ( k in unique( ct )){
			s = names( ct )[ ct == k ]
			if ( length( s ) > 1 ){
				u = s[ which.max(sts[s])[1]  ]
				drop <- c( drop , setdiff( s, u ) )
			}
		}
		if ( length( drop  ) > 1 ){
			keep <- setdiff( rownames(D), drop )
			D <- D[keep, keep]
			sts <- sts[ rownames(D)]
		}
		cat( paste( 'Removed', length( drop ), 'identical sequences\n' ) )
	}
	tr <- keep.tip( tr0, keep )#nj( D )
	tr$edge.length <- pmax( tr$edge.length , minEdge)
	# treedater designed to be fast rather than accurate for outlier detection: 
	td <- suppressWarnings( dater( unroot( tr ), sts, s = 29e3, omega0=.001, numStartConditions=1, maxit=1,meanRateLimits=c(.00099,.00101) , ncpu = ncpu )  )
	ot0 = outlierTips(td) 
	outs0 = as.character( ot0$taxon )[ ot0$q < q_threshold ]
	tr1 <- unroot( drop.tip( tr, outs0 ) )
	#td1 <- suppressWarnings( dater( unroot( tr1 ), sts, s = 29e3, omega0 = .001,numStartConditions=1, maxit=1,meanRateLimits=c(.00099,.00101),  ...) )
	
	d1 = d[ tr1$tip.label, ]
	
	n1 <- Ntip( tr1 )
	cat( paste('Removed', n0-n1, 'sequences.', n1, 'remaining. Aligment saved to', path_to_save, '\n') )
	
	rv = list( alignment = d1, time_tree = td , fasttree = tr0 , sts = sts )
	if ( !is.null( path_to_save ))
		saveRDS( rv, file = path_to_save )
		#write.dna( d1, file = path_to_save, format = 'fasta' )
	
	rv
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
	library( ape ) 
	library( treedater )
	library( lubridate )

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
#' AND including close distance matches to a subset of sequences from a particular region
#' All sequences from given location will be included 
#'
#' Select samples based on quantile of sample time distribution. Requires date to be at and of sequence label. Alternatively can do a simple random sample within a region
#'
#' @param region_regex Sample names matching this regular expression will be retained and closest matches also retained
#" @param path_to_align A DNAbin alignment *or* a system path (type character) where original alignment can be found, such as /gisaid/gisaid_cov2020_sequences_March14_aligned.fas
#' @param path_to_save Where to store (as fasta) the filtered alignment
#' @param n sample size from outside region 
#' @param nregion sample size within region; if null will include everything in region
#' @param q_threshold Clock outlier threshold 
#' @param minEdge minimum branch length (substitutions per site) to stabilize clock inference 
#' @param time_stratify_region If TRUE (default) will perform a time stratified sample within region, otherwise will do a simple random sample 
#' 
#' @return A DNAbin alignment. Will also save to path_to_save 
#' @export 
region_time_stratified_sample <- function(region_regex, n, path_to_align, nregion = NULL,  path_to_save = NULL, time_stratify_region=TRUE ) {
	library( ape ) 
	library( treedater )
	library( lubridate )

	if ( inherits( path_to_align, 'DNAbin' ) )
		d = path_to_align
	else
		d = read.dna( path_to_align, format = 'fasta')
	
	dr = d[ grepl(pattern= region_regex, x = rownames(d) ) , ]
	dnr = d[ setdiff( rownames(d), rownames(dr)), ]
	
	stsdr <- sapply( strsplit( rownames(dr), '\\|' ) , function(x){
		decimal_date( ymd( tail(x,1)))
	})
	names( stsdr ) <- rownames( dr )
	mrsid = names( stsdr )[ which.max( stsdr ) ] # most recent in region 
	
	nr <- nrow(dr ) 
	if (!is.null( nregion )){
		if ( nregion < nr ){
			if ( time_stratify_region){
				#dr <- dr[ sample( rownames(dr), replace=FALSE, size = nregion), ]
				dr <- time_stratified_sample(nregion, dr, path_to_save = NULL )$alignment
			} else{ # do a simple random sample in region 
				include = sample( setdiff( rownames(dr), mrsid ), size = nregion-1, replace=FALSE )
				include <- c( include, mrsid )
				dr <- dr[ include , ] 
			}
		}
	}
	
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

#' Change the name of sequences to recognize numeric time of sampling and deme 
#'
#' @param path_to_algn A DNAbin alignment *or* a system path (type character) where original alignment can be found
#' @param deme The name of the deme to add to the end of each label 
#' @param regexprs A vector of regular expressions that can be used to match tip labels and categorize in to demes. These must include *all* sequcence IDs, and each ID should match at most one regular expression
#' @param invert_regexpr A logical vector specifying if the i'th regular expression should be an inverse match 
#' @param demes A character vector of deme names to be appended to corresponding sequence IDs 
#'
#" @return DNAbin alignment 
#' @export 
#~ prep_tip_labels_phydyn <- function( path_to_align, path_to_save = NULL, deme = 'Il'  ){
prep_tip_labels_phydyn <- function( path_to_align, path_to_save = NULL
 , regexprs = c( '.*/Netherlands/.*', '.*/Netherlands/.*' ) 
 , invert_regexpr = c( FALSE, TRUE )
 , demes = c( 'Il'  , 'exog'  )
 ){
	library( ape ) 
	library( treedater )
	library( lubridate )
	
	if ( inherits( path_to_align, 'DNAbin' ) )
		d = path_to_align
	else
		d = read.dna( path_to_align, format = 'fasta')
	
	sids = rownames(d) 
	if ( length( regexprs ) != length( demes ))
		stop('Must provide equal numbers of regex and deme names ') 
	
	demegroups = lapply( 1:length(demes), function(k) {
		x = regexprs[k]
		if ( invert_regexpr[k] ) {
			return( sids[ !grepl( pattern = x , sids ) ] )
		}else{
			return( sids[ grepl( pattern = x , sids ) ] )
		}
	})
	int <- do.call( intersect, demegroups )
	uni = do.call( c, demegroups )
	if ( length( int ) > 0 )
		stop( 'Intersection of deme groups is non-empty. Each regex must match a unique set.' )
	if ( length( uni ) < length(sids) ){
		print( setdiff( sids, uni ))
		stop( 'There were some sequence IDs that did not match a regex. ' )
	}
	
	deme <- setNames( rep(demes[1], nrow(d)), sids )
	for ( k in 1:length( demes )){
		deme[ demegroups[[k]] ] <- demes[k]
	}
	
	sts <- sapply( strsplit( rownames(d), '\\|' ) , function(x){
		decimal_date( ymd( tail(x,1)))
	})
	rownames(d) <- paste(sep='|', rownames(d), sts, paste0('_', deme) )
	rownames(d) <- gsub( rownames(d), pattern = '\\s' , replacement = '_')
	if ( !is.null( path_to_save ))
		write.dna( d, file = path_to_save, format = 'fasta' )
	d
}
	
# debug 
if (FALSE){
	p0 = '/home/erikvolz/git/sarscov2-phylodynamics/gisaid/gisaid_cov2020_sequences_aligned_March14_noGaps.fasta'
	o = filter_quality0( p0, NULL, ncpu = 6 )
#~ 	o1 = time_stratified_sample( 50, o$alignment,  '/home/erikvolz/git/sarscov2-phylodynamics/gisaid/gisaid_cov2020_sequences_March14_aligned_filter0_tempsamp50.fas')
	o1 = region_time_stratified_sample( '.*/Netherlands/.*' , 25, o$alignment,  '/home/erikvolz/git/sarscov2-phylodynamics/gisaid/gisaid_cov2020_sequences_aligned_March14_noGaps_filter0_netherlands.fas')
	o1 = region_time_stratified_sample( '.*/USA/WA.*' , 25, o$alignment,  '/home/erikvolz/git/sarscov2-phylodynamics/gisaid/gisaid_cov2020_sequences_aligned_March14_noGaps_filter0_WA.fas')
	
		library( ape ) 
	library( treedater )
	library( lubridate )
	o = prep_tip_labels_phydyn( '../gisaid_cov2020_sequences_aligned_March14_noGaps_filter0_netherlands.fas'
	 , regexprs = c( '.*/Netherlands/.*', '.*/Netherlands/.*' ) 
	 , demes = c( 'Il'  , 'exog'  )
	 , path_to_save = 'tmp.fasta' 
	)

}
