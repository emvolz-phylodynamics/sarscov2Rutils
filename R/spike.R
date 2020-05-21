#' Determine the spike 614 genotype for the given fasta alignment
#"
#' NOTE: requires mafft which will be used to realign input data. Requires seqinr package 
#' 
#' @param algnfn path to fasta with nucleotide alignment. This will be re-aligned, so no need for standard coordinates
#' @param algnfn2 path to which realigned sequences will be saved 
#' @return A data frame with a column giving the genotype 
#' @export 
compute_spike614genotype <- function( algnfn, algnfn2 = 'standard_coords.fasta' )
{
	reffn = tempfile() 
	spikereffn = tempfile() 
	file.copy( system.file( 'extdata/ref.fasta' , package = 'sarscov2'), reffn )
	file.copy( system.file( 'extdata/spikeref.fasta' , package = 'sarscov2'), spikereffn )
	
	system( paste( 'mafft --thread -1 --keeplength --add', algnfn, reffn , ' > ' , algnfn2) )
		
	library( seqinr )
	spike_coords <- 21563:25384
	a0 = read.fasta( algnfn2 )
	b1 = lapply( a0, function(x){
		seqinr::as.SeqFastaAA( translate( x[spike_coords] )
		  , name=attr(a0[[1]], 'name' )
		  , Annot= attr(a0[[1]], 'Annot' ))
	})
	resultdf = data.frame( s614= sapply( b1, '[', 614 ) )
	saveRDS( resultdf , file = 's614df.rds' )
	cat( 'Done. Check the dataframe so lineages have the appropriate genotype (D or G for spike 614)\n' )
	resultdf 
}
#~ o = compute_spike614genotype( 'algn3.fasta', 'algn4.fasta' )

#' Run skygrowth and treestructure on clades defined by s614 polymorphism 
#'
#' @param tds List or multiPhylo of treedater trees or ape::phylo
#' @param s614 A data frame produced by compute_spike614genotype
#' @param ... Additional arguments passed to sarscov2::skygrowth1
#' @return A list with skygrowth for G and D, and a list of treestructure tests for each tree 
#' @examples
#' \dontrun{
#' library( sarscov2 ) 
#'  spikedf = compute_spike614genotype( 'algn3.fasta', 'algn4.fasta' )
#'  # inspect spikedf to make sure it a assigned D or G 
#'  # some precomputed treedater trees corersponding to algn3.fasta: 
#'  tds = readRDS('tds.rds' )  
#'  result = s614_phylodynamics( tds, spikedf)
#'  print( result$treestructureTests  )
#'  # make a plot 
#'  add_ggR_sarscov2skygrowth( result$sgD, result$sgG )
#' }
#' @export 
s614_phylodynamics <- function( tds, s614 , ... )
{
	library( sarscov2 )
	library( treestructure ) 
	library( ape )
	
	tres = lapply( tds, function(x) {class(x) = 'phylo'; x } )
	
	s614$taxon <- rownames( s614 ) 

	s614.region = s614[ grepl(s614$taxon, patt = '_Il') , ]
	Dtips = s614.region$taxon[ s614.region$s614=='D' ]
	Gtips = s614.region$taxon[ s614.region$s614=='G' ]
	trs.region = lapply( tres, function(x) keep.tip( x, s614.region$taxon ) )
	trs.D <- lapply( tres, function(x) keep.tip( x, Dtips ))
	trs.G <- lapply( tres, function(x) keep.tip( x, Gtips ))
	
	# test structure 
	tests = lapply( trs.region  , function(tre) treestructure.test( tre, Gtips, Dtips, nsim = 1e4 ))
	
	sgD = skygrowth1( trs.D , ... ) 
	sgG = skygrowth1( trs.G , ... ) 
	
	list( sgG = sgG , sgD = sgD, treestructureTests = tests )
}

#~ s614 = readRDS('s614df.rds' ) 
#~ tds = readRDS('tds.rds' ) 
#~ o = s614_phylodynamics( tds, s614 )

