
.date_trait <- function(sts){
	paste( collapse=',' , paste( sep = '=', names(sts), unname(sts)) )
}

.seq_format <- function( d ){
	y = '<sequence id="NAME" spec="Sequence" taxon="NAME" totalcount="4" value="SEQUENCE"/>'
	seqd = lapply(  rownames(d) , function(sid){
		y1 = gsub( pattern = 'NAME', replace=sid, y )
		seq = paste( collapse='', as.character( d[sid, ] )[1,]  ) 
		y2 = gsub( pattern ='SEQUENCE', replace= seq, y1 )
		y2 
	})
	paste( seqd, collapse = '\n' )
}
#~ y = .seq_format ( d )

#' Generates a runnable XML by inserting sequence data, a starting tree, tip dates, and a model start time to a XML skeleton
#'
#' The names of the sequence data and trees must match and must be in the standard form for inferring tip dates. The XML skeleton must have these entries to be overwritten: 
#' SEQUENCES
#' START_TREE
#' DATE_TRAIT
#' SEIJR_START
#' 
#' A XML file will be generated for each starting tree provided 
#' 
#' @param xmlfn The file name of the XML skeleton
#' @param fastafn The file name of the sequence data in fasta format 
#' @param treefn The file name of the newick format starting trees. 
#' @param start The numeric date when the model dynamics (SEIJR) will initiate 
#' @return Character string of runnable xml. Individual XMLs for each starting tree will be written to disk 
#' @export 
format_xml0 <- function( xmlfn , fastafn, treefn, start = 2020.085 ){
	
	cat( paste( 'Setting start time of seijr dynamics to be', start , '\n'))
	
	#~ 	Attributes to overwrite: 
	#~ SEQUENCES
	#~ START_TREE
	#~ DATE_TRAIT
	#~ SEIJR_START
	
	tres = read.tree( treefn ) 
	ntres <- length( tres )
	
	d = read.dna( fastafn, format = 'fasta' )
	seqdata = .seq_format( d )
	
	sids = rownames(d)
	sts <- sapply( strsplit( sids, '\\|' ), function(x){
		as.numeric( tail(x,2)[1] )
	})
	names(sts) <- sids
	datedata = .date_trait( sts ) 
	
	
	x = readLines( xmlfn ) 
	xmlofn = gsub( xmlfn, pattern='skeleton', replacement='' )
	
	for ( k in 1:ntres ){
		xk0 = gsub( x , pattern = 'START_TREE', replacement = write.tree( tres[[k]] )  )
		xk1 = gsub( xk0, pattern = 'DATE_TRAIT', replacement = datedata )
		xk2  = gsub( xk1, pattern='SEIJR_START', replacement= as.character(start) )
		xk3 = gsub( xk2, pattern='SEQUENCES', replacement = seqdata )
		
		if ( !grepl( pattern = '\\.xml$', xmlofn )  )
			writeLines( xk3, con =  paste0( xmlofn, '.', k, '.xml' )  )
		else 
			writeLines( xk3, con =  gsub( pattern = '\\.xml$', replacement = paste0('\\.',k,'\\.xml'), xmlofn ) )
	}
	invisible( xk3 )
}

# DEBUG 
if (FALSE){
	library( sarscov2 )
	library( treedater )
	library( ape ) 
	y = format_xml0( 'seijr0.1.0_skeleton.xml', 'england2.fasta', 'trees.nwk' )
}
