
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
#' SUSC_SIZE 
#' 
#' A XML file will be generated for each starting tree provided 
#' 
#' @param xmlfn The file name of the XML skeleton
#' @param fastafn The file name of the sequence data in fasta format 
#' @param treefn The file name of the newick format starting trees. 
#' @param start The numeric date when the model dynamics (SEIJR) will initiate 
#' @param susc_size The initial number of susceptible, replaces SUSC_SIZE in xml, also sets scale for exponential size prior; supported in seijr0.1.3_skeleton.xml
#' @return Character string of runnable xml. Individual XMLs for each starting tree will be written to disk 
#' @export 
format_xml0 <- function( xmlfn , fastafn, treefn, start = 2020.085, susc_size = NULL ){
	library( ape )
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
		if (!is.null(susc_size))
			xk4 = gsub( xk3, pattern='SUSC_SIZE', replacement = susc_size )
		
		if ( !grepl( pattern = '\\.xml$', xmlofn )  )
			writeLines( xk4, con =  paste0( xmlofn, '.', k, '.xml' )  )
		else 
			writeLines( xk4, con =  gsub( pattern = '\\.xml$', replacement = paste0('\\.',k,'\\.xml'), xmlofn ) )
	}
	invisible( xk4 )
}



#' Assemble a XML fragment that will define beta(t) aas a linear spline with a given number of change points 
gen_gmrf_betablock <- function(   tend, tstart = 2020.038, numchange = 3){
	#tstart = 2020.085
	#tend = 2020.25
	#numchange = 3
	library( glue )
	
	tknots = seq( tstart, tend, length.out = numchange+2 )
	ctknots <- as.list( as.character( tknots ))
	n <- length( ctknots )
	names(ctknots) <- paste0( 't', 1:n)
	
	betalb = '<definition id="Definition.betalb" spec="phydyn.model.Definition" value="betalb = if ( t !>= 2020.035 ) then 0'
	betalb = paste(sep='\n', betalb
		, with( ctknots, glue( 'else if (t  !>= {t2} ) then beta1' ) )
		)
	for ( k in 3:(n-1)){
		dlogbetasum = paste(collapse=' + ', sapply( 0:(k-3), function(j) with ( list(y=as.character(j)), glue('dlogbeta[{y}]') ) ) ) 
		tindex = ctknots[ paste0( 't', k ) ]
		betalb = paste(sep='\n', betalb
		, glue( 'else if (t  !>= {tindex} ) then beta1 * exp( {dlogbetasum} )' )
		)
	}
	dlogbetasum = paste(collapse=' + ', sapply( 0:(n-3), function(j) with ( list(y=as.character(j)), glue('dlogbeta[{y}]') ) ) ) 
	betalb = paste(sep='\n', betalb
		, glue( 'else ( beta1 * exp( {dlogbetasum} ) )' )
		, '"/>'
		)
	betalb = glue( betalb )
	
	betaub = '<definition id="Definition.betaub" spec="phydyn.model.Definition" value="betaub = if ( t !>= 2020.035 ) then 0'
	for ( k in 2:(n-1)){
		dlogbetasum = paste(collapse=' + ', sapply( 0:(k-2), function(j) with ( list(y=as.character(j)), glue('dlogbeta[{y}]') ) ) ) 
		tindex = ctknots[ paste0( 't', k ) ]
		betaub = paste(sep='\n', betaub
		, glue( 'else if (t  !>= {tindex} ) then beta1 * exp( {dlogbetasum} )' )
		)
	}
	dlogbetasum = paste(collapse=' + ', sapply( 0:(n-3), function(j) with ( list(y=as.character(j)), glue('dlogbeta[{y}]') ) ) ) 
	betaub = paste(sep='\n', betaub
		, glue( 'else ( beta1 * exp( {dlogbetasum} ) )' )
		, '"/>'
		)
	betaub = glue( betaub )

	tlb =  with( ctknots, glue( '<definition id="Definition.tlb" spec="phydyn.model.Definition" value=" tlb  = if (t  !>=  2020.035 ) then {t1} ' ) )
	for ( k in 3:(n)){
		tindex0 = ctknots[ paste0( 't', k-1 ) ]
		tindex1 = ctknots[ paste0( 't', k ) ]
		tlb = paste(sep='\n', tlb
		, glue( 'else if (t  !>= {tindex1} ) then {tindex0}' )
		)
	}
	tindex = ctknots[ paste0( 't', n ) ]
	tlb = glue( paste(sep='\n', tlb,   'else({tindex}) "/>' ) )

	difftime = as.character( tknots[2] - tknots[1] )
	propub = '<definition id="Definition.propub" spec="phydyn.model.Definition" value=" propub = max(0, min(1, (t-tlb)/ ({difftime})  ) ) "/>'
	
	beta = '<definition id="Definition.beta" spec="phydyn.model.Definition" value=" beta = max(0, propub * betaub + (1-propub)*betalb) "/>'
	
	glue( paste( sep = '\n', betalb, betaub, tlb, propub, beta  ) )
}

#' Generates a runnable XML by inserting sequence data, a starting tree, tip dates, and number of gmrf change points into a XML skeleton
#'
#' The names of the sequence data and trees must match and must be in the standard form for inferring tip dates. The XML skeleton must have these entries to be overwritten: 
#' SEQUENCES
#' START_TREE
#' DATE_TRAIT
#' BETADEFS
#' DLOGBETADIM
#' 
#' A XML file will be generated for each starting tree provided 
#' 
#' @param xmlfn The file name of the XML skeleton
#' @param fastafn The file name of the sequence data in fasta format 
#' @param treefn The file name of the newick format starting trees. 
#' @param numchange integer number of change points for beta(t) 
#' @return Character string of runnable xml. Individual XMLs for each starting tree will be written to disk 
#' @export 
format_xml_gmrf_exogsir <- function( xmlfn , fastafn, treefn, numchange =3 ){
	library( ape )
	cat( paste( 'Setting number of change points in beta(t)', numchange , '\n'))
	
	#~ Attributes to overwrite: 
	#~ SEQUENCES
	#~ START_TREE
	#~ DATE_TRAIT
	#~ DLOGBETADIM
	#~ BETADEFS
	
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
	
	betablock = gen_gmrf_betablock(tend=max(sts), tstart = 2020.038, numchange = numchange)
	
	x = readLines( xmlfn ) 
	xmlofn = gsub( xmlfn, pattern='skeleton', replacement='' )
	
	for ( k in 1:ntres ){
		xk0 = gsub( x , pattern = 'START_TREE', replacement = write.tree( tres[[k]] )  )
		xk1 = gsub( xk0, pattern = 'DATE_TRAIT', replacement = datedata )
		xk2  = gsub( xk1, pattern='DLOGBETADIM', replacement= as.character(numchange) )
		xk3 = gsub( xk2, pattern = 'BETADEFS', replacement = betablock )
		xk4 = gsub( xk3, pattern='SEQUENCES', replacement = seqdata )
		
		if ( !grepl( pattern = '\\.xml$', xmlofn )  )
			writeLines( xk4, con =  paste0( xmlofn, '.', k, '.xml' )  )
		else 
			writeLines( xk4, con =  gsub( pattern = '\\.xml$', replacement = paste0('\\.',k,'\\.xml'), xmlofn ) )
	}
	invisible( xk4 )
}





# DEBUG 
if (FALSE){
	library( sarscov2 )
	library( treedater )
	library( ape ) 
	y = format_xml0( 'seijr0.1.0_skeleton.xml', 'england2.fasta', 'trees.nwk' )
}
