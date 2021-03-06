
#~ library( sarscov2)
#~ library( ape ) 
#~ library( treedater )
#~ library( phangorn )
#~ library( ggtree )
#~ library( phangorn )
#~ library( Hmisc)

#' Extact all clades in given tree descended from nodes with high confidece of being within region (parsimony)
#'
#' @export
region_clades<- function(td, type = c('ACCTRAN', 'MPR')){
	library (ape )
	library( treedater )
	library( phangorn )
	tr2 <- td ; class(tr2) = 'phylo'
	tr2$tip.label <- unname( tr2$tip.label ) 
	
	
	if (! ( 'sts' %in% names(tr2))){
		sts = as.numeric( sapply( strsplit( tr2$tip.label, '\\|'), function(x) tail(x, 2)[1] ) ) 
		tr2$sts = sts 
	}
	if (! ('Ti' %in% names(tr2) ) ){
		nel = node.depth.edgelength( tr2 ) 
		nel <- nel - max(nel) + max( tr2$sts )
		tr2$Ti <-  nel [ (Ntip(tr2)+1) : (Ntip(tr2)+Nnode(tr2) ) ]
	}
	
	region <- sapply( strsplit( tr2$tip.label, '_|\\|'), function(x) tail(x,1))
	names(region) <- tr2$tip.label
	
	region_levels <- c('Il', 'exog')
	region_pd <- as.phyDat(region, type ='USER' , levels = region_levels )
	
	ancestral.pars(tr2, region_pd , type = type[1] ) -> region_ap
	#plotAnc(tr2,  region_ap, cex = 0) 

	nodeStates = ns <- t( sapply(  subset(region_ap, select = 1, site.pattern = TRUE), as.numeric ) )
	iedge = rev( postorder( tr2 ) )
	
	ismrca_region = rep(FALSE, Nnode(tr2) + Ntip( tr2 ))
	rootnode = tr2$edge[iedge[1],1]
	
	region_ancestor = rep( NA, Nnode(tr2) + Ntip( tr2 ))
	
	for ( k in 1:length( iedge )){ #preorder traversal 
		i <- iedge[k] 
		uv = tr2$edge[i, ]
		u <- uv[1]
		v = uv[2]
		isexog_u <- nodeStates[u,2]==1
		isexog_v <- nodeStates[v,2]==1
		isregion_u <- ( nodeStates[u,1]==1)
		isregion_v <- ( nodeStates[v,1]==1)
		
		if ( !is.na( region_ancestor[u] ))
		{
			region_ancestor[v] <- region_ancestor[u] 
		} else if ( isregion_v ){
			region_ancestor[v] <- v
		}
		
	}
	region_ancestors = unique( na.omit( region_ancestor ))
	region_ancestors <- setdiff( region_ancestors, 1:Ntip(tr2))
#~ browser()
	region_clades = lapply( region_ancestors, function(u){
		tr= extract.clade(tr2, u ) 
		drop.tip( tr, tr$tip.label[ grepl(tr$tip.label, patt = 'exog$') ] )
	})
	
	Ti = c( tr2$sts, tr2$Ti )
	
	list( tres = region_clades, ancestors = region_ancestors, tmrcas = Ti[ region_ancestors] )
}

#' Extact all clades in given tree descended from nodes with high confidece of being within region (parsimony)
#' 
#' Can specify parsimony method (acctrans by default)
#'
#' @export
region_clades2 <- function(td, type = c('ACCTRAN', 'MPR')){
	library (ape )
	library( treedater )
	library( phangorn )
	tr2 <- td ; class(tr2) = 'phylo'
	tr2$tip.label <- unname( tr2$tip.label ) 
	
	
	if (! ( 'sts' %in% names(tr2))){
		sts = as.numeric( sapply( strsplit( tr2$tip.label, '\\|'), function(x) tail(x, 2)[1] ) ) 
		tr2$sts = sts 
	}
	if (! ('Ti' %in% names(tr2) ) ){
		nel = node.depth.edgelength( tr2 ) 
		nel <- nel - max(nel) + max( tr2$sts )
		tr2$Ti <-  nel [ (Ntip(tr2)+1) : (Ntip(tr2)+Nnode(tr2) ) ]
	}
	
	region <- sapply( strsplit( tr2$tip.label, '_|\\|'), function(x) tail(x,1))
	names(region) <- tr2$tip.label
	
	region_levels <- c('Il', 'exog')
	region_pd <- as.phyDat(region, type ='USER' , levels = region_levels )
	
	ancestral.pars(tr2, region_pd , type = type[1] ) -> region_ap
	#plotAnc(tr2,  region_ap, cex = 0) 

	nodeStates = ns <- t( sapply(  subset(region_ap, select = 1, site.pattern = TRUE), as.numeric ) )
	iedge =  postorder( tr2 )
	
	ismrca_region = rep(FALSE, Nnode(tr2) + Ntip( tr2 ))
	rootnode = tr2$edge[tail(iedge,1),1]
	
	#region_ancestor 
	clade_mrca = rep( FALSE, Nnode(tr2) + Ntip( tr2 ))
	
	if ( nodeStates[rootnode,1]==1 )
		clade_mrca[rootnode] <- TRUE 
	
	#postorder traversal 
	for ( k in 1:length( iedge )){ 
		i <- iedge[k] 
		uv = tr2$edge[i, ]
		u <- uv[1]
		v = uv[2]
		isexog_u <- nodeStates[u,2]==1
		isexog_v <- nodeStates[v,2]==1
		isregion_u <- ( nodeStates[u,1]==1)
		isregion_v <- ( nodeStates[v,1]==1)
		
		if ( (!isregion_u) & isregion_v & (v>Ntip(tr2)))
			clade_mrca[v] <- TRUE 
		
	}
	clade_mrcas <- which( clade_mrca ) 
	region_clades = lapply( clade_mrcas, function(u){
		tr= extract.clade(tr2, u ) 
		drop.tip( tr, tr$tip.label[ grepl(tr$tip.label, patt = 'exog$') ] )
	})
	# deduplicate overlapping clades 
	region_clades2 <- list() 
	for ( k in 1:length(region_clades) ){
		tr <- region_clades[[k]]
		for (l in (1:length(region_clades))[-k] ){
			tr1 <- region_clades[[l]] 
			if ( Ntip ( tr1 ) < Ntip( region_clades[[k]] )){
				todrop <- intersect( tr$tip.label, tr1$tip.label )
				tr <- drop.tip( tr, todrop )
			}
		}
		region_clades2[[k]] <- tr 
	}
	
	Ti = c( tr2$sts, tr2$Ti )
	list( tres = region_clades2, ancestors = clade_mrcas, tmrcas = Ti[ clade_mrcas ] )
}




#' @export
.compute_timports <- function(td){
	library (ape )
	library( treedater )
	library( phangorn )
	tr2 <- td ; class(tr2) = 'phylo'
	tr2$tip.label <- unname( tr2$tip.label ) 
	
	region <- sapply( strsplit( tr2$tip.label, '_' ), function(x) tail(x,1))
	names(region) <- tr2$tip.label
	isregion <- setNames( grepl( tr2$tip.label, patt = 'Il$' ) ,  tr2$tip.label )
	
	region_levels <- c('Il', 'exog')
	region_pd <- as.phyDat(region, type ='USER' , levels = region_levels )
	
	ancestral.pars(tr2, region_pd ) -> region_ap
	#plotAnc(tr2,  region_ap, cex = 0) 

	nodeStates = ns <- t( sapply(  subset(region_ap, select = 1, site.pattern = TRUE), as.numeric ) )
	iedge = rev( postorder( tr2 ) )
	exog_ancestor = rep( NA, Nnode(tr2) + Ntip( tr2 ))
	ismrca_monoregion = rep(FALSE, Nnode(tr2) + Ntip( tr2 ))
	rootnode = tr2$edge[iedge[1],1]
	exog_ancestor[ rootnode ] <- rootnode 
	for ( k in 1:length( iedge )){
		i <- iedge[k] 
		uv = tr2$edge[i, ]
		u <- uv[1]
		v = uv[2]
		isexog_u <- nodeStates[u,2]==1
		isexog_v <- nodeStates[v,2]==1
		isregion_u <- ( nodeStates[u,1]==1)
		isregion_v <- ( nodeStates[v,1]==1)
			
		if ( isexog_u ){
			exog_ancestor[u] = exog_ancestor[v] <- u 
		}else {
			exog_ancestor[v] <- exog_ancestor[u]
		}
		
		if (  isregion_v & (!isregion_u) ){
			ismrca_monoregion[v] <- TRUE
		}
	}
	Ti = c( tr2$sts, tr2$Ti )
	tops = Ti [ ismrca_monoregion ]
	bottoms = Ti[ exog_ancestor[ ismrca_monoregion ] ]
	mids = (tops + bottoms)/2
	w <- sapply( which( ismrca_monoregion ), function(u){
		if ( u <= Ntip (tr2))
			return(1)
		tr3 = extract.clade( tr2, u )
		sum( isregion[ tr3$tip.label ]  )
	})
	pldf = data.frame( yplus = tops ,  yminus = bottoms, y = mids , w = w)
	pldf 
}

.compute_timports_pb <- function( td, numpb , ncpu = 6, overrideSearchRoot=FALSE, excludeBefore = 2020)
{
	library( treedater )
	if ( numpb > 1 ){
		pb = parboot( td, nreps = numpb , overrideTempConstraint = FALSE, overrideSearchRoot=overrideSearchRoot, ncpu = ncpu )
		tis = lapply( pb$trees, .compute_timports )
	} else{
		tis = list( .compute_timports( td ) )
		if ( nrow ( tis[[1]] ) == 0 )
			return( list( pldf = tis[[1]], tis = tis ) )
	}
	
	yp = do.call( cbind, lapply( tis, '[[', 'yplus' ))
	ym = do.call( cbind, lapply( tis, '[[', 'yminus' ))
	mi = do.call( cbind, lapply( tis, '[[', 'y' ))
	w = do.call( cbind, lapply( tis, '[[', 'w' ))
	
	pldf = data.frame( y = apply( mi, MAR=1, FUN = median )
	 , yminus = apply( ym, MAR=1, FUN = function(x) unname( quantile(x, .025)) )
	 , yplus = apply( yp, MAR=1, FUN = function(x) unname( quantile(x, .975)) )
	 , w = apply( w, MAR=1, FUN= median )
	)
	pldf <- pldf [ pldf$y > excludeBefore , ]
	list(pldf =pldf , tis = tis )
}
#~ pldf = .compute_timports_pb( tds[[1]], numpb=100 )


#' Plot a errbar of importation times given a data frame from compute_timports
#'
#' @param pldf A data frame made by compute_timports 
#' @param showDates Will plot dates instead of numeric if TRUE 
#' @param ... passed to errbar 
#' @export
plot_importations <- function(pldf, showDates=TRUE , ...)
{
	library( lubridate )
	pldf <- pldf[ order( pldf$y ) , ]
	pldf$x = 1:nrow(pldf )
	if (showDates )
		with( pldf, Hmisc::errbar ( x, date_decimal( y ), date_decimal( yplus) , date_decimal( yminus)  , xlab = 'Cluster index', ylab = 'Cluster origin time', ... ))
	else
		with( pldf, Hmisc::errbar ( x, y, yplus, yminus , xlab = 'Cluster index', ylab = 'Cluster origin time' , ...))
	invisible(pldf )
}
#~ plot.importations( pldf , bty = 'n')

#' Plot a histogram of importation times given a data frame from compute_timports
#'
#' @param pldf A data frame made by compute_timports 
#' @param showDates Will plot dates instead of numeric if TRUE 
#' @param breaks 
#' @param xlab 
#' @export
plot_importationsHist <- function(pldf, showDates = TRUE, breaks=10 , xlab='Cluster origin',main = '', ... ){
	if (showDates)
		hist( as.Date( date_decimal( pldf$y )  ) , main = main, breaks = breaks , xlab = xlab, ylab = '', ...)
	else
		hist( ( pldf$y )  , main = '', breaks = breaks, xlab = xlab, ylab ='', ...)
	invisible( pldf )
}
#~ X11(); plot.importationsHist( pldf , xlab ='Netherlands cluster time of origin')


#' Computes the time of importation events using a one or more treedater time trees and maximum parsimony
#'
#' Time of importations is taken as the midpoing on branches separating an exogenous parent and a regional daughter; the states of nodes are determined by parsimony 
#'
#" @param tds A list of treedater trees 
#' @param numpb A number of parametric bootstrap replicates for each treedater tree. 
#' @param ncpu Number of cpus to use
#' @param excludeBefore Removes any dates before this point (some outlier sequences merge before then (guangdong))
#' @return A data frame that can be used for plotting etc. *weighted_seed_date* shows the mean and CI of the seed date weighted by number of region tips descended from that node. This gives much less weight to later introductions with few descendents. 
#' @examples
#' \dontrun{ 
#' PATT = 'Netherlands'
#' o1 = region_time_stratified_sample( PATT , n=50, path_to_align=algn, D = D, path_to_save='temp.fasta', nregion = Inf, time_stratify_region=TRUE)
#' prep_tip_labels_phydyn( 'temp.fasta'
#'  , regexprs = c(  PATT, PATT  ) 
#'  , invert_regexpr = c( FALSE, TRUE )
#'  , demes = c( 'Il'  , 'exog'  )
#'  , path_to_save = 'temp.fasta'
#' )
#' tds = make_starting_trees ( 'temp.fasta', ntres = 20, ncpu = 6, plotout=NULL )
#' pldf <- compute_timports( tds )
#' plot.importations( pldf )
#' X11(); plot.importationsHist( pldf )
#' }
#' @export 
compute_timports <- function(tds, numpb = 100, ncpu = 6, excludeBefore = 2020 ){
	if (class( tds ) == 'treedater' )
		tds <- list( tds )
	ntds <- length( tds )
	results <- lapply( tds, function(td) .compute_timports_pb( td, numpb, ncpu , excludeBefore = excludeBefore) )
	pldfs <- lapply( results, '[[', 'pldf' )
	tis <- lapply( results, '[[', 'tis' )
	pldf <- do.call(rbind, pldfs )
	ws = apply( pldf[, c('y', 'yminus','yplus')], MAR=2, FUN=function(x) weighted.mean(x, pldf$w ))
	list( pldf = pldf, tis = tis , weighted_seed_date = lubridate::date_decimal(ws))
}

