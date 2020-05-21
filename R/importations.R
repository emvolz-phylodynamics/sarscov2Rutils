
#~ library( sarscov2)
#~ library( ape ) 
#~ library( treedater )
#~ library( phangorn )
#~ library( ggtree )
#~ library( phangorn )
#~ library( Hmisc)


#' @export
.compute_timports <- function(td){
	library (ape )
	library( treedater )
	library( phangorn )
	tr2 <- td ; class(tr2) = 'phylo'
	tr2$tip.label <- unname( tr2$tip.label ) 
	
	region <- sapply( strsplit( tr2$tip.label, '_' ), function(x) tail(x,1))
	names(region) <- tr2$tip.label
	region_pd <- as.phyDat(region, type ='USER' , levels = unique( region ) )
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
	pldf = data.frame( yplus = tops ,  yminus = bottoms, y = mids )
	pldf 
}

.compute_timports_pb <- function( td, numpb , ncpu = 6, overrideSearchRoot=FALSE, excludeBefore = 2020)
{
	library( treedater )
	pb = parboot( td, nreps = numpb , overrideTempConstraint = FALSE, overrideSearchRoot=overrideSearchRoot, ncpu = ncpu )
	tis = lapply( pb$trees, .compute_timports )
	yp = do.call( cbind, lapply( tis, '[[', 'yplus' ))
	ym = do.call( cbind, lapply( tis, '[[', 'yminus' ))
	mi = do.call( cbind, lapply( tis, '[[', 'y' ))
	
	pldf = data.frame( y = apply( mi, MAR=1, FUN = median )
	 , yminus = apply( ym, MAR=1, FUN = function(x) unname( quantile(x, .025)) )
	 , yplus = apply( yp, MAR=1, FUN = function(x) unname( quantile(x, .975)) )
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
#' @return A data frame that can be used for plotting etc 
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
	list( pldf = pldf, tis = tis )
}


