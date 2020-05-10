
#' Compute skygrowth Ne(t) and R(t) estimates for a given list of time scaled trees
#' 
#' If tstart is provided than you can supply a list of ape::phylo trees instead of treedater trees 
#' 
#' @param tds A list or multiPhylo, containing ape::phylo or treedater trees
#' @param tstart A numeric time to consider as the beginning of the epidemic in region; for best results, select a time when exponential growth in region is well underway. Can be NULL in which case the sarscov2::timports function is used to guess an appropriate date
#' @param tau0 Precision parameter for skygrowth. Default value corresponds to 1 per cent change in growth per week 
#' @param see skygrowth, if omitted will guess a good value
#' @param  numpb Number of parametric bootstrap to use if using timports
#' @param ncpu Number CPUs to use 
#' @param gamma1 Death/recovery rates used for translating growth rates into R(t) 
#' @param ... Additional arguments passed to skygrowth (eg mhsteps)
#' @return A list with dataframes for Ne, growth, and R(t)
#' @examples
#' \dontrun{
#' library( ape ) ; library( sarscov2 )
#' tres = read.tree( 'startTrees.nwk' ) 
#' o = skygrowth0 ( tres, start = decimal_date( as.Date('2020-02-15') ) )
#' }
#' @export 
skygrowth0 <- function( tds , tstart = decimal_date( as.Date('2020-02-15') ), tau0 = 1/365 / 36.5^2, res = NULL,  numpb = 10, ncpu = 6, gamma1 = (-log(.5)) * 365 / 6.5, ...){
	stopifnot( require( skygrowth )  )
	library( lubridate )
	
	if ( is.null( tstart )){
		tis = compute_timports( tds , numpb = numpb , ncpu = ncpu )
		timport_quantile = .25
		tstart <- quantile( tis[, 'y' ] , timport_quantile )
	}
	if ( !is.numeric( tstart ))
		tstart <- decimal_date( tstart )
	
	tres = lapply( tds, function(x) { class(x)='phylo'; x } )
	
	tr <- tres[[1]] 
	sts = sapply( strsplit( tr$tip, '\\|') , function(x) as.numeric( tail(x,2)[1] ) )
	msts = max( sts ) 
	
	.skygrowth <- function( tr ) {
		tr = drop.tip( tr, tr$tip[ grepl( tr$tip, patt='_exog' ) ]  ) 
		.res <- res
		if ( is.null( res )) # guess a good res based on sample size 
			.res <-  floor( 10 + 30*max(0, min( (Ntip(tr) - 50)/(300-50), 1 ) ) )
		
		a = capture.output( { 
				sg = skygrowth.mcmc( tr , res = .res , tau0 = tau0, ... )# , mhsteps = 1e6)
		})
		sg
	}
	
	
	sgs = parallel::mclapply( tres, .skygrowth, mc.cores = ncpu  )
	taxis = seq( tstart, msts , length =  floor((msts - tstart)*365) )
	dates = date_decimal( taxis )
	
	# Ne
	Ne <- do.call( cbind, lapply( sgs , function(sg){
		apply( sg$ne, MAR=1, FUN= function(x) approx( sg$time + msts, x , rule=2, xout = taxis)$y )
	}))
	NeCI = t( apply( Ne, MAR=1, FUN = function(x) quantile(x, c(.025, .5, .975 )) ) )
	colnames(NeCI)= c( 'pc2.5', 'pc50', 'pc97.5' ) 
	NeCI = data.frame( time = dates, NeCI )
	
	# growth rate 
	GR <- do.call( cbind, lapply( sgs , function(sg){
		apply( sg$growthrate, MAR=1, FUN= function(x) approx( sg$time+msts, x , rule=2, xout = taxis)$y )
	}))
	GRCI = t( apply( GR, MAR=1, FUN = function(x) quantile(x, c(.025, .5, .975 )) ) )
	colnames(GRCI)= c( 'pc2.5', 'pc50', 'pc97.5' ) 
	GRCI = data.frame( time = dates, GRCI )
	
	Tg = 1/gamma1
	# growthrate = (R0-1)/Tg
	Rmat <-  GR*Tg + 1
	Rmat[ Rmat < 0 ] <- 0#NA 
	RCI =  t( apply( Rmat, MAR=1, FUN = function(x) quantile(na.omit(x), c(.025, .5, .975 )) ) )
	colnames(RCI) = c( 'pc2.5', 'pc50', 'pc97.5' ) 
	RCI = data.frame( time = dates, RCI )
	
	structure(
	list(
	  taxis = taxis
	  , date = dates
	  , Ne = NeCI
	  , growth= GRCI
	  , R =RCI
	  , GRmat = GR 
	  , Rmat = Rmat 
	)
	, class = 'sarscov2skygrowth' 
	)
}

#todo 
.phylodyn <- function( tr ) {
	tr = drop.tip( tr, tr$tip[ grepl( tr$tip, patt='_exog' ) ]  ) 
	.res <- res
	if ( is.null( res )) # guess a good res based on sample size 
		.res <- 10 + 30*max(0, min( (Ntip(tr) - 50)/(300-50), 1 ) )
	
	phylodyn::BNPR( tr, lengthout = .res )
}

#~ library( ape ) 
#~ tres = read.tree( 'startTrees.nwk' ) 
#~ o = skygrowth0 ( tres[1:5], start = decimal_date( as.Date('2020-02-15') ) )

#' @export
plot.sarscov2skygrowth <- function(x, ...)
{
	stopifnot( inherits( x, 'sarscov2skygrowth' ))
	y = x$Ne 
	taxis = as.Date( y$time )
	plot( taxis, y$pc50 , col='black', lty = 1, lwd = 2, ylim = c( min( na.omit( y$pc2.5)) , max( na.omit(y$pc97.5)) ), type = 'l', xlab='', ylab='Ne(t)', ...)
	lines( taxis, y$pc2.5, col='black' , lty = 3 )
	lines( taxis,y$pc97.5 , col='black', lty = 3)
}

#' @export
plotR.sarscov2skygrowth <- function(x, ...)
{
	stopifnot( inherits( x, 'sarscov2skygrowth' ))
	y = x$R 
	taxis = as.Date( y$time )
	plot( taxis, y$pc50 , col='black', lty = 1, lwd = 2, ylim = c( min( na.omit( y$pc2.5)) , max( na.omit(y$pc97.5)) ), type = 'l', ylab = 'R(t)', xlab = '',  ...)
	lines( taxis, y$pc2.5, col='black' , lty = 3 )
	lines( taxis,y$pc97.5 , col='black', lty = 3)
	abline( h = 1, col = 'red')
}
#~ plot.sarscov2skygrowthR( sg0 )

#' @export 
ggNe_sarscov2skygrowth <- function(x,  date_limits = c( as.Date( '2020-02-01'), NA ) ,... )
{
	require(lubridate)
	require(ggplot2)
	stopifnot( inherits( x, 'sarscov2skygrowth' ))
	y = x$Ne 
	taxis = as.Date( y$time )
	
	if ( is.na( date_limits[2]) )
		date_limits[2] <- as.Date( date_decimal( max(taxis)  ) )
	#qs <- c( .5, .025, .975 )
	
	pldf <- data.frame( Date = taxis , reported=FALSE )
	pldf$Ne = y$pc50
	pldf$`2.5%` = y$pc2.5
	pldf$`97.5%` = y$pc97.5
	
	pldf <- pldf[ with( pldf, Date > date_limits[1] & Date <= date_limits[2] ) , ]
	pl = ggplot( pldf ) + 
	  geom_path( aes(x = Date, y = Ne ), lwd=1.25) + 
	  geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`) , alpha = .25 ) 
	
	pl <- pl + theme_minimal()  + xlab('') + 
	 ylab ('Effective population size ' ) 
	pl
}
#~ ggplot.sarscov2skygrowth( sg0 )

#' @export 
ggR_sarscov2skygrowth <- function(x,  date_limits = c( as.Date( '2020-02-01'), NA ) ,... )
{
	require(ggplot2)
	require(lubridate)
	stopifnot( inherits( x, 'sarscov2skygrowth' ))
	y = x$R
	taxis = as.Date( y$time )
	
	if ( is.na( date_limits[2]) )
		date_limits[2] <- as.Date( date_decimal( max(taxis)  ) )
	#qs <- c( .5, .025, .975 )
	
	pldf <- data.frame( Date = taxis , reported=FALSE )
	pldf$R = y$pc50
	pldf$`2.5%` = y$pc2.5
	pldf$`97.5%` = y$pc97.5
	
	pldf <- pldf[ with( pldf, Date > date_limits[1] & Date <= date_limits[2] ) , ]
	pl = ggplot( pldf ) + 
	  geom_path( aes(x = Date, y = R ), lwd=1.25) + 
	  geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`) , alpha = .25 ) 
	
	pl <- pl + geom_hline( aes(yintercept = 1 ), colour = 'red' )
	pl <- pl + theme_minimal()  + xlab('') + 
	 ylab ('Effective reproduction number' ) 
	pl
}
#~ ggRplot.sarscov2skygrowth( sg0 )

#~ y = 4.125 * sg$ne_ci[,2] / (6.5/365)
#~ taxis2 = seq( -.2, 0, length = 1e3 )
#~ Yt = approx( sg$time, y , rule=2, xout = taxis2 )$y
#~ gamma1 = (-log(.5)) * 365 / 6.5
#~ p = skygrowth::R.plot( sg , gamma = gamma1 )
#~ Rt = pmax(0, approx( p$data$t, p$data$med , rule = 2, xout = taxis2 )$y )
#~ dtaxis2 <- diff( taxis2 )[1] 
#~ sum( Yt * (Rt * gamma1) * dtaxis2 ) 
