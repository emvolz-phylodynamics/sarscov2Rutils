

#' Combine BEAST log files after removing burnin
#'
#" Additionally, if some BEAST runs converged to a different tree space with different posterior, will exclude log files from runs with significantly lower posterior (analysis of variance)
#'
#" @param logfns Vector of log files 
#' @param bunrProportion Will remove this proportion from the beginning of logs
#' @param ofn If not NULL will save the combined log to this file name 
#' @return Data frame with combined logs 
#' @export
combine_logs <- function(logfns, burnProportion = .5 , ofn = NULL ){
	
	Xs <- lapply( logfns, function(fn){
		d = read.table( fn, header=TRUE, stringsAsFactors=FALSE) 
		i <- floor( burnProportion * nrow(d))
		tail( d, nrow(d) - i )
	})

	if ( length( Xs ) > 1 ){
		medlogpos <- sapply( Xs, function(X) median (X$posterior ))
		Xs <- Xs[ order( medlogpos, decreasing=TRUE ) ]
		
		xdf = data.frame(logpo = NULL, logfn = NULL ) 
		k <- 0
		for (X in Xs ){
			k <- k + 1
			xdf <- rbind ( xdf ,  data.frame( logpo = X$posterior , logfn = logfns[k] ) )
		}

		m <- aov(  logpo ~ logfn , data = xdf )
		a <- TukeyHSD( m )
		ps <- sapply( 1:(length(Xs)-1), function(k) {
			a$logfn[ k, 4 ]
		})
		keep <- c( 1, 1 + which ( ps > .05 ) )
		Xs1 <- Xs[ keep ]
		
		cat( "These logs were retained:\n" )
		print( names(Xs1) )

		X <- do.call( rbind, Xs1)
	} else{
		X <- Xs[[1]]
	}
	
	if ( !is.null( ofn ))
		saveRDS(X , file = ofn )
	
	X
}

#' Combine PhyDyn trajectory files after removing burnin
#'
#" @param trajfns Vector of trajectory files. Don't pass low quality BEAST fits; see the combine_logs function for filtering by quality 
#' @param bunrProportion Will remove this proportion from the beginning of logs
#' @param ntraj This integer number of trajectories will be sampled from each trajectory file 
#' @param ofn If not NULL will save the combined trajectory log to this file name 
#' @return Data frame with combined trajectories 
#' @export
combine_traj <- function(trajfns, burnProportion = .50, ntraj = 100, ofn = NULL )
{
	cat( 'NOTE: only pass trajectory fns that pass the quality filter in the combine_logs function\n' )
	Xs <- lapply( trajfns, function(fn){
		d = read.table( fn, header=TRUE, stringsAsFactors=FALSE)
		ids = sort (unique( d$Sample ))
		keep <- sample( tail(ids, floor(length(ids)-burnProportion*length(ids)))
		 , size = ntraj, replace=FALSE)
		d = d[ d$Sample %in% keep , ]
		d$Sample <- paste0( as.character(d$Sample) , '.', fn)
		d
	})
	if ( length( Xs ) == 1)
		X <- Xs[[1]]
	else 
		X <- do.call( rbind, Xs )
	if ( !is.null ( ofn ))
		saveRDS( X, file = ofn )
	
	X
}



#' Compute R0, growth rate and doubling time for the SEIJR.0.0 model 
#'
#' Also prints to the screen a markdown table with the results. This can be copied into reports. 
#'
#" @param X a data frame with the posterior trace, can be produced by 'combine_logs' function
#' @param gamma1 Rate of recovery (per capita per year)
#' @param tau Transm rate increase in J 
#' @param p_h Proportion in J 
#' @return Data frame with tabulated results and CI 
#' @export 
SEIJR_reproduction_number <- function( X, gamma0 = 89.042, gamma1 = 96, tau = 74, p_h = 0.20 , precision =3 ) {
	cat( 'Double check that you have provided the correct gamma0 and gamma1 parameters\n' )
	
	Rs = ((1-p_h)*X$seir.b/gamma1 + tau*p_h*X$seir.b/gamma1) 
	qR = signif( quantile( Rs, c(.5, .025, .975 )), precision )
	
	# growth rates 
	beta = (1-p_h)*X$seir.b + p_h*tau*X$seir.b
	r = (-(gamma0 + gamma1) + sqrt( (gamma0-gamma1)^2 + 4*gamma0*beta )) / 2
	qr =  signif( quantile( r/365, c( .5, .025, .975 )), precision ) # growth rate per day 
	
	# double times 
	dbl = 1 / (r /365)
	qdbl  = signif( quantile( dbl, c( .5, .025, .975 )), precision ) # days 
	
	#cat ( 'R0\n')
	#print( qR )
	
	O = data.frame( 
	  `Reproduction number` = qR
	  , `Growth rate (per day)` = qr
	  , `Doubling time (days)` = qdbl
	 )
	print( knitr::kable(O) )
	O 
}
#~ SEIJR_reproduction_number( d) 


#' Plot the cumulative infections through time from a SEIJR trajectory sample
#' 
#" Also computes CIs of various dynamic variables 
#'
#" @param trajdf Either a dataframe or a path to rds containing a data frame with a posterior sample of trajectories (see combine_traj)
#' @param case_data An optional dataframe containing reported/confirmed cases to be plotted alongside estimates. *Must* contain columns 'Date' and 'Cumulative'. Ensure Date is not a factor or character (see as.Date )
#' @param date_limits  a 2-vector containing bounds for the plotting window. If the upper bound is missing, will use the maximum time in the trajectories
#' @param path_to_save Will save a png here 
#' @return a list with derived outputs from the trajectories. The first element is a ggplot object if you want to further customize the figure 
#' @export 
SEIJR_plot_size <- function(trajdf
  , case_data = NULL
  , date_limits = c( as.Date( '2020-02-01'), NA ) 
  , path_to_save='size.png'
  , ...
) {
	library( ggplot2 ) 
	library( lubridate )
	
	if (!is.data.frame(trajdf)){
		# assume this is path to a rds 
		readRDS(trajdf) -> trajdf 
	}
	
	#~ 	r <- read.csv( '../seir21.0/weifangReported.csv', header=TRUE , stringsAsFactors=FALSE)
	#~ 	r$Date <- as.Date( r$Date )
	#~ 	r$reported = TRUE
	#~ 	r$`Cumulative confirmed` = r$Cumulative.confirmed.cases

	dfs <- split( trajdf, trajdf$Sample )
	taxis <- dfs[[1]]$t 
	
	if ( is.na( date_limits[2]) )
		date_limits[2] <- as.Date( date_decimal( max(taxis)  ) )

	qs <- c( .5, .025, .975 )

	# infectious
	Il = do.call( cbind, lapply( dfs, '[[', 'Il' ))
	Ih = do.call( cbind, lapply( dfs, '[[', 'Ih' ))
	I = Il + Ih 
	t(apply( I, MAR=1, FUN= function(x) quantile(x, qs ))) -> Ici 

	# cases 
	cases <- do.call( cbind, lapply( dfs, '[[', 'infections' ))
	t(apply( cases, MAR=1, FUN=function(x) quantile(x,qs))) -> casesci 

	#exog 
	exog <- do.call( cbind, lapply( dfs, '[[', 'exog' ))
	t(apply( exog, MAR=1, FUN=function(x) quantile(x, qs )	)) -> exogci 

	#E
	E <- do.call( cbind, lapply( dfs, '[[', 'E' ))
	t(apply( E, MAR=1, FUN=function(x) quantile(x, qs )	)) -> Eci 


	#~ ------------
	pldf <- data.frame( Date = ( date_decimal( taxis ) ) , reported=FALSE )
	pldf$`Cumulative infections` = casesci[,1]
	pldf$`2.5%` = casesci[,2]
	pldf$`97.5%` = casesci[,3] 
	
	if ( !is.null( case_data )){
		case_data$reported = TRUE
		pldf <- merge( pldf, case_data , all = TRUE ) 
	}
	
	pldf <- pldf[ with( pldf, Date > date_limits[1] & Date <= date_limits[2] ) , ]
	pl = ggplot( pldf ) + 
	  geom_path( aes(x = Date, y = `Cumulative infections` , group = !reported), lwd=1.25) + 
	  geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`, group = !reported) , alpha = .25 ) 
	
	if ( !is.null(case_data) ) {
		if( !(class(case_data$Date)=='Date') ){
			stop('case_data Date variable must be class *Date*, not character, integer, or POSIXct. ') 
		}
		pl <- pl + geom_point( aes( x = Date, y = Cumulative ) ) 
	}
	pl <- pl + theme_minimal()  + xlab('') + 
	 ylab ('Cumulative estimated infections (ribbon) 
	 Cumulative confirmed (points)' )  #+  scale_y_log10()
	
	if (!is.null(path_to_save))
		ggsave(pl, file = path_to_save)
	
	list(
	 pl
	  , taxis = taxis 
	  , Il = Il
	  , Ih = Ih
	  , E = E 
	  , I = I 
	  , cases = cases 
	  , exog = exog 
	  , pldf =pldf 
	  , case_data = case_data 
	)
}





# debug 
if (FALSE)
{
	X = combine_logs( logfns = 'nl5.log', burnProportion = .5 , ofn = 'tmp.rds' )
	Y = combine_traj( trajfns = 'seir.nl5.traj', burnProportion = .50, ntraj = 200, ofn = 'tmptraj.rds' )
	Z = SEIJR_reproduction_number( X, gamma1 = 96, tau = 74, p_h = 0.20  ) 
	
	O = SEIJR_plot_size(trajdf = 'tmptraj.rds'
	  , case_data = NULL
	  , date_limits = c( as.Date( '2020-02-01'), NA ) 
	  , path_to_save='size.png')

}
