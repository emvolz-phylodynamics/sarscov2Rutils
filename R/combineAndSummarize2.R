
.comp_beta2 <- function( taxis , theta,  tknots)
{
	numchange = length(tknots) - 2
	dlogbetanames = paste0('seir.dlogbeta.t',1:numchange )
	betas = theta$seir.b * c(1,  exp( cumsum(unlist(theta[dlogbetanames])  )) )
	approx( tknots[1:(numchange+1)], betas, rule = 2 , xout = taxis )$y
}


#' Plot the cumulative infections through time from a SEIJR trajectory sample
#' 
#" Also computes CIs of various dynamic variables 
#'
#" @param trajdf Either a dataframe or a path to rds containing a data frame with a posterior sample of trajectories (see combine_traj)
#' @param case_data An optional dataframe containing reported/confirmed cases to be plotted alongside estimates. *Must* contain columns 'Date' and 'Cumulative'. Ensure Date is not a factor or character (see as.Date )
#' @param date_limits  a 2-vector containing bounds for the plotting window. If the upper bound is missing, will use the maximum time in the trajectories
#' @param path_to_save Will save a png here 
#' @param log_y_axis 
#' @param k Overdispersion parameter (default 16 per cent, Lloyd Smith 2005)
#' @return a list with derived outputs from the trajectories. The first element is a ggplot object if you want to further customize the figure 
#' @export 
gmrf_exogsir_plot_size <- function(trajdf
  , case_data = NULL
  , date_limits = c( as.Date( '2020-02-01'), NA ) 
  , path_to_save='size.png'
  , log_y_axis = F
  , k = 0.16 # overdispersion parameter
  , ...
) {
	inf_adj <- (2*k + 1 ) / (2*k)
	
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
	I = Il = inf_adj * do.call( cbind, lapply( dfs, '[[', 'Il' ))
	t(apply( I, MAR=1, FUN= function(x) quantile(x, qs ))) -> Ici 

	# cases 
	cases <- inf_adj * do.call( cbind, lapply( dfs, '[[', 'infections' ))
	t(apply( cases, MAR=1, FUN=function(x) quantile(x,qs))) -> casesci 

	#exog 
	exog <- do.call( cbind, lapply( dfs, '[[', 'exog' ))
	t(apply( exog, MAR=1, FUN=function(x) quantile(x, qs )	)) -> exogci 

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
	
	if(log_y_axis == T)
	  pl <- pl +  scale_y_log10()
	
	if (!is.null(path_to_save))
		ggsave(pl, file = path_to_save)
	
	list(
	 pl = pl
	  , taxis = taxis 
	  , I = I 
	  , cases = cases 
	  , exog = exog 
	  , pldf =pldf 
	  , case_data = case_data 
	)
}
#~ p0 = gmrf_exogsir_plot_size( 'traj.rds' )



#' Plot the daily new infections through time from a GMRF SEIJR trajectory sample
#' 
#" Also computes CIs of various dynamic variables 
#'
#" @param trajdf Either a dataframe or a path to rds containing a data frame with a posterior sample of trajectories (see combine_traj)
#' @param logdf Either a dataframe or a path to rds containing a data frame with posterior logs
#' @param case_data An optional dataframe containing reported/confirmed cases to be plotted alongside estimates. *Must* contain columns 'Date' and 'Confirmed'. Ensure Date is not a factor or character (see as.Date )
#' @param date_limits a 2-vector containing bounds for the plotting window. If the upper bound is missing, will use the maximum time in the trajectories
#' @param path_to_save Will save a png here 
#' @param k Overdispersion parameter (default 16 per cent, Lloyd Smith 2005)
#' @return a list with derived outputs from the trajectories. The first element is a ggplot object if you want to further customize the figure 
#' @export 
gmrf_exogsir_plot_daily_inf <- function(trajdf
  , logdf
  , case_data = NULL
  , date_limits = c( as.Date( '2020-02-01'), NA ) 
  , path_to_save='daily.png'
  , log_y_axis = F
  , k = 0.16
  , ...
) {
	
	inf_adj <- (2*k + 1) / (2*k)
	
	library( ggplot2 ) 
	library( lubridate )
	
	if (!is.data.frame(trajdf)){
		# assume this is path to a rds 
		readRDS(trajdf) -> trajdf 
	}
	
	if (!is.data.frame( logdf ))
		X <- readRDS(logdf ) else X <- logdf
		
	dfs <- split( trajdf, trajdf$Sample )
	taxis <- dfs[[1]]$t 
	
	pnames = colnames(X)
	dlbnames = pnames [ grepl(pnames, pattern='^seir\\.dlogbeta\\.t[0-9]+$') ] 
	numchange = length( dlbnames )
	tend = max(taxis)
	tknots = seq( 2020.038, tend, length.out = numchange+2 )
	
	s <- sapply( dfs, function(df) df$Sample[1] )
	X1 <- X[match( s, X$Sample ), ]
	
	if ( is.na( date_limits[2]) )
		date_limits[2] <- as.Date( date_decimal( max(taxis)  ) )

	qs <- c( .5, .025, .975 )

	# infectious
	Il = do.call( cbind, lapply( dfs, '[[', 'Il' ))
	I = Il 
	t(apply( I, MAR=1, FUN= function(x) quantile(x, qs ))) -> Ici 
	
	# daily new inf 
	Y = lapply( 1:length(dfs) , function(k) {
		x = X1[k, ]
		theta = ( X1[k, ] )
		Il = dfs[[k]]$Il
		b = .comp_beta2 (taxis, theta, tknots)
		b *(1/365)*Il*inf_adj
	})
	Y = do.call( cbind, Y )
	Yci = t(apply( Y, MAR=1, FUN= function(x) quantile(na.omit(x), qs ))) 

	#~ ------------
	pldf <- data.frame( Date = ( date_decimal( taxis ) ) , reported=FALSE )
	pldf$`New infections` = Yci[,1]
	pldf$`2.5%` = Yci[,2]
	pldf$`97.5%` = Yci[,3] 
	
	if ( !is.null( case_data )){
		case_data$reported = TRUE
		pldf <- merge( pldf, case_data , all = TRUE ) 
	}
	
	pldf <- pldf[ with( pldf, Date > date_limits[1] & Date <= date_limits[2] ) , ]
	pl = ggplot( pldf ) + 
	  geom_path( aes(x = Date, y = `New infections` , group = !reported), lwd=1.25) + 
	  geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`, group = !reported) , alpha = .25 ) 
	
	if ( !is.null(case_data) ) {
		if( !(class(case_data$Date)=='Date') ){
			stop('case_data Date variable must be class *Date*, not character, integer, or POSIXct. ') 
		}
		pl <- pl + geom_point( aes( x = Date, y = Confirmed ) ) 
	}
	pl <- pl + theme_minimal()  + xlab('') + 
	 ylab ('Estimated daily new infections (ribbon) 
	 Daily confirmed cases (points)' )  #+  scale_y_log10()
	
	if(log_y_axis == T)
	  pl <- pl +  scale_y_log10()
	
	if (!is.null(path_to_save))
		ggsave(pl, file = path_to_save)
	
	list(
	 pl = pl
	  , taxis = taxis 
	  , Il = Il
	  , I = I
	  , Y = Y 
	  , pldf =pldf 
	  , case_data = case_data 
	)
}
#~ p1 = gmrf_exogsir_plot_daily_inf( 'traj.rds' , 'logs.rds')


#' Plot reproduction number through time from a SEIJR trajectory sample
#'
#" @param trajdf Either a dataframe or a path to rds containing a data frame with a posterior sample of trajectories (see combine_traj)
#' @param logdf Either a dataframe or a path to rds containing a data frame with posterior logs
#' @param date_limits  a 2-vector containing bounds for the plotting window. If the upper bound is missing, will use the maximum time in the trajectories
#' @param gamma1 rate of recovery once infectious
#' @param path_to_save Will save a png here 
#' @param changePoints the time where beta changes
#' @return a list with derived outputs from the trajectories. The first element is a ggplot object if you want to further customize the figure 
#' @export 
gmrf_exogsir_plot_Rt <- function(trajdf
  , logdf
  , gamma1 = 56.15385
  , date_limits = c( as.Date( '2020-02-01'), NA ) 
  , path_to_save='Rt.png'
  , ...
) {
	
	library( ggplot2 ) 
	library( lubridate )
	
	if (!is.data.frame(trajdf)){
		# assume this is path to a rds 
		readRDS(trajdf) -> trajdf 
	}
	
	if (!is.data.frame( logdf ))
		X <- readRDS(logdf ) else X <- logdf
		
	dfs <- split( trajdf, trajdf$Sample )
	taxis <- dfs[[1]]$t 
	
	pnames = colnames(X)
	dlbnames = pnames [ grepl(pnames, pattern='^seir\\.dlogbeta\\.t[0-9]+$') ] 
	numchange = length( dlbnames )
	tend = max(taxis)
	tknots = seq( 2020.038, tend, length.out = numchange+2 )
	
	s <- sapply( dfs, function(df) df$Sample[1] )
	X1 <- X[match( s, X$Sample ), ]
	
	if ( is.na( date_limits[2]) )
		date_limits[2] <- as.Date( date_decimal( max(taxis)  ) )
	
	qs <- c( .5, .025, .975 )
	
	Y = lapply( 1:length(dfs) , function(k) {
		x = X1[k, ]
		theta = ( X1[k, ] )
		b = .comp_beta2 (taxis, theta, tknots)
		b / gamma1
	})
	Y = do.call( cbind, Y )
	Yci = t(apply( Y, MAR=1, FUN= function(x) quantile(na.omit(x), qs ))) 

	#~ ------------
	pldf <- data.frame( Date = ( date_decimal( taxis ) ) , reported=FALSE )
	pldf$`R(t)` = Yci[,1]
	pldf$`2.5%` = Yci[,2]
	pldf$`97.5%` = Yci[,3] 
	
	
	pldf <- pldf[ with( pldf, Date > date_limits[1] & Date <= date_limits[2] ) , ]
	pl = ggplot( pldf ) + 
	  geom_path( aes(x = Date, y = `R(t)` , group = !reported), lwd=1.25) + 
	  geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`, group = !reported) , alpha = .25 ) 
	
	pl <- pl + theme_minimal()  + xlab('') + 
	 ylab ('Eff reproduction number through time' )  #+  scale_y_log10()
	
	if (!is.null(path_to_save))
		ggsave(pl, file = path_to_save)
	
	list(
	 plot = pl
	  , taxis = taxis 
	  , pldf =pldf 
	)
}
#~ trajdf = 'traj.rds' 
#~ logdf = 'logs.rds'
#~ p2 = gmrf_exogsir_plot_Rt( trajdf, logdf )
