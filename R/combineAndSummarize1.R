# summary functions for gmrf_seijr model 

invisible('
<definition id="Definition.betalb" spec="phydyn.model.Definition" value="betalb = if ( t !>= SEIJR_START ) then 0
	else if (t  !>= 2020.123 ) then beta1 
	else if (t  !>= 2020.161 ) then beta1 * exp( dlogbeta[0] )
	else if ( t !>= 2020.199 ) then beta1 * exp( dlogbeta[0] + dlogbeta[1] )
	else (beta1 * exp( dlogbeta[0] + dlogbeta[1] + dlogbeta[2]) ) "/>
<definition id="Definition.betaub" spec="phydyn.model.Definition" value="betaub = if ( t  !>= SEIJR_START) then 0
	else if (t  !>= 2020.123 ) then beta1 * exp(  dlogbeta[0] )
	else if (t  !>= 2020.161 ) then beta1 * exp( dlogbeta[0] + dlogbeta[1] )
	else if ( t !>= 2020.199 ) then beta1 * exp( dlogbeta[0] + dlogbeta[1] + dlogbeta[2] )
	else (beta1 * exp( dlogbeta[0] + dlogbeta[1] + dlogbeta[2] )) "/>
<definition id="Definition.tlb" spec="phydyn.model.Definition" value=" tlb  = if (t  !>=  SEIJR_START ) then 2020.085 
	else if ( t !>= 2020.161) then 2020.123
	else if ( t  !>= 2020.199) then 2020.161 
	else if ( t  !>=  2020.238) then 2020.199
	else ( 2020.238 ) "/>
<definition id="Definition.propub" spec="phydyn.model.Definition" value=" propub = max(0, min(1, (t-tlb)/ (2020.238 - 2020.199)  ) ) "/>
<definition id="Definition.beta" spec="phydyn.model.Definition" value=" beta = max(0, propub * betaub + (1-propub)*betalb) "/>
') 

.comp_beta <- function( taxis, theta, changePoints = c(2020.123, 2020.161, 2020.199, 2020.238) ) {
	dlogbetanames = paste0('seir.dlogbeta.t',1:(length( changePoints )-1) )
	betas = theta$seir.b * c(1,  exp( cumsum(unlist(theta[dlogbetanames])  )) )
	betaslb = rep(NA, length(taxis ))
	betasub = rep(NA, length(taxis ))
	betas_t = rep(NA, length(taxis ))
	betaslb[ taxis < changePoints[1] ] =  betasub[ taxis < changePoints[1] ] <- betas[1]  #first interval 
	betaslb[ taxis >= tail(changePoints,1) ] = betasub[ taxis >= tail(changePoints,1) ] <- tail( betas,1) # last interval 
	
	betas_t[ taxis < changePoints[1] ] <- betas[1]
	betas_t [  taxis >= tail(changePoints,1)  ] <-  tail( betas,1)
	for ( i in 1:(length(changePoints)-1)){
		j <- which( (taxis >= changePoints[i] ) & (taxis < changePoints[i+1] ) )
		betaslb [ j ] <- betas[i]
		betasub [ j ] <- betas[i+1]
		propub = (taxis[j] - changePoints[i]) / abs(changePoints[i] - changePoints[i+1])
		betas_t[j] <- propub * betasub[j] + (1-propub)*betaslb[j]
	}
#~ browser()
	betas_t 
}

#debug: 
#~ library( sarscov2 )
#~ library( ggtree ) 
#~ library( ggplot2 )
#~ logfn = 'gmrf.SHG3.2.xml.log'
#~ trajfn = 'gmrf.SHG3.2.xml.traj'
#~ X = read.table( logfn , header=TRUE )
#~ J = read.table( trajfn, header = TRUE )




#' Plot the daily new infections through time from a GMRF SEIJR trajectory sample
#' 
#" Also computes CIs of various dynamic variables 
#'
#" @param trajdf Either a dataframe or a path to rds containing a data frame with a posterior sample of trajectories (see combine_traj)
#' @param logdf Either a dataframe or a path to rds containing a data frame with posterior logs
#' @param case_data An optional dataframe containing reported/confirmed cases to be plotted alongside estimates. *Must* contain columns 'Date' and 'Confirmed'. Ensure Date is not a factor or character (see as.Date )
#' @param date_limits  a 2-vector containing bounds for the plotting window. If the upper bound is missing, will use the maximum time in the trajectories
#' @param path_to_save Will save a png here 
#' @param changePoints the time where beta changes
#' @return a list with derived outputs from the trajectories. The first element is a ggplot object if you want to further customize the figure 
#' @export 
GMRF_SEIJR_plot_daily_inf <- function(trajdf
  , logdf
  , case_data = NULL
  , date_limits = c( as.Date( '2020-02-01'), NA ) 
  , path_to_save='daily.png'
  , log_y_axis = F
  , changePoints =  c(2020.123, 2020.161, 2020.199, 2020.238)
  , ...
) {
cat( 'Double check changePoints for your model\n' ) 
print( changePoints )

	library( ggplot2 ) 
	library( lubridate )
	
	if (!is.data.frame(trajdf)){
		# assume this is path to a rds 
		readRDS(trajdf) -> trajdf 
	}
	
	if (!is.data.frame( logdf ))
		X <- readRDS(logdf ) else X <- logdf
	
		if(is.null(X$seir.tau))
		  X$seir.tau <- 74; 
		
		if(is.null(X$seir.p_h))
		  X$seir.p_h <- .2
		
	dfs <- split( trajdf, trajdf$Sample )
	taxis <- dfs[[1]]$t 
	
	s <- sapply( dfs, function(df) df$Sample[1] )
	X1 <- X[match( s, X$Sample ), ]
	
	if ( is.na( date_limits[2]) )
		date_limits[2] <- as.Date( date_decimal( max(taxis)  ) )

	qs <- c( .5, .025, .975 )

	# infectious
	Il = do.call( cbind, lapply( dfs, '[[', 'Il' ))
	Ih = do.call( cbind, lapply( dfs, '[[', 'Ih' ))
	I = Il + Ih 
	t(apply( I, MAR=1, FUN= function(x) quantile(x, qs ))) -> Ici 
	
	# daily new inf 
	Y = lapply( 1:length(dfs) , function(k) {
		x = X1[k, ]
		Il = dfs[[k]]$Il
		Ih = dfs[[k]]$Ih
		E = dfs[[k]]$E
		tau = X1$seir.tau[k]
		p_h = X1$seir.p_h[k] 
		#b = X1$seir.b[k]
		theta = ( X1[k, ] )
		b = .comp_beta (dfs[[k]]$t, theta, changePoints )
		y = (1/365) * ( b * Il + b * tau * Ih   )
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
	  , Ih = Ih
	  , Y
	  , pldf =pldf 
	  , case_data = case_data 
	)
}
#~ p2 = GMRF_SEIJR_plot_daily_inf( J
#~   , X
#~   , case_data = NULL
#~   , date_limits = c( as.Date( '2020-02-01'), NA ) 
#~   , path_to_save='daily.png'
#~   , log_y_axis = F
#~   , changePoints =  c(2020.123, 2020.161, 2020.199, 2020.238)
#~ )



#' Plot reproduction number through time from a SEIJR trajectory sample
#'
#" @param trajdf Either a dataframe or a path to rds containing a data frame with a posterior sample of trajectories (see combine_traj)
#' @param logdf Either a dataframe or a path to rds containing a data frame with posterior logs
#' @param date_limits  a 2-vector containing bounds for the plotting window. If the upper bound is missing, will use the maximum time in the trajectories
#' @param gamma0 rate of becoming infectious during incubation period
#' @param gamma1 rate of recovery once infectious
#' @param path_to_save Will save a png here 
#' @param changePoints the time where beta changes
#' @return a list with derived outputs from the trajectories. The first element is a ggplot object if you want to further customize the figure 
#' @export 
GMRF_SEIJR_plot_Rt <- function(trajdf
  , logdf
  , gamma0 = 73, gamma1 = 121.667
  , date_limits = c( as.Date( '2020-02-01'), NA ) 
  , path_to_save='Rt.png'
  , changePoints = c(2020.123, 2020.161, 2020.199, 2020.238)
  , ...
) {
cat( 'Double check changePoints for your model\n' ) 
print( changePoints )
	library( ggplot2 ) 
	library( lubridate )
	
	if (!is.data.frame(trajdf)){
		# assume this is path to a rds 
		readRDS(trajdf) -> trajdf 
	}
	
	if (!is.data.frame( logdf ))
		X <- readRDS(logdf ) else X <- logdf
		
		if(is.null(X$seir.tau))
		  X$seir.tau <- 74; 
		
		if(is.null(X$seir.p_h))
		  X$seir.p_h <- .2
	
	dfs <- split( trajdf, trajdf$Sample )
	taxis <- dfs[[1]]$t 
	
	s <- sapply( dfs, function(df) df$Sample[1] )
	X1 <- X[match( s, X$Sample ), ]
	
	if ( is.na( date_limits[2]) )
		date_limits[2] <- as.Date( date_decimal( max(taxis)  ) )

	qs <- c( .5, .025, .975 )
	
	Y = lapply( 1:length(dfs) , function(k) {
		x = X1[k, ]
		tau = X1$seir.tau[k]
		p_h = X1$seir.p_h[k] 
		theta = ( X1[k, ] )
		b = .comp_beta (dfs[[k]]$t, theta, changePoints )
		((1-p_h)*b/gamma1 + tau*p_h*b/gamma1)
	})
	Y = do.call( cbind, Y )
	Yci = t(apply( Y, MAR=1, FUN= function(x) quantile(na.omit(x), qs ))) 

#~ browser()
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
#~ p3 = GMRF_SEIJR_plot_Rt(J
#~   , X
#~   , gamma0 = 73, gamma1 = 121.667
#~   , date_limits = c( as.Date( '2020-02-01'), NA ) 
#~   , path_to_save='Rt.png'
#~ )



