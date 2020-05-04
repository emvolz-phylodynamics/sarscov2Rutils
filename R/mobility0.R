
#' Compute time series cross correlations and plots using google mobility and R(t) estimates
#'
#' This function requires the pika package https://github.com/mrc-ide/pika/
#' NOTE this function computes the rolling correlation for a range of lag days, but does not currently use this information
#' The lag data is returned. 
#' If R(t) is out of sync with mobility, better correlations can be found by shifting the dates of the R(t) data 
#'
#' @param mobdf A data frame with google mobility data for the selected region
#' @param J Combined trajectory file from PhyDyn analysis 
#' @param X Combined log file from PhyDyn analysis 
#' @param region A title for the plots describing the region of analysis
#' @param regionshort A short title that will go in file names 
#' @param The output of a function to compute R(t), such as SEIJR_plot_Rt or gmrf_exogsir_plot_Rt
#' \dontrun{
#' # load the mobility data:
#' mobdf0 = read.csv( system.file( 'extdata/googmob-1may2020.csv', package = 'sarscov2' ) , stringsAs=FALSE)
#' # extract new york: 
#' mobdf <- mobdf0[ mobdf0$sub_region_1=='New York' & mobdf0$sub_region_2=='New York County', ]
#' # make a nice header for the plot
#' region = 'New York City, USA'
#' # load the trajectory and log files:
#' J = readRDS('traj.rds' ) 
#' X = readRDS('logs.rds' )
#' # run it: 
#' o = googmobility( mobdf, J, X, region = 'New York City, USA', regionshort='newyork' , rt = gmrf_exogsir_plot_Rt(J,  X )  )
#' }
#' @return A list with the tabulated time series and outputs of pika::rolling_corr
#' @export
googmobility <- function(mobdf, J, X, region = '', regionshort='', rt = NULL ) 
{
	stopifnot( require( pika )  )
	stopifnot( require( dplyr )  )
	
	
	if (is.null(rt))
		rt = SEIJR_plot_Rt(J,  X ) 
	
	df = rt$pldf 
	df$Date = as.Date( df$Date ) 
	colnames( df) <- c( 'date', 'reported', 'R', 'lb','ub' ) 
	
	mobdf$date = as.Date( as.character( mobdf$date ))
	df1 = merge( df, mobdf, by.y='date' ) 
	colnames(df1) = sapply( strsplit( colnames(df1 ), '_' ) , '[', 1 )
	df1$grp = 1 
	colnames(df1)[6] = 'country1'
	colnames(df1)[8] = 'sub1'
	
	
	mobnames = c( 'grocery', 'parks',  'transit', 'workplaces',  'residential' )
	#~ scale( df1[, mobnames] )
	df1[, mobnames ] <- scale( df1[, mobnames] )
	df1[, 'parks_workplaces_grocery'] <- rowMeans( df1[, c('parks', 'workplaces', 'grocery') ] )
	
	lags <- cross_corr(dat = df1,
	   date_var = "date",
	   x_var = "R",
	   y_var = "transit",
	   max_lag = 10, 
	   grp_var = 'grp'
	)
	
	#~ plot( with( df1, TTR::runCor( R, transit ) ) )
	
	.pl_scale_mob <- function(){
		par(bty= 'n')
		matplot( lubridate::decimal_date(df1$date), cbind( scale(df1$R), df1[, c('transit', 'residential', 'parks_workplaces_grocery')] )
		  , type = 'l'
		  , ylab = ''
		  , xlab = ''
		  , main = region
		  , lty = 1 
		  , lwd= c( 4, 2, 2, 2)
		) 
		legend( x = 'topleft', y = NULL,   c('R', 'transit', 'residential', 'other') , col = 1:4, lty = 1)
	}
	
	.corplot <- function(yvar = 'transit')
	{
		df2 <- df1 %>% select(date, grp, R, lb, ub)
		df2[[yvar]] <- df1[[yvar]] 
		data_corr <- rolling_corr(dat = df2,
								  grp_var = "grp",
								  x_var = "R",
								  y_var = yvar,
								  n = 14)
		
		my_legend = c("Correlation", "Reproduction Number", "Movement")
		p = plot_corr(dat = data_corr,
				  date_var = "date",
				  grp_var = "grp",
				  x_var = "R",
				  y_var = yvar,
				  x_var_lower = "lb",
				  x_var_upper = "ub",
				  facet_labels = c('1' = 'grp'),
				  legend_labels = my_legend
		)
		ggsave(p, file= paste0('corplot-', yvar, '.pdf') )
		ggsave(p, file= paste0('corplot-', yvar, '.png') )
		ggsave(p, file= paste0('corplot-', yvar, '.svg') )
		data_corr 
	}
	
	
	dctransit = .corplot( 'transit' )
	dcresidence = .corplot( 'residential' )
	
	pdf(paste0(regionshort, '-scalemob.pdf')) ; .pl_scale_mob() ; dev.off() 
	png(paste0(regionshort, '-scalemob.png')) ; .pl_scale_mob() ; dev.off() 
	svg(paste0(regionshort, '-scalemob.pdf')) ; .pl_scale_mob() ; dev.off() 
	
	list( data = df1 , transit = dctransit, residential = dcresidence , lags = lags )
}
if (FALSE){
# load the mobility data:
mobdf0 = read.csv( system.file( 'extdata/googmob-1may2020.csv', package = 'sarscov2' ) , stringsAs=FALSE)
# extract new york: 
mobdf <- mobdf0[ mobdf0$sub_region_1=='New York' & mobdf0$sub_region_2=='New York County', ]
# make a nice header for the plot
region = 'New York City, USA'
# load the trajectory and log files:
J = readRDS('traj.rds' ) 
X = readRDS('logs.rds' )
# run it: 
o = googmobility( mobdf, J, X, region = 'New York City, USA', regionshort='newyork' , rt = gmrf_exogsir_plot_Rt(J,  X )  )
}
