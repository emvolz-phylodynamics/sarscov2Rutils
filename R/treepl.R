library( ape ) 
library( ggtree ) 
library ( treeio )



#' Wrapper for logcombiner 
#' 
#' @param burnin Integer number of MCMC iters to treat as burnin 
#' @param fns List of tree log files. If not provided, will recursively search the directory for logs 
#' @export 
tree_combiner_helper <- function( burnin , fns = NULL, ofn = 'combined.trees'){
	if ( is.null( fns ))
		fns = list.files( pattern='trees$', recursive=TRUE)
	cat( 'NOTE: these tree logs are being combined. Double check that these are the files you want to combine\n' )
	print( fns )
	
	command = paste( 'logcombiner', '-trees', '-burnin', format(burnin, scientific=FALSE), paste(collapse=' ', fns ) , ofn )	

		if(Sys.info()["sysname"] == "Windows") {
	  
		  # logcombiner in Windows seems to take different commands, so it's rewritten here.
		  # you have to paste in your filepath to logcombiner.bat -- someone can generalise this if needs be
		  # note that burnin is required as a percentage (eg 50) instead of an absolute number of trees
		  # this command does not use the burnin arg
		  
		  # the line below runs without downsampling; this made 60000 trees in the combined tree, and I didn't have enough memory to annotate it
		  # command = paste( 'C:/Users/lilyl/Downloads/BEAST_with_JRE.v2.6.2.Windows/BEAST/bat/logcombiner.bat', '-log', paste(collapse=' ', fns ), '-b', format(50, scientific=FALSE) ,'-o', ofn ) 
		  # below is the same but downsampling even further, in this case to every 10000 trees
		  command = paste( 'C:/Users/lilyl/Downloads/BEAST_with_JRE.v2.6.2.Windows/BEAST/bat/logcombiner.bat', '-log', paste(collapse=' ', fns ), '-b', format(50, scientific=FALSE) ,'-o', ofn, '-resample', 10000 )
		  
	}
	  
	
	print ( command ) 
	system ( command ) 
	TRUE
}
#~ tree_combiner_helper( 5000000 )

#' Tree annotator in Windows
#' 
#' @param burnin Integer number of MCMC iters to treat as burnin 
#' @param fns List of tree log files. If not provided, will recursively search the directory for logs 
#' @export 
tree_annotator_windows <- function( inputfile = "combined.trees", outputfile = "mcc.nex" ) {
  command = paste("C:/Users/lilyl/Downloads/BEAST_with_JRE.v2.6.2.Windows/BEAST/bat/treeannotator.bat -limit 0.5 -burnin 0", inputfile, outputfile)
  print ( command ) 
  system ( command ) 
  TRUE
}
#~ tree_annotator_windows()

#' Plot the maximum clade credibility tree and showing HPD node heights for nodes with >50 per cent node support. Tips sampled from within the specified demes will be coloured red.
#'
#' You will need to combine your tree logs like this:
#' logcombiner -trees -burnin <integer corresponding to half of sample, eg 5000000> <treelog1> ... <treelog_n>  <combined.trees>
#' You can use the *tree_combiner_helper* function. See ?tree_combiner_helper
#' NOTE make sure that tree logs included here only correspond to those beast runs that passed checks in *combine_logs_and_traj* 
#" The MCC nexus file can be made using treeannotator2 like this: 
#' treeannotator2 -limit 0.5 -burnin 0 <combined>.trees <outputFileName>.nex
#'
#' @param nexfn Path to nexus file containing annotated MCC tree (output of treeannotator2)
#' @param mostRecentSampleDate A character string containign the date of the most recent sample in the form 2020-03-17
#' @param regionDemes A character vector of deme names which will be coloured red 
#' @param ofn Output file name of figure 
#' @return A ggtree object which can be customized further
#' @export 
mcc_tree_plot <- function( nexfn, mostRecentSampleDate, regionDemes = c( 'Il', 'Ih', 'E' ), ofn ='mcc.png' )
{
	tr = treeio::read.beast ( nexfn )
	trd = treeio::get.tree( tr ) 
	btr = ggtree(tr, mrsd=mostRecentSampleDate, ladderize=TRUE) + geom_range(range='height_0.95_HPD', color='steelblue', alpha=.4, size=2) + theme_tree2() 
	tipdeme <- sapply( strsplit( trd$tip.label, '_' ), function(x) tail(x,1))
	tipdata <- data.frame( 
	  taxa =trd$tip.label
	  , region =  tipdeme %in% regionDemes
	)
	tipdata$size <- .25
	tipdata$size[ !tipdata$region ] <- 0
	tipdata$region[ !tipdata$region ] <- NA
	btr <- btr %<+% tipdata 
	btr = btr + geom_tippoint( aes(color = region, size = size), na.rm=TRUE, show.legend=FALSE, size =1.25) + theme_tree2( legend.position = "none" )
	#~ > decimal_date( as.Date( '2020-01-10' ) )
	#~ [1] 2020.025

	ggsave( btr, file= ofn , width = 4, height=7)

	btr
}
# mcc_tree_plot( 'mcc.nex', mostRecentSampleDate = '2020-03-17' )


#' Plot a rooted tree (may be a time tree) with tips highlighting a matching regex
#'
#' @param td A rooted tree, may be a treedater output 
#' @param region A regular expression matching the tips that you want to highlight 
#' @param maxdate If a time tree, you can optionally give the date of the most recent sample
#' @return A ggtree plot
#' @export 
quick_region_treeplot = function( td, region , maxdate = NULL)
{ #date_decimal( max(tr$sts))
	require(ggtree )
	require(ggplot2)
	tr = td 
	class( tr ) = 'phylo'
	btr = ggtree(tr, mrsd= maxdate, ladderize=TRUE)  + theme_tree2() 
	tipdeme <-  grepl( tr$tip.label, pat = region ) 
	tipdata <- data.frame( 
	  taxa =tr$tip.label
	  , region =  tipdeme 
	)
	tipdata$size <- .75
	tipdata$size[ !tipdata$region ] <- 0
	tipdata$region[ !tipdata$region ] <- NA
	btr <- btr %<+% tipdata 
	btr = btr + geom_tippoint( aes(color = region, size = size), na.rm=TRUE, show.legend=FALSE, size =1.25) + theme_tree2( legend.position = "none" )

	btr + ggplot2::ggtitle( region )
}

#~ --------------------------
#~ ml and treedater 
if(FALSE)
{
library( treedater )
tr = unroot( read.tree( '../algn.21.2.fasta.treefile' ) )
sts <- sapply( strsplit( tr$tip.label , '_' ), function(x) as.numeric( tail(x,2)[1] ))
names(sts) <- tr$tip.label 
td = dater( tr, sts, s = 29e3, omega0 = .001 )
trtd = td; class( trtd ) <- 'phylo' 
library( ggtree ) 
btr = ggtree(ladderize(trtd), mrsd="2020-02-10") + theme_tree2() 
tipdata <- data.frame( taxa = trtd$tip.label, weifang = grepl('WFCDC', trtd$tip.label) )
tipdata$size <- .25
tipdata$size[ !tipdata$weifang ] <- 0
tipdata$weifang[ !tipdata$weifang ] <- NA
btr <- btr %<+% tipdata 
btr = btr + geom_tippoint( aes(color = weifang, size = size), na.rm=TRUE, show.legend=FALSE, size =1.25) + theme_tree2( legend.position = "none" )

ggsave( btr, file='treedatertree.pdf', width = 4, height=7)
ggsave( btr, file='treedatertree.svg', width = 4, height=7)
}

# rtt 
#~ > rootToTipRegressionPlot( td ) 
#~ Root-to-tip mean rate: 0.00095787150868136 
#~ Root-to-tip p value: 3.93364980244933e-06 
#~ Root-to-tip R squared (variance explained): 0.3439734573654 
if (FALSE){
rttpl <- function()
{
	par(mai = c(.6, .91, .05, .15 ))
	rootToTipRegressionPlot( td, show.tip.labels=F, pch = 20, cex = 1,bty='n' ) 
}
pdf( 'rtt.pdf', width=3.25, height=2.8); rttpl(); dev.off() 
svg( 'rtt.svg', width=3.25, height=2.8); rttpl(); dev.off() 
}

