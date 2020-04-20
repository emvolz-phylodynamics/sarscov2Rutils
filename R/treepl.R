library( ape ) 
library( ggtree ) 
library ( treeio )



#' Wrapper for logcombiner 
#' 
#' @param burnin Unix: Integer number of MCMC iters to treat as burnin; In windows this is the percentage burnin as a number
#' @param fns List of tree log files. If not provided, will recursively search the directory for logs 
#' @param resample set this to downsample logs further when combining
#' @export 
tree_combiner_helper <- function( burnin , fns = NULL, ofn = 'combined.trees', resample = NULL){
  if ( is.null( fns ))
    fns = list.files( pattern='trees$', recursive=TRUE)
  cat( 'NOTE: these tree logs are being combined. Double check that these are the files you want to combine\n' )
  print( fns )
  
  command = ifelse(is.null(resample),
                   paste( 'logcombiner', '-trees', '-burnin', format(burnin, scientific=FALSE), paste(collapse=' ', fns ) , ofn ),
                   paste( 'logcombiner', '-trees', '-burnin', format(burnin, scientific=FALSE), paste(collapse=' ', fns ) , ofn, '-resample', resample))
  
  if(Sys.info()["sysname"] == "Windows") {
    if(system("logcombiner -help >NUL 2>NUL")==0){
      command = ifelse(is.null(resample),
                       paste( 'logcombiner', '-log', paste(collapse=' ', fns ), '-b', format(burnin, scientific=FALSE) ,'-o', ofn),
                       paste( 'logcombiner', '-log', paste(collapse=' ', fns ), '-b', format(burnin, scientific=FALSE) ,'-o', ofn, '-resample', resample))
    } else {
    message("
            logcombiner is not in your path.
            Either add the BEAST/bat folder to path
            Or run the command with the full filepath to logcombiner.bat
            command = ", paste( 'C:/Users/lilyl/Downloads/BEAST_with_JRE.v2.6.2.Windows/BEAST/bat/logcombiner.bat', '-log', paste(collapse=' ', fns ), '-b', format(50, scientific=FALSE) ,'-o', ofn ) 
            )
      command = NULL
    }
      # you have to paste in your filepath to logcombiner.bat -- someone can generalise this if needs be
    # note that burnin is required as a percentage (eg 50) instead of an absolute number of trees

    
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
#' @param lowMem Set to TRUE to run treeanotator in low memory mode.
#' @export 
tree_annotator_windows <- function( inputfile = "combined.trees", outputfile = "mcc.nex", lowMem=FALSE ) {
  if(system("treeannotator -help >NUL 2>NUL")==0){
    command = ifelse(lowMem,
    paste( 'treeannotator', '-lowMem -limit 0.5 -burnin 0', inputfile, outputfile),
    paste( 'treeannotator', '-limit 0.5 -burnin 0', inputfile, outputfile))
    
  } else {
    message("
            treeannotator is not in your path.
            Either add the BEAST/bat folder to path
            Or run the command with the full filepath to treeannotator.bat
            command = paste('C:/Users/lilyl/Downloads/BEAST_with_JRE.v2.6.2.Windows/BEAST/bat/treeannotator.bat -limit 0.5 -burnin 0', inputfile, outputfile)
            ")
    command = NULL
  }
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

#' plotting tree with background sequences labelled
#' 
#' Edit of mcc_tree_plot allowing non-region tips to be coloured by continent.
#' requires nexus input tree, most recent sample date and dictionary converting second element of 
#' sequence names to continent.
#' 
#' @param nexfn Path to nexus file containing annotated MCC tree (output of treeannotator2)
#' @param mostRecentSampleDate A character string containign the date of the most recent sample in the form 2020-03-17
#' @param country_dictionary data dictionary to convert country names to world bank continents
#' @param internal A character string giving the name of the region of interest
#' @param style plot style from 1,2 or 3
#' @param regionDemes A character vector of deme names which will be coloured red 
#' @param ofn Output file name of figure 
#' @return A ggtree object which can be customized further
#' @export 
mcc_col_tree_plot <- function(nexfn
                              , mostRecentSampleDate
                              , country_dictionary
                              , internal = "Region"
                              , style = 1
                              , regionDemes = c( 'Il', 'Ih', 'E' )
                              , ofn = 'mcc.png'
  ){
  
  country_dict <- utils::read.table(country_dictionary, header = TRUE)
  tr = treeio::read.beast ( nexfn )
  trd = treeio::get.tree( tr )
 
  btr = ggtree::ggtree(tr, mrsd=mostRecentSampleDate, ladderize=TRUE) + ggtree::geom_range(range='height_0.95_HPD', color='steelblue', branch.length="height", alpha=.4, size=2) + ggtree::theme_tree2() 
  
  
  tipdeme <- sapply( strsplit( trd$tip.label, '_' ), function(x) tail(x,1))
  tipcols <- sapply( strsplit( trd$tip.label, '/' ), function(x) x[[2]])
  
  ## recode tipcols w/ dictionary
  tipcols <- country_dict$continent[match(tipcols, country_dict$name)]
  
  continent_names <- c("Africa", "Americas", "Asia", "Europe", "Oceania", internal)
  levels(tipcols) <- continent_names
  if (style == 1){
  col_pal <- c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#999999", "#f50000")
  }else{
  col_pal <- c("#cbf4ae", "#ffd0a8", "#ffb1b1", "#d9d1ff", "#b7efff", "#f50000")
  }
  names(col_pal) <- continent_names
  shape <- ifelse(tipdeme %in% regionDemes, 19,1)
  size <- ifelse(tipdeme %in% regionDemes, 2,1.5)
  tipdata <- data.frame( 
    taxa =trd$tip.label
    , region =  tipdeme %in% regionDemes
    , Location = tipcols
    , shape = shape
    , size = size
    , stringsAsFactors = F
  )
  tipdata[tipdata$region==TRUE,]$Location <- internal
  
  btr <- btr %<+% tipdata
   
    if(style == 1){
      btr = btr +
        ggtree::geom_tippoint( aes(color = Location, shape = shape), size = 2) +
        ggplot2::scale_shape_identity() +
        ggtree::theme_tree2( ) +
        ggplot2::scale_color_manual(values = col_pal, 
                                    guide = ggplot2::guide_legend(title.position="top",
                                                                  title.hjust = 0.5,
                                                                  nrow=2,
                                                                  byrow=TRUE ))+
        ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_text(face="bold", size=10))
        
        }
    if(style == 2){
      btr = btr +
        ggtree::geom_tippoint( aes(color = Location), size = 2) + 
        ggtree::theme_tree2( ) +
        ggplot2::scale_color_manual(values = col_pal, 
                                    guide = ggplot2::guide_legend(title.position="top",
                                                                  title.hjust = 0.5,
                                                                  nrow=2,
                                                                  byrow=TRUE ))+
        ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_text(face="bold", size=10))
    }
    if(style == 3){
      btr = btr +
        ggtree::geom_tippoint( aes(color = Location, size = size)) +
        ggplot2::scale_size_identity()+
        ggtree::theme_tree2( ) +
        ggplot2::scale_color_manual(values = col_pal, 
                                    guide = ggplot2::guide_legend(title.position="top",
                                                                  title.hjust = 0.5,
                                                                  nrow=2,
                                                                  byrow=TRUE))+
        ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_text(face="bold", size=10))
        
    }

  ggsave( btr, file= "mcc2.png" , width = 4, height=7)
  return(btr)
}
# mcc_col_tree_plot("./mcc.nex", "2020-03-29", "country_dict.txt", internal = "Reykjavik", style = 3)
		    
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

