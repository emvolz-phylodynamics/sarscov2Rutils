

.instantiate_sir7 <- function()
{
	invisible('
	 E , I 
	Rt logistic 
	')

	library( phydynR ) 
	library( lubridate ) 
	
	parms <- list( 
		t0 = 2020
		# transmission 
		, R0 = 3.3 
		, min_R = 0.5 
		, logistic_location = lubridate::decimal_date( as.Date('2020-03-21') ) 
		, logistic_rate = 48
		#init cond
		, I0 = 1e-4
		# natural history 
		, gamma0 = 365 / 5.1
		, gamma1 =  365 / 2.5
		, treeIndex = 1 
	)
	
	# phydyn model 
	
	demes = c('E', 'I' )
	nd = length(demes)
	nondemes = c('cinf')
	
	default_x0 <- c( E = parms$I0, I = parms$I0 ,  cinf = 0  )
	
	
	## Rt
	parms$.Rt <- function( t, p )
	{
		tau = t-p$logistic_location #location
		y = exp( -p$logistic_rate * tau ) * ( 1 + exp( -p$logistic_rate * tau )  )^(-1) #rate
		(Rt = p$min_R  +  (p$R0 - p$min_R)*y ) #scale 
	}

	## eqns 
	births <- matrix('0', nrow=nd, ncol=nd)
	rownames(births)=colnames(births) <- demes

	migs <- matrix('0', nrow=nd, ncol=nd)
	rownames(migs)=colnames(migs) <- demes

	deaths <- rep('0',nd)
	names(deaths) <- demes 

	nonDemeDynamics <- c() 

	### transmission
	births['I', 'E'] <- 'with(parms, I * (.Rt(t,parms)*gamma1)  )' 
	### death 
	deaths['I'] <- 'I * parms$gamma1' 
	###progression
	migs['E', 'I'] <- 'E * parms$gamma0'
	###ndd 
	nonDemeDynamics['cinf'] <-  births['I', 'E']

	dm <- build.demographic.process( births = births
	 , migrations = migs 
	 , deaths = deaths 
	 , nonDemeDynamics = nonDemeDynamics
	 , parameterNames = names(parms)
	 , sde = FALSE
	 , rcpp = FALSE
	)
	
	list( dm = dm, parms = parms, default_x0 = default_x0, demes = demes, nondemes = nondemes, nd = nd )
}


fit_model7 <- function(
  tds
  
  , continueFits = NULL 
  , MCMCITER = 2000*20 
  , NCPU = 20
  
  , MAXHEIGHT_TIME = lubridate::decimal_date(as.Date('2020-02-01') ) 
  , TIME0 = 2020.00
  
  , LIKSWITCH='PL2'
  , STEP_SIZE_RES = 5
  
  , outfn = NULL 
)
{
	library( phydynR ) 
	library( treedater ) 
	library( ape ) 
	library( lubridate ) 
	library( BayesianTools ) 
	
	NTREES = length(tds)
	sts = tds[[1]]$sts
	MAXHEIGHT <- max(sts) - MAXHEIGHT_TIME
	ssts1 <- matrix( c(0, 1), nrow = Ntip(tds[[1]]), byrow=T, ncol =2)
	colnames(ssts1) <- c( 'E', 'I')
	rownames(ssts1) <- tds[[1]]$tip.label
	
	attach( .instatiate_sir7())
	
	bdts = lapply( tds, function(x) { 
		tr = x; class(tr) <- 'phylo'
		DatedTree(tr, x$sts[tr$tip.label ], ssts1[tr$tip.label, ], tol  = .01)  
	})
	
	
	comp.tri <- function( x ){
		max(1 ,floor( (NTREES * ( (1e6 * sqrt(abs(x))  %%1 ) %%1 )  ) ) )
	}
	
	ESTIMATE <- c('R0',  'min_R', 'logistic_location', 'logistic_rate', 'I0', 'treeIndex',  'mu' )
	ELOWER <- c( 
	  R0 = 1
	  , min_R = .1
	  , logistic_location = 2020.1
	  , logistic_rate = 12
	  , I0= 1e-12
	  , treeIndex = -Inf 
	  , mu = .01
	)
	EUPPER <- c( 
	  R0 = 20
	  , min_R = 1
	  , logistic_location = 2020.35
	  , logistic_rate = 200
	  , I0 = 1
	  , treeIndex = Inf
	  , mu = 24
	)


	.of <- function(x, likswitch = 'PL1', ssr = STEP_SIZE_RES , res = 400 ){
		names(x) = ESTIMATE 
		p = modifyList( parms, as.list( x ))
		x0 = default_x0 
		x0['E'] <- unname( p[['I0']] )
		x0['I'] <- unname( p[['I0']] )
		
		tri <- comp.tri ( p[['treeIndex']] )
		bdt <- bdts[[ tri ]]
		ll = colik( bdt, p, dm, x0, TIME0, res = res, AgtY_penalty=0 , maxHeight = MAXHEIGHT , likelihood = likswitch, step_size_res = ssr ) 
		if (is.na(ll))
			ll <- -Inf 
		
		cat( '-------------------------\n' )
		print( Sys.time() )
		print( x )
		print( ll ) 
		ll 
	}
	
	
	.sampler <- function( n=1 )
	{
		s = list( 
		  R0 = rlnorm( n, log( 3.3 ), 1 )
		  , min_R = rlnorm( n, log( .1 ), 1 )
		  , logistic_location = rnorm( n,2020.219, 0.025)
		  , logistic_rate = runif( n, 12, 200 )

		  , I0 = rexp( n, rate = 100 )
		  , treeIndex = runif( n )
		) 
		
		if ( n == 1 ) {
			s <- pmax( ELOWER, unlist( s ) )
			s <- pmin( EUPPER, s )
			names(s) <- ESTIMATE 
		}else{
			s <- do.call( cbind, s ) 
			for ( i in 1:n){
				s[i,] <- pmax( ELOWER, s[i,] )
				s[i,] <- pmin( EUPPER, s[i,] )
			}
		}
		s
	}
	#~ .sampler() 


	.priordensity <- function(x)
	{
		names(x) <- ESTIMATE 
		if ( is.vector( x)){
			y = with( as.list( x ), {
				dlnorm( R0, log(3.3), 1, log=TRUE ) +
				
				dlnorm( min_R, log(.1), 1, log=TRUE) + 
				dnorm( logistic_location,2020.219, 0.025, log=TRUE) +
				
				dexp( I0, rate = 100 , log = TRUE )  + 
				dlnorm( mu, log(1), .5 , log = TRUE)
			})
		} else{
			y = with( as.data.frame( x ), {
				dlnorm( R0, log(3.3), 1, log=TRUE ) +
				
				dlnorm( min_R, log(.25), 1, log=TRUE) + 
				dnorm( logistic_location,2020.202, 0.025, log=TRUE) +
				
				dexp( I0, rate =  1 / 10, log = TRUE )  + 
				dlnorm( mu, log(1), .5 , log = TRUE)
			})
		}		
		unname( y  )
	}
	#~ .priordensity( .sampler( 3 ))


	prior <- createPrior(density = .priordensity,
						 sampler = .sampler,
						 lower = ELOWER ,
						 upper = EUPPER 
	)
	
	
	cl <- parallel::makeCluster(NCPU )
	.of1 <- function(x, likswitch = LIKSWITCH){
		if ( is.vector(x))
			return( .of( x, likswitch ))
		
		parallel::parApply(cl = cl, X = x, MARGIN = 1, FUN = .of, likswitch = likswitch)
	}
	
	parallel::clusterEvalQ(cl, library(BayesianTools))
	parallel::clusterEvalQ(cl, library(phydynR))
	parallel::clusterEvalQ(cl, library(parallel))
	parallel::clusterExport(cl=cl, varlist = c('.instantiate_sir7', '.of', 'bdts', 'ESTIMATE', 'MAXHEIGHT', 'NCPU', 'NTREES', '.sampler', 'ELOWER', 'EUPPER', 'TIME0' , 'comp.tri', 'STEP_SIZE_RES'))
	parallel::clusterEvalQ(cl, attach( .instantiate_sir7() ) )
	
	st0 = Sys.time() 
	if (!is.null(continueFits)){
		f <- runMCMC( continueFits 
		  , sampler = "DEzs"
		  , settings =  list(iterations = MCMCITER, startValue = .sampler(NCPU))
		)
		
	} else{
		f <- runMCMC(
		  bayesianSetup = createBayesianSetup(likelihood = .of1 , prior = prior, names = ESTIMATE , parallel = 'external')
		  , sampler = "DEzs"
		  , settings =  list(iterations = MCMCITER, nrChains=1, startValue = .sampler(NCPU ))
		)
		saveRDS( f0, paste0( ofnpre, '-f0.rds' ) )
	}
	parallel::stopCluster(cl)
	st1 = Sys.time() 
	
	if ( is.null( outfn ))
		saveRDS( f, paste0( 'fit_model7-', Sys.getpid(), '.rds' ) )
	else 
		saveRDS(f, outfn )
	
	f
}
