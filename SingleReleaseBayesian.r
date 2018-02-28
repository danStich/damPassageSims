# Front-end code -----
# Software needs
  # The most recent version of JAGS needs to be installed on
  # the machine. This can be accessed through the files
  # page hosted through SourceForge:
  # https://sourceforge.net/projects/mcmc-jags/files/latest/download?source=files

# Package install and load
  #install.packages('R2jags')
  library(R2jags)

# Set working directory
  setwd("/media/stich/Windows/Users/STICHDS/Desktop/damPassageSim")

# Function definitions
  # All functions (and most comments) in this section are from
  # Kery and Schaub (2011) Bayesian Population Analysis

	# Define function to simulate a capture-history (CH) matrix
		simul.cjs <- function(PHI, P, marked){
			n.occasions <- dim(PHI)[2] + 1
			CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
			# Define a vector with the occasion of marking
				mark.occ <- rep(1:length(marked), marked[1:length(marked)])
			# Fill the CH matrix
				for (i in 1:sum(marked)){
					CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
					if (mark.occ[i]==n.occasions) next
					for (t in (mark.occ[i]+1):n.occasions){
						# Bernoulli trial: does individual survive occasion?
							sur <- rbinom(1, 1, PHI[i,t-1])
							if (sur==0) break		# If dead, move to next individual 
							# Bernoulli trial: is individual recaptured? 
								rp <- rbinom(1, 1, P[i,t-1])
								if (rp==1) CH[i,t] <- 1
					} #t
				} #i
			return(CH)
		 }

	# Initial values for the latent state z. At all occasions when an individual
  # was observed, its state is z = 1 for sure. In addition, if an individual 
  # was not observed at an occasion, but was alive for sure, because it was
  # observed before and thereafter (i.e. has a capture history of e.g. {101}
  # or {10001}), then we know that the individual was alive at all of these
  # occasions, and thus z = 1. Therefore, we should provide initial values of
  # z = 1 at these positions as well. The following function provides such
  # initial values from the observed capture histories:
	known.state.cjs <- function(ch){
		 state <- ch
		 for (i in 1:dim(ch)[1]){
				n1 <- min(which(ch[i,]==1))
				n2 <- max(which(ch[i,]==1))
				state[i,n1:n2] <- 1
				state[i,n1] <- NA
				}
		 state[state==0] <- NA
		 return(state)
		 }

	# Function to create a matrix of initial values
	# for latent state z
	  cjs.init.z <- function(ch,f){
		  for (i in 1:dim(ch)[1]){
				if (sum(ch[i,])==1) next
				n2 <- max(which(ch[i,]==1))
				ch[i,f[i]:n2] <- NA
				}
		  for (i in 1:dim(ch)[1]){
			  ch[i,1:f[i]] <- NA
		  }
		  return(ch)
		}

# Simulation settings -----
  # Sample sizes
    N <- seq(50, 300, 50)
	  
	# Number of runs
    total = 500

  # Number of dams
    ndams <- 3
    
  # Empty vector to hold bias and standard error 
  # estimates from posteriors. These are based on the
  # most downstream dam in the system because any 
  # error inflation or systematic bias should be most
  # evident in those reaches (we could do them separate).
    Bias = vector(mode='double', length=total)
    SE = vector(mode='double', length=total)
  # Empty vectors to hold mean parameter estimates for
  # each simulation run.
    Nm = vector(mode='double', length=total)
    Pd = vector(mode='double', length=total)
    rPhi1 = vector(mode='double', length=total)
    dPhi1 = vector(mode='double', length=total)
    dPhi2 = vector(mode='double', length=total)
    dPhi3 = vector(mode='double', length=total)
    rPhi2 = vector(mode='double', length=total)

# Simulation -----
	for(i in 1:total){
	# Parameter definitions
		n <- sample(N, 1, replace=FALSE)                                       # Sample size
		r <- 1.00                                       # Proportion released downstream
		riverPhi <- 0.98                                # In-river survival
		F <- .92/riverPhi                               # Proportional reduction in survival through dam
		damPhi <- F*riverPhi                            # Reach survival that contains background mortality
		p <- round(runif(1, 0.86, 1), 2)                 # Probability of detection   
		marked <- c(n*(r), 0, 0, 0)                     # Number released above dam
		phi <- c(riverPhi, 
		         rep(damPhi*riverPhi, ndams),
		         riverPhi)                              # Reach-specific survival rates for fish released above dam
		n.occasions = length(phi)+1
		det <- rep(p, n.occasions-1)                    # Static probabilty of detection

	# Define matrices with survival and recapture probabilities
		PHI <- matrix(phi, ncol = (n.occasions-1), nrow = sum(marked), byrow=TRUE)
		P <- matrix(det, ncol = (n.occasions-1), nrow = sum(marked))

	# Execute function to simulate individual capture histories and combine them
	# into a single matrix
		CH <- simul.cjs(PHI, P, marked)

	# Concatenate the capture histories from both groups
	# Create vector with occasion of marking
		get.first <- function(x) min(which(x!=0))
		f <- apply(CH, 1, get.first)

	# Specify model in BUGS language
		modelFileName <- 'singleRelease.txt'
		cat("
		model {

		# Derived quantities
			bias <- mean.phi[4] - 0.92
	 
		# Priors and constraints
			# Input parameters for likelihood
				for (i in 1:nind){
					for (t in 1:(n.occasions-1)){
						phi[i, t] <- mean.phi[t]
						p[i, t] <- mean.p
						} #t
				 } #i
				
			# Prior for detection, static across groups and reaches
				mean.p ~ dunif(0, 1)  

      # Prior for survival, separate for each reach
				for(t in 1:(n.occasions-1)){
					mean.phi[t] ~ dunif(0, 1)
				}         

		# Likelihood 
			for (i in 1:nind){
				 # Define latent state at first capture
				 z[i,f[i]] <- 1
				 for (t in (f[i]+1):n.occasions){
						# State process
						z[i,t] ~ dbern(mu1[i,t])
						mu1[i,t] <- phi[i,t-1] * z[i,t-1]
						# Observation process
						y[i,t] ~ dbern(mu2[i,t])
						mu2[i,t] <- p[i,t-1] * z[i,t]
						} #t
				 } #i
			}
			",fill = TRUE, file=modelFileName)

	# Bundle data
		jags.data <- list(
			y = CH,
			f = f,
			nind = dim(CH)[1],
			n.occasions = dim(CH)[2],
			z = known.state.cjs(CH)
			)

	# Declare initial values
		inits <- function(){
			list(
				mean.phi = runif(5, 0, 1),
				mean.p = runif(1, 0, 1),
				z = cjs.init.z(CH, f)
			)
		}
	 
	# Parameters monitored
		parameters <- c("mean.phi",
		                "mean.p",
		                "bias"
		                )

	# MCMC settings
		ni <- 3600  # Number of iterations, including burn in
		nt <- 6     # Thinning parameter
		nb <- 600   # Number of burn-in runs
		nc <- 3     # Number of chains

	# Call JAGS from R
		sRel <- jags(jags.data,
		             inits,
		             parameters,
		             'singleRelease.txt',
		             n.chains = nc,
		             n.thin = nt,
		             n.iter = ni,
		             n.burnin = nb,
			           working.directory = getwd()
			           )

		Nm[i] = n
    Pd[i] = p
		rPhi1[i] = median(sRel$BUGSoutput$sims.list$mean.phi[,1])
		dPhi1[i] = median(sRel$BUGSoutput$sims.list$mean.phi[,2])
    dPhi2[i] = median(sRel$BUGSoutput$sims.list$mean.phi[,3])
    dPhi3[i] = median(sRel$BUGSoutput$sims.list$mean.phi[,4])
    rPhi2[i] = median(sRel$BUGSoutput$sims.list$mean.phi[,5])
		Bias[i] = mean(sRel$BUGSoutput$sims.list$bias)
    SE[i] = sd(sRel$BUGSoutput$sims.list$mean.phi[,4])

		# View the results of the simulation as they come out
			hist(SE[SE!=0], main='', xlab='Standard error', col='gray87')
			mtext(side=3, sprintf('%.4f', i/length(SE)), adj=1, cex=2)
      mtext(side=3, paste(n), adj=.5, cex=2)
			mtext(side=3, sprintf('%.6f', mean(SE[SE!=0])), adj=0, cex=2)
      
	# Summarize posteriors
		#print(sRel, digits = 3)
	}

# Collect the inputs and outputs
	resS = data.frame(Nm, 
	                  Pd, 
	                  F, 
	                  r, 
	                  rPhi1, 
	                  dPhi1, 
	                  dPhi2, 
	                  dPhi3, 
	                  rPhi2,
	                  Bias,
	                  SE) 

# Write the results to a text file
  write.table(resS, 
              paste('./BayesResultsSingle',
                    runif(1,0,100),
                    '.csv',
                    sep=''
              ),
              sep=',', 
              row.names=FALSE,
              col.names=TRUE,
              quote=FALSE,
              append=FALSE)


