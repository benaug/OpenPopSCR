#' Run MCMC algorithm for fixed parameter or primary period-specific parameter open population SCR model.
#' @param data a list produced by simOpenSCR or in the same format.  See description.
#' @param niter an integer indicating the number of MCMC iterations to run
#' @param  nburn an integer indicating the number of MCMC iterations to discard as burn in
#' @param nthin an integer indicating the MCMC thinning parameter. Record output on every nthin iterations.  nthin=1 corresponds to no thinning
#' @param M an integer indicating the size of the augmented superpopulation.  M should be much larger than then
#' number of individuals alive across all primary periods.
#' @param inits a list of user-supplied starting values with elements for lam0, sigma, gamma, phi, psi,
#' and sigma_t (if using a structure movement model).
#' @param proppars a list of tuning parameters for the proposal distributions with elements for lam0, sigma, gamma,
#' s1x, s2y (for fixed and metamu models), s2x, s2y (for independent, markov, and metamu models), sex, and sigma_t (if using a structured movement model).
#' @param storeLatent a logical indicating whether or not to store and return the posteriors for z, s1, and/or s2
#' @param jointZ a logical indicating whether you want to use the sequential or joint z update.
#' @param Rcpp a logical indicating whether or not to use Rcpp
#' @param ACtype A character indicating the type of activity centers.  "fixed" activity centers do not move between primary periods, "metamu" assumes there is a
#' meta activity center around which the primary period activity centers distributed following a bivariate normal distribution
#' with spatial scale parameter sigma_t.  Primary period activity centers are required to stay inside the state space.
#' "metamu2" is identical to "metamu" except primary period activity centers are allowed to leave the state space.
#'  markov" assumes activity centers in primary period 1 are distributed uniformly across the state space and the
#'  activy centers in primary period l+1 are a bivariate normal draw centered around the activity center in
#'  primary period l (but must stay within the state space).  "markov2" is the same as "markov" except each
#'  dispersal considers the available distances to disperse to in a use vs. availability framework. Discrete
#'  state space is required.  Finally, "independent" assumes animals randomly mix between primary periods.
#'  (spatial uniformity in all primary periods with independence between primary periods).
#'
#'  The "metamu" and "markov" models are truncated on rectangular state spaces.  This is required to estimate
#'  sigma_t and to a lesser degree, the Ns without bias.  If using a polygon state space, the "metamu2" model
#'  will provide unbiased estimates.  If using a polygon state space with a Markov model, switching to a discrete
#'  state space and ACype="markov2" will provide unbiased estimates.
#' @param obstype a character indicating the observation model "bernoulli" or "poisson"
#' @param primary an optional vector of length t with entries 1 if data was recorded in that primary period and 0 if not. This allows
#' population dynamics to occur at equal interval primary periods even if data was not recorded at each primary period. If not
#' entered, the population is assumed to have been sampled in all primary periods.
#' @param dSS an optional (N_SS x 2) matrix for a discrete state space locations that overrules the buff or vertices objects in "data".  A matrix with columns for x and y locations
#' @return  a list with the posteriors for the open population SCR parameters (out), s (s1xout,s1yout,s2xout,s2yout with
#' s1 being meta ACs and s2 being primary period ACs), and z.  s1x and yout are of dimension niter x M and s2x and yout and z are
#' of dimension niter x M x T
#' @param dualACup a logical to be used with discrete state spaces indicating whether a second, interpatch activity
#' center proposal is used.  Currently under development--probably should not use.
#' @author Ben Augustine, Richard Chandler
#' @description This function runs an MCMC algorithm for an open population SCR model.  The data list should have the following elements:
#' 1.  y, a n x maxJ x T capture history where maxJ is the maximum number of traps across primary periods, T is the number of primary periods, and n
#' is the number of animals captured.  If population was not observed in primary period l, enter all zeros
#' 2.  X,  a list with elements that consists of J[l] x 2 matrices housing the X and Y trap locations for each primary period, l.
#' If the traps do not vary across primary periods, just repeat the same traps in each list element. If population was not
#' observed in primary period l, enter NA instead of a matrix
#' 3.  K, a vector of size T indicating how many capture occasions there were in each primary period. If population
#' was not observed in primary period l, enter NA
#' 4.  J, a vector of size T indicating how many traps there were in each primary period.  If population was not
#' observed in primary period l, enter NA
#' 5. either buff or vertices.  buff is the fixed buffer for the traps to produce the state space.  It is applied to the minimum and maximum
#' X and Y trap locations across primary periods, producing a square or rectangular state space.  vertices is a *list* of matrices with the X and Y coordinates
#' of a polygonal state space with one polygon in each list element.  If there is just one polygon, the list is of length 1.
#' If there are many polygons separated by large distances, you should think about the implications of activity centers
#' possibly being stuck inside the polygons and check the activity center posteriors to see if this is happening.
#' 6. tf is an optional list of vectors of length T containing the trap operation information.  Each vector has one element for each trap
#' and indicates how many occasions each trap was operational in that primary period.
#' 7. primary is an optional a vector of length T with entries 1 if the population is to be observed in primary period l and 0 otherwise.
#'
#' inits sets the initial values and determines if parameters are fixed or primary period-specific. It must have elements "lam0"
#' "sigma", "gamma", "phi", and "psi".  If there is an element "sigma_t", the parameters of a "meta mu" or Markov mobile
#' activity center model will be estimated.  If length(lam0)=1, a single lam0 will be estimated while if length(lam0)=t,
#' lam0 will be primary period-specific.  This goes for sigma, gamma, and phi as well, except for gamma and phi length needs to be t-1
#' for between primary period-specific parameters.  A decent starting value for psi is (hypothesized) N/M.
#'
#' proppars is a list containing the tuning parameters for the parameters that use a Metropolis-Hastings update.
#' It must have elements "lam0", "sigma", "gamma", "s1x", "s1y", s2x", "s2y", and "propz" (if jointZ=FALSE). If a parameter is primary period-specific, it needs
#' the appropriate number of proppars. propz is the number of data augmentation z's to update in primary periods 2,...,t, so it should
#' be of length t-1.  Increasing propz improves mixing (up to a point) but increases computation time. If you set
#' an initial value for sigma_t, you need to provide proppars for sigma_t".  Finally, proppars$s1x
#' and proppars$s1y tune fixed activity centers or the meta mus and proppars$s2x and proppars$s2y tune primary period
#' activity centers.
#'
#' A note on the z samplers.  jointZ=TRUE will update all the z's for each individual at the same time while jointZ=FALSE
#' will update them sequentially.  For T=3-6ish, you get a greater effective sample size per unit time with the joint update than with the sequential update.
#' The joint update always mixes better, but takes longer as T increases.  It must calculate the likelihood for all possible
#' z histories for a given T on each update.
#'
#' A note on Rcpp.  I've seen speed gains of 3-7x so far using Rcpp, but this improvement appears to be much reduced
#' when switching to discrete state spaces where R has a fast vectorized operation for the majority of the calculations
#' used in this update.
#'
#' A note on discrete state spaces.  Currently, only contiguous discrete state spaces are supported, or state spaces with only
#' small gaps between patches.  If using a complicated discrete state space, the activity centers may get stuck in patches
#' in which they were initilized.  To see if the current algorithms work, you should run several chains starting at different
#' values of sigma_t and see if they converge.  A better algorithm is in the works.
#'
#' A final note on splitting data sets to get group-specific parameters (e.g. sex).  This will alter the interpretation
#' of per capita recruitment.  For example, if you run sexes separatly, you are estimating recruitment per number
#' of males or females in the population, which probably does not make sense.  See the mcmc.OpenSCR() function for
#' sex-specific parameters.
#'
#' @examples
#' \dontrun{
#' library(coda)
#' #You'll need to run all these examples longer for convergence and inference.
#' #2 primary periods of data,  all parameters fixed, fixed activity centers (default setting)
#' t=2
#' N=100
#' p0=0.5
#' lam0=-log(1-p0) #because we're using the hazard half-normal detection function
#' sigma=0.75
#' phi=0.8
#' gamma=0.3
#' buff=3 #trap buffer to define state space
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9)) #Trap XY can vary across primary periods
#' K=c(10,10) #Number of occasions can vary across primary periods as well
#' M=250
#' #Simulate some data
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,obstype="bernoulli")
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M) #Single initial values for a parameter will estimate a single parameter across primary periods
#' niter=1000
#' nburn=0
#' nthin=1
#' proppars=list(lam0=0.025,sigma=0.025,gamma=0.1,phi=0.1,s1x=0.2,s1y=0.2)#No proppars$propz because using default jointZ=TRUE
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars)
#' plot(mcmc(out$out))
#' summary(mcmc(out$out))
#' #Note, realized N for primary periods>1 will vary if simulating using gamma.
#' data$N
#'
#' #same but t=3
#' t=3
#' N=100
#' p0=0.5
#' lam0=-log(1-p0)
#' sigma=0.750
#' phi=0.7
#' gamma=0.3
#' buff=3
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9)) #Note we need 3 trap objects stuffed into the trap list
#' K=c(10,10,10) #and 3 numbers of occasions within primary period
#' M=425
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,obstype="bernoulli")
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M)
#' niter=1000
#' nburn=0
#' nthin=1
#' proppars=list(lam0=0.025,sigma=0.025,gamma=0.1,phi=0.1,s1x=0.2,s1y=0.2)
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars)
#' plot(mcmc(out$out))
#' summary(mcmc(out$out))
#' data$N
#'
#'#Now detection function varies by primary period
#' t=3
#' N=100
#' p0=c(0.4,0.5,0.6) #need 3 of these
#' lam0=-log(1-p0)
#' sigma=c(0.75,0.65,0.55) #these, too
#' phi=0.7
#' gamma=0.3
#' buff=3
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9))
#' K=c(10,10,10)
#' M=425
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,obstype="bernoulli")
#' #Starting at simulated values gives 3 lam0 and sigma starting values
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M) #Note now 3 initial values are specivied here
#' niter=1000
#' nburn=0
#' nthin=1
#' #Now we need 3 proppars for lam0 and sigma
#' proppars=list(lam0=rep(0.025,3),sigma=rep(0.025,3),gamma=0.1,phi=0.1,s1x=0.2,s1y=0.2)
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars)
#' plot(mcmc(out$out))
#' summary(mcmc(out$out))
#' data$N
#'
#' #Let's set lam0 and sigma back to fixed and vary gamma and phi
#' t=3
#' N=100
#' p0=0.5
#' lam0=-log(1-p0)
#' sigma=0.750
#' phi=c(0.7,0.8) #3 primary periods allows for 2 periods for population dynamics
#' gamma=c(0.2,0.3)
#' buff=2
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9))
#' K=c(10,10,10)
#' M=425 #Not sure if this is still enough...
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,obstype="bernoulli")
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M) #now phi and gamma are getting 2 inits
#' niter=1000
#' nburn=0
#' nthin=1
#' #phi and gamma get 2 proppars
#' proppars=list(lam0=0.025,sigma=0.025,gamma=c(0.1,0.1),phi=c(0.1,0.1),s1x=0.2,s1y=0.2)
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars)
#' plot(mcmc(out$out))
#' summary(mcmc(out$out))
#' data$N
#'
#' #Let's set everything back to fixed and do mobile activity centers between primary periods
#' #We're using a bivariate normal "metamu" model which is truncated because we are using
#' a rectangular state space.
#' t=3
#' N=100
#' p0=0.5
#' lam0=-log(1-p0)
#' sigma=0.750
#' sigma_t=0.7
#' #We'll use this to simulate mobile activity centers
#' ACtype="metamu"
#' phi=0.7
#' gamma=0.3
#' buff=4 #increase buffer to not constrain movement as much
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9))
#' K=c(10,10,10)
#' M=425
#' #Note, we're putting sigma_t in now
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,sigma_t=sigma_t,K=K,X=X,t=t,M=M,buff=buff,
#' ACtype=ACtype)
#' #Need an init for sigma_t
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M,sigma_t=sigma_t)
#' niter=1000
#' nburn=0
#' nthin=1
#' #Note there is a proppar for sigma_t, the meta mus (s1) and primary period activity centers (s2)
#' proppars=list(lam0=0.055,sigma=0.025,gamma=0.15,phi=0.1,s1x=0.5,s1y=0.5,s2x=0.4,s2y=0.4,sigma_t=0.04)
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,ACtype=ACtype)
#' plot(mcmc(out$out))
#' summary(mcmc(out$out))
#' data$N
#'
#' #Let's look at some activity centers
#' par(mfrow=c(1,1))
#' par(ask=FALSE)
#' idx=1 #which guy to plot
#' plot(out$s1xout[,idx],out$s1yout[,idx],xlim=c(0,12),ylim=c(0,12)) #meta mu posterior
#' points(out$s2xout[,idx,1],out$s2yout[,idx,1],col="red")
#' points(out$s2xout[,idx,2],out$s2yout[,idx,2],col="green")
#' points(out$s2xout[,idx,3],out$s2yout[,idx,3],col="blue")
#' #Plot true values
#' points(data$mu[idx,1],data$mu[idx,2],pch=4,col="yellow",cex=3,lwd=5)
#' points(data$s[idx,1,1],data$s[idx,1,2],pch=4,col="yellow",lwd=5,cex=3)
#' points(data$s[idx,2,1],data$s[idx,2,2],pch=4,col="yellow",lwd=5,cex=3)
#' points(data$s[idx,3,1],data$s[idx,3,2],pch=4,col="yellow",lwd=5,cex=3)
#'
#' #Now, let's see how well our chains mixed, etc. using the coda package (must fit this last model first)
#' library(coda)
#' #First, let's look at the effective sample size.  For credible intervals, we should shoot for at least 400
#' #although for N, I've seen that you can get accurate 95% credible intervals with much less than that.
#' #This might have something to do with it being discrete.
#' effectiveSize(mcmc(out$out))
#' #We see that we haven't run the chains long enough, or perhaps the mixing can be improved or maybe both
#' #Let's look at the mixing.  We want an acceptance rate around 23%
#' 1-rejectionRate(mcmc(out$out))
#' #We can't adjust the sampler for N other than switching from the joint to the sequential sampler
#' #and since we use a Gibbs update for phi, we can't adjust that.  But we can adjust the proppars
#' #for lam0, sigma, gamma, and sigma_t.  Here, everything looks good because I've already adjusted the
#' tuning parameters.
#'
#' But what about the activity center acceptance?  Here, we look at the acceptance rates for the primary period
#' #ACs in primary period 1.  Note, they will be identical for the x and y dimension
#' 1-rejectionRate(mcmc(out$s2xout[,,1]))
#' #You'll want to focus on the individuals that were actually captured in this primary period
#' (generally towards the beginning). You can also look at the other primary periods. I think you should aim for
#' no acceptance rates less than 0.2. Inevitabely, the primary period activity centers for individuals not captured
#' in a given year will have high acceptance rates.
#'
#' Now for the meta mus.
#' 1-rejectionRate(mcmc(out$s1xout))
#' #These look OK, but we could probably raise the tuning parameters a bit.
#'
#' #So we can do short runs (say of length 500) to assess and adjust the tuning parameters until we're seeing
#' #acceptance rates in the 20-40% range and then run a chain long enough to get an acceptable effective
#' #sample size for all parameters (but usually the N's are the limiting factor or sigma_t if in the model).
#' #If at any point a parameter is not updating, the tuning parameter is probably too large.  Ideally, you will
#' run multiple chains with some variation in starting values, especially sigma_t if it's in the model.
#'
#' #So which point estimates and credible intervals should we use?  I've found the posterior mode
#' #to be the closest to unbiased point estimator, with the posterior median being pretty good, too.
#' #I would not use the posterior mean, especially for N. Here is a function to get the mode.
#' library(mcmcGLMM)
#' posterior.mode(mcmc(out$out))
#' #You can get the posterior means, modes, and quantile credible intervals using coda's summary
#' summary(mcmc(out$out))
#' #Note the 2.5% and 97.5% quantiles can be used as the credible interval limits.  However, I've found
#' #the highest posteior density intervals to have slightly better frequentist coverage. coda will do that for you
#' HPDinterval(mcmc(out$out))
#'
#' #Some more model types:
#' #Here is a markov example
#' t=3
#' N=100
#' p0=0.5
#' lam0=-log(1-p0)
#' sigma=0.750
#' sigma_t=0.7
#' #We'll use this to simulate mobile activity centers
#' ACtype="markov"
#' phi=0.7
#' gamma=0.3
#' buff=4 #increase buffer to not constrain movement as much
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9))
#' K=c(10,10,10)
#' M=425
#' #Note, we're putting sigma_t in now
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,sigma_t=sigma_t,K=K,X=X,t=t,M=M,buff=buff,
#' ACtype=ACtype)
#' #Need an init for sigma_t
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M,sigma_t=sigma_t)
#' niter=1000
#' nburn=0
#' nthin=1
#' #Note there is only proppars for primary period activity centers (s2; no meta mus).
#' proppars=list(lam0=0.055,sigma=0.025,gamma=0.15,phi=0.1,s2x=0.2,s2y=0.2,sigma_t=0.04)
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,ACtype=ACtype)
#' plot(mcmc(out$out))
#' summary(mcmc(out$out))
#' data$N
#'
#' #Here is a discrete state space example with ACtype="markov2"
#' t=3
#' N=50
#' p0=0.3
#' lam0=-log(1-p0)
#' sigma=1
#' sigma_t=1 #dispersal sigma
#' phi=0.7
#' gamma=0.3
#' buff=2
#' X=list(expand.grid(4:11,4:11),expand.grid(4:11,4:11),expand.grid(4:11,4:11))
#' K=c(10,10,10)
#' M=125
#' ACtype="markov2"
#' obstype="poisson"
#' dSS=expand.grid(seq(2,13,0.25),seq(2,13,0.25)) #here is the discrete state space
#' plot(dSS)
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,sigma_t=sigma_t,K=K,X=X,t=t,M=M,buff=buff,
#'                ACtype=ACtype,obstype=obstype,dSS=dSS)
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M,sigma_t=sigma_t)
#' niter=1000
#' nburn=0
#' nthin=1
#' proppars=list(lam0=0.01,sigma=0.02,gamma=0.1,s2x=0.2,s2y=0.2,sigma_t=0.09)
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,
#'                 proppars=proppars,ACtype=ACtype,obstype=obstype,dSS=dSS)
#' plot(mcmc(out$out))
#'
#' #Here is the second example with t=3, but the population is not observed in year 2
#' t=3
#' N=100
#' p0=0.5
#' lam0=-log(1-p0)
#' sigma=0.750
#' phi=0.7
#' gamma=0.3
#' buff=3
#' #Note, there are no traps or K on occasion 2
#' X=list(expand.grid(4:9,4:9),NA,expand.grid(4:9,4:9))
#' K=c(10,NA,10)
#' M=425
#' primary=c(1,0,1) #primary period 2 not observed
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,
#' obstype="bernoulli",primary=primary)
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M)
#' niter=1000
#' nburn=0
#' nthin=1
#' proppars=list(lam0=0.025,sigma=0.025,gamma=0.1,phi=0.1,s1x=0.2,s1y=0.2) #Need 1 more propz
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars)
#' plot(mcmc(out$out))
#' #Note that if lam0 and/or sigma are primary period-specific, you still need to enter starting values
#' and tuning parameters, but these parameters are not updated. You can ignore them in the MCMC output.
#'}
#'@export

mcmc.OpenSCR <-
  function(data,niter=1000,nburn=0, nthin=1, K=NA,M = NA, inits=NA,proppars=NA,jointZ=TRUE,storeLatent=TRUE,Rcpp=TRUE,ACtype="fixed",obstype="bernoulli",dSS=NA,dualACup=FALSE){
    if(Rcpp==TRUE){ #Do we use Rcpp?
      out2=SCRmcmcOpenRcpp(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,
                           jointZ=jointZ,ACtype=ACtype,obstype=obstype,dSS=dSS,dualACup=dualACup,storeLatent=storeLatent)
    }else{#Don't use Rcpp
      out2=SCRmcmcOpen(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,
                       jointZ=jointZ,ACtype=ACtype,obstype=obstype,dSS=dSS,dualACup=dualACup,storeLatent=storeLatent)
    }
    if(storeLatent==TRUE){
      if(ACtype%in%c("markov","markov2","independent")){
        list(out=out2$out, s2xout=out2$s2xout, s2yout=out2$s2yout, zout=out2$zout,dSS=dSS)
      }else if(ACtype%in%c("metamu","metamu2")){
        list(out=out2$out, s1xout=out2$s1xout, s1yout=out2$s1yout,s2xout=out2$s2xout, s2yout=out2$s2yout, zout=out2$zout,dSS=dSS)
      }else{
        list(out=out2$out, s1xout=out2$s1xout, s1yout=out2$s1yout, zout=out2$zout,dSS=dSS)
      }
    }else{
      list(out=out2$out,dSS=dSS)
    }
  }


