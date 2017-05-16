#' Run MCMC algorithm for Open population SCR model.
#' @param data a list produced by simOpenSCR or in the same format
#' @param niter number of MCMC iterations to run
#' @param  nburn number of MCMC iterations to discard as burn in
#' @param nthin MCMC thinning parameter. Record output on every nthin iterations.  nthin=1 corresponds to no thinning
#' @param M The size of the augmented superpopulation
#' @param inits a list of user-supplied initial values.
#' @param proppars a list of tuning parameters for the proposal distributions
#' @param keepACs a logical indicating whether or not to keep the posteriors for z, s
#' @param jointZ a logical indicating whether you want to use the sequential or joint z update.
#' @param Rcpp a logical indicating whether or not to use Rcpp
#' @param ACtype a character indicating the model for between year activity center dependence. 'fixed' assumes activity centers are uniformly
#' distributed across the state space and do not move between years.  'independent" assumes activity centers for each year are uniformly
#' distributed across the state space, implying no dependence of activity centers between years. "markov" assumes activity centers in the
#' first year are uniformly distributed across the state space and activity centers in year l is a bivariate normal draw from N(s_{i,l-1},sigma_t).
#' "metamu" assumes animals have meta activity centers which are distributed uniformly across the landscape and the realized yearly activity
#' centers are a bivariate normal draw from N(mu_i,sigma_t), but must stay in the state space.  "metamu2" enforces the metamus to stay in
#' the state space, but the yearly ACs may leave.
#' @param obstype a character indicating the observation model "bernoulli" or "poisson"
#' @param dSS a discrete state space that overrules the buff or vertices objects in "data".  A matrix with columns for x and y locations
#'
#' @return  a list with the posteriors for the open population SCR parameters (out), s (s1xout,s1yout,s2xout,s2yout with
#' s1 being meta ACs and s2 being yearly ACs), and z.  s1x and yout are of dimension niter x M and s2x and yout and z are
#' of dimension niter x M x T
#' @author Ben Augustine, Richard Chandler
#' @description This function runs the MCMC algorithm for an open population SCR model.  The data list should have the following elements:
#' 1.  y, a n x J x T capture history where J is the maximum number of traps across years, T is the number of years, and n
#' is the number of animals captured
#' 2.  X,  a list with elements that consists of a matrix with the X and Y trap locations in the first two columns for each year.
#' If the traps do not vary across years, just repeat the same traps in each year
#' 3.  K, a vector of size T indicating how many capture occasions there were in each year
#' 4.  J, a vector of size T indicating how many traps there were in each year
#' 5. either buff or vertices.  buff is the fixed buffer for the traps to produce the state space.  It is applied to the minimum and maximum
#' X and Y trap locations across years, producing a square or rectangular state space.  vertices is a *list* of matrices with the X and Y coordinates
#' of a polygonal state space with one polygon in each list element.  If there is just one polygon, the list is of length 1.
#' If there are many polygons separated by large distances, you should think about the implications of activity centers
#' possibly being stuck inside the polygons.
#' 6. tf is an optional list of vectors of length T containing the trap operation information.  Each vector has one element for each trap
#' and indicates how many occasions each trap was operational.
#'
#' inits sets the initial values and determines if parameters are fixed or year-specific. It must have elements "lam0"
#' "sigma", "gamma", "phi", and "psi".  If there is an element "sigma_t", the parameters of a bivariate normal or Markov mobile
#' activity center model will be estimated.  If length(lam0)=1, a single lam0 will be estimated while if length(lam0)=t,
#' lam0 will be year-specific.  This goes for sigma, gamma, and phi as well, except for gamma and phi length needs to be t-1
#' or year-specific paramters.  A decent starting value for psi is (hypothesized) N/M.
#'
#' proppars is a list containing the tuning parameters for the parameters that use a Metropolis-Hastings update.
#' It must have elemnts "lam0", "sigma", "gamma", "s2x", "s2y", and "propz" (if jointZ=FALSE). If a parameter is year-specific, it needs
#' the appropriate number of proppars. propz is the number of data augmentation z's to update in years 2,...,t, so it should
#' be of length t-1.  Increasing propz improves mixing (up to a point) but increases computation time. Finally, if you set
#' an initial value for sigma_t, you need to provide proppars for "s1x", "s1y", and "sigma_t".
#'
#' A note on the z samplers.  jointZ=TRUE will update all the z's for each individual at the same time while jointZ=FALSE
#' will update them sequentially.  For T=3-6ish, you get a greater effective sample size per unit time with the joint update than with the sequential update.
#' The joint update always mixes better, but takes longer as t increases.
#'
#' A note on the activity center models. I think little is known about data requirements for the more complex models at this point or
#' how consequential it is to hold them fixed when they do indeed move between years.  Anectodally, the "metamu" and "markov" models
#' underestimate sigma_t because they force yearly activity centers to remain in the state space.  "metamu2" lets the yearly
#' ACs leave the state space, but they are still associated with the metamus that stay in the state space.  One could let the ACs
#' leave the state space in the Markov model, but then density would decrease through time, which seems problematic.  The
#' metamu2 model seems most sensible to me.  Email me if you have any other sensible models.
#'
#' A final note on splitting data sets to get group-specific parameters (e.g. sex).  This will alter the interpretation
#' of per capita recruitment.  For example, if you run sexes separetly, you are estimating recruitment per number
#' of males or females in the population, which probably does not make sense.  See the mcmc.OpenSCR() function for
#' sex-specific parameters.
#'
#' @examples
#' \dontrun{
#' library(coda)
#' #2 years of data,  all parameters fixed, stationary activity centers
#' t=2
#' N=100
#' p0=0.5
#' lam0=-log(1-p0)
#' sigma=0.750
#' phi=0.8
#' gamma=0.3
#' buff=2
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9)) #Trap XY can vary across years
#' K=c(10,10) #Number of occasions can vary across years as well
#' M=250
#' #Simulate some data
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,obstype="bernoulli")
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M) #Single initial values for a parameter will estimate a single parameter across years
#' niter=1000
#' nburn=0
#' nthin=1
#' proppars=list(lam0=0.025,sigma=0.025,gamma=0.1,phi=0.1,s2x=0.2,s2y=0.2,propz=c(30)) #The number of proppars must match the number of initial values
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,Rcpp=TRUE)
#' plot(mcmc(out$out))
#'
#' #same but t=3
#' t=3
#' N=100
#' p0=0.5
#' lam0=-log(1-p0)
#' sigma=0.750
#' phi=0.7
#' gamma=0.3
#' buff=2
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9)) #Note we need 3 trap objects stuffed into the trap list
#' K=c(10,10,10) #and 3 numbers of occasions within year
#' M=425
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,obstype="bernoulli")
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M)
#' niter=1000
#' nburn=0
#' nthin=1
#' proppars=list(lam0=0.025,sigma=0.025,gamma=0.1,phi=0.1,s2x=0.2,s2y=0.2,propz=c(30,30)) #Need 1 more propz
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,Rcpp=TRUE)
#' plot(mcmc(out$out))
#' summary(mcmc(out$out))
#'
#'#Now detection function varies by year
#' t=3
#' N=100
#' p0=c(0.4,0.5,0.6) #need 3 of these
#' lam0=-log(1-p0)
#' sigma=c(0.75,0.65,0.55) #these, too
#' phi=0.7
#' gamma=0.3
#' buff=2
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9)) #Note we need 3 trap objects stuffed into the trap list
#' K=c(10,10,10) #and 3 numbers of occasions within year
#' M=425
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,obstype="bernoulli")
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M) #Note now 3 initial values are specivied here
#' niter=1000
#' nburn=0
#' nthin=1
#' proppars=list(lam0=rep(0.025,3),sigma=rep(0.025,3),gamma=0.1,phi=0.1,s2x=0.2,s2y=0.2,propz=c(30,30)) #Note we need 3 proppars for lam0 and sigma
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,Rcpp=TRUE)
#' plot(mcmc(out$out))
#' summary(mcmc(out$out))
#'
#' #Let's set lam0 and sigma back to fixed and vary gamma and phi
#' t=3
#' N=100
#' p0=0.5
#' lam0=-log(1-p0)
#' sigma=0.750
#' phi=c(0.7,0.8)
#' gamma=c(0.2,0.3)
#' buff=2
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9)) #Note we need 3 trap objects stuffed into the trap list
#' K=c(10,10,10) #and 3 numbers of occasions within year
#' M=425 #Not sure if this is still enough...
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,obstype="bernoulli")
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M) #now phi and gamma are getting 2 inits
#' niter=1000
#' nburn=0
#' nthin=1
#' proppars=list(lam0=0.025,sigma=0.025,gamma=c(0.1,0.1),phi=c(0.1,0.1),s2x=0.2,s2y=0.2,propz=c(30,30)) #and phi and gamma get 2 proppars
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,Rcpp=TRUE)
#' plot(mcmc(out$out))
#' summary(mcmc(out$out))
#'
#' #Let's set everything back to fixed and do mobile activty centers between years
#' #We're using a bivariate normal "metamu2" model.
#' t=3
#' N=100
#' p0=0.5
#' lam0=-log(1-p0)
#' sigma=0.750
#' sigma_t=0.7
#' #We'll use this to simulate mobile activity centers
#' ACtype="metamu2"
#' phi=0.7
#' gamma=0.3
#' buff=4 #increase buffer to not constrain movement as much
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9))
#' K=c(10,10,10) #and 3 numbers of occasions within year
#' M=425
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,sigma_t=sigma_t,K=K,X=X,t=t,M=M,buff=buff,ACtype=ACtype)#Note, we're putting sigma_t in now
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M,sigma_t=sigma_t)#Need an init for sigma_t
#' niter=1000
#' nburn=0
#' nthin=1
#' proppars=list(lam0=0.025,sigma=0.025,gamma=0.1,phi=0.1,s1x=0.05,s1y=0.05,s2x=0.2,s2y=0.2,propz=c(30,30),sigma_t=0.025) #Note there is a proppar for sigma_t and the meta mus
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,ACtype=ACtype)
#' plot(mcmc(out$out))
#' summary(mcmc(out$out))
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
#' #for lam0, sigma, gamma, and sigma_t.  Here, we're accepting too much, so we should make the proppars
#' #larger.  But what about the activity centers?  Here, we look at the acceptance rates for the yearly
#' #ACs in year 1.  Note, they will be identical for the x and y dimension
#' 1-rejectionRate(mcmc(out$s2xout[,,1]))
#' #It looks like we need to raise the proppars for the yearly ACs. Now for the meta mus.
#' 1-rejectionRate(mcmc(out$s1xout))
#' #We definitely need to raise their proppars! So we can do short runs (say of length 500) starting at
#' #inits close to the converged values to assess and adjust the tuning parameters until we're seeing
#' #acceptance rates near 23% and then run a chain long enough to get an acceptable effective
#' #sample size for all parameters 9(but usually the N's are the limiting factor).  If at any point a
#' #parameter is not updating, the tuning parameter is probably too large.
#' summary(mcmc(out$out))
#'
#' #So which point estimates and credible intervals should we use?  I've found the posterior mode
#' #to be the closest to unbiased point estimator, with the posterior median being pretty good, too.
#' #I would not use the posterior mean, especially for N. Here is a function to get the mode.
#' library(mcmcGLMM)
#' #posterior.mode(mcmc(out$out))
#' #You can get the posterior means, modes, and quantile credible intervals using coda's summary
#' summary(mcmc(out$out))
#' #Note the 2.5% and 97.5% quantiles can be used as the credible interval limits.  However, I've found
#' #the highest posteior density intervals to have slightly better frequentist coverage. coda will do that for you
#' HPDinterval(mcmc(out$out))
#'}
#'@export

mcmc.OpenSCR <-
  function(data,niter=1000,nburn=0, nthin=1, K=NA,M = NA, inits=NA,proppars=NA,jointZ=TRUE,keepACs=TRUE,Rcpp=TRUE,ACtype="fixed",obstype="bernoulli",dSS=NA){
    if(Rcpp==TRUE){ #Do we use Rcpp?
      out2=SCRmcmcOpenRcpp(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,jointZ=jointZ,ACtype=ACtype,obstype=obstype,dSS=dSS)
    }else{#Don't use Rcpp
      out2=SCRmcmcOpen(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,jointZ=jointZ,ACtype=ACtype,obstype=obstype,dSS=dSS)
    }
    if(keepACs==TRUE){
      if("s2xout"%in%names(out2)){
        list(out=out2$out, s1xout=out2$s1xout, s1yout=out2$s1yout,s2xout=out2$s2xout, s2yout=out2$s2yout, zout=out2$zout,dSS=dSS)
      }else{
        list(out=out2$out, s1xout=out2$s1xout, s1yout=out2$s1yout, zout=out2$zout,dSS=dSS)
      }
    }else{
      list(out=out2$out)
    }
  }


