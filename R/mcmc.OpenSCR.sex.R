#' Run MCMC algorithm for fixed or sex-specific parameter open population SCR model.
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
#' @param ACtype a character string indicating the activity center model.  ACtypes can be divided into two types, those that operate on both
#' continuous and discrete state space models and those that operate only on
#' discrete state spaces.  First, those operating on both.  1.  "fixed": activity centers do not move between primary periods.
#' 2.  "independent": activity centers for each individual are independent between primary periods (e.g. random population mixing).
#' Second, continuous state space models:  3. "metamu": primary period activity centers follow a bivariate normal distribution
#' around a meta activity center with sigma_t determining the spread of primary period activity centers around the meta activity center.
#' Primary period activity centers are required to stay inside the state space  4. "metamu2": is the same as "metamu" except only the
#' meta activity centers are required to stay inside the state space.  The primary period ACs float around nearby the state space,
#' linked to their meta activity center.  5. "markov": activity centers in primary period 1 are spatially uniform.  Activity centers
#' in primary period l+1 are a bivariate normal draw centered around the activity center in primary period l
#'  (but must stay within the state space).
#'
#'  The "metamu" and "markov" models are truncated on rectangular state spaces.  This is required to estimate
#'  sigma_t and to a lesser degree, the Ns without bias.  If using a polygon state space, the "metamu2" model
#'  will provide unbiased estimates.  If using a non-rectangular state space with a Markov model, see the discrete
#'  state space model below.
#'
#' ACtypes operating only on discrete state spaces are currently limited to
#' "markov2":  Activity centers in primary period 1 are spatially uniform.  Activity centers in primary period
#'  l+1 follow an exponential dispersal kernel centered at the
#' activity center in primary period l; however the availability of dispersal distances in the state space are factored into the
#' dispersal decision.  This will be formally described elsewhere, but uses use vs. availability ideas to correct for
#' restricted availabilities of dispersal distances.
#' @param obstype a character indicating the observation model "bernoulli" or "poisson"
#' @param dualACup a logical to be used with discrete state spaces indicating whether a second, interpatch activity
#' center proposal is used.  Currently under development--probably should not use.
#' @param dSS an optional (N_SS x 2) matrix for a discrete state space locations that overrules the buff or vertices objects in "data".  A matrix with columns for x and y locations
#' @param primary an optional vector of length t with entries 1 if data was recorded in that primary period and 0 if not. This allows
#' population dynamics to occur at equal interval primary periods even if data was not recorded at each primary period.
#' If not entered, the population is assumed to have been sampled in all primary periods.
#' @return  a list with the posteriors for the open population SCR parameters (out), s (s1xout,s1yout,s2xout,s2yout with
#' s1 being meta ACs and s2 being primary period ACs), and z.  s1x and yout are of dimension niter x M and s2x and yout and z are
#' of dimension niter x M x T.  Posteriors for the sex of each individual could be returned--email me if you want this.
#' The "psex" parameter in out is the realized sex ratio in primary period 1.  You can compute realized and expected primary period sex
#' ratios by postprocessing the output.  Email me if you want to know more.
#' @author Ben Augustine
#' @description This function runs an MCMC algorithm for a sex-specific open population SCR model.  The data list should have the following elements:
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
#' 8. sex, a vector with 1, 2, and NA indicating male, female, and missing sex, respectively
#'
#' inits sets the initial values and determines if parameters are fixed or sex-specific. It must have elements "lam0"
#' "sigma", "gamma", "phi", and "psi".  If there is an element "sigma_t", the parameters of a "meta mu" or Markov mobile
#' activity center model will be estimated.  If length(lam0)=1, a single lam0 will be estimated while if length(lam0)=2,
#' lam0 will be sex-specific.  This goes for sigma, gamma, and phi as well.
#'
#' proppars is a list containing the tuning parameters for the parameters that use a Metropolis-Hastings update.
#' It must have elemnts "lam0", "sigma", "gamma", "s2x", "s2y", "propz" (if jointZ=FALSE), and "sex". If a parameter is sex-specific, it needs
#' two proppars. propz is the number of data augmentation z's to update in primary periods 2,...,t, so it should
#' be of length t-1.  Increasing propz improves mixing (up to a point) but increases computation time.  The sex element
#' determines how many latent sexes are updated on each iteration.  If you set an initial value for
#' sigma_t, you need to provide proppars for "sigma_t".   Finally, proppars$s1x
#' and proppars$s1y tune fixed activity centers or the meta mus and proppars$s2x and proppars$s2y tune primary period
#' activity centers.
#'
#' A note on the z samplers.  jointZ=TRUE will update all the z's for each individual at the same time while jointZ=FALSE
#' will update them sequentially.  For T=3-6ish, you get a greater effective sample size per unit time with the joint update than with the sequential update.
#' The joint update always mixes better, but takes longer as T increases.
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
#' A note on sex-specific gamma.  This algorithm estimates the number of males recruited per N and the number of females
#' recruited per N and both will be smaller than what is estimated by mcmc.OpenSCR(), the number of individuals
#' recruited per N.  An alternative model would be the number of each sex recruited per females in the population.
#'
#' Finally, a note on the speed of this sampler.  It is substantially slower than the year-specific sampler
#' because of the sex update that (potentially) depends on sex-specific lam0, sigma, gamma, phi, and sigma_t.
#' All of these parameters in combination with where and how many times an individual was captured inform
#' their sex.  The sex-specific per capita recruitment also adds some computational cost.
#'
#' @examples
#' \dontrun{
#' #Here is an example with all sex-specific parameters. To fix a parameter, just use 1 init and proppar
#' library(coda)
#' t=3
#' N=c(30,30) #order is M,F
#' p0=c(0.25,0.5)
#' lam0=-log(1-p0)
#' sigma=c(0.75,0.5)
#' sigma_t=c(0.5) #single sigma_t
#' phi=c(0.4,0.8)
#' gamma=c(0.4,0.1)
#' buff=3
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9))
#' K=c(10,10,10)
#' M=225
#' pIDsex=0.95 #probability you will ID a sex
#' obstype="bernoulli"
#' ACtype="metamu"
#' data=simOpenSCRsex(N=N,gamma=gamma,phi=phi,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,
#'                    ACtype=ACtype,obstype="bernoulli",pIDsex=pIDsex)
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=(sum(N)/M),psex=0.5,sigma_t=sigma_t)
#' niter=2000
#' nburn=0
#' nthin=1
#' proppars=list(lam0=c(0.075,0.115),sigma=c(0.055,0.045),gamma=c(0.115,0.085),s1x=0.4,s1y=0.4,
#' s2x=0.5,s2y=0.5,sigma_t=c(0.04),sex=100)
#' out=mcmc.OpenSCR.sex(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,ACtype=ACtype,
#'                    storeLatent=TRUE,jointZ=TRUE)
#' plot(mcmc(out$out))
#' #need to run longer for proper inference
#' see the help file of mcmc.OpenSCR for examples with other ACtypes and state space options that work the
#' same in this sampler and for some MCMC tips.
#'}
#'@export

mcmc.OpenSCR.sex <-
  function(data,niter=1000,nburn=0, nthin=1,M = NA, inits=NA,proppars=NA,jointZ=TRUE,storeLatent=TRUE,
           Rcpp=FALSE,ACtype="fixed",obstype="bernoulli",dSS=NA,dualACup=FALSE){
    if(Rcpp==TRUE){ #Do we use Rcpp?
      out2=SCRmcmcOpensexRcpp(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,
                              jointZ=jointZ,ACtype=ACtype,obstype=obstype,dSS=dSS,dualACup=dualACup,storeLatent=storeLatent)
    }else{#Don't use Rcpp
      out2=SCRmcmcOpensex(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,
                          jointZ=jointZ,ACtype=ACtype,obstype=obstype,dSS=dSS,dualACup=dualACup,storeLatent=storeLatent)
    }
    if(storeLatent==TRUE){
      if(ACtype%in%c("markov","markov2","independent")){
        list(out=out2$out, s2xout=out2$s2xout, s2yout=out2$s2yout, zout=out2$zout,sexout=out2$sexout,dSS=dSS)
      }else if(ACtype%in%c("metamu","metamu2")){
        list(out=out2$out, s1xout=out2$s1xout, s1yout=out2$s1yout,s2xout=out2$s2xout, s2yout=out2$s2yout,
             zout=out2$zout,sexout=out2$sexout,dSS=dSS)
      }else{
        list(out=out2$out, s1xout=out2$s1xout, s1yout=out2$s1yout, zout=out2$zout,sexout=out2$sexout,dSS=dSS)
      }
    }else{
      list(out=out2$out,dSS=dSS)
    }
  }


