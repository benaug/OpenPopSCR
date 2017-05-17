#' Run MCMC algorithm for sex-specific Open population SCR model.
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
#' @return  a list with the posteriors for the open population SCR parameters (out), s (s1xout,s1yout,s2xout,s2yout with
#' s1 being meta ACs and s2 being yearly ACs), and z.  s1x and yout are of dimension niter x M and s2x and yout and z are
#' of dimension niter x M x T.  Posteriors for the sex of each individual could be returned--email me if you want this.
#' The "psex" parameter in out is the realized sex ratio in year 1.  You can compute realized and expected yearly sex
#' ratios by postprocessing the output.  Email me if you want to know more.
#' @author Ben Augustine
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
#' 7. sex, a vector with 1, 2, and NA indicating male, female, and missing sex, respectively
#'
#' inits sets the initial values and determines if parameters are fixed or sex-specific. It must have elements "lam0"
#' "sigma", "gamma", "phi", and "psi".  If there is an element "sigma_t", the parameters of a bivariate normal or Markov mobile
#' activity center model will be estimated.  If length(lam0)=1, a single lam0 will be estimated while if length(lam0)=2,
#' lam0 will be sex-specific.  This goes for sigma, gamma, and phi as well.
#'
#' proppars is a list containing the tuning parameters for the parameters that use a Metropolis-Hastings update.
#' It must have elemnts "lam0", "sigma", "gamma", "s2x", "s2y", "propz" (if jointZ=FALSE), and "sex". If a parameter is sex-specific, it needs
#' two proppars. propz is the number of data augmentation z's to update in years 2,...,t, so it should
#' be of length t-1.  Increasing propz improves mixing (up to a point) but increases computation time.  The sex element
#' determines how many latent sexes are updated on each iteration.   Finally, if you set an initial value for
#' sigma_t, you need to provide proppars for "s1x", "s1y", and "sigma_t".
#'
#' A note on the z samplers.  jointZ=TRUE will update all the z's for each individual at the same time while jointZ=FALSE
#' will update them sequentially.  For T=3-6ish, you get a greater effective sample size per unit time with the joint update than with the sequential update.
#' The joint update always mixes better, but takes longer as t increases.
#'
#' A note on sex-specific gamma.  This algorithm estimates the number of males recruited per N and the number of females
#' recruited per N and both will be smaller than what is estimated by mcmc.OpenSCR(), the number of individuals
#' recruited per N.  An alternative model would be the number of each sex recruited per females in the population.
#'
#' @examples
#' \dontrun{
#' #Here is an example with all sex-specific parameters. To fix a parameter, just use 1 init and proppar
#' library(coda)
#' t=3
#' N=c(40,40) #order is M,F
#' p0=c(0.25,0.5)
#' lam0=-log(1-p0)
#' sigma=c(0.75,0.5)
#' sigma_t=c(0.5,0.25)
#' phi=c(0.4,0.8)
#' gamma=c(0.4,0.1)
#' buff=3
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9))
#' K=c(10,10,10)
#' M=225
#' psex=0.5 #this is used in the simulation, but does not determine the sex ratio.
#' pIDsex=0.95 #probability you will ID a sex
#' obstype="bernoulli"
#' ACtype="metamu2"
#' data=simOpenSCRsex(N=N,gamma=gamma,phi=phi,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,
#'                    ACtype=ACtype,obstype="bernoulli",pIDsex=pIDsex)
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=(sum(N)/M),psex=0.5,sigma_t=sigma_t)
#' niter=500
#' nburn=0
#' nthin=1
#' proppars=list(lam0=c(0.075,0.115),sigma=c(0.055,0.045),gamma=c(0.115,0.085),s1x=0.4,s1y=0.4,s2x=0.5,s2y=0.5,sigma_t=c(0.04,0.03),propz=c(20,20),sex=100)
#' out=mcmc.OpenSCR.sex(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,ACtype=ACtype,
#'                    keepACs=TRUE,jointZ=TRUE)
#' plot(mcmc(out$out))
#' see the help file of mcmc.OpenSCR for examples with other ACtypes and for some MCMC tips.
#'}
#'@export

mcmc.OpenSCR.sex <-
  function(data,niter=1000,nburn=0, nthin=1, K=NA,M = NA, inits=NA,proppars=NA,jointZ=TRUE,keepACs=TRUE,Rcpp=TRUE,ACtype="fixed",obstype="bernoulli",dSS=NA){
    if(Rcpp==TRUE){ #Do we use Rcpp?
      stop("Sorry, Rcpp isn't done, yet! =)")
      # out2=SCRmcmcOpenRcpp(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,jointZ=jointZ,ACtype=ACtype,obstype=obstype,dSS=dSS)
    }else{#Don't use Rcpp
      out2=SCRmcmcOpensex(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,jointZ=jointZ,ACtype=ACtype,obstype=obstype,dSS=dSS)
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

