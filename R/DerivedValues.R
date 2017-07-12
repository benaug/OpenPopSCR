#' Function to post-process MCMC results of the year-specific sampler
#' @param N matrix containing the posteriors for abundance
#' @param gamma matrix or vector containing the posteriors for gamma(s)
#' @param phi matrix or vector containing the posteriors for phi(s)
#' @param psi vector containing the posterior for psi
#' @param M an integer indicating the i dimension of the augmented data
#' @return a list containing realized and expected population growth rates and population sizes,
#' @description This function derives realized and expected population growth rates and population sizes
#' from the inputted MCMC chains.
#' @author Ben Augustine
#' @examples
#' \dontrun{
#' #Run a model
#' t=3
#' N=50
#' p0=0.5
#' lam0=-log(1-p0)
#' sigma=0.750
#' phi=0.7
#' gamma=0.3
#' buff=3
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9)) #Note we need 3 trap objects stuffed into the trap list
#' K=c(10,10,10) #and 3 numbers of occasions within primary period
#' M=225
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,obstype="bernoulli")
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M)
#' niter=1000
#' nburn=0
#' nthin=1
#' proppars=list(lam0=0.025,sigma=0.025,gamma=0.1,phi=0.1,s1x=0.2,s1y=0.2,propz=c(30,30)) #Need 1 more propz
#' out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,Rcpp=TRUE)
#' DV=DerivedValues(out$out[,5:7],out$out[,3],out$out[,4],out$out[,8],M)
#' str(DV)
#' library(coda)
#' plot(mcmc(DV$EN))
#' plot(mcmc(DV$lambda))
#' plot(mcmc(DV$Elambda))
#' #should discard burnin
#' }
#' @export

DerivedValues=function(N,gamma,phi,psi,M){
  t=ncol(N)
  if(!is.matrix(gamma)){
    gamma=matrix(rep(gamma,t-1),ncol=(t-1))
  }
  if(!is.matrix(phi)){
    phi=matrix(rep(phi,t-1),ncol=(t-1))
  }
  iters=nrow(N)
  EN=matrix(NA,ncol=t,nrow=iters)
  Elambda=matrix(NA,ncol=t-1,nrow=iters)
  lambda=matrix(NA,ncol=t-1,nrow=iters)
  EN[,1]=M*psi
  for(l in 2:t){
    EN[,l]=EN[,l-1]*(phi[,l-1]+gamma[,l-1])
    Elambda[,l-1]=EN[,l]/EN[,l-1]
    lambda[,l-1]=N[,l]/N[,l-1]
  }
  return(out=list(EN=EN,lambda=lambda,Elambda=Elambda))
}

