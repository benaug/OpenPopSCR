#' Function to post-process MCMC results for the sex-specific sampler
#' @param Nm matrix containing the posteriors for the number of males
#' @param Nf matrix containing the posteriors for the number of females
#' @param gamma matrix or vector (if sex-specific) containing the posteriors for gamma(s)
#' @param phi matrix or vector (if sex-specific) containing the posteriors for phi(s)
#' @param psex vector containing the posterior for psex
#' @param psi vector containing the posterior for psi
#' @param M an integer indicating the i dimension of the augmented data
#' @param extrap an optional integer indicating how many primary periods to extrapolate beyond the
#' observed primary periods for EN, ENm, and ENf. These are abundance projections from the population growth model
#' @return a list containing realized and expected population growth rates, population sizes,
#' and sex ratios
#' @description This function derives realized and expected population growth rates, population sizes, and sex ratios
#' from the inputted MCMC chains.
#' @author Ben Augustine
#' @examples
#' \dontrun{
#' #Run a sex-specific model
#' t=3
#' N=c(40,40)
#' p0=c(0.25,0.5)
#' lam0=-log(1-p0)
#' sigma=c(0.75,0.5)
#' phi=c(0.4,0.8)
#' gamma=c(0.4,0.1)
#' buff=3
#' X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9))
#' K=c(10,10,10)
#' M=250
#' pIDsex=0.75
#' obstype="bernoulli"
#' ACtype="metamu"
#' data=simOpenSCRsex(N=N,gamma=gamma,phi=phi,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,
#'                    ACtype=ACtype,obstype="bernoulli",pIDsex=pIDsex,sigma_t=sigma_t)
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=(sum(N)/M),psex=0.5)
#' niter=1000
#' nburn=0
#' nthin=1
#' proppars=list(lam0=c(0.075,0.115),sigma=c(0.055,0.045),gamma=c(0.115,0.085),s1x=0.4,s1y=0.4,
#' sigma_t=c(0.04,0.03),sex=100)
#' out=mcmc.OpenSCR.sex(data,niter=niter,nburn=nburn, nthin=nthin, M=M, inits=inits,proppars=proppars,
#' ACtype=ACtype)
#' DV=DerivedValuesSex(out$out[,12:14],out$out[,15:17],out$out[,5:6],out$out[,7:8],out$out[,18],out$out[,19],M)
#' str(DV)
#' library(coda)
#' plot(mcmc(DV$lambdaM))
#' plot(mcmc(DV$lambdaF))
#' #etc.
#' #should discard burnin
#' }
#' @export

DerivedValuesSex=function(Nm,Nf,gamma,phi,psex,psi,M,extrap=0){
  if(!is.matrix(gamma)){
    gamma=cbind(gamma,gamma)
  }
  if(!is.matrix(phi)){
    phi=cbind(phi,phi)
  }
  t=ncol(Nm)
  iters=nrow(Nm)
  ENm=ENf=EN=matrix(NA,ncol=t,nrow=iters)
  Elambda=ElambdaM=ElambdaF=matrix(NA,ncol=t-1,nrow=iters)
  lambda=lambdaM=lambdaF=matrix(NA,ncol=t-1,nrow=iters)
  Esex=sex=matrix(NA,ncol=t,nrow=iters)
  ENm[,1]=M*psi*(1-psex)
  ENf[,1]=M*psi*psex
  EN[,1]=ENm[,1]+ENf[,1]
  N=Nm+Nf
  Esex[,1]=psex
  sex[,1]=Nf[,1]/N[,1]
  for(l in 2:t){
    ENm[,l]=ENm[,l-1]*(phi[,1]+gamma[,1])+ENf[,l-1]*gamma[,1]
    ENf[,l]=ENf[,l-1]*(phi[,2]+gamma[,2])+ENm[,l-1]*gamma[,2]
    EN[,l]=ENm[,l]+ENf[,l]
    ElambdaM[,l-1]=ENm[,l]/ENm[,l-1]
    ElambdaF[,l-1]=ENf[,l]/ENf[,l-1]
    Elambda[,l-1]=EN[,l]/EN[,l-1]
    lambdaM[,l-1]=Nm[,l]/Nm[,l-1]
    lambdaF[,l-1]=Nf[,l]/Nf[,l-1]
    lambda[,l-1]=N[,l]/N[,l-1]
    Esex[,l]=ENf[,l]/EN[,l]
    sex[,l]=Nf[,l]/N[,l]
  }
  if(extrap>0){
    ENm2=ENf2=matrix(NA,ncol=extrap,nrow=iters)
    #extrap period 1
    ENm2[,1]=ENm[,t]*(phi[,1]+gamma[,1])+ENf[,t]*gamma[,1]
    ENf2[,1]=ENf[,t]*(phi[,2]+gamma[,2])+ENm[,t]*gamma[,2]
    if(extrap>1){
      for(l in 2:extrap){
        ENm2[,l]=ENm2[,l-1]*(phi[,1]+gamma[,1])+ENf2[,l-1]*gamma[,1]
        ENf2[,l]=ENf2[,l-1]*(phi[,2]+gamma[,2])+ENm2[,l-1]*gamma[,2]
      }
    }
    EN2=ENm2+ENf2
    EN=cbind(EN,EN2)
    ENm=cbind(ENm,ENm2)
    ENf=cbind(ENf,ENf2)
  }
    return(out=list(EN=EN,ENm=ENm,ENf=ENf,lambda=lambda,lambdaM=lambdaM,lambdaF=lambdaF,
                    Elambda=Elambda,ElambdaM=ElambdaM,ElambdaF=ElambdaF,sex=sex,Esex=Esex))

}

