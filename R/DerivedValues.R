#' Function to post-process MCMC results. Need to extend for model without sex.
#' @param Nm matrix containing the posteriors for the number of males
#' @param Nf matrix containing the posteriors for the number of females
#' @param gamma matrix of vector containing the posteriors for gamma(s)
#' @param phi matrix of vector containing the posteriors for phi(s)
#' @param psex vector containing the posterior for psex
#' @return a list containing realized and expected population growth rates, population sizes,
#' and sex ratios
#' @description This function derives realized and expected population growth rates, population sizes, and sex ratios
#' from the inputted MCMC chains.
#' @author Ben Augustine
#' @export

DerivedValues=function(Nm,Nf,gamma,phi,psex,psi,M){
  if(!is.matrix(gamma)){
    gamma=cbind(gamma,gamma)
  }
  if(!is.matrix(phi)){
    gamma=cbind(phi,phi)
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
  return(out=list(EN=EN,ENm=ENm,ENf=ENf,lambda=lambda,lambdaM=lambdaM,lambdaF=lambdaF,
                   Elambda=Elambda,ElambdaM=ElambdaM,ElambdaF=ElambdaF,sex=sex,Esex=Esex))
}

