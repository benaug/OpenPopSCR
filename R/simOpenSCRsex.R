#' Simulate data from a Open population SCR study with sex (or other group of size 2) differences.
#' @param N a 1 x 2 vector or t x 2 matrix indicating the number of individuals to simulate. Column 1 is males, column 2 is females.
#' If only the first year is entered, provide a gamma to determine N in subsequent years.  Otherwise, size is t, the number of years.
#' @param lam0 a vector containing the detection function expected number of captures at distance 0. If size 1, no sex difference is assumed.
#' Otherwise, sexes differ.
#' @param sigma a vector containing the spatial scale parameters. If size 1, no sex difference is assumed.
#' Otherwise, sexes differ.  Do not enter a gamma if N for all years specified.
#' @param phi a vector containing the survival rates. If size 1, no sex difference is assumed.
#' Otherwise, sexes differ.
#' @param gamma a vector containing the per capita recruitment rates for each sex.  If size 1, gamma is the same for each sex.
#' Note, this is the number of males recruited per N and number of females recruited per N and each will be smaller than
#' the number of individuals recruited per N.
#' @param K  a vector containing the number of capture occasions in each year
#' @param X a list of trap locations in each year.  Each list element is a J[l] x 2 matrix of trap locations, with J[l] being
#' the number of traps in each year.
#' @param psex the probability a row of z will be male.  This does not determine the sex ratio, but may be needed to
#' simulate data from very sex-skewed populations without increasing M.
#' @param pIDsex the probability that you can ID the sex of a captured individual
#' @param buff the distance to buffer the trapping array in the X and Y dimensions to produce the state space
#' @param obstype observation type, either "bernoulli" or "poisson"
#' @param ACtype Type of activity centers.  "fixed" don't move between years, "metamu" assume a bivariate normal distribution
#' with a meta mu and sigma_t with yearly activity centers required to stay inside the state space , "metamu2", is
#' the same as "metamu" except only meta mus are required to stay inside the state space.  markov" assumes activity
#' centers in year t+1 is a bivariate normal draw centered around the activity center in year t (but must stay within
#' the state space), and "independent" assumes animals randomly mix between years.
#' @param sigma_t a numeric indicating the between year spatial scale parameter for ACtypes "metamu" "metamu2", and "markov".  If size 1, no sex difference is assumed.
#' Otherwise, sexes differ.
#' @param M an integer indicating the level of data augmentation to use during simulation.
#' @return a list containing the capture history, activity centers, trap object, and several other data objects and summaries.
#' @description This function simulates data from an open population SCR model.
#' @author Ben Augustine
#' @export

simOpenSCRsex <-
  function(N=cbind(c(30,30),c(35,45,55)),psex=0.5,pIDsex=0.75,gamma=NULL,gamma.sex=TRUE,phi=rep(0.8,2),lam0=rep(0.2,2),sigma=rep(0.50,2),K=rep(10,3),X=X,t=3,M=M,sigma_t=NULL,buff=3,
           obstype="bernoulli",ACtype="fixed",vertices=NA,maxprop=10000){
    #Check for user errors
    if(!is.matrix(N)==1&is.null(gamma)){
      stop("Must provide gamma if N is a vector")
    }
    if(is.matrix(N)){
      if(nrow(N)==t&!is.null(gamma)){
        stop("Do not provide gamma if nrow(N)=t")
      }
    }
    if(length(K)!=t){
      stop("Must supply a K for each year")
    }
    if(length(X)!=t){
      stop("Must supply a X (trap locations) for each year")
    }
    if((ACtype=="metamu"|ACtype=="metamu2"|ACtype=="markov")&is.null(sigma_t)){
      stop("If ACtype is metamu, metamu2, or markov, must specify sigma_t")
    }
    if(!(ACtype=="metamu"|ACtype=="metamu2"|ACtype=="markov")&!is.null(sigma_t)){
      stop("ACtype must be metamu, metamu2, or markov when inputting a sigma_t")
    }
    storeparms=list(N=N,gamma=gamma,lam0=lam0,sigma=sigma,phi=phi)

    J=unlist(lapply(X,nrow))
    maxJ=max(J)
    #Get state space extent such that buffer is at least buff on all grids
    #minmax=rbind(apply(X[[1]],2,min),apply(X[[1]],2,max))
    useverts=!is.na(vertices)[1]
    if(useverts){
      if(!is.list(vertices)){
        stop("vertices must be a list")
      }
    }
    if(!useverts){
      minmax=array(NA,dim=c(length(X),2,2))
      for(i in 1:length(X)){
        minmax[i,,]=rbind(apply(X[[i]],2,min),apply(X[[i]],2,max))
      }
      xylim=rbind(apply(minmax,3,min),apply(minmax,3,max))
      xlim=xylim[,1]
      ylim=xylim[,2]
    }else{
      xlim=c(min(unlist(lapply(vertices,function(x){min(x[,1])}))),max(unlist(lapply(vertices,function(x){max(x[,1])}))))
      ylim=c(min(unlist(lapply(vertices,function(x){min(x[,2])}))),max(unlist(lapply(vertices,function(x){max(x[,2])}))))
    }
    s=array(NA,dim=c(M,t,2))
    D=lamd=array(NA,dim=c(M,maxJ,t))
    sex=rbinom(M,1,psex)+1
    if(length(lam0)==1){
      lam0=rep(lam0,2)
    }
    if(length(sigma)==1){
      sigma=rep(sigma,2)
    }
    if(length(sigma_t)==1){
      sigma_t=rep(sigma_t,2)
    }

    if(ACtype=="fixed"){#no movement
      if(useverts){
        mu<- cbind(runif(M, xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
        for(i in 1:M){
          inside=any(unlist(lapply(vertices,function(x){inout(mu[i,],x)})))
          while(inside==FALSE){
            mu[i,]=c(runif(1, xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside=any(unlist(lapply(vertices,function(x){inout(mu[i,],x)})))
          }
        }
      }else{
        mu<- cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff))
      }

      for(i in 1:t){
        s[,i,]=mu
        D[,1:nrow(X[[i]]),i]=e2dist(s[,i,],X[[i]])
        lamd[,,i]=lam0[sex]*exp(-D[,,i]^2/(2*sigma[sex]*sigma[sex]))
      }
    }else if(ACtype=="metamu"){
      if(useverts){
        mu<- cbind(runif(M, xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
        for(i in 1:M){
          inside=any(unlist(lapply(vertices,function(x){inout(mu[i,],x)})))
          while(inside==FALSE){
            mu[i,]=c(runif(1, xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside=any(unlist(lapply(vertices,function(x){inout(mu[i,],x)})))
          }
        }
      }else{
        mu<- cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff))
      }
      for(i in 1:t){#meta mu movement
        countout=0
        for(j in 1:M){
          out=1
          while(out==1){
            s[j,i,1]=rnorm(1,mu[j,1],sigma_t[sex[j]])
            s[j,i,2]=rnorm(1,mu[j,2],sigma_t[sex[j]])
            if(useverts){
              inside=any(unlist(lapply(vertices,function(x){inout(s[j,i,],x)})))
            }else{
              inside=s[j,i,1] < xlim[2]+buff & s[j,i,1] > xlim[1]-buff & s[j,i,2] < ylim[2]+buff & s[j,i,2] > ylim[1]-buff
            }
            if(inside){
              out=0
            }
            countout=countout+1
            if(countout==maxprop){#So you don't end up in an infinite loop if you set sigma_t too large relative to state space size
              stop("ACs proposed outside of state space exceeds maxprop. Should probably reconsider settings.")
            }
          }
        }
        D[,1:nrow(X[[i]]),i]=e2dist(s[,i,],X[[i]])
        lamd[,,i]=lam0[sex]*exp(-D[,,i]^2/(2*sigma[sex]*sigma[sex]))
      }
    }else if(ACtype=="metamu2"){
      if(useverts){
        mu<- cbind(runif(M, xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
        for(i in 1:M){
          inside=any(unlist(lapply(vertices,function(x){inout(mu[i,],x)})))
          while(inside==FALSE){
            mu[i,]=c(runif(1, xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside=any(unlist(lapply(vertices,function(x){inout(mu[i,],x)})))
          }
        }
      }else{
        mu<- cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff))
      }
      for(i in 1:t){#meta mu movement, only meta mus stay in SS
        for(j in 1:M){
          s[j,i,1]=rnorm(1,mu[j,1],sigma_t[sex[j]])
          s[j,i,2]=rnorm(1,mu[j,2],sigma_t[sex[j]])
        }
        D[,1:nrow(X[[i]]),i]=e2dist(s[,i,],X[[i]])

        lamd[,,i]=lam0[sex]*exp(-D[,,i]^2/(2*sigma[sex]*sigma[sex]))

      }

    }else if(ACtype=="markov"){
      if(useverts){
        s[,1,]=cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff)) #initial locations
        for(i in 1:M){
          inside=any(unlist(lapply(vertices,function(x){inout(mu[i,],x)})))
          while(inside==FALSE){
            s[,i,]=c(runif(1, xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside=any(unlist(lapply(vertices,function(x){inout(s[,i,],x)})))
          }
        }
      }else{
        s[,1,]=cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff)) #initial locations
      }
      D[,1:nrow(X[[1]]),1]=e2dist(s[,1,],X[[1]])
      lamd[,,1]=lam0[sex]*exp(-D[,,1]^2/(2*sigma[sex]*sigma[sex]))
      for(i in 2:t){
        countout=0
        for(j in 1:M){
          out=1
          while(out==1){
            s[j,i,1]=rnorm(1,s[j,i-1,1],sigma_t[sex[j]])
            s[j,i,2]=rnorm(1,s[j,i-1,2],sigma_t[sex[j]])
            if(useverts){
              inside=any(unlist(lapply(vertices,function(x){inout(s[j,i,],x)})))
            }else{
              inside=s[j,i,1] < xlim[2]+buff & s[j,i,1] > xlim[1]-buff & s[j,i,2] < ylim[2]+buff & s[j,i,2] > ylim[1]-buff
            }
            if(inside){
              out=0
            }
            countout=countout+1
            if(countout==maxprop){#So you don't end up in an infinite loop if you set sigma_t too large relative to state space size
              stop("ACs proposed outside of state space exceeds maxprop. Should probably reconsider settings.")
            }
          }
        }
        D[,1:nrow(X[[i]]),i]=e2dist(s[,i,],X[[i]])
        lamd[,,i]=lam0[sex[i]]*exp(-D[,,i]^2/(2*sigma[sex[i]]*sigma[sex[i]]))
      }
    }else if(ACtype=="independent"){
      for(i in 1:t){
        s[,i,]=cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff))
        if(useverts){
          countout=0
          for(j in 1:M){
            inside=any(unlist(lapply(vertices,function(x){inout(s[j,i,],x)})))

            while(inside==FALSE){
              s[j,i,]=c(runif(1, xlim[1]-buff,xlim[2]+buff), runif(1,ylim[1]-buff,ylim[2]+buff))
              inside=any(unlist(lapply(vertices,function(x){inout(s[j,i,],x)})))
              countout=countout+1
              if(countout==maxprop){#So you don't end up in an infinite loop if you set sigma_t too large relative to state space size
                stop("ACs proposed outside of state space exceeds maxprop. Should probably reconsider settings.")
              }
            }
          }
        }
        D[,1:nrow(X[[i]]),i]=e2dist(s[,i,],X[[i]])
        lamd[,,i]=lam0[sex]*exp(-D[,,i]^2/(2*sigma[sex]*sigma[sex]))
      }

    }else{
      stop("ACtype not recognized")
    }
    psi=N[1]/M
    #Calculate (per capita) gamma or N
    if(length(phi)==1){
      phi=rep(phi,2)
    }

    if(is.null(gamma)){
      gamma=matrix(NA,nrow=t-1,ncol=2)
      for(i in 1:2){
        for(l in 2:t){
          gamma[l-1,i]=(N[l,i]-phi[i]*N[l-1,i])/(N[l-1,i])
        }
      }
      storeparms$gamma=gamma
    }else{
      if(length(gamma)==1){
        gamma=rep(gamma,2)
      }
      N=rbind(N,matrix(NA,ncol=2,nrow=t-1))
      for(i in 1:2){
        for(l in 2:t){
          N[l,i]=round(N[l-1,i]*phi[i]+N[l-1,i]*gamma[i])
        }
      }
      gamma=cbind(rep(gamma[1],t-1),rep(gamma[2],t-1))
    }
    #####Population Dynamics#############  sex specific!
    z=a=matrix(0,M,t)
    z[sex==1&cumsum(sex==1)<=N[1,1],1]=1
    z[sex==2&cumsum(sex==2)<=N[1,2],1]=1
    a[,1]= 1-z[,1] # Available to be recruited?
    gamma.prime.M=gamma.prime.F=rep(NA,t-1)
    for(i in 2:t) {
      ER.M <- sum(z[,i-1]==1)*gamma[i-1,1] # Expected number of recruits, not sex-specific
      ER.F <- sum(z[,i-1]==1)*gamma[i-1,2]
      A.M <- sum(a[,i-1]==1&sex==1) # nAvailable to be recruited
      A.F <- sum(a[,i-1]==1&sex==2)
      if(ER.M>A.M){
        stop("M is too low. There aren't any males left to be recruited")
      }
      if(ER.F>A.F){
        stop("M is too low. There aren't any females left to be recruited")
      }
      gamma.prime.M[i-1] <- ER.M/A.M # individual-level recruitment *probability*
      gamma.prime.F[i-1] <- ER.F/A.F
      Ez.M <- 1*(z[,i-1]==1&sex==1)*phi[1] + 1*(a[,i-1]==1&sex==1)*gamma.prime.M[i-1]
      Ez.F <- 1*(z[,i-1]==1&sex==2)*phi[2] + 1*(a[,i-1]==1&sex==2)*gamma.prime.F[i-1]
      z[,i] <- rbinom(M, 1, Ez.M)+rbinom(M, 1, Ez.F)
      a[,i] <- apply(z[,1:i]<1, 1, all)
    }
    RN=cbind(colSums(z==1&sex==1),colSums(z==1&sex==2))

    #######Capture process######################
    # Simulate encounter history
    y <-array(0,dim=c(M,maxJ,t))
    if(obstype=="bernoulli"){
      pd=1-exp(-lamd)
      for(l in 1:t){
        for(i in 1:M){
          for(j in 1:J[l]){
            y[i,j,l]=rbinom(1,K[l],pd[i,j,l]*z[i,l])
          }
        }
      }
    }else if(obstype=="poisson"){
      for(l in 1:t){
        for(i in 1:M){
          for(j in 1:J[l]){
            y[i,j,l]=rpois(1,K[l]*lamd[i,j,l]*z[i,l])
          }
        }
      }
    }else{
      stop("observation model not recognized")
    }
    sfull=s
    yfull=y
    if(ACtype%in%c("metamu","metamu2")){
      mufull=mu
    }
    caps=apply(y,1,sum)
    idx=order(caps,decreasing=TRUE)
    y=y[idx,,]
    s=s[idx,,]
    if(ACtype%in%c("metamu","metamu2")){
      mu=mu[idx,]
    }
    sex=sex[idx]
    keep=which(rowSums(y)>0)
    y=y[keep,,]
    s=s[keep,,]
    if(ACtype%in%c("metamu","metamu2")){
      mu=mu[keep,]
    }
    sex=sex[keep]
    sexID=sex
    sexID[rbinom(length(sex),1,pIDsex)==0]=NA
    n=sum(caps>0)
    caps2d=apply(y,c(1,3),sum)
    n2d=colSums(caps2d>0)
    if(length(storeparms$gamma)==1){
      gamma=gamma[1]
    }else{
      gamma=gamma
    }

    if(!ACtype%in%c("metamu","metamu2")){
      if(!missing(vertices)){
        out<-list(y=y,s=s,yfull=yfull,sfull=sfull,X=X,K=K,n=n,n2d=n2d,buff=buff,J=J,EN=N,N=RN,z=z,gamma=gamma,
                  phi=storeparms$phi,obstype=obstype,ACtype=ACtype,vertices=vertices,sex=sex,sexID=sexID)
      }else{
        out<-list(y=y,s=s,yfull=yfull,sfull=sfull,X=X,K=K,n=n,n2d=n2d,buff=buff,J=J,EN=N,N=RN,z=z,gamma=gamma,
                  phi=storeparms$phi,obstype=obstype,ACtype=ACtype,sexfull=sex,sex=sexID)
      }
    }else{
      if(!missing(vertices)){
        out<-list(y=y,mu=mu,s=s,yfull=yfull,sfull=sfull,mufull=mufull,X=X,K=K,n=n,n2d=n2d,buff=buff,J=J,EN=N,N=RN,z=z,gamma=gamma,
                  phi=storeparms$phi,obstype=obstype,ACtype=ACtype,vertices=vertices,sexfull=sex,sex=sexID)
      }else{
        out<-list(y=y,mu=mu,s=s,yfull=yfull,sfull=sfull,mufull=mufull,X=X,K=K,n=n,n2d=n2d,buff=buff,J=J,EN=N,N=RN,z=z,gamma=gamma,
                  phi=storeparms$phi,obstype=obstype,ACtype=ACtype,sexfull=sex,sex=sexID)
      }
    }
    return(out)
  }
