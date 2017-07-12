e2dist<-function (x, y)
{
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

#' Simulate data from an open population SCR study with fixed or primary period-specific parameters.
#' @param N a vector indicating the number of individuals to simulate. If size 1, provide one or multiple gammas (see gamma below) to determine N in subsequent primary periods.
#' Otherwise, size is t, the number of primary periods.
#' @param lam0 a vector containing the detection function expected number of captures at distance 0. If size 1, a constant rate across primary periods is assumed.
#' Otherwise, size is t, the number of primary periods.
#' @param sigma a vector containing the detection function spatial scale parameter.  If size 1, sigma is fixed across primary periods.
#' #' Otherwise, size is t, the number of primary periods.
#' @param gamma a vector containing the per capita recruitment rates between primary periods.  If size 1, gamma is fixed across primary periods.
#' Otherwise, size is t-1.  Do not enter a gamma if N for all primary periods specified.
#' @param phi a vector containing the survival rates between primary periods. If size 1, phi is fixed across primary periods.
#' Otherwise, size is t-1.
#' @param K  a vector containing the number of capture occasions in each of the t primary periods
#' @param X a list of trap locations in each primary period.  Each list element is a J[l] x 2 matrix of trap locations, with J[l] being
#' the number of traps in primary period l.
#' @param buff the distance to buffer the trapping array in the X and Y dimensions to produce the state space in which the population lives
#' @param obstype observation type, either "bernoulli" or "poisson"
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
#' @param sigma_t a numeric indicating the between primary period spatial scale parameter for ACtypes "metamu" "metamu2", "markov" and "markov2"
#' This is the only parameter that is not primary period-specific. Very rich data sets would be required to estimate a year-specific sigma_t.
#' @param M an integer indicating the level of data augmentation to use during simulation.  This should be
#' larger than the total number of individuals ever alive.
#' @param vertices an optional list of polygon vertices to use for the state space.  Each list element should be a matrix with 2
#' columns, corresponding to the X and Y coordinates of polygon vertices.  The vertices must close the polygon (e.g. the first
#' and last vertices should be the same).  Cannot enter both vertices and dSS.
#' @param dSS an optional (N_SS x 2) matrix of discrete state space locations.  The matrix should have 2 columns corresponding to
#' the X and Y coordinates of each state space element.  Cannot enter both vertices and dSS.
#' @param primary a vector of length T with entries 1 if the population is to be observed in primary period l and 0 otherwise.
#' @return a list containing the capture history, activity centers, trap object, and several other data objects and summaries.
#' @description This function simulates data from an open population SCR model.
#' @author Ben Augustine
#' @export

simOpenSCR <-
  function(N=c(40,60,80),gamma=NULL,phi=rep(0.8,2),lam0=rep(0.2,3),sigma=rep(0.50,3),K=rep(10,3),X=X,t=3,M=M,sigma_t=NULL,buff=3,
           obstype="bernoulli",ACtype="fixed",vertices=NA,maxprop=10000,dSS=NA,primary=NA){
    #Check for user errors
    if(length(N)==1&is.null(gamma)){
      stop("Must provide gamma if length(N)==1")
    }
    if(length(N)==t&!is.null(gamma)){
      stop("Do not provide gamma if length(N)=t")
    }
    if(length(K)!=t){
      stop("Must supply a K for each year")
    }
    if(length(X)!=t){
      stop("Must supply a X (trap locations) for each year")
    }
    if((ACtype%in%c("metamu","metamu2","markov","markov2"))&is.null(sigma_t)){
      stop("If ACtype is metamu, metamu2, markov or markov2, must specify sigma_t")
    }
    if(!(ACtype%in%c("metamu","metamu2","markov","markov2"))&!is.null(sigma_t)){
      warning("Ignoring sigma_t because ACtype is not metamu, metamu2, markov, or markov2")
    }
    if(is.data.frame(dSS)){
      dSS=as.matrix(dSS)
    }
    if(ACtype=="markov2"){
      if(is.na(dSS[1])){
        stop("If ACtype is markov2, must input a discrete state space")
      }
      if(!is.na(vertices)){
        warning("Ignoring vertices since ACtype is markov2")
      }
    }
    if(is.na(primary[1])){
      primary=rep(1,t)
    }else{
      for(l in 1:t){
        if(primary[l]==0){
          X[[l]]=matrix(nrow=0,ncol=0)
        }
      }
      if(primary[1]==0){
        stop("first primary period must be a 1 (data recorded)")
      }
    }
    storeparms=list(N=N,gamma=gamma,lam0=lam0,sigma=sigma,phi=phi)
    J=unlist(lapply(X,nrow))
    maxJ=max(J)
    #Get state space extent such that buffer is at least buff on all grids
    #minmax=rbind(apply(X[[1]],2,min),apply(X[[1]],2,max))
    useverts=!is.na(vertices)[1]
    usedSS=!is.na(dSS)[1]
    if(useverts&usedSS){
      stop("Cannot input both vertices and dSS")
    }
    if(useverts){
      if(!is.list(vertices)){
        stop("vertices must be a list")
      }
      if(any(!unlist(lapply(vertices,is.matrix)))){
        stop("not all vertices list elements are matrices")
      }
    }
    if(!useverts){
      minmax=array(NA,dim=c(sum(primary==1),2,2))
      idx=1
      for(i in 1:length(X)){
        if(primary[i]==1){
          minmax[idx,,]=rbind(apply(X[[i]],2,min),apply(X[[i]],2,max))
          idx=idx+1
        }
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
    if(length(lam0)==1){
      lam0=rep(lam0,t)
    }
    if(length(sigma)==1){
      sigma=rep(sigma,t)
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
      }else if(usedSS){
        NdSS=nrow(dSS)
        mu=dSS[sample(1:NdSS,M,replace=TRUE),1:2]
      }else{
        mu<- cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff))
      }
      for(i in 1:t){
        s[,i,]=mu
        if(primary[i]){
          D[,1:nrow(X[[i]]),i]=e2dist(s[,i,],X[[i]])
          lamd[,,i]=lam0[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
        }
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
      }else if(usedSS){
        stop("Ben hasn't coded this, yet.  Try using vertices in continuous space")
      }else{
        mu<- cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff))
      }
      for(i in 1:t){#meta mu movement
        countout=0
        for(j in 1:M){
          out=1
          while(out==1){
            s[j,i,1]=rnorm(1,mu[j,1],sigma_t)
            s[j,i,2]=rnorm(1,mu[j,2],sigma_t)
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
        if(primary[i]){
          D[,1:nrow(X[[i]]),i]=e2dist(s[,i,],X[[i]])
          lamd[,,i]=lam0[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
        }
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
      }else if(usedSS){
        stop("Ben hasn't coded this, yet.  Try using vertices in continuous space")
      }else{
        mu<- cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff))
      }
      for(l in 1:t){#meta mu movement, only meta mus stay in SS
        for(j in 1:M){
          s[j,l,1]=rnorm(1,mu[j,1],sigma_t)
          s[j,l,2]=rnorm(1,mu[j,2],sigma_t)
        }
        if(primary[l]){
          D[,1:nrow(X[[l]]),l]=e2dist(s[,l,],X[[l]])
          lamd[,,l]=lam0[l]*exp(-D[,,l]^2/(2*sigma[l]*sigma[l]))
        }
      }

    }else if(ACtype=="markov"){
      if(useverts){
        s[,1,]=cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff)) #initial locations
        for(i in 1:M){
          inside=any(unlist(lapply(vertices,function(x){inout(s[i,1,],x)})))
          while(inside==FALSE){
            s[i,1,]=c(runif(1, xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside=any(unlist(lapply(vertices,function(x){inout(s[i,1,],x)})))
          }
        }
      }else if(usedSS){
        stop("Ben hasn't coded this, yet.  Try using markov2 with dSS or vertices in continuous space with markov")
      }else{
        s[,1,]=cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff)) #initial locations
      }
      D[,1:nrow(X[[1]]),1]=e2dist(s[,1,],X[[1]])
      lamd[,,1]=lam0[1]*exp(-D[,,1]^2/(2*sigma[1]*sigma[1]))
      for(i in 2:t){
        countout=0
        for(j in 1:M){
          out=1
          while(out==1){
            s[j,i,1]=rnorm(1,s[j,i-1,1],sigma_t)
            s[j,i,2]=rnorm(1,s[j,i-1,2],sigma_t)
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
        if(primary[i]){
          D[,1:nrow(X[[i]]),i]=e2dist(s[,i,],X[[i]])
          lamd[,,i]=lam0[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
        }
      }
    }else if(ACtype=="markov2"){
      #initial locs
      if(usedSS==FALSE){
        stop("must enter dSS with ACtype=markov2")
      }
      NdSS=nrow(dSS)
      s[,1,]=dSS[sample(1:NdSS,M,replace=TRUE),1:2] #initial locations
      for(l in 2:t){
        dists=e2dist(s[,l-1,],dSS[,1:2])
        for(i in 1:M){
          probs=dexp(dists[i,],1/sigma_t)
          probs=probs/sum(probs)
          s[i,l,]=dSS[sample(1:NdSS,1,prob=probs),1:2]
        }
      }
      for(l in 1:t){
        if(primary[l]){
          D[,1:nrow(X[[l]]),l]=e2dist(s[,l,],X[[l]])
          lamd[,,l]=lam0[l]*exp(-D[,,l]^2/(2*sigma[l]*sigma[l]))
        }
      }
    }else if(ACtype=="independent"){
      for(i in 1:t){
        if(useverts){
          s[,i,]=cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff))
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
        }else if(usedSS){
          NdSS=nrow(dSS)
          for(l in 1:t){
            s[,l,]=dSS[sample(1:NdSS,M,replace=TRUE),1:2]
          }
        }else{
          s[,i,]=cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff))
        }
        for(i in 1:t){
          if(primary[i]){
            D[,1:nrow(X[[i]]),i]=e2dist(s[,i,],X[[i]])
            lamd[,,i]=lam0[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
          }
        }
      }
    }else{
      stop("ACtype not recognized")
    }
    psi=N[1]/M
    #Calculate (per capita) gamma or N
    if(length(phi)==1){
      phi=rep(phi,t-1)
    }
    if(is.null(gamma)){
      gamma=rep(NA,t-1)
      for(l in 2:t){
        gamma[l-1]=(N[l]-phi[t-1]*N[l-1])/(N[l-1])
      }
      storeparms$gamma=gamma
    }else{
      if(length(gamma)==1){
        gamma=rep(gamma,t-1)
      }
      N=c(N,rep(NA,t-1))
      for(l in 2:t){
        N[l]=round(N[l-1]*phi[t-1]+N[l-1]*gamma[l-1])
      }
    }
    #####Population Dynamics############# Stolen from Richard Chandler
    z=a=matrix(NA,M,t)
    # z[,1] <- rbinom(M, 1, EN1/M)
    z[1:N[1],1]=1
    z[(N[1]+1):M]=0
    a[,1]= 1-z[,1] # Available to be recruited?
    gamma.prime=rep(NA,t-1)
    for(i in 2:t) {
      ER <- sum(z[,i-1])*gamma[i-1] # Expected number of recruits
      A <- sum(a[,i-1]) # nAvailable to be recruited
      if(ER>A){
        stop("M is too low. There aren't any individuals left to be recruited")
      }
      gamma.prime[i-1] <- ER/A # individual-level recruitment *probability*
      Ez <- z[,i-1]*phi[i-1] + a[,i-1]*gamma.prime[i-1]
      z[,i] <- rbinom(M, 1, Ez)
      a[,i] <- apply(z[,1:i]<1, 1, all)
    }

    #######Capture process######################
    # Simulate encounter history
    y <-array(0,dim=c(M,maxJ,t))
    if(obstype=="bernoulli"){
      pd=1-exp(-lamd)
      for(l in 1:t){
        if(primary[l]){
          for(i in 1:M){
            for(j in 1:J[l]){
              y[i,j,l]=rbinom(1,K[l],pd[i,j,l]*z[i,l])
            }
          }
        }
      }
    }else if(obstype=="poisson"){
      for(l in 1:t){
        if(primary[l]){
          for(i in 1:M){
            for(j in 1:J[l]){
              y[i,j,l]=rpois(1,K[l]*lamd[i,j,l]*z[i,l])
            }
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
    keep=which(rowSums(y)>0)
    y=y[keep,,]
    s=s[keep,,]
    if(ACtype%in%c("metamu","metamu2")){
      mu=mu[keep,]
    }
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
        out<-list(y=y,s=s,yfull=yfull,sfull=sfull,X=X,K=K,n=n,n2d=n2d,buff=buff,J=J,EN=N,N=colSums(z),
                  z=z,gamma=gamma,phi=storeparms$phi,obstype=obstype,ACtype=ACtype,vertices=vertices,
                  primary=primary)
      }else{
        out<-list(y=y,s=s,yfull=yfull,sfull=sfull,X=X,K=K,n=n,n2d=n2d,buff=buff,J=J,EN=N,N=colSums(z),
                  z=z,gamma=gamma,phi=storeparms$phi,obstype=obstype,ACtype=ACtype,primary=primary)
      }
    }else{
      if(!missing(vertices)){
        out<-list(y=y,mu=mu,s=s,yfull=yfull,sfull=sfull,mufull=mufull,X=X,K=K,n=n,n2d=n2d,buff=buff,
                  J=J,EN=N,N=colSums(z),z=z,gamma=gamma,phi=storeparms$phi,obstype=obstype,
                  ACtype=ACtype,vertices=vertices,primary=primary)
      }else{
        out<-list(y=y,mu=mu,s=s,yfull=yfull,sfull=sfull,mufull=mufull,X=X,K=K,n=n,n2d=n2d,buff=buff,
                  J=J,EN=N,N=colSums(z),z=z,gamma=gamma,phi=storeparms$phi,obstype=obstype,
                  ACtype=ACtype,primary=primary)
      }
    }
    return(out)
  }
