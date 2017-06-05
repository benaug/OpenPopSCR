SCRmcmcOpen <-
  function(data,niter=2400,nburn=1200, nthin=5,M = 200, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),
           jointZ=TRUE,keepACs=TRUE,ACtype="fixed",obstype="bernoulli",dSS=NA){
    library(abind)
    t=dim(data$y)[3]
    y<-data$y
    X<-data$X
    #make sure X list elements are matrices
    for(i in 1:length(X)){
      X[[i]]=as.matrix(X[[i]])
    }
    if(!is.na(dSS[1])&"vertices"%in%names(data)){
      rem=which(names(data)=="vertices")
      data[[rem]]=NULL
      warning("Discarding vertices since dSS supplied")
    }
    if(length(X)!=t){
      stop("must input traps for each year")
    }

    J<-data$J
    maxJ=max(J)
    K<-data$K
    n=dim(data$y)[1]
    n2d<-colSums(apply(data$y,c(3),rowSums)>0)
    ####Error checks
    if(length(K)!=t){
      stop("Must supply a K for each year")
    }
    #If using polygon state space
    if("vertices"%in%names(data)){
      vertices=data$vertices
      useverts=TRUE
      if(!is.list(vertices)){
        stop("vertices must be a list")
      }
      xlim=c(min(unlist(lapply(vertices,function(x){min(x[,1])}))),max(unlist(lapply(vertices,function(x){max(x[,1])}))))
      ylim=c(min(unlist(lapply(vertices,function(x){min(x[,2])}))),max(unlist(lapply(vertices,function(x){max(x[,2])}))))
    }else if("buff"%in%names(data)){
      buff<- data$buff
      xlim<- c(min(unlist(lapply(X,function(x){min(x[,1])}))),min(unlist(lapply(X,function(x){max(x[,1])}))))+c(-buff, buff)
      ylim<- c(min(unlist(lapply(X,function(x){min(x[,2])}))),min(unlist(lapply(X,function(x){max(x[,2])}))))+c(-buff, buff)
      vertices=list(rbind(c(xlim[1],ylim[1]),c(xlim[1],ylim[2]),c(xlim[2],ylim[2]),c(xlim[2],ylim[1])))
      useverts=FALSE
    }else{
      stop("user must supply either 'buff' or 'vertices' in data object")
    }
    if("tf"%in%names(data)){
      tf=data$tf
      if(length(tf)!=length(X)){
        stop("If using a trap operation file, must input one for each year")
      }
      for(i in 1:length(X)){
        if(is.matrix(tf[[i]])){
          stop("If using trap operation file, must enter vector of number of occasions operational, not matrix of trap by occasion operation")
        }
        if(nrow(X[[i]])!=length(tf[[i]])){
          stop("If using trap operation file, must enter operation for every trap")
        }
      }
    }else{
      tf=vector("list",t)
      for(l in 1:t){
        tf[[l]]=rep(K[l],nrow(X[[i]]))
      }
    }
    #make tf into matrix
    for(l in 1:t){
      tf[[l]]=matrix(rep(tf[[l]],M),ncol=J[l],nrow=M,byrow=TRUE)
    }


    ##pull out initial values
    lam0<- inits$lam0
    sigma<- inits$sigma
    sigma_t=inits$sigma_t
    gamma=inits$gamma
    phi=inits$phi
    psi=inits$psi
    if(!length(lam0)%in%c(1,t)){
      stop("Input either 1 or t initial values for lam0")
    }
    if(!length(sigma)%in%c(1,t)){
      stop("Input either 1 or t initial values for sigma")
    }
    if(!length(gamma)%in%c(1,t-1)){
      stop("Input either 1 or t-1 initial values for gamma (psi and fixed gamma)")
    }
    if(!length(phi)%in%c(1,t-1)){
      stop("Input either 1 or t-1 initial values for phi")
    }
    #Check proppars
    #Check proppars
    if(jointZ==FALSE&(length(proppars$propz)!=(t-1))){
      stop("must supply t-1 proppars for propz when using sequential sampler")
    }
    if(jointZ==TRUE&("propz"%in%names(proppars))){
      warning("ignoring propz because z joint z sampler specified")
    }
    if(length(lam0)!=length(proppars$lam0)){
      stop("Must supply a tuning parameter for each lam0")
    }
    if(length(sigma)!=length(proppars$sigma)){
      stop("Must supply a tuning parameter for each sigma")
    }
    if(length(gamma)!=length(proppars$gamma)){
      stop("Must supply a tuning parameter for each gamma")
    }
    if(!ACtype%in%c("fixed","independent","metamu","metamu2","markov")){
      stop("ACtype must be 'fixed','independent','metamu', 'metamu2', or 'markov'")
    }
    if(ACtype%in%c("metamu","markov")){
      if(!"sigma_t"%in%names(proppars)){
        stop("must supply proppars$sigma_t if ACtype is metamu or markov")
      }
      if(!"s1x"%in%names(proppars)|!"s1y"%in%names(proppars)){
        stop("must supply proppars$s1x and proppars$s1y if ACtype is metamu or markov")
      }
      if(is.null(sigma_t)){
        stop("must supply inits$sigma_t if ACtype is metamu or markov")
      }
    }
    #augment data
    y<- abind(y,array(0, dim=c( M-dim(y)[1],maxJ, t)), along=1)
    known.vector=c(rep(1,data$n),rep(0,M-data$n))

    #Initialize z, r, and a consistent with y
    known.matrix=1*(apply(y,c(1,3),sum)>0)#make z consistent with y
    if(t>2){
      for(l in 2:(t-1)){#Turn on zeros with 1's on either side
        known.matrix[known.matrix[,l]==0&known.matrix[,l-1]==1&rowSums(matrix(known.matrix[,(l+1):t],nrow=M))>0,l]=1
      }
    }
    z=known.matrix
    # r=array(0,dim=dim(z))
    #turn on z's with caps on either side so we know they were in pop
    for(i in 1:n){
      idx=which(z[i,]==1)
      z[i,min(idx):max(idx)]=1
    }
    #turn on some augmented guys for z1
    # z[(n+1):M,1]=rbinom(M-n,1,psi)#add augmented guys to t=1 with psi.. not enough
    z1deal=psi*M-sum(z[,1])
    if(z1deal<0){
      stop("initial psi is too small given M to turn on any uncaptured z[,1] guys ")
    }
    z[sample((n+1):M,z1deal),1]=1
    a=matrix(1,nrow=M,ncol=t) #a is available to be recruited
    a[which(z[,1]==1),]=0#turn off guys caught year 1
    if(length(gamma)==(t-1)){
      gammasim=gamma
    }else{
      gammasim=rep(gamma,t-1)
    }
    if(length(phi)==(t-1)){
      phisim=phi
    }else{
      phisim=rep(phi,t-1)
    }

    for(i in 2:t){
      #Recruitment
      nrecruit=round(sum(z[,i-1])*gammasim[i-1])#how many should we recruit?
      recruits=which(z[1:n,i-1]==0&z[1:n,i]==1)#who is already recruited based on y constraints?
      a[which(z[,i]==1),i:t]=0 #anyone alive in time t can't be recruited
      nleft=nrecruit-length(recruits)
      if(nleft>0){
        cands=setdiff(which(a[,i-1]==1),recruits)#Who is available to recruit?
        cands=cands[!cands%in%1:n]#Don't mess with known guys
        if(nleft>length(cands)){
          nleft=length(cands)
        }
        if(length(cands)>1){
          pick=sample(cands,nleft)
        }else if(length(cands)==1&nleft==1){
          pick=cands
        }else{
          next
        }
        z[pick,i]=1
        a[pick,i:t]=0 #no longer available for recruit on any occasion
      }else{
        warning("Can't initialize all E[recruits] given initial value for gamma. Should probably raise M?")
      }
      #Survival
      nlive=rbinom(1,sum(z[,i-1]),phisim[i-1])#How many to live
      livers=which(z[1:n,i-1]==1&z[1:n,i]==1)
      nleft=nlive-length(livers)
      a[livers,i]=0
      if(nleft>0){
        cands=which(apply(matrix(z[,i:t],nrow=M),1,sum)==0&z[,i-1]==1)
        cands=cands[!cands%in%1:n]
        if(nleft>length(cands)){
          nleft=length(cands)
        }
        if(length(cands)>1){
          pick=sample(cands,nleft)
        }else if(length(cands)==1&nleft==1){
          pick=cands
        }else{
          next
        }
        pick=sample(cands,nleft)
        z[pick,i]=1
        a[pick,i:t]=0
      }
    }
    N=colSums(z)

    ll.z=matrix(0, M, t)
    Ez=Ez.cand=matrix(NA, M, t-1)
    ll.z[,1]=dbinom(z[,1], 1, psi, log=TRUE)
    gamma.prime <- numeric(t-1)
    if(length(gamma)==1){
      gammause=rep(gamma,t-1)
    }else{
      gammause=gamma
    }
    if(length(phi)==1){
      phiuse=rep(phi,t-1)
    }else{
      phiuse=phi
    }
    for(l in 2:t) {
      gamma.prime[l-1]=N[l-1]*gammause[l-1] / sum(a[,l-1])
      Ez[,l-1]=z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime[l-1]
      ll.z[,l]=dbinom(z[,l], 1, Ez[,l-1], log=TRUE)
    }
    ll.z.cand=ll.z
    gamma.prime.cand=gamma.prime
    #Optimize starting locations given where they are trapped. Initalizing s1 and s2 at same locs

    s1<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    idx=which(known.vector==1) #switch for those actually caught
    for(i in idx){
      trps=matrix(0,nrow=0,ncol=2)
      for(j in 1:t){ #loop over t to get all cap locs
        trps<- rbind(trps,X[[j]][which(y[i,,j]>0),1:2])
      }
      s1[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }
    if(length(dSS)>1){
      usedSS=TRUE
    }else{
      usedSS=FALSE
    }
    if(!usedSS){
      #check to make sure everyone is in polygon
      if(useverts==TRUE){
        inside=rep(NA,nrow(s1))
        for(i in 1:nrow(s1)){
          # inside[i]=inout(s1[i,],vertices)
          inside[i]=any(unlist(lapply(vertices,function(x){inout(s1[i,],x)})))
        }
        idx2=which(inside==FALSE)
        if(any(idx2%in%idx)){
          warning("Vertices too complicated for this AC initialization algorithm to provide good starting values.
                  Either hassle Ben to fix it or use a discrete state space")
        }
        if(length(idx2)>0){
          for(i in 1:length(idx)){
            while(inside[idx2[i]]==FALSE){
              s1[idx2[i],]=c(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
              # inside[idx[i]]=inout(s1[idx[i],],vertices)
              inside[idx2[i]]=any(unlist(lapply(vertices,function(x){inout(s1[idx2[i],],x)})))
            }
          }
        }
      }
    }else{#discrete SS
      dSS=as.matrix(dSS)
      #snap everyone back to closest place in state space
      for(i in idx){
        dists=sqrt((s1[i,1]-dSS[,1])^2+(s1[i,2]-dSS[,2])^2)
        s1[i,]=dSS[which(dists==min(dists))[1],]
      }
      #randomly assign uncaptured guys to cells
      idx2=setdiff(1:M,idx)
      s1[idx2,]=dSS[sample(1:nrow(dSS),length(idx2)),]
    }
    #Initialize s2
    #Assign s2 to be s1 for all occasions
    s2=array(NA,dim=c(M,t,2))
    for(l in 1:t){
      s2[,l,]=s1
    }
    if(ACtype%in%c("metamu","metamu2","markov")){
      #update s2s for guys captured each year and add noise for uncaptured guys. More consistent with sigma_t>0
      #should be OK for markov and independent
      for(l in 1:t){
        idx=which(rowSums(y[,,l])>0) #switch for those actually caught
        for(i in 1:M){
          if(i%in%idx){
            trps<- X[[l]][which(y[i,,l]>0),1:2]
            if(is.matrix(trps)){
              s2[i,l,]<- c(mean(trps[,1]),mean(trps[,2]))
            }else{
              s2[i,l,]=trps
            }
          }else{
            inside=FALSE
            while(inside==FALSE){
              s2[i,l,]=c(rnorm(1,s1[i,1],sigma_t),rnorm(1,s1[i,2],sigma_t))
              # inside=inout(s2[i,l,],vertices)
              inside=any(unlist(lapply(vertices,function(x){inout(s2[i,l,],x)})))
            }
          }
        }
      }
      if(ACtype%in%c("metamu","metamu2")){
        ll.s2=(dnorm(s2[,,1],s1[,1],sigma_t,log=TRUE)+dnorm(s2[,,2],s1[,2],sigma_t,log=TRUE))
        ll.s2.cand=ll.s2
      }else if(ACtype=="markov"){
        ll.s2=matrix(NA,nrow=M,ncol=t-1)
        for(l in 2:t){
          ll.s2[,l-1]=(dnorm(s2[,l,1],s2[,l-1,1],sigma_t,log=TRUE)+dnorm(s2[,l,2],s2[,l-1,2],sigma_t,log=TRUE))
        }
        ll.s2.cand=ll.s2
      }
    }
    if(ACtype=="independent"){
      #update s2s for guys captured each year and put uncaptured guys OFF the GRID so they can get turned on.
      #general not as good alternative, make sure they're not within buff/2 of a trap
      for(l in 1:t){
        idx=which(rowSums(y[,,l])>0) #switch for those actually caught
        for(i in 1:M){
          if(i%in%idx){
            trps<- X[[l]][which(y[i,,l]>0),1:2]
            if(is.matrix(trps)){
              s2[i,l,]<- c(mean(trps[,1]),mean(trps[,2]))
            }else{
              s2[i,l,]=trps
            }
          }else{
            inside=FALSE
            while(inside==FALSE){
              s2[i,l,]=cbind(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
              # dists=sqrt((s2[i,l,1]-X[[l]][,1])^2+(s2[i,l,2]-X[[l]][,2])^2)
              # inside=inout(s2[i,l,],vertices)#&(!any(dists<buff/2))
              inside=any(unlist(lapply(vertices,function(x){inout(s2[i,l,],x)})))

            }
          }
        }
      }
    }
    if(usedSS){
      #snap everyone back to closest place in state space
      for(l in 1:t){
        for(i in 1:M){
          dists=sqrt((s2[i,l,1]-dSS[,1])^2+(s2[i,l,2]-dSS[,2])^2)
          s2[i,l,]=dSS[which(dists==min(dists))[1],]
        }
      }
    }
    # some objects to hold the MCMC simulation output
    if(niter<(nburn)){
      stop("niter is smaller than nburn")
    }
    nstore=(niter-nburn)/nthin
    if((nburn)%%nthin!=0){
      nstore=nstore+1
    }
    if(length(lam0)==t){
      lam0names=paste("lam0",1:t,sep="")
    }else{
      lam0names="lam0"
    }
    if(length(sigma)==t){
      sigmanames=paste("sigma",1:t,sep="")
    }else{
      sigmanames="sigma"
    }
    if(length(gamma)==(t-1)){
      gammanames=paste("gamma",1:(t-1),sep="")
    }else{
      gammanames="gamma"
    }
    if(length(phi)==(t-1)){
      phinames=paste("phi",1:(t-1),sep="")
    }else{
      phinames="phi"
    }
    Nnames=paste("N",1:t,sep="")
    if(ACtype%in%c("metamu","metamu2","markov")){
      out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t+1)
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,"sigma_t")
      s1xout<- s1yout<- matrix(NA,nrow=nstore,ncol=M)
      zout<-array(NA,dim=c(nstore,M,t))
      s2xout<- s2yout<-array(NA,dim=c(nstore,M,t))
    }else if(ACtype=="independent"){
      out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t)
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames)
      s1xout<- s1yout<- matrix(NA,nrow=nstore,ncol=M)
      zout<-array(NA,dim=c(nstore,M,t))
      s2xout<- s2yout<-array(NA,dim=c(nstore,M,t))
    }else{
      out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t)
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames)
      s1xout<- s1yout<- matrix(NA,nrow=nstore,ncol=M)
      zout<-array(NA,dim=c(nstore,M,t))
    }
    idx=1 #for storing output not recorded every iteration

    D=lamd=ll.y=ll.y.cand=array(NA,dim=c(M,maxJ,t))
    D[is.na(D)]=Inf  #hack to allow years with different J and K to fit in one array
    for(l in 1:t){
      D[,1:J[l],l]=e2dist(s2[,l,],X[[l]])
      if(length(lam0)==t&length(sigma)==t){
        lamd[,1:J[l],l]=lam0[l]*exp(-D[,1:J[l],l]^2/(2*sigma[l]*sigma[l]))
      }else if(length(lam0)==1&length(sigma)==t){
        lamd[,1:J[l],l]=lam0*exp(-D[,1:J[l],l]^2/(2*sigma[l]*sigma[l]))
      }else if(length(lam0)==t&length(sigma)==1){
        lamd[,1:J[l],l]=lam0[l]*exp(-D[,1:J[l],l]^2/(2*sigma*sigma))
      }else{
        lamd[,1:J[l],l]=lam0*exp(-D[,1:J[l],l]^2/(2*sigma*sigma))
      }
    }
    #Calculate ll for observation model
    if(obstype=="bernoulli"){
      pd=pd.cand=1-exp(-lamd)
      for(l in 1:t){
        ll.y[,1:J[l],l]= dbinom(y[,1:J[l],l],tf[[l]],pd[,1:J[l],l]*z[,l],log=TRUE)
      }
    }else if(obstype=="poisson"){
      for(l in 1:t){
        ll.y[,1:J[l],l]= dpois(y[,1:J[l],l],tf[[l]]*lamd[,1:J[l],l]*z[,l],log=TRUE)
      }
    }else{
      stop("obstype must be 'bernoulli' or 'poisson'")
    }
    ll.y.cand=ll.y
    ll.y.t.sum=ll.y.cand.t.sum=apply(ll.y,3,sum) #ll summed for each year
    ll.y.sum=sum(ll.y.t.sum) #full ll sum
    lamd.cand=lamd
    #Check ll.y.sum for inf
    if(!is.finite(ll.y.sum)){
      stop("Detection function starting values produce -Inf log likelihood values. Try increasing sigma and/or lam0")
    }
    if(jointZ==TRUE){
      #Figure out all possible z histories
      zpossible=cbind(c(1,1,0),c(1,0,1))
      if (t > 2) {#4 lines from From Rcapture histpost()
        for (i in (3:t)) {
          zpossible=cbind(c(rep(1,2^(i-1)),rep(0,((2^(i-1))-1))), rbind(zpossible, rep(0, (i-1)), zpossible))
        }
        #remove zombie histories
        illegal=rep(FALSE,nrow(zpossible))
        for(l in 1:(t-2)){
          latecaps=rep(0,nrow(zpossible))
          for(l2 in (l+2):(t)){
            latecaps=latecaps+zpossible[,l2]
          }
          illegal=illegal|(zpossible[,l]==1&zpossible[,l+1]==0&latecaps>0)
        }
        zpossible=zpossible[which(illegal==FALSE),]
      }

      zpossible=rbind(zpossible,rep(0,t))#add on all zero history
      nzpossible=nrow(zpossible)
      apossible=matrix(1,nrow=nzpossible,ncol=t)
      apossible[zpossible[,1]==1,1]=0
      for(l in 2:t){
        apossible[apossible[,l-1]==1&zpossible[,l]==1,l]=0
        apossible[apossible[,l-1]==0,l]=0
      }
      Ezpossible=matrix(NA,nrow=nzpossible,ncol=t-1)
      ll_z_possible=matrix(NA,nrow=nzpossible,ncol=t)
      #Can always update anyone who wasn't captured on every occasion
      upz3=which(rowSums(known.matrix)!=t)
      #Zero out known matrix years for both swapped
      cancel=matrix(1,nrow=M,ncol=nzpossible)
      for(i in 1:M){
        fixed=which(known.matrix[i,]==1)
        for(i2 in 1:nzpossible){
          if(!all(fixed%in%which(zpossible[i2,]==1))){
            cancel[i,i2]=0
          }
        }
      }
    }

    for(iter in 1:niter){
      ll.y.t.sum=apply(ll.y,3,sum) #only needed for detection parameters, changes in z and AC updates
      ll.y.sum=sum(ll.y.t.sum)
      # Update lam0
      if(length(lam0)==t){ #if lam0 is year-specific
        for(l in 1:t){
          lam0.cand<- rnorm(1,lam0[l],proppars$lam0[l])
          if(lam0.cand > 0){
            if(length(sigma)==t){#if sigma is year specific
              lamd.cand[,1:J[l],l]<- lam0.cand*exp(-D[,1:J[l],l]^2/(2*sigma[l]*sigma[l]))
            }else{#fixed sigma
              lamd.cand[,1:J[l],l]<- lam0.cand*exp(-D[,1:J[l],l]^2/(2*sigma*sigma))
            }
            if(obstype=="bernoulli"){
              pd.cand[,1:J[l],l]=1-exp(-lamd.cand[,1:J[l],l])
              ll.y.cand[,1:J[l],l]= dbinom(y[,1:J[l],l],tf[[l]],pd.cand[,1:J[l],l]*z[,l],log=TRUE) #only need to update this year
            }else{
              ll.y.cand[,1:J[l],l]= dpois(y[,1:J[l],l],tf[[l]]*lamd.cand[,1:J[l],l]*z[,l],log=TRUE)
            }
            ll.y.cand.t.sum[l]=sum(ll.y.cand[,1:J[l],l])#just 1 year
            if(runif(1) < exp(ll.y.cand.t.sum[l] -ll.y.t.sum[l])){
              lam0[l]<- lam0.cand
              lamd[,1:J[l],l]=lamd.cand[,1:J[l],l]
              if(obstype=="bernoulli"){
                pd[,1:J[l],l]=pd.cand[,1:J[l],l]
              }
              ll.y[,1:J[l],l]=ll.y.cand[,1:J[l],l]
              ll.y.t.sum[l]=ll.y.cand.t.sum[l]
            }
          }
        }
        ll.y.sum=sum(ll.y.t.sum)
      }else{#fixed lam0
        lam0.cand<- rnorm(1,lam0,proppars$lam0)
        if(lam0.cand > 0){
          if(length(sigma)==t){#if sigma is year specific
            for(l in 1:t){
              lamd.cand[,1:J[l],l]<- lam0.cand*exp(-D[,1:J[l],l]^2/(2*sigma[l]*sigma[l]))
            }
          }else{#fixed sigma
            for(l in 1:t){
              lamd.cand[,1:J[l],l]<- lam0.cand*exp(-D[,1:J[l],l]^2/(2*sigma*sigma))
            }
          }
          if(obstype=="bernoulli"){
            pd.cand=1-exp(-lamd.cand)
            for(l in 1:t){
              ll.y.cand[,1:J[l],l]= dbinom(y[,1:J[l],l],tf[[l]],pd.cand[,1:J[l],l]*z[,l],log=TRUE)
            }
          }else{
            for(l in 1:t){
              ll.y.cand[,1:J[l],l]= dpois(y[,1:J[l],l],tf[[l]]*lamd.cand[,1:J[l],l]*z[,l],log=TRUE)
            }
          }
          ll.y.cand.t.sum=apply(ll.y.cand,3,sum)
          ll.y.cand.sum=sum(ll.y.cand.t.sum)
          if(runif(1) < exp(ll.y.cand.sum - ll.y.sum)){
            lam0= lam0.cand
            lamd=lamd.cand
            if(obstype=="bernoulli"){
              pd=pd.cand
            }
            ll.y=ll.y.cand
            ll.y.sum=ll.y.cand.sum
            ll.y.t.sum=ll.y.cand.t.sum
          }
        }
      }
      #Update sigma
      if(length(sigma)==t){ #if sigma is year-specific
        for(l in 1:t){
          sigma.cand<- rnorm(1,sigma[l],proppars$sigma[l])
          if(sigma.cand > 0){
            if(length(lam0)==t){#if lam0 is year specific
              lamd.cand[,1:J[l],l]<- lam0[l]*exp(-D[,1:J[l],l]^2/(2*sigma.cand*sigma.cand))
            }else{#fixed lam0
              lamd.cand[,1:J[l],l]<- lam0*exp(-D[,1:J[l],l]^2/(2*sigma.cand*sigma.cand))
            }
            if(obstype=="bernoulli"){
              pd.cand[,1:J[l],l]=1-exp(-lamd.cand[,1:J[l],l])
              ll.y.cand[,1:J[l],l]= dbinom(y[,1:J[l],l],tf[[l]],pd.cand[,1:J[l],l]*z[,l],log=TRUE) #only need to update this year
            }else{
              ll.y.cand[,1:J[l],l]= dpois(y[,1:J[l],l],tf[[l]]*lamd.cand[,1:J[l],l]*z[,l],log=TRUE) #only need to update this year
            }
            ll.y.cand.t.sum[l]=sum(ll.y.cand[,1:J[l],l])#just 1 year
            if(runif(1) < exp(ll.y.cand.t.sum[l] -ll.y.t.sum[l])){
              sigma[l]<- sigma.cand
              lamd[,1:J[l],l]=lamd.cand[,1:J[l],l]
              if(obstype=="bernoulli"){
                pd[,1:J[l],l]=pd.cand[,1:J[l],l]
              }
              ll.y[,1:J[l],l]=ll.y.cand[,1:J[l],l]
              ll.y.t.sum[l]=ll.y.cand.t.sum[l]
            }
          }
        }
        ll.y.sum=sum(ll.y.t.sum)
      }else{#fixed sigma
        sigma.cand<- rnorm(1,sigma,proppars$sigma)
        if(sigma.cand > 0){
          if(length(lam0)==t){#if lam0 is year specific
            for(l in 1:t){
              lamd.cand[,1:J[l],l]<- lam0[l]*exp(-D[,1:J[l],l]^2/(2*sigma.cand*sigma.cand))
            }
          }else{#fixed lam0
            for(l in 1:t){
              lamd.cand[,1:J[l],l]<- lam0*exp(-D[,1:J[l],l]^2/(2*sigma.cand*sigma.cand))
            }          }
          if(obstype=="bernoulli"){
            pd.cand=1-exp(-lamd.cand)
            for(l in 1:t){
              ll.y.cand[,1:J[l],l]= dbinom(y[,1:J[l],l],tf[[l]],pd.cand[,1:J[l],l]*z[,l],log=TRUE)
            }
          }else{
            for(l in 1:t){
              ll.y.cand[,1:J[l],l]= dpois(y[,1:J[l],l],tf[[l]]*lamd.cand[,1:J[l],l]*z[,l],log=TRUE)
            }
          }
          ll.y.cand.t.sum=apply(ll.y.cand,3,sum)
          ll.y.cand.sum=sum(ll.y.cand.t.sum)
          if(runif(1) < exp(ll.y.cand.sum - ll.y.sum)){
            sigma<- sigma.cand
            lamd=lamd.cand
            if(obstype=="bernoulli"){
              pd=pd.cand
            }
            ll.y=ll.y.cand
            ll.y.sum=ll.y.cand.sum#dont really need these anymore
            ll.y.t.sum=ll.y.cand.t.sum
          }
        }
      }
      # Update z[,1]
      if(jointZ==FALSE){
        if(t>3){
          upz=which(!(z[,1]==0&z[,2]==0&rowSums(z[,3:t]>0))&known.matrix[,1]==0)
        }else if(t==3){
          upz=which(!(z[,1]==0&z[,2]==0&z[,3]>0)&known.matrix[,1]==0)
        }else{
          upz=which(known.matrix[,1]==0)
        }
        for(i in upz){
          gamma.prime.cand <- gamma.prime
          #Do we need to modify a in more than one year.
          z.cand <- z #use full z to calculate correct proposed Ez.cand
          z.cand[i,1] <- 1-z[i,1]
          if(obstype=="bernoulli"){
            ll.y.cand[i,1:J[1],1]=dbinom(y[i,1:J[1],1],tf[[1]][i,],pd[i,1:J[1],1]*z.cand[i,1],log=TRUE)
          }else{
            ll.y.cand[i,1:J[1],1]=dpois(y[i,1:J[1],1],tf[[1]][i,]*lamd[i,1:J[1],1]*z.cand[i,1],log=TRUE)
          }
          if(((z.cand[i,1]==1&sum(z[i,])==0)|(sum(z[i,])==1&z.cand[i,1]==0&z[i,1]==1))&(t>2)){#Are we turning on a guy that was never on before? or turning off a guy that was only on on z1?
            a.cand <- a
            #only a option is all on or all off
            if(z.cand[i,1]==1&sum(a[i,])==t){
              a.cand[i,]=0 #all off
            }else{
              a.cand[i,]=1 #all on
            }
            #Calculate gamma.prime, Ez, and ll.z candidates
            ll.z.cand[i,1] <- dbinom(z.cand[i,1], 1, psi, log=TRUE)
            #Make object to put N1_prop in to but use current values of the other Ns
            Ntmp=N
            Ntmp[1]=sum(z.cand[,1])
            for(l in 2:t){
              gamma.prime.cand[l-1]=(Ntmp[l-1]*gammause[l-1]) / sum(a.cand[,l-1])
              if(gamma.prime.cand[l-1] > 1) { # E(Recruits) must be < nAvailable
                warning("Rejected z due to low M")
                next
              }
              Ez.cand[,l-1]=z.cand[,l-1]*phiuse[l-1] + a.cand[,l-1]*gamma.prime.cand[l-1]
              ll.z.cand[,l]=dbinom(z.cand[,l], 1, Ez.cand[,l-1], log=TRUE)
            }
            if(runif(1) < exp((sum(ll.y.cand[i,1:J[1],1])+ ll.z.cand[i,1]+sum(ll.z.cand[,-1]))-(sum(ll.y[i,1:J[1],1])+ll.z[i,1]+sum(ll.z[,-1])) )) {
              ll.y[i,,1] = ll.y.cand[i,,1]
              ll.z[i,1] = ll.z.cand[i,1]
              ll.z[,-1] = ll.z.cand[,-1]
              Ez = Ez.cand
              z[i,1] = z.cand[i,1]
              gamma.prime = gamma.prime.cand
              a=a.cand
            }
          }else{#Don't need to modify more than one year
            z1.cand <- z[,1]
            z1.cand[i] <- 1-z[i,1]
            a1.cand <- a[,1]
            a1.cand[i] <- 1-z1.cand[i]
            gamma.prime.cand[1] <- sum(z1.cand)*gamma[1] / sum(a1.cand)
            if(gamma.prime.cand[1] > 1) { # E(Recruits) must be < nAvailable
              warning("Rejected z due to low M")
              next
            }
            ll.z.cand[i,1] <- dbinom(z1.cand[i], 1, psi, log=TRUE)
            Ez.cand[,1]=z1.cand*phi[1] + a1.cand*gamma.prime.cand[1]
            ll.z.cand[,2] <- dbinom(z[,2], 1, Ez.cand[,1], log=TRUE)
            if(runif(1) < exp((sum(ll.y.cand[i,1:J[1],1])+ ll.z.cand[i,1]+sum(ll.z.cand[,2]))-(sum(ll.y[i,1:J[1],1])+ll.z[i,1]+sum(ll.z[,2])) )) {#z1 and z2 matter
              ll.y[i,1:J[1],1] = ll.y.cand[i,1:J[1],1]
              ll.z[i,1] = ll.z.cand[i,1]
              ll.z[,2] = ll.z.cand[,2]
              Ez[,1] = Ez.cand[,1]
              z[i,1] = z1.cand[i]
              gamma.prime[1] = gamma.prime.cand[1]
              a[i,1]=a1.cand[i]
            }
          }
        }
        N[1]=sum(z[,1])

        for(l in 2:t){
          #Determine who can be updated
          upz=which(known.matrix[,l]==0)#guys not caught on or on either side of this occ
          #Figure out illegal moves.  Depends on t.
          #Could use apply function, but need something to convert to Rcpp
          if(t>3){
            if(l==2&l==(t-2)){#only happens if t=4
              rem=which(z[,1]>0&rowSums(z[,(l+1):t])>0)#guys in pop before and after
              rem2=which(z[,l]==0&z[,l+1]==0&z[,t]>0)#guys with 0 0 and subsequent 1 can't be turned on
            }else if(l==2){
              rem=which(z[,1]>0&rowSums(z[,(l+1):t])>0)#guys in pop before and after
              rem2=which(z[,l]==0&z[,l+1]==0&rowSums(z[,(l+2):t])>0) #guys with 0 0 and subsequent 1 can't be turned on
            }else if(l==(t-2)){
              rem=which(rowSums(z[,1:(l-1)])>0&rowSums(z[,(l+1):t])>0)#guys in pop before and after
              rem2=which(z[,l]==0&z[,l+1]==0&z[,t]>0)#guys with 0 0 and subsequent 1 can't be turned on
            }else if(l==(t-1)){
              rem=which(rowSums(z[,1:(l-1)])>0&z[,t]>0)#guys in pop before and after
              rem2=integer()
            }else if (l==t){
              rem=integer()
              rem2=integer()
            }else{
              rem=which(rowSums(z[,1:(l-1)])>0&rowSums(z[,(l+1):t])>0)#guys in pop before and after
              rem2=which(z[,l]==0&z[,l+1]==0&rowSums(z[,(l+2):t])>0)#guys with 0 0 and subsequent 1 can't be turned on
            }
          }else if(t==3){
            if(l==2){
              rem=which(z[,1]>0&z[,t]>0)#guys in pop before and after
              rem2=integer()
            }else{#l==3
              rem=integer()
              rem2=integer()
            }
          }else{#t==2
            rem=integer()
            rem2=integer()
          }
          rem3=which( (z[,l-1]==0) & (a[,l-1]==0)) #dead guys.
          remall=c(rem,rem2,rem3)
          upz=upz[!upz%in%remall]#We can change these guys
          navail=length(upz)
          if(navail<1){
            next
          }
          if(navail < proppars$propz[l-1]) {
            propz=navail
            warning("M isn't big enough to propose all propz")
          }else{
            propz=proppars$propz[l-1]
          }
          swapz=upz[sample.int(navail, propz)]
          #Update swapz one at a time
          for(i in swapz){
            zt.cand=z[,l]
            ####Try proposing not based on Ez. Just swap it. Take out prop and back probs. Should help mixing
            #Normal stuff
            zt.cand[i]=1-z[i,l]
            at.cand=1*(a[,l-1]==1&zt.cand==0) #who was available on last occasion and not proposed to be captured?
            if(obstype=="bernoulli"){
              ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,],pd[i,1:J[l],l]*zt.cand[i],log=TRUE)
            }else{
              ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd[i,1:J[l],l]*zt.cand[i],log=TRUE)
            }
            # ll.z.cand[,l] <- dbinom(zt.cand, 1, Ez[,l-1], log=TRUE) ## Don't subset z
            ll.z.cand[i,l] <- dbinom(zt.cand[i], 1, Ez[i,l-1], log=TRUE) ## why not?
            # prior.z <- sum(ll.z[i,l])
            # prior.z.cand <- sum(ll.z.cand[i,l])
            prior.z=ll.z[i,l]
            prior.z.cand=ll.z.cand[i,l]
            #New stuff
            fix1=zt.cand[i]==1&sum(z[i,])==0 #guys never in pop proposed to be turned on
            fix2=sum(z[i,])==1&zt.cand[i]==0&z[i,l]==1 #guys in pop only once and proposed to be turned off
            if((fix1|fix2)&(t>3)&(l<t)){
              a.cand <- a
              a.cand[,l]=at.cand
              z.cand=z
              z.cand[,l]=zt.cand
              if(fix1){
                a.cand[i,l:t]=0 #all off
              }
              if(fix2){
                a.cand[i,l:t]=1 #all on
              }
              #Calculate gamma.prime, Ez, and ll.z for l:t
              reject=FALSE
              Ntmp=N
              Ntmp[l]=sum(zt.cand)
              for(l2 in l:(t-1)){
                gamma.prime.cand[l2]=(Ntmp[l2]*gammause[l2]) / sum(a.cand[,l2])
                if(gamma.prime.cand[l2] > 1){
                  reject=TRUE
                }
                Ez.cand[,l2]=z.cand[,l2]*phiuse[l2] + a.cand[,l2]*gamma.prime.cand[l2]
                ll.z.cand[,l2+1]=dbinom(z.cand[,l2+1], 1, Ez.cand[,l2], log=TRUE)
                prior.z <- prior.z + sum(ll.z[,l2+1])
                prior.z.cand <- prior.z.cand + sum(ll.z.cand[,l2+1])
              }
              if(reject){
                warning("Rejected z due to low M")
                next
              }
            }else{
              #Calculate gamma.prime, Ez, and ll.z for l
              if(l<t){ ## NOTE: Don't subset with swapz
                gamma.prime.cand[l] <- sum(zt.cand)*gammause[l] / sum(at.cand)
                if(gamma.prime.cand[l] > 1){
                  warning("Rejected z due to low M")
                  next
                }
                Ez.cand[,l] <- zt.cand*phiuse[l] + at.cand*gamma.prime.cand[l]
                ll.z.cand[,l+1] <- dbinom(z[,l+1], 1, Ez.cand[,l], log=TRUE)
                prior.z <- prior.z + sum(ll.z[,l+1])
                prior.z.cand <- prior.z.cand + sum(ll.z.cand[,l+1])
              }
            }
            if(runif(1) < exp((sum(ll.y.cand[i,1:J[l],l]) + prior.z.cand) - (sum(ll.y[i,1:J[l],l]) +prior.z) )) {
              ll.y[i,1:J[l],l] <- ll.y.cand[i,1:J[l],l]
              ll.z[i,l] <- ll.z.cand[i,l]
              if((fix1|fix2)&(t>3)&(l<t)){
                z=z.cand
                a=a.cand
                ll.z[,l:t]=ll.z.cand[,l:t]
                Ez[,l:(t-1)]=Ez.cand[,l:(t-1)]
                gamma.prime[l:(t-1)]= gamma.prime.cand[l:(t-1)]
              }else{
                z[,l] <- zt.cand
                a[,l] <- at.cand
                if(l < t) {
                  ll.z[,l+1] <- ll.z.cand[,l+1]
                  Ez[,l] <- Ez.cand[,l]
                  gamma.prime[l] <- gamma.prime.cand[l]
                }
              }
              N[l] <- sum(z[,l])
            }
          }
        }
      }else{
        #jointZ update
        #ll.z[,1] won't change across i
        ll_z_possible[,1]=dbinom(zpossible[,1], 1, psi,log=TRUE)
        #Get likelihood for all possible z histories. Need to update when accepted
        for(l in 2:t){
          Ezpossible[,l-1]=zpossible[,l-1]*phiuse[l-1] + apossible[,l-1]*gamma.prime[l-1]
          ll_z_possible[,l]=dbinom(zpossible[,l], 1, Ezpossible[,l-1],log=TRUE)
        }
        for(i in upz3){
          #new z stuff
          propto1=rowSums(exp(ll_z_possible))*(1*(cancel[i,]==1))
          propto=propto1/sum(propto1)
          zchoose=sample(1:nzpossible,1,prob=propto)
          zprop=zpossible[zchoose,]
          if(all(zprop==z[i,])) next #abort if proposal is current history
          prop.prob=propto[zchoose]
          #old z stuff
          currz=which(apply(zpossible,1,function(x){all(x==z[i,])}))
          back.prob=propto[currz]

          #Because a and z changes, must update gamma.prime and Ez
          #Don't need to update all years every time, but not figuring that out for now
          aprop=apossible[zchoose,]

          #Calculate gamma.prime, Ez, and ll.z candidates
          z.cand=z
          z.cand[i,]=zprop
          a.cand=a
          a.cand[i,]=aprop
          Ntmp=colSums(z.cand)
          ll.z.cand[i,1] <- dbinom(z.cand[i,1], 1, psi, log=TRUE)
          for(l in 2:t){
            gamma.prime.cand[l-1]=(Ntmp[l-1]*gammause[l-1]) / sum(a.cand[,l-1])
            if(gamma.prime.cand[l-1] > 1) { # E(Recruits) must be < nAvailable
              warning("Rejected z due to low M")
              next
            }
            Ez.cand[,l-1]=z.cand[,l-1]*phiuse[l-1] + a.cand[,l-1]*gamma.prime.cand[l-1]
            ll.z.cand[,l]=dbinom(z.cand[,l], 1, Ez.cand[,l-1], log=TRUE)
          }
          #update ll.y
          if(obstype=="bernoulli"){
            for(l in 1:t){
              ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,],pd[i,1:J[l],l]*z.cand[i,l],log=TRUE)
            }
          }else{
            for(l in 1:t){
              ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd[i,1:J[l],l]*z.cand[i,l],log=TRUE)
            }
          }

          #MH step
          if(runif(1)<exp((sum(ll.y.cand[i,1:J[l],])+ll.z.cand[i,1]+sum(ll.z.cand[,2:t]))-(sum(ll.y[i,1:J[l],])+ll.z[i,1]+sum(ll.z[,2:t])))*(back.prob/prop.prob)){
            z[i,]=zprop
            a[i,]=aprop
            gamma.prime=gamma.prime.cand
            N=Ntmp
            Ez=Ez.cand
            ll.z[i,1]=ll.z.cand[i,1]
            ll.z[,2:t]=ll.z.cand[,2:t]
            ll.y[i,1:J[l],]=ll.y.cand[i,1:J[l],]
            #Update likelihood for all possible z histories
            for(l in 2:t){
              Ezpossible[,l-1]=zpossible[,l-1]*phiuse[l-1] + apossible[,l-1]*gamma.prime[l-1]
              ll_z_possible[,l]=dbinom(zpossible[,l], 1, Ezpossible[,l-1],log=TRUE)
            }
          }
        }
      }
      # if(t>3){
      #   #a for last year isn't updated so it does not matter if it comes back on.
      #   sanity=any(rowSums(z)==0&rowSums(a[,-t])!=(t-1))
      #   for(l2 in 2:(t-2)){
      #     sanity=c(sanity,any(a[,l2]==0&a[,l2-1]==1&a[,l2+1]==1))
      #   }
      #   if(any(sanity)){stop("insanity")}
      # }
      #update psi
      psi <- rbeta(1, 1+N[1], 1+M-N[1])
      ll.z[,1] <- ll.z.cand[,1] <- dbinom(z[,1], 1, psi, log=TRUE)

      #Update phi
      if(length(phi)==(t-1)){#if time-specific survival
        for(l in 2:t){
          survive=sum(z[,l-1]==1&z[,l]==1)
          dead=sum(z[,l-1]==1&z[,l]==0)
          phi[l-1]=rbeta(1, 1+survive, 1+dead)
        }
      }else{
        survive=sum(z[,-t]==1&z[,-1]==1)
        dead=sum(z[,-t]==1&z[,-1]==0)
        phi=rbeta(1, 1+survive, 1+dead)
      }

      ## Update gamma
      ## NOTE: Must update ll.z, Ez, etc...
      if(length(phi)==1){
        phiuse=rep(phi,t-1)
      }else{
        phiuse=phi
      }
      if(length(gamma)==1){
        gamma.cand <- rnorm(1, gamma, proppars$gamma)
        gamma.cand.ok <- TRUE
        for(l in 2:t) {
          gamma.prime.cand[l-1] <- (N[l-1]*gamma.cand) / sum(a[,l-1])
          if(gamma.prime.cand[l-1] > 1){ ## Note don't break loop b/c ll.z needs updating because phi changed
            gamma.cand.ok=FALSE
          }
          Ez[,l-1] <- z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime[l-1]
          ll.z[,l] <- dbinom(z[,l], 1, Ez[,l-1], log=TRUE)
        }
        if(gamma.cand>0 & gamma.cand.ok) {
          #Only update ll.z for a=1 cases. originally. Changed from Chandler. changed back.
          for(l in 2:t) {
            Ez.cand[,l-1] <- z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime.cand[l-1]
            # ll.z.cand[a[,l-1]==1,l] <- dbinom(z[a[,l-1]==1,l], 1, Ez.cand[a[,l-1]==1,l-1], log=TRUE)
            ll.z.cand[,l] <- dbinom(z[,l], 1, Ez.cand[,l-1], log=TRUE)
          }
          if(runif(1) < exp(sum(ll.z.cand[,-1]) - sum(ll.z[,-1]))) {
            gamma <- gamma.cand
            gamma.prime <- gamma.prime.cand
            Ez <- Ez.cand
            ll.z[,-1] <- ll.z.cand[,-1]
          }
        }
      }else{
        for(l in 2:t){
          gamma.cand <- rnorm(1, gamma[l-1], proppars$gamma[l-1])
          gamma.cand.ok <- TRUE
          gamma.prime.cand[l-1] <- N[l-1]*gamma.cand / sum(a[,l-1])
          if(gamma.prime.cand[l-1] > 1){ ## Note don't break loop b/c ll.z needs updating because phi changed
            gamma.cand.ok=!gamma.cand.ok
          }
          Ez[,l-1] <- z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime[l-1]
          ll.z[,l] <- dbinom(z[,l], 1, Ez[,l-1], log=TRUE)
          if(gamma.cand>0 & gamma.cand.ok) {
            Ez.cand[,l-1] <- z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime.cand[l-1]
            # ll.z.cand[a[,l-1]==1,l] <- dbinom(z[a[,l-1]==1,l], 1, Ez.cand[a[,l-1]==1,l-1], log=TRUE)
            ll.z.cand[,l] <- dbinom(z[,l], 1, Ez.cand[,l-1], log=TRUE)
            if(runif(1) < exp(sum(ll.z.cand[,l]) - sum(ll.z[,l]))) {
              gamma[l-1] <- gamma.cand
              gamma.prime[l-1] <- gamma.prime.cand[l-1]
              Ez[,l-1] <- Ez.cand[,l-1]
              ll.z[,l] <- ll.z.cand[,l]
            }
          }
        }
      }
      #Update gamma use
      if(length(gamma)==1){
        gammause=rep(gamma,t-1)
      }else{
        gammause=gamma
      }
      ## Now we have to update the activity centers
      if(ACtype%in%c("metamu","metamu2")){
        #Update within year ACs
        for (i in 1:M){
          for(l in 1:t){
            Scand=c(rnorm(1, s2[i,l,1], proppars$s2x), rnorm(1, s2[i,l,2], proppars$s2y))
            if(ACtype=="metamu"){
              if(usedSS){
                dists=sqrt((Scand[1]-dSS[,1])^2+(Scand[2]-dSS[,2])^2)
                Scand=dSS[which(dists==min(dists)),]
              }
              if(useverts==FALSE){
                inbox=Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
              }else{
                # inbox=inout(Scand,vertices)
                inbox=any(unlist(lapply(vertices,function(x){inout(Scand,x)})))
              }
            }else{#don't force them to stay in
              inbox=TRUE
            }
            if(inbox){
              dtmp=sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
              if(length(lam0)==1&length(sigma==1)){
                lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
              }else if(length(lam0)==t&length(sigma)==1){
                lamd.cand[i,1:J[l],l]<- lam0[l]*exp(-dtmp*dtmp/(2*sigma*sigma))
              }else{
                lamd.cand[i,1:J[l],l]<- lam0[l]*exp(-dtmp*dtmp/(2*sigma[l]*sigma[l]))
              }
              if(obstype=="bernoulli"){
                pd.cand[i,1:J[l],l]=1-exp(-lamd.cand[i,1:J[l],l])
                ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,], pd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
              }else{
                ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
              }
              ll.s2.cand[i,l]<- dnorm(Scand[1],s1[i,1],sigma_t,log=TRUE)+dnorm(Scand[2],s1[i,2],sigma_t,log=TRUE)
              if(runif(1) < exp((sum(ll.y.cand[i,1:J[l],l])+ll.s2.cand[i,l]) -(sum(ll.y[i,1:J[l],l])+ll.s2[i,l]))){
                s2[i,l,] <- Scand
                D[i,1:nrow(X[[l]]),l] <- dtmp
                lamd[i,1:J[l],l] <- lamd.cand[i,1:J[l],l]
                if(obstype=="bernoulli"){
                  pd[i,1:J[l],l]=pd.cand[i,1:J[l],l]
                }
                ll.y[i,1:J[l],l]=ll.y.cand[i,1:J[l],l]
                ll.s2[i,l]=ll.s2.cand[i,l]
              }
            }
          }
        }
        #Update meta mus
        for (i in 1:M){
          Scand <- c(rnorm(1,s1[i,1],proppars$s1x), rnorm(1,s1[i,2],proppars$s1y))
          if(usedSS){
            dists=sqrt((Scand[1]-dSS[,1])^2+(Scand[2]-dSS[,2])^2)
            Scand=dSS[which(dists==min(dists)),]
          }
          if(useverts==FALSE){
            inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
          }else{
            # inbox=inout(Scand,vertices)
            inbox=any(unlist(lapply(vertices,function(x){inout(Scand,x)})))
          }
          if(inbox){
            #Count z==0 guys
            ll.s2.cand[i,]<-dnorm(s2[i,,1],Scand[1],sigma_t,log=TRUE)+dnorm(s2[i,,2],Scand[2],sigma_t,log=TRUE)
            if (runif(1) < exp(sum(ll.s2.cand[i,]) - sum(ll.s2[i,]))) {
              s1[i, ]=Scand
              ll.s2[i,]=ll.s2.cand[i,]
            }
          }
        }
        #Update sigma_t
        sigma_t.cand <- rnorm(1,sigma_t,proppars$sigma_t)
        if(sigma_t.cand > 0){
          ll.s2.cand=dnorm(s2[,,1],s1[,1],sigma_t.cand,log=TRUE)+dnorm(s2[,,2],s1[,2],sigma_t.cand,log=TRUE)
          if (runif(1) < exp(sum(ll.s2.cand) - sum(ll.s2))) {
            sigma_t=sigma_t.cand
            ll.s2=ll.s2.cand
          }
        }
      }else if(ACtype=="fixed"){#Stationary ACs
        for (i in 1:M) {
          Scand <- c(rnorm(1, s1[i, 1], proppars$s2x), rnorm(1, s1[i, 2], proppars$s2y))
          if(usedSS){
            dists=sqrt((Scand[1]-dSS[,1])^2+(Scand[2]-dSS[,2])^2)
            Scand=dSS[which(dists==min(dists)),]
          }
          if(useverts==FALSE){
            inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
          }else{
            # inbox=inout(Scand,vertices)
            inbox=any(unlist(lapply(vertices,function(x){inout(Scand,x)})))
          }
          if (inbox) {
            dtmp=matrix(Inf,maxJ,t)
            for(l in 1:t){
              dtmp[1:J[l],l] <- sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
            }
            if(length(lam0)==1&length(sigma==1)){
              for(l in 1:t){
                lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp[1:J[l],l]*dtmp[1:J[l],l]/(2*sigma*sigma))
              }
            }else if(length(lam0)==t&length(sigma)==1){
              for(l in 1:t){
                lamd.cand[i,1:J[l],l]<- lam0[l]*exp(-dtmp[1:J[l],l]*dtmp[1:J[l],l]/(2*sigma*sigma))
              }
            }else{
              for(l in 1:t){
                lamd.cand[i,1:J[l],l]<- lam0[l]*exp(-dtmp[1:J[l],l]*dtmp[1:J[l],l]/(2*sigma[l]*sigma[l]))
              }
            }
            if(obstype=="bernoulli"){
              pd.cand[i,1:J[l],]=1-exp(-lamd.cand[i,1:J[l],])
            }
            ll.y.cand[i,1:J[l],l]=ll.y[i,1:J[l],l]
            for(l in 1:t) {
              if(z[i,l]==0)
                next
              if(obstype=="bernoulli"){
                ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,], pd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
              }else{
                ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
              }
            }
            if(runif(1) < exp(sum(ll.y.cand[i,1:J[l],]) -sum(ll.y[i,1:J[l],]))){
              s1[i, ] <- Scand
              D[i,1:J[l], ] <- dtmp
              lamd[i,1:J[l], ] <- lamd.cand[i,1:J[l],]
              if(obstype=="bernoulli"){
                pd[i,1:J[l],]=pd.cand[i,1:J[l],]
              }
              ll.y[i,1:J[l],]=ll.y.cand[i,1:J[l],]
            }
          }
        }

      }else if(ACtype=="markov"){#Markov ACs
        for (i in 1:M){
          for(l in 1:t){
            Scand=c(rnorm(1, s2[i,l,1], proppars$s2x), rnorm(1, s2[i,l,2], proppars$s2y))
            if(usedSS){
              dists=sqrt((Scand[1]-dSS[,1])^2+(Scand[2]-dSS[,2])^2)
              Scand=dSS[which(dists==min(dists)),]
            }
            if(useverts==FALSE){
              inbox=Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
            }else{
              # inbox=inout(Scand,vertices)
              inbox=any(unlist(lapply(vertices,function(x){inout(Scand,x)})))
            }
            if(inbox) {
              dtmp=sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
              if(length(lam0)==1&length(sigma==1)){
                lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
              }else if(length(lam0)==t&length(sigma)==1){
                lamd.cand[i,1:J[l],l]<- lam0[l]*exp(-dtmp*dtmp/(2*sigma*sigma))
              }else{
                lamd.cand[i,1:J[l],l]<- lam0[l]*exp(-dtmp*dtmp/(2*sigma[l]*sigma[l]))
              }
              if(obstype=="bernoulli"){
                pd.cand[i,1:J[l],l]=1-exp(-lamd.cand[i,1:J[l],l])
                ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,], pd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
              }else{
                ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
              }
              ll.s2.cand[i,]=ll.s2[i,]
              if(l==1){#only ll.s2[i,1] matters
                #time 1 to 2
                ll.s2.cand[i,1]=dnorm(s2[i,2,1],Scand[1],sigma_t,log=TRUE)+dnorm(s2[i,2,2],Scand[2],sigma_t,log=TRUE)
              }else if(l>1&l<t){#ll.s2[i,l-1] and ll.s2[i,l] matter
                #time l-1 to time l
                ll.s2.cand[i,l-1]=dnorm(Scand[1],s2[i,l-1,1],sigma_t,log=TRUE)+dnorm(Scand[2],s2[i,l-1,2],sigma_t,log=TRUE)
                #time l to l+1
                ll.s2.cand[i,l]=dnorm(s2[i,l+1,1],Scand[1],sigma_t,log=TRUE)+dnorm(s2[i,l+1,2],Scand[2],sigma_t,log=TRUE)
              }else{#only ll.s2[i,t-1] matters
                #time t-1 to t
                ll.s2.cand[i,t-1]=dnorm(Scand[1],s2[i,t-1,1],sigma_t,log=TRUE)+dnorm(Scand[2],s2[i,t-1,2],sigma_t,log=TRUE)
              }
              if(runif(1) < exp((sum(ll.y.cand[i,1:J[l],l])+sum(ll.s2.cand[i,])) -(sum(ll.y[i,1:J[l],l])+sum(ll.s2[i,])))){
                s2[i,l,] <- Scand
                D[i,1:nrow(X[[l]]),l] <- dtmp
                lamd[i,,l] <- lamd.cand[i,,l]
                if(obstype=="bernoulli"){
                  pd[i,,l]=pd.cand[i,,l]
                }
                ll.y[i,1:J[l],l]=ll.y.cand[i,1:J[l],l]
                ll.s2[i,]=ll.s2.cand[i,]
              }
            }
          }
        }
        #Update sigma_t
        sigma_t.cand <- rnorm(1,sigma_t,proppars$sigma_t)
        if(sigma_t.cand > 0){
          for(l in 2:t){
            for(i in 1:M){
              ll.s2.cand[i,l-1]=dnorm(s2[i,l,1],s2[i,l-1,1],sigma_t.cand,log=TRUE)+dnorm(s2[i,l,2],s2[i,l-1,2],sigma_t.cand,log=TRUE)
            }
          }
          if (runif(1) < exp(sum(ll.s2.cand) - sum(ll.s2))) {
            sigma_t=sigma_t.cand
            ll.s2=ll.s2.cand
          }
        }

      }else{#independent ACS
        for(l in 1:t){
          for (i in 1:M) {
            # if(z[i,l]==0) next
            Scand <- c(rnorm(1, s2[i,l,1], proppars$s2x), rnorm(1, s2[i,l,2], proppars$s2y))
            if(usedSS){
              dists=sqrt((Scand[1]-dSS[,1])^2+(Scand[2]-dSS[,2])^2)
              Scand=dSS[which(dists==min(dists)),]
            }
            if(useverts==FALSE){
              inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
            }else{
              # inbox=inout(Scand,vertices)
              inbox=any(unlist(lapply(vertices,function(x){inout(Scand,x)})))
            }
            if (inbox) {
              dtmp<- sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
              if(length(lam0)==1&length(sigma==1)){
                lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
              }else if(length(lam0)==t&length(sigma)==1){
                lamd.cand[i,1:J[l],l]<- lam0[l]*exp(-dtmp*dtmp/(2*sigma*sigma))
              }else{
                lamd.cand[i,1:J[l],l]<- lam0[l]*exp(-dtmp*dtmp/(2*sigma[l]*sigma[l]))
              }
              if(obstype=="bernoulli"){
                pd.cand[i,1:J[l],l]=1-exp(-lamd.cand[i,,l])
                ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,], pd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
              }else{
                ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
              }
              if(runif(1) < exp(sum(ll.y.cand[i,1:J[l],l]) -sum(ll.y[i,1:J[l],l]))){
                s2[i,l, ] <- Scand
                D[i,1:J[l],l] <- dtmp
                lamd[i,1:J[l], l] <- lamd.cand[i,1:J[l],l]
                if(obstype=="bernoulli"){
                  pd[i,1:J[l],l]=pd.cand[i,1:J[l],l]
                }
                ll.y[i,1:J[l],l]=ll.y.cand[i,1:J[l],l]
              }
            }
          }
        }
      }

      #Do we record output on this iteration?
      if(iter>(nburn)&iter%%nthin==0){
        s1xout[idx,]<- s1[,1]
        s1yout[idx,]<- s1[,2]
        zout[idx,,]<- z
        if(ACtype%in%c("metamu","metamu2","markov")){
          out[idx,]<- c(lam0,sigma ,gamma,phi,N,sigma_t)
          s2xout[idx,,]<- s2[,,1]
          s2yout[idx,,]<- s2[,,2]
        }else if (ACtype=="independent"){
          out[idx,]<- c(lam0,sigma ,gamma,phi,N)
          s2xout[idx,,]<- s2[,,1]
          s2yout[idx,,]<- s2[,,2]
        }else{
          out[idx,]<- c(lam0,sigma ,gamma,phi,N)
        }
        idx=idx+1
      }
    }  # end of MCMC algorithm

    if(keepACs==TRUE){
      if(ACtype%in%c("metamu2","metamu","markov","independent")){
        list(out=out, s1xout=s1xout, s1yout=s1yout,s2xout=s2xout, s2yout=s2yout, zout=zout)
      }else{
        list(out=out, s1xout=s1xout, s1yout=s1yout, zout=zout)
      }
    }else{
      list(out=out)
    }
  }

