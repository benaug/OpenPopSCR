SCRmcmcOpensexRcpp <-
  function(data,niter=2400,nburn=1200, nthin=5,M = 200, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),
           jointZ=TRUE,keepACs=TRUE,ACtype="fixed",obstype="bernoulli",dSS){
    library(abind)
    t=dim(data$y)[3]
    y<-data$y
    X<-data$X
    #make sure X list elements are matrices
    for(i in 1:length(X)){
      X[[i]]=as.matrix(X[[i]])
    }
    if(length(X)!=t){
      stop("must input traps for each year")
    }
    if(!is.na(dSS[1])&"vertices"%in%names(data)){
      rem=which(names(data)=="vertices")
      data[[rem]]=NULL
      warning("Discarding vertices since dSS supplied")
    }
    sex=data$sex

    J<-data$J
    maxJ=max(J)
    K<-data$K
    n=dim(data$y)[1]
    n2d<-colSums(apply(data$y,c(3),rowSums)>0)
    if(length(sex)!=n){
      stop("sex must have an entry for each captured individual")
    }
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
    tf=vector("list",t)
    for(l in 1:t){
      tf[[l]]=rep(K[l],nrow(X[[l]]))
    }

    ##pull out initial values
    lam0<- inits$lam0
    sigma<- inits$sigma
    sigma_t=inits$sigma_t
    gamma=inits$gamma
    phi=inits$phi
    psi=inits$psi
    psex=inits$psex
    if(is.null(psex)){
      stop("must supply initial value for psex")
    }

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
      if(is.null(sigma_t)){
        stop("must supply inits$sigma_t if ACtype is metamu or markov")
      }
    }
    if(is.null(proppars$sex)){
      stop("must enter a proppar for sex update")
    }
    sexparms=vector("list")
    if(length(lam0)==1){
      sexparms$lam0="fixed"
    }else{
      sexparms$lam0="sex"
    }
    if(length(sigma)==1){
      sexparms$sigma="fixed"
    }else{
      sexparms$sigma="sex"
    }
    if(length(phi)==1){
      sexparms$phi="fixed"
      # phi=rep(phi,2)
    }else{
      sexparms$phi="sex"
    }
    if(length(gamma)==1){
      sexparms$gamma="fixed"
      # gamma=rep(gamma,2)
    }else{
      sexparms$gamma="sex"
    }
    if(length(sigma_t)==1){
      sexparms$sigma_t="fixed"
      sigma_t=rep(sigma_t,2)
      proppars$sigma_t=rep(proppars$sigma_t,2)
    }else{
      sexparms$sigma_t="sex"
    }


    #augment data
    y<- abind(y,array(0, dim=c( M-dim(y)[1],maxJ, t)), along=1)
    known.vector=c(rep(1,data$n),rep(0,M-data$n))
    sex=c(sex,rep(NA,M-data$n))
    known.sex=1*(!is.na(sex))
    sex[is.na(sex)]=rbinom(sum(is.na(sex)),1,0.5)+1
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
    if(length(gamma)==2){
      gammasim=gamma
    }else{
      gammasim=rep(gamma,2)
    }
    if(sexparms$phi=="sex"){
      phiMuse=rep(phi[1],t-1)
      phiFuse=rep(phi[2],t-1)
    }else{
      phiMuse=rep(phi,t-1)
      phiFuse=rep(phi,t-1)
    }


    for(i in 2:t){
      #Recruitment males first
      nrecruitM=round(sum(z[,i-1])*gammasim[1])#how many males hould we recruit?
      recruitsM=which(z[1:n,i-1]==0&z[1:n,i]==1&sex[1:n]==1)#who is already recruited based on y constraints?
      a[which(z[,i]==1),i:t]=0 #anyone alive in time t can't be recruited
      nleftM=nrecruitM-length(recruitsM)
      if(nleftM>0){
        candsM=setdiff(which(a[,i-1]==1&sex==1),recruitsM)#Who is available to recruit?
        candsM=candsM[!candsM%in%1:n]#Don't mess with known guys
        if(nleftM>length(candsM)){
          nleftM=length(candsM)
        }
        if(length(candsM)>1){
          pick=sample(candsM,nleftM)
        }else if(length(candsM)==1&nleftM==1){
          pick=candsM
        }else{
          next
        }
        z[pick,i]=1
        a[pick,i:t]=0 #no longer available for recruit on any occasion
      }else{
        warning("Can't initialize all E[recruits] given initial value for gamma. Should probably raise M?")
      }
      #recruit females
      nrecruitF=round(sum(z[,i-1])*gammasim[2])#how many males hould we recruit?
      recruitsF=which(z[1:n,i-1]==0&z[1:n,i]==1&sex[1:n]==2)#who is already recruited based on y constraints?
      # a[which(z[,i]==1),i:t]=0 #anyone alive in time t can't be recruited
      nleftF=nrecruitF-length(recruitsF)
      if(nleftF>0){
        candsF=setdiff(which(a[,i-1]==1&sex==2),recruitsF)#Who is available to recruit?
        candsF=candsF[!candsF%in%1:n]#Don't mess with known guys
        if(nleftF>length(candsF)){
          nleftF=length(candsF)
        }
        if(length(candsF)>1){
          pick=sample(candsF,nleftF)
        }else if(length(candsF)==1&nleftF==1){
          pick=candsF
        }else{
          next
        }
        z[pick,i]=1
        a[pick,i:t]=0 #no longer available for recruit on any occasion
      }else{
        warning("Can't initialize all E[recruits] given initial value for gamma. Should probably raise M?")
      }
      #Survival
      nliveM=rbinom(1,sum(z[,i-1]==1&sex==1),phiMuse[i-1])
      nliveF=rbinom(1,sum(z[,i-1]==1&sex==2),phiFuse[i-1])
      liversM=which(z[1:n,i-1]==1&z[1:n,i]==1&sex[1:n]==1)
      liversF=which(z[1:n,i-1]==1&z[1:n,i]==1&sex[1:n]==2)
      nleftM=nliveM-length(liversM)
      nleftF=nliveF-length(liversF)
      a[liversM,i]=0
      a[liversF,i]=0
      if(nleftM>0){
        candsM=which(apply(matrix(z[,i:t],nrow=M),1,sum)==0&z[,i-1]==1&sex==1)
        candsM=candsM[!candsM%in%1:n]
        if(nleftM>length(candsM)){
          nleftM=length(candsM)
        }
        if(length(candsM)>1){
          pick=sample(candsM,nleftM)
        }else if(length(candsM)==1&nleftM==1){
          pick=candsM
        }else{
          next
        }
        z[pick,i]=1
        a[pick,i:t]=0
      }
      if(nleftF>0){
        candsF=which(apply(matrix(z[,i:t],nrow=M),1,sum)==0&z[,i-1]==1&sex==2)
        candsF=candsF[!candsF%in%1:n]
        if(nleftF>length(candsF)){
          nleftF=length(candsF)
        }
        if(length(candsF)>1){
          pick=sample(candsF,nleftF)
        }else if(length(candsF)==1&nleftF==1){
          pick=candsF
        }else{
          next
        }
        z[pick,i]=1
        a[pick,i:t]=0
      }
    }
    N=colSums(z)
    Nm=colSums(z==1&sex==1)
    Nf=colSums(z==1&sex==2)

    ll.z=matrix(0, M, t)
    Ez=Ez.cand=matrix(NA, M, t-1)
    ll.z[,1]=dbinom(z[,1], 1, psi, log=TRUE)
    gamma.primeM=gamma.primeF=numeric(t-1)
    if(sexparms$gamma=="sex"){
      gammaMuse=rep(gamma[1],t-1)
      gammaFuse=rep(gamma[2],t-1)
    }else{
      gammaMuse=gammaFuse=rep(gamma,t-1)
    }
    for(l in 2:t) {
      gamma.primeM[l-1]=N[l-1]*gammaMuse[l-1] / sum(a[sex==1,l-1])
      gamma.primeF[l-1]=N[l-1]*gammaFuse[l-1] / sum(a[sex==2,l-1])
      Ez[sex==1,l-1]=z[sex==1,l-1]*phiMuse[l-1] + a[sex==1,l-1]*gamma.primeM[l-1]
      Ez[sex==2,l-1]=z[sex==2,l-1]*phiFuse[l-1] + a[sex==2,l-1]*gamma.primeF[l-1]
      ll.z[,l]=dbinom(z[,l], 1, Ez[,l-1], log=TRUE)
    }
    if(any(c(gamma.primeM,gamma.primeF)>1)){
      stop("raise M or lower gamma inits")
    }
    ll.z.cand=ll.z
    gamma.primeM.cand=gamma.primeM
    gamma.primeF.cand=gamma.primeF
    #Optimize starting locations given where they are trapped. Initalizing s1 and s2 at same locs
    s1<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    idx=which(known.vector==1) #switch for those actually caught
    for(i in idx){
      trps=matrix(0,nrow=0,ncol=2)
      for(j in 1:t){ #loop over t to get all cap locs
        trps<- rbind(trps,X[[j]][y[i,,j]>0,1:2])
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
            trps<- X[[l]][y[i,,l]>0,1:2]
            if(is.matrix(trps)){
              s2[i,l,]<- c(mean(trps[,1]),mean(trps[,2]))
            }else{
              s2[i,l,]=trps
            }
          }else{
            inside=FALSE
            while(inside==FALSE){
              s2[i,l,]=c(rnorm(1,s1[i,1],sigma_t[sex[i]]),rnorm(1,s1[i,2],sigma_t[sex[i]]))
              # inside=inout(s2[i,l,],vertices)
              inside=any(unlist(lapply(vertices,function(x){inout(s2[i,l,],x)})))
            }
          }
        }
      }
      if(ACtype%in%c("metamu","metamu2")){
        ll.s2=(dnorm(s2[,,1],s1[,1],sigma_t[sex],log=TRUE)+dnorm(s2[,,2],s1[,2],sigma_t[sex],log=TRUE))
        ll.s2.cand=ll.s2
      }else if(ACtype=="markov"){
        ll.s2=matrix(NA,nrow=M,ncol=t-1)
        for(l in 2:t){
          ll.s2[,l-1]=(dnorm(s2[,l,1],s2[,l-1,1],sigma_t[sex],log=TRUE)+dnorm(s2[,l,2],s2[,l-1,2],sigma_t[sex],log=TRUE))
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
            trps<- X[[l]][y[i,,l]>0,1:2]
            if(is.matrix(trps)){
              s2[i,l,]<- c(mean(trps[,1]),mean(trps[,2]))
            }else{
              s2[i,l,]=trps
            }
          }else{
            inside=FALSE
            while(inside==FALSE){
              s2[i,l,]=cbind(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
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
    if(length(lam0)==2){
      lam0names=c("lam0M","lam0F")
    }else{
      lam0names="lam0"
    }
    if(length(sigma)==2){
      sigmanames=c("sigmaM","sigmaF")
    }else{
      sigmanames="sigma"
    }
    if(length(gamma)==2){
      gammanames=c("gammaM","gammaF")
    }else{
      gammanames="gamma"
    }
    if(length(phi)==2){
      phinames=c("phiM","phiF")
    }else{
      phinames="phi"
    }
    Nnames=paste("N",1:t,sep="")
    NMnames=paste("Nm",1:t,sep="")
    NFnames=paste("Nf",1:t,sep="")
    if(ACtype%in%c("metamu","metamu2","markov")){
      if(sexparms$sigma_t=="fixed"){
        out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t*3+2)
        colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,NMnames,NFnames,"sigma_t","psex")
      }else{
        out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t*3+3)
        colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,NMnames,NFnames,"sigma_t_M","sigma_t_F","psex")
      }
      s1xout<- s1yout<- matrix(NA,nrow=nstore,ncol=M)
      zout<-array(NA,dim=c(nstore,M,t))
      s2xout<- s2yout<-array(NA,dim=c(nstore,M,t))
    }else if(ACtype=="independent"){
      out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t*3+1)
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,NMnames,NFnames,"psex")
      s1xout<- s1yout<- matrix(NA,nrow=nstore,ncol=M)
      zout<-array(NA,dim=c(nstore,M,t))
      s2xout<- s2yout<-array(NA,dim=c(nstore,M,t))
    }else{
      out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t*3+1)
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,NMnames,NFnames,"psex")
      s1xout<- s1yout<- matrix(NA,nrow=nstore,ncol=M)
      zout<-array(NA,dim=c(nstore,M,t))
    }
    idx=1 #for storing output not recorded every iteration

    D=lamd=ll.y=ll.y.cand=array(NA,dim=c(M,maxJ,t))
    D[is.na(D)]=Inf  #hack to allow years with different J and K to fit in one array
    for(l in 1:t){
      D[,1:nrow(X[[l]]),l]=e2dist(s2[,l,],X[[l]])
      if(length(lam0)==2&length(sigma)==2){
        lamd[,,l]=lam0[sex]*exp(-D[,,l]^2/(2*sigma[sex]*sigma[sex]))
      }else if(length(lam0)==1&length(sigma)==2){
        lamd[,,l]=lam0*exp(-D[,,l]^2/(2*sigma[sex]*sigma[sex]))
      }else if(length(lam0)==2&length(sigma)==1){
        lamd[,,l]=lam0[sex]*exp(-D[,,l]^2/(2*sigma*sigma))
      }else{
        lamd[,,l]=lam0*exp(-D[,,l]^2/(2*sigma*sigma))
      }
    }
    #Calculate ll for observation model
    if(obstype=="bernoulli"){
      pd=pd.cand=1-exp(-lamd)
      for(l in 1:t){
        ll.y[,,l]= dbinom(y[,,l],tf[[l]],pd[,,l]*z[,l],log=TRUE)
      }
    }else if(obstype=="poisson"){
      for(l in 1:t){
        ll.y[,,l]= dpois(y[,,l],tf[[l]]*lamd[,,l]*z[,l],log=TRUE)
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
      EzpossibleM=EzpossibleF=matrix(NA,nrow=nzpossible,ncol=t-1)
      ll_z_possibleM=ll_z_possibleF=matrix(NA,nrow=nzpossible,ncol=t)
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
    #sex ll
    nmale=sum(sex[known.sex==1&z[,1]==1]==1)
    nsex=length(sex[known.sex==1&z[,1]==1])
    ll.sex=rep(NA,M)
    for(i in 1:M){
      ll.sex[i]=dbinom(sex[i]-1,1,psex,log=TRUE)
    }
    ll.sex.cand=ll.sex

    Xidx=unlist(lapply(X,dim))[seq(1,2*length(X)-1,2)]
    Xcpp=array(NA,dim=c(t,max(Xidx),2))
    for(l in 1:t){
      for(j in 1:Xidx[l]){
        Xcpp[l,j,]=as.numeric(X[[l]][j,])
      }
    }
    each=unlist(lapply(inits,length))[1:4]
    npar=sum(each)+3*t+1#added one for psex
    if(ACtype%in%c("metamu","metamu2","markov")){
      if(length(sigma_t)==1){
        npar=npar+1
      }else{
        npar=npar+2
      }
    }
    if(length(lam0)==1){
      lam0=rep(lam0,2)
    }
    if(length(sigma)==1){
      sigma=rep(sigma,2)
    }
    if(length(phi)==1){
      phi=rep(phi,2)
    }
    if(length(gamma)==1){
      gamma=rep(gamma,2)
    }
    if(length(sigma_t)==1){
      sigma_t=rep(sigma_t,2)
    }
    #So these aren't modified by Rcpp
    lam0in=lam0
    sigmain=sigma
    gammain=gamma
    phiin=phi
    psexin=psex
    if(ACtype=="fixed"){
      ACtype=1
    }else if(ACtype=="metamu"){
      ACtype=2
    }else if(ACtype=="markov"){
      ACtype=3
    }else if(ACtype=="independent"){
      ACtype=4#independent
    }else{#metamu2
      ACtype=5
    }
    if(obstype=="bernoulli"){
      obstype2=1
    }else{
      obstype2=2
    }
    tf2=matrix(NA,nrow=maxJ,ncol=t)
    for(l in 1:t){
      tf2[1:J[l],l]=tf[[l]]
    }
    if(!usedSS){
      dSS=matrix(0.5,nrow=2,ncol=2)#dummy to fool rcpp
    }
    if(is.null(sigma_t)){
      sigma_t=0.5 #dummy for Rcpp
    }
    sexparmsin=rep(FALSE,5)
    for(i in 1:5){
      if(sexparms[i]=="sex"){
        sexparmsin[i]=TRUE
      }
    }
    store=mcmc_Open_sex(lam0in,  sigmain,  gammain, gamma.primeM,gamma.primeF,  phiin, psexin,D,lamd, y, z, a,s1,s2,
                    ACtype, useverts, vertices, xlim, ylim, sex,known.matrix,Xidx, Xcpp, K, Ez,  psi,
                    N, proppars$lam0, proppars$sigma, proppars$propz,  proppars$gamma, proppars$s1x,  proppars$s1y,
                    proppars$s2x,proppars$s2y,proppars$sigma_t,proppars$sex,sigma_t,niter,nburn,nthin,npar,each,jointZ,
                    zpossible,apossible,cancel,obstype2,tf2,dSS,usedSS,sexparmsin,which(known.sex==0)-1)

    out=store[[1]]
    s1xout=store[[2]]
    s1yout=store[[3]]
    s2xout=store[[4]]
    s2yout=store[[5]]
    zout=store[[6]]
    warn=store[[7]]
    # aout=store[[8]]
    # llzout=store[[8]]
    # storeupz=store[[14]]
    # storeswapz=store[[15]]
    if(warn>0){
      warning(paste("gamma proposal led to gamma prime >1",warn,"times. May want to raise M"))
    }


    # name out
    if(sexparmsin[1]){
      lam0names=c("lam0M","lam0F")
    }else{
      lam0names="lam0"
    }
    if(sexparmsin[2]){
      sigmanames=c("sigmaM","sigmaF")
    }else{
      sigmanames="sigma"
    }
    if(sexparmsin[3]){
      gammanames=c("gammaM","gammaF")
    }else{
      gammanames="gamma"
    }
    if(sexparmsin[4]){
      phinames=c("phiM","phiF")
    }else{
      phinames="phi"
    }
    if(sexparmsin[5]){
      sigmatnames=c("sigma_tM","sigma_tF")
    }else{
      sigmatnames="sigma_t"
    }
    Nnames=paste("N",1:t,sep="")
    Nmnames=paste("Nm",1:t,sep="")
    Nfnames=paste("Nf",1:t,sep="")
    if(ACtype%in%c(2,3,5)){
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,Nmnames,Nfnames,sigmatnames,"psex")
    }else{
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,Nmnames,Nfnames,"psex")
    }
    if(keepACs==TRUE){
      if(ACtype%in%c(2,3,4,5)){
        list(out=out, s1xout=s1xout, s1yout=s1yout,s2xout=s2xout, s2yout=s2yout, zout=zout)
      }else{
        list(out=out, s1xout=s1xout, s1yout=s1yout, zout=zout)
      }
    }else{
      list(out=out,zout=zout)
    }
  }

