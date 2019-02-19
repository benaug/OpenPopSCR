SCRmcmcOpensex <-
  function(data,niter=2400,nburn=1200, nthin=5,M = 200, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),
           jointZ=TRUE,storeLatent=TRUE,ACtype="fixed",obstype="bernoulli",dSS=NA,dualACup=FALSE){

   library(abind)
    t=dim(data$y)[3]
    y<-data$y
    X<-data$X
    if("primary"%in%names(data)){
      primary=data$primary
      if(primary[1]==0){
        stop("first primary period must be a 1 (data recorded)")
      }
    }else{
      primary=rep(1,t)
    }

    #make sure X list elements are matrices
    for(l in 1:length(X)){
      if(primary[l]==1){
        X[[l]]=as.matrix(X[[l]])
      }else{
        for(l in 1:t){
          if(primary[l]==0){
            X[[l]]=matrix(nrow=0,ncol=0)
          }
        }
      }
    }
    dSS=as.matrix(dSS)
    if(!is.na(dSS[1])){
      usedSS=TRUE
      if("vertices"%in%names(data)){
        rem=which(names(data)=="vertices")
        data[[rem]]=NULL
        warning("Discarding vertices since dSS supplied")
      }
      NdSS=nrow(dSS)
      useverts=FALSE
    }else{
      if(ACtype=="markov2"){
        stop("Must enter dSS for markov2")
      }
      usedSS=FALSE
    }
    sex=data$sex
    J<-data$J
    maxJ=max(J,na.rm=TRUE)
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
      minmax=array(NA,dim=c(sum(primary==1),2,2))
      idx=1
      for(l in 1:length(X)){
        if(primary[l]==1){
          minmax[idx,,]=rbind(apply(X[[l]],2,min),apply(X[[l]],2,max))
          idx=idx+1
        }
      }
      xylim=rbind(apply(minmax,3,min),apply(minmax,3,max))
      xlim=xylim[,1]+c(-buff,buff)
      ylim=xylim[,2]+c(-buff,buff)
      vertices=list(rbind(c(xlim[1],ylim[1]),c(xlim[1],ylim[2]),
                          c(xlim[2],ylim[2]),c(xlim[2],ylim[1]),c(xlim[1],ylim[1])))
      useverts=FALSE
    }else{
      stop("user must supply either 'buff' or 'vertices' in data object")
    }
    if(usedSS){
      if(max(dSS[,1])>xlim[2]|min(dSS[,1])<xlim[1]|max(dSS[,2])>ylim[2]|min(dSS[,2])<ylim[1]){
        stop("dSS dimensions exceed xlim or ylim. Change dSS or buff")
      }
    }
    if("tf"%in%names(data)){
      tf=data$tf
      if(length(tf)!=length(X)){
        stop("If using a trap operation file, must input one for each year")
      }
      for(i in 1:length(X)){
        if(primary[l]==1){
          if(is.matrix(tf[[l]])){
            stop("If using trap operation file, must enter vector of number of occasions operational, not matrix of trap by occasion operation")
          }
          if(nrow(X[[l]])!=length(tf[[l]])){
            stop("If using trap operation file, must enter operation for every trap")
          }
        }
      }
    }else{
      tf=vector("list",t)
      for(l in 1:t){
        if(primary[l]==1){
          tf[[l]]=rep(K[l],nrow(X[[l]]))
        }
      }
    }
    #make tf into matrix
    for(l in 1:t){
      if(primary[l]==1){
        tf[[l]]=matrix(rep(tf[[l]],M),ncol=J[l],nrow=M,byrow=TRUE)
      }
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
    if(!ACtype%in%c("fixed","independent","metamu","metamu2","markov","markov2")){
      stop("ACtype must be 'fixed','independent','metamu', 'metamu2','markov', or 'markov2'")
    }
    if(ACtype%in%c("metamu","markov","metamu2","markov2")){
      if(!"sigma_t"%in%names(proppars)){
        stop("must supply proppars$sigma_t if ACtype is metamu or markov")
      }
      if(is.null(sigma_t)){
        stop("must supply inits$sigma_t if ACtype is metamu or markov")
      }
    }
    if(ACtype%in%c("metamu","metamu2","fixed")){
      if(!"s1x"%in%names(proppars)|!"s1y"%in%names(proppars)){
        stop("must supply proppars$s1x and proppars$s1y if ACtype is metamu, metamu2, or fixed")
      }
    }
    if(dualACup){
      if(is.null(proppars$dualAC)){
        stop("If dualACup=TRUE, must specify proppars$dualAC")
      }
    }

    if(is.null(proppars$sex)){
      stop("must enter a proppar for sex update")
    }
    sexparms=vector("list")
    if(length(lam0)==1){
      sexparms$lam0="fixed"
      # lam0=rep(lam0,2)
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
    if(!is.null(sigma_t)){
      if(length(sigma_t)==1){
        sexparms$sigma_t="fixed"
        sigma_t=rep(sigma_t,2)
        proppars$sigma_t=rep(proppars$sigma_t,2)
      }else{
        sexparms$sigma_t="sex"
      }
    }

    #augment data
    y<- abind(y,array(0, dim=c( M-dim(y)[1],maxJ, t)), along=1)
    known.vector=c(rep(1,n),rep(0,M-n))
    sex=c(sex,rep(NA,M-n))
    known.sex=1*(!is.na(sex))
    sex[is.na(sex)]=rbinom(sum(is.na(sex)),1,0.5)+1
    #Initialize z, r, and a consistent with y
    known.matrix=1*(apply(y,c(1,3),sum)>0)#make z consistent with y
    if(t>2){
      for(l in 2:(t-1)){#Turn on zeros with 1's on either side
        known.matrix[known.matrix[,l]==0&known.matrix[,l-1]==1&rowSums(matrix(known.matrix[,(l+1):t],nrow=M))>0,l]=1
      }
    }
    if(sum(known.sex==0)<proppars$sex){
      warning("Fewer unknown sexes than number to update in proppars$sex. Updating all unknown sexes.")
      proppars$sex=sum(known.sex==0)
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
      stop("Raise M, lower gamma, phi and/or psi inits. Ran out of inds to recruit during init.")
    }
    ll.z.cand=ll.z
    gamma.primeM.cand=gamma.primeM
    gamma.primeF.cand=gamma.primeF

    #Initialize ACs
    if(!usedSS){
      #make sure all traps are inside a polygon
      for(l in 1:t){
        if(primary[l]==1){
          for(j in 1:J[l]){
            if(!any(unlist(lapply(vertices,function(x){inout(X[[l]][j,],x)})))){
              stop(paste("Trap",j,"in year",l,"is not in the state space!"))
            }
          }
        }
      }
      #initialize s1
      s1<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
      idx=which(known.vector==1) #switch for those actually caught
      for(i in idx){
        trps=matrix(0,nrow=0,ncol=2)
        for(l in 1:t){ #loop over t to get all cap locs
          if(primary[l]==1){
            trps<- rbind(trps,X[[l]][which(y[i,,l]>0),1:2])
          }
        }
        trps=as.matrix(trps,ncol=2)
        if(nrow(trps)>1){
          s1[i,]<- trps[1,]
        }else{
          s1[i,]=trps
        }
        # s1[i,]<- c(mean(trps[,1]),mean(trps[,2]))
      }
      #check to make sure everyone is in polygon
      inside=rep(NA,nrow(s1))
      for(i in 1:nrow(s1)){
        # inside[i]=inout(s1[i,],vertices)
        inside[i]=any(unlist(lapply(vertices,function(x){inout(s1[i,],x)})))
      }
      idx2=which(inside==FALSE)
      if(any(idx2%in%idx)){#shouldn't happen now
        warning("Vertices too complicated for this AC initialization algorithm to provide good starting values.
                  Either hassle Ben to fix it or use a discrete state space")
      }
      if(length(idx2)>0){
        for(i in 1:length(idx2)){
          while(inside[idx2[i]]==FALSE){
            s1[idx2[i],]=c(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            # inside[idx[i]]=inout(s1[idx[i],],vertices)
            inside[idx2[i]]=any(unlist(lapply(vertices,function(x){inout(s1[idx2[i],],x)})))
          }
        }
      }
      #initialize s2
      s2=array(NA,dim=c(M,t,2))
      if(ACtype%in%c("metamu","metamu2","markov")){
        #update s2s for guys captured each year and add noise for uncaptured guys. More consistent with sigma_t>0
        for(l in 1:t){
          idx=which(rowSums(y[,,l])>0) #switch for those actually caught
          for(i in 1:M){
            if(i%in%idx){
              trps<- X[[l]][which(y[i,,l]>0),1:2]
              if(!is.matrix(trps)){
                trps=matrix(trps,ncol=2,nrow=1)
              }
              if(nrow(trps)>1){
                s2[i,l,]<- c(mean(trps[,1]),mean(trps[,2]))
              }else{
                s2[i,l,]=trps
              }
              inside=any(unlist(lapply(vertices,function(x){inout(s2[i,l,],x)})))
              if(!inside){
                if(nrow(trps)>1){
                  s2[i,l,]<- trps[1,]
                }else{
                  s2[i,l,]=trps
                }
              }

            }else{
              if(ACtype%in%c("metamu","metamu2")){#this doesn't work well for markov with larger sigma
                inside=FALSE
                while(inside==FALSE){
                  s2[i,l,]=c(rnorm(1,s1[i,1],sigma_t[sex[i]]),rnorm(1,s1[i,2],sigma_t[sex[i]]))
                  inside=any(unlist(lapply(vertices,function(x){inout(s2[i,l,],x)})))
                }
              }
            }
          }
        }
        if(ACtype=="markov"){#smarter starting values for markov mvmt
          for(i in 1:M){
            if(is.na(s2[i,1,1])){
              s2[i,1,]=s1[i,]
            }
            for(l in 2:t){
              if(is.na(s2[i,l,1])){
                inside=FALSE
                while(inside==FALSE){
                  s2[i,l,1]=rnorm(1,s2[i,l-1,1],sigma_t[sex[i]])
                  s2[i,l,2]=rnorm(1,s2[i,l-1,2],sigma_t[sex[i]])
                  inside=any(unlist(lapply(vertices,function(x){inout(s2[i,l,],x)})))
                }
              }
            }
          }
        }
        #calculate likelihoods
        if(ACtype%in%c("metamu","metamu2")){
          if(ACtype=="metamu2"|(ACtype=="metamu"&useverts)){
            ll.s2=(dnorm(s2[,,1],s1[,1],sigma_t[sex],log=TRUE)+dnorm(s2[,,2],s1[,2],sigma_t[sex],log=TRUE))
          }else{
            ll.s2=log(dnorm(s2[,,1],s1[,1],sigma_t[sex])/
                        (pnorm(xlim[2],s1[,1],sigma_t[sex])-pnorm(xlim[1],s1[,1],sigma_t[sex])))
            ll.s2=ll.s2+log(dnorm(s2[,,2],s1[,2],sigma_t[sex])/
                              (pnorm(ylim[2],s1[,2],sigma_t[sex])-pnorm(ylim[1],s1[,2],sigma_t[sex])))
          }
          ll.s2.cand=ll.s2
        }else if(ACtype%in%c("markov")){
          ll.s2=matrix(NA,nrow=M,ncol=t-1)
          for(l in 2:t){
            if(!useverts){#truncate
              ll.s2[,l-1]=log(dnorm(s2[,l,1],s2[,l-1,1],sigma_t[sex])/
                                (pnorm(xlim[2],s2[,l-1,1],sigma_t[sex])-pnorm(xlim[1],s2[,l-1,1],sigma_t[sex])))
              ll.s2[,l-1]=ll.s2[,l-1]+log(dnorm(s2[,l,2],s2[,l-1,2],sigma_t[sex])/
                                            (pnorm(ylim[2],s2[,l-1,2],sigma_t[sex])-pnorm(ylim[1],s2[,l-1,2],sigma_t[sex])))
            }else{
              ll.s2[,l-1]=(dnorm(s2[,l,1],s2[,l-1,1],sigma_t[sex],log=TRUE)+dnorm(s2[,l,2],s2[,l-1,2],sigma_t[sex],log=TRUE))
            }
          }
          ll.s2.cand=ll.s2
        }
      }else if(ACtype=="independent"){
        for(l in 1:t){
          s2[,l,]=s1
          idx=which(rowSums(y[,,l])>0) #switch for those actually caught
          for(i in 1:M){
            if(i%in%idx){
              trps<- X[[l]][which(y[i,,l]>0),1:2]
              if(!is.matrix(trps)){
                trps=matrix(trps,ncol=2,nrow=1)
              }
              if(nrow(trps)>1){
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
      }else{#fixed
        for(l in 1:t){
          s2[,l,]=s1
        }
      }
    }else{#use dSS
      #not optimal but should work ok for all options
      s1=dSS[sample(1:NdSS,M,replace=TRUE),1:2] #initial locations
      #update if caught at all
      idx=which(rowSums(y)>0) #switch for those actually caught to first trap
      for(i in 1:M){
        if(i%in%idx){
          trps=matrix(0,nrow=0,ncol=2)
          for(l in 1:t){ #loop over t to get all cap locs
            if(primary[l]==1){
              trps<- rbind(trps,X[[l]][which(y[i,,l]>0),1:2])
            }
          }
          if(!is.matrix(trps)){
            trps=matrix(trps,ncol=2,nrow=1)
          }
          if(nrow(trps)>1){#nearest trap
            dist=sqrt((trps[1,1]-dSS[,1])^2+(trps[1,2]-dSS[,2])^2)
          }else{
            dist=sqrt((trps[1]-dSS[,1])^2+(trps[2]-dSS[,2])^2)
          }
          s1[i,]=dSS[which(dist==min(dist))[1],1:2]
        }
      }
      #record s1 dSS
      s1.cell=rep(NA,M)
      for(i in 1:M){
        s1.cell[i]=which(dSS[,1]==s1[i,1]&dSS[,2]==s1[i,2])
      }
      #update s2 in each year if caught in that year
      s2=array(NA,dim=c(M,t,2))
      s2.cell=matrix(NA,nrow=M,ncol=t)
      if(ACtype%in%c("independent","metamu","metamu2","markov","markov2")){
        for(l in 1:t){
          # s2[,l,]=s1
          if(primary[l]==1){
            idx=which(rowSums(y[,,l])>0)
            for(i in 1:M){
              if(i%in%idx){
                trps<- X[[l]][which(y[i,,l]>0),1:2]
                if(!is.matrix(trps)){
                  trps=matrix(trps,ncol=2,nrow=1)
                }
                if(nrow(trps)>1){
                  dist=sqrt((trps[1,1]-dSS[,1])^2+(trps[1,2]-dSS[,2])^2)
                }else{
                  dist=sqrt((trps[1]-dSS[,1])^2+(trps[2]-dSS[,2])^2)
                }
                s2[i,l,]=dSS[which(dist==min(dist))[1],1:2]
              }
            }
          }
        }
      }
      #Get smart values for ids/years not captured and calculate movement likelihoods
      if(ACtype%in%c("metamu","metamu2")){
        for(l in 1:t){
          for(i in 1:M){
            if(is.na(s2[i,l,1])){
              inside=FALSE
              while(inside==FALSE){
                s2[i,l,]=c(rnorm(1,s1[i,1],sigma_t[sex[i]]),rnorm(1,s1[i,2],sigma_t[sex[i]]))
                inside=any(unlist(lapply(vertices,function(x){inout(s2[i,l,],x)})))
              }
            }
          }
        }
        if(ACtype=="metamu2"|(ACtype=="metamu"&useverts)){
          ll.s2=(dnorm(s2[,,1],s1[,1],sigma_t[sex],log=TRUE)+dnorm(s2[,,2],s1[,2],sigma_t[sex],log=TRUE))
        }else{
          ll.s2=log(dnorm(s2[,,1],s1[,1],sigma_t[sex])/
                      (pnorm(xlim[2],s1[,1],sigma_t[sex])-pnorm(xlim[1],s1[,1],sigma_t[sex])))
          ll.s2=ll.s2+log(dnorm(s2[,,2],s1[,2],sigma_t[sex])/
                            (pnorm(ylim[2],s1[,2],sigma_t[sex])-pnorm(ylim[1],s1[,2],sigma_t[sex])))
        }
        ll.s2.cand=ll.s2
      }else if(ACtype%in%c("markov")){
        ll.s2=matrix(NA,nrow=M,ncol=t-1)
        for(i in 1:M){
          if(is.na(s2[i,1,1])){
            s2[i,1,]=s1[i,]
          }
          for(l in 2:t){
            if(is.na(s2[i,l,1])){
              inside=FALSE
              while(inside==FALSE){
                s2[i,l,1]=rnorm(1,s2[i,l-1,1],sigma_t[sex[i]])
                s2[i,l,2]=rnorm(1,s2[i,l-1,2],sigma_t[sex[i]])
                inside=any(unlist(lapply(vertices,function(x){inout(s2[i,l,],x)})))
              }
            }
          }
        }
        for(l in 2:t){
          if(!useverts){#truncate
            ll.s2[,l-1]=log(dnorm(s2[,l,1],s2[,l-1,1],sigma_t[sex])/
                              (pnorm(xlim[2],s2[,l-1,1],sigma_t[sex])-pnorm(xlim[1],s2[,l-1,1],sigma_t[sex])))
            ll.s2[,l-1]=ll.s2[,l-1]+log(dnorm(s2[,l,2],s2[,l-1,2],sigma_t[sex])/
                                          (pnorm(ylim[2],s2[,l-1,2],sigma_t[sex])-pnorm(ylim[1],s2[,l-1,2],sigma_t[sex])))
          }else{
            ll.s2[,l-1]=(dnorm(s2[,l,1],s2[,l-1,1],sigma_t[sex],log=TRUE)+dnorm(s2[,l,2],s2[,l-1,2],sigma_t[sex],log=TRUE))
          }
        }
        ll.s2.cand=ll.s2
      }else if(ACtype%in%"markov2"){#get smarter s2 starts
        ll.s2=matrix(NA,nrow=M,ncol=t-1)
        for(i in 1:M){
          if(is.na(s2[i,1,1])){
            s2[i,1,]=s1[i,]
          }
          #snap
          dists2=sqrt((s2[i,1,1]-dSS[,1])^2+(s2[i,1,2]-dSS[,2])^2)
          #record cell
          s2.cell[i,1]=which(dists2==min(dists2))
        }
        for(l in 2:t){
          dists=e2dist(s2[,l-1,],dSS[,1:2])#cells each ind can choose
          for(i in 1:M){
            probs=dexp(dists[i,],1/sigma_t[sex[i]])
            probs=probs/sum(probs)
            if(is.na(s2[i,l,1])){#only move if not captured
              s2.cell[i,l]=sample(1:NdSS,1,prob=probs)#overwrite above to be consistent with sigma_t
              s2[i,l,]=dSS[s2.cell[i,l],1:2]
            }else{
              #snap
              dists2=sqrt((s2[i,l,1]-dSS[,1])^2+(s2[i,l,2]-dSS[,2])^2)
              #record cell
              s2.cell[i,l]=which(dists2==min(dists2))
              s2[i,l,]=dSS[s2.cell[i,l],1:2]
            }
            pick=rep(0,NdSS)
            pick[s2.cell[i,l]]=1#you just picked this one
            ll.s2[i,l-1]=dmultinom(pick,1,probs,log=TRUE)
          }
        }
        ll.s2.cand=ll.s2

      }else if(ACtype%in%c("fixed")){
        for(l in 1:t){
          s2[,l,]=s1
        }
      }else if(ACtype%in%c("independent")){
        cap=1*(!is.na(s2[,1,1]))
        for(i in 1:n){
          for(l in 1:t){
            if(is.na(s2[i,l,1])&cap[i]==1){
              s2[i,l,]=s2[i,l-1,]
            }else if(is.na(s2[i,l,1])&cap[i]==0){
              s2[i,l,]=s2[i,which(!is.na(s2[i,,1]))[1],]
              cap[i]=1
            }
          }
        }
        for(i in (n+1):M){
          for(l in 1:t){
            s2[i,l,]=s1[i,]
          }
        }
      }
      if(ACtype%in%c("fixed","independent","metamu","metamu2","markov")){
        #record s2 cells
        for(l in 1:t){
          for(i in 1:M){
            #snap
            dists=sqrt((s2[i,l,1]-dSS[,1])^2+(s2[i,l,2]-dSS[,2])^2)
            #record cell
            s2.cell[i,l]=which(dists==min(dists))
            s2[i,l,]=dSS[s2.cell[i,l],1:2]
          }
        }
      }
      ##precalculate distance matrix
      distances=matrix(NA,nrow=NdSS,ncol=NdSS)
      for(i in 1:NdSS){
        distances[i,]=sqrt((dSS[i,1]-dSS[,1])^2+(dSS[i,2]-dSS[,2])^2)
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
    if(ACtype%in%c("markov","markov2")){
      if(sexparms$sigma_t=="fixed"){
        out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t*3+2+1)
        colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,NMnames,NFnames,"sigma_t","psex","psi")
      }else{
        out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t*3+3+1)
        colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,NMnames,NFnames,"sigma_t_M","sigma_t_F","psex","psi")
      }
      if(storeLatent){
        zout<-array(NA,dim=c(nstore,M,t))
        s2xout<- s2yout<-array(NA,dim=c(nstore,M,t))
        sexout=matrix(NA,nrow=nstore,ncol=M)
      }
    }else if(ACtype%in%c("metamu","metamu2")){
      if(sexparms$sigma_t=="fixed"){
        out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t*3+2+1)
        colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,NMnames,NFnames,"sigma_t","psex","psi")
      }else{
        out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t*3+3+1)
        colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,NMnames,NFnames,"sigma_t_M","sigma_t_F","psex","psi")
      }
      if(storeLatent){
        s1xout<- s1yout<- matrix(NA,nrow=nstore,ncol=M)
        zout<-array(NA,dim=c(nstore,M,t))
        s2xout<- s2yout<-array(NA,dim=c(nstore,M,t))
        sexout=matrix(NA,nrow=nstore,ncol=M)
      }
    }else if(ACtype=="independent"){
      out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t*3+2)
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,NMnames,NFnames,"psex","psi")
      if(storeLatent){
        zout<-array(NA,dim=c(nstore,M,t))
        s2xout<- s2yout<-array(NA,dim=c(nstore,M,t))
        sexout=matrix(NA,nrow=nstore,ncol=M)
      }
    }else{
      out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t*3+2)
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,NMnames,NFnames,"psex","psi")
      if(storeLatent){
        s1xout<- s1yout<- matrix(NA,nrow=nstore,ncol=M)
        zout<-array(NA,dim=c(nstore,M,t))
        sexout=matrix(NA,nrow=nstore,ncol=M)
      }
    }
    idx=1 #for storing output not recorded every iteration

    D=lamd=ll.y=ll.y.cand=array(0,dim=c(M,maxJ,t))
    # D[is.na(D)]=Inf  #hack to allow years with different J and K to fit in one array
    for(l in 1:t){
      if(primary[l]==1){
        D[,1:J[l],l]=e2dist(s2[,l,],X[[l]])
        if(length(lam0)==2&length(sigma)==2){
          lamd[,1:J[l],l]=lam0[sex]*exp(-D[,1:J[l],l]^2/(2*sigma[sex]*sigma[sex]))
        }else if(length(lam0)==1&length(sigma)==2){
          lamd[,1:J[l],l]=lam0*exp(-D[,1:J[l],l]^2/(2*sigma[sex]*sigma[sex]))
        }else if(length(lam0)==2&length(sigma)==1){
          lamd[,1:J[l],l]=lam0[sex]*exp(-D[,1:J[l],l]^2/(2*sigma*sigma))
        }else{
          lamd[,1:J[l],l]=lam0*exp(-D[,1:J[l],l]^2/(2*sigma*sigma))
        }
      }
    }
    #Calculate ll for observation model
    if(obstype=="bernoulli"){
      pd=pd.cand=1-exp(-lamd)
      for(l in 1:t){
        if(primary[l]==1){
          ll.y[,1:J[l],l]= dbinom(y[,1:J[l],l],tf[[l]],pd[,1:J[l],l]*z[,l],log=TRUE)
        }
      }
    }else if(obstype=="poisson"){
      for(l in 1:t){
        if(primary[l]==1){
          ll.y[,1:J[l],l]= dpois(y[,1:J[l],l],tf[[l]]*lamd[,1:J[l],l]*z[,l],log=TRUE)
        }
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

    ######################################################################################
    for(iter in 1:niter){
      ll.y.sum=sum(ll.y)
      # Update lam0
      if(sexparms$lam0=="sex"){#if lam0 is sex-specific
        for(i in 1:2){
          lam0.cand<- rnorm(1,lam0[i],proppars$lam0[i])
          if(lam0.cand > 0){
            lamd.cand=lamd
            for(l in 1:t){
              if(primary[l]==1){
                if(length(sigma)==2){#if sigma is sex specific
                  lamd.cand[sex==i,1:J[l],l]<- lam0.cand*exp(-D[sex==i,1:J[l],l]^2/(2*sigma[i]*sigma[i]))
                }else{#fixed sigma
                  lamd.cand[sex==i,1:J[l],l]<- lam0.cand*exp(-D[sex==i,1:J[l],l]^2/(2*sigma*sigma))
                }
                if(obstype=="bernoulli"){
                  pd.cand[,1:J[l],l]=1-exp(-lamd.cand[,1:J[l],l])
                  ll.y.cand[,1:J[l],l]= dbinom(y[,1:J[l],l],tf[[l]],pd.cand[,1:J[l],l]*z[,l],log=TRUE)
                }else{
                  ll.y.cand[,1:J[l],l]= dpois(y[,1:J[l],l],tf[[l]]*lamd.cand[,1:J[l],l]*z[,l],log=TRUE)
                }
              }
            }
            ll.y.cand.sum=sum(ll.y.cand)
            if(runif(1) < exp(ll.y.cand.sum-ll.y.sum)){
              lam0[i]<- lam0.cand
              lamd=lamd.cand
              if(obstype=="bernoulli"){
                pd=pd.cand
              }
              ll.y=ll.y.cand
              ll.y.sum=ll.y.cand.sum
            }
          }
        }
      }else{#fixed lam0
        lam0.cand<- rnorm(1,lam0,proppars$lam0)
        if(lam0.cand > 0){
          lamd.cand=lamd
          for(l in 1:t){
            if(primary[l]==1){
              if(length(sigma)==t){#if sigma is sex specific
                lamd.cand[sex==1,1:J[l],l]<- lam0.cand*exp(-D[sex==1,1:J[l],l]^2/(2*sigma[1]*sigma[1]))
                lamd.cand[sex==2,1:J[l],l]<- lam0.cand*exp(-D[sex==2,1:J[l],l]^2/(2*sigma[2]*sigma[2]))
              }else{#fixed sigma
                lamd.cand[,1:J[l],l]<- lam0.cand*exp(-D[,1:J[l],l]^2/(2*sigma*sigma))
              }
              if(obstype=="bernoulli"){
                pd.cand[,1:J[l],l]=1-exp(-lamd.cand[,1:J[l],l])
                ll.y.cand[,1:J[l],l]= dbinom(y[,1:J[l],l],tf[[l]],pd.cand[,1:J[l],l]*z[,l],log=TRUE)
              }else{
                ll.y.cand[,1:J[l],l]= dpois(y[,1:J[l],l],tf[[l]]*lamd.cand[,1:J[l],l]*z[,l],log=TRUE)
              }
            }
          }
          ll.y.cand.sum=sum(ll.y.cand)
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
      if(sexparms$sigma=="sex"){ #if sigma is sex-specific
        for(i in 1:2){
          sigma.cand<- rnorm(1,sigma[i],proppars$sigma[i])
          if(sigma.cand > 0){
            lamd.cand=lamd
            for(l in 1:t){
              if(primary[l]==1){
                if(length(lam0)==2){#if lam0 is year specific
                  lamd.cand[sex==i,1:J[l],l]<- lam0[i]*exp(-D[sex==i,1:J[l],l]^2/(2*sigma.cand*sigma.cand))
                }else{#fixed lam0
                  lamd.cand[sex==i,1:J[l],l]<- lam0*exp(-D[sex==i,1:J[l],l]^2/(2*sigma.cand*sigma.cand))
                }
                if(obstype=="bernoulli"){
                  pd.cand[,1:J[l],l]=1-exp(-lamd.cand[,1:J[l],l])
                  ll.y.cand[,1:J[l],l]= dbinom(y[,1:J[l],l],tf[[l]],pd.cand[,1:J[l],l]*z[,l],log=TRUE) #only need to update this year
                }else{
                  ll.y.cand[,1:J[l],l]= dpois(y[,1:J[l],l],tf[[l]]*lamd.cand[,1:J[l],l]*z[,l],log=TRUE) #only need to update this year
                }
              }
            }
            ll.y.cand.sum=sum(ll.y.cand)
            if(runif(1) < exp(ll.y.cand.sum -ll.y.sum)){
              sigma[i]<- sigma.cand
              lamd=lamd.cand
              if(obstype=="bernoulli"){
                pd=pd.cand
              }
              ll.y=ll.y.cand
              ll.y.sum=ll.y.cand.sum
            }
          }
        }
      }else{#fixed sigma
        sigma.cand<- rnorm(1,sigma,proppars$sigma)
        if(sigma.cand > 0){
          for(l in 1:t){
            if(primary[l]==1){
              if(length(lam0)==2){#if lam0 is sex specific
                lamd.cand[sex==1,1:J[l],l]<- lam0[1]*exp(-D[sex==1,1:J[l],l]^2/(2*sigma.cand*sigma.cand))
                lamd.cand[sex==2,1:J[l],l]<- lam0[2]*exp(-D[sex==2,1:J[l],l]^2/(2*sigma.cand*sigma.cand))
              }else{#fixed lam0
                lamd.cand[,1:J[l],l]<- lam0*exp(-D[,1:J[l],l]^2/(2*sigma.cand*sigma.cand))
              }
              if(obstype=="bernoulli"){
                pd.cand[,1:J[l],l]=1-exp(-lamd.cand[,1:J[l],l])
                for(l in 1:t){
                  ll.y.cand[,1:J[l],l]= dbinom(y[,1:J[l],l],tf[[l]],pd.cand[,1:J[l],l]*z[,l],log=TRUE)
                }
              }else{
                for(l in 1:t){
                  ll.y.cand[,1:J[l],l]= dpois(y[,1:J[l],l],tf[[l]]*lamd.cand[,1:J[l],l]*z[,l],log=TRUE)
                }
              }
            }
          }
          ll.y.cand.sum=sum(ll.y.cand)
          if(runif(1) < exp(ll.y.cand.sum - ll.y.sum)){
            sigma<- sigma.cand
            lamd=lamd.cand
            if(obstype=="bernoulli"){
              pd=pd.cand
            }
            ll.y=ll.y.cand
            ll.y.sum=ll.y.cand.sum#dont really need these anymore
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
          gamma.primeM.cand <- gamma.primeM
          gamma.primeF.cand <- gamma.primeF
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
            ll.z.cand[i,1] <- dbinom(z.cand[i,1], 1, psi, log=TRUE)
            #Make object to put N1_prop in to but use current values of the other Ns
            Ntmp=N
            Ntmp[1]=sum(z.cand[,1])
            reject=FALSE
            for(l in 2:t){
              gamma.primeM.cand[l-1]=(Ntmp[l-1]*gammaMuse[l-1]) / sum(a.cand[sex==1,l-1])
              gamma.primeF.cand[l-1]=(Ntmp[l-1]*gammaFuse[l-1]) / sum(a.cand[sex==2,l-1])
              if(gamma.primeM.cand[l-1] > 1 | gamma.primeF.cand[l-1] > 1) { # E(Recruits) must be < nAvailable
                warning("Rejected z due to low M")
                reject=TRUE
                next
              }
              Ez.cand[sex==1,l-1]=z.cand[sex==1,l-1]*phiMuse[l-1] + a.cand[sex==1,l-1]*gamma.primeM.cand[l-1]
              Ez.cand[sex==2,l-1]=z.cand[sex==2,l-1]*phiFuse[l-1] + a.cand[sex==2,l-1]*gamma.primeF.cand[l-1]
              ll.z.cand[,l]=dbinom(z.cand[,l], 1, Ez.cand[,l-1], log=TRUE)
            }
            if(reject==FALSE){
              if(runif(1) < exp((sum(ll.y.cand[i,1:J[1],1])+ll.z.cand[i,1]+sum(ll.z.cand[,2:t]))-(sum(ll.y[i,1:J[1],1])+ll.z[i,1]+sum(ll.z[,2:t])) )) {
                ll.y[i,1:J[1],1] = ll.y.cand[i,1:J[1],1]
                ll.z[i,1]=ll.z.cand[i,1]
                ll.z[,2:t] = ll.z.cand[,2:t]
                Ez = Ez.cand
                gamma.primeM = gamma.primeM.cand
                gamma.primeF = gamma.primeF.cand
                z[i,1] = z.cand[i,1]
                a=a.cand
              }
            }
          }else{#Don't need to modify more than one year
            z1.cand <- z[,1]
            z1.cand[i] <- 1-z[i,1]
            a1.cand <- a[,1]
            a1.cand[i] <- 1-z1.cand[i]
            gamma.primeM.cand[1] <- sum(z1.cand)*gammaMuse[1] / sum(a1.cand[sex==1])
            gamma.primeF.cand[1] <- sum(z1.cand)*gammaFuse[1] / sum(a1.cand[sex==2])
            reject=FALSE
            if(gamma.primeM.cand[1] > 1| gamma.primeF.cand[1] > 1) { # E(Recruits) must be < nAvailable
              warning("Rejected z due to low M")
              reject=TRUE
            }
            ll.z.cand[i,1] <- dbinom(z1.cand[i], 1, psi, log=TRUE)
            Ez.cand[sex==1,1]=z1.cand[sex==1]*phiMuse[1] + a1.cand[sex==1]*gamma.primeM.cand[1]
            Ez.cand[sex==2,1]=z1.cand[sex==2]*phiFuse[1] + a1.cand[sex==2]*gamma.primeF.cand[1]
            ll.z.cand[,2] <- dbinom(z[,2], 1, Ez.cand[,1], log=TRUE)
            if(reject==FALSE){
              if(runif(1) < exp((sum(ll.y.cand[i,1:J[1],1])+ll.z.cand[i,1]+sum(ll.z.cand[,2:t]))-(sum(ll.y[i,1:J[1],1])+ll.z[i,1]+sum(ll.z[,2:t])) )) {
                ll.y[i,1:J[1],1] = ll.y.cand[i,1:J[1],1]
                ll.z[i,1]=ll.z.cand[i,1]
                ll.z[,2] = ll.z.cand[,2]
                Ez[,1] = Ez.cand[,1]
                gamma.primeM[1]= gamma.primeM.cand[1]
                gamma.primeF[1]= gamma.primeF.cand[1]
                z[i,1] = z1.cand[i]
                a[i,1]=a1.cand[i]
              }
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
            if(primary[l]==1){
              if(obstype=="bernoulli"){
                ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,],pd[i,1:J[l],l]*zt.cand[i],log=TRUE)
              }else{
                ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd[i,1:J[l],l]*zt.cand[i],log=TRUE)
              }
            }
            ll.z.cand[i,l] <- dbinom(zt.cand[i], 1, Ez[i,l-1], log=TRUE)
            prior.z=ll.z[i,l]
            prior.z.cand=ll.z.cand[i,l]
            fix1=zt.cand[i]==1&sum(z[i,])==0 #guys never in pop proposed to be turned on
            fix2=sum(z[i,])==1&zt.cand[i]==0&z[i,l]==1 #guys in pop only once and proposed to be turned off
            if((t>3)&(l<t)&(fix1|fix2)){
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
                gamma.primeM.cand[l2]=(Ntmp[l2]*gammaMuse[l2]) / sum(a.cand[sex==1,l2])
                gamma.primeF.cand[l2]=(Ntmp[l2]*gammaFuse[l2]) / sum(a.cand[sex==2,l2])
                if(gamma.primeM.cand[l2] > 1 | gamma.primeF.cand[l2] > 1){
                  reject=TRUE
                  warning("Rejected z due to low M")
                  next
                }
                Ez.cand[sex==1,l2]=z.cand[sex==1,l2]*phiMuse[l2] + a.cand[sex==1,l2]*gamma.primeM.cand[l2]
                Ez.cand[sex==2,l2]=z.cand[sex==2,l2]*phiFuse[l2] + a.cand[sex==2,l2]*gamma.primeF.cand[l2]
                ll.z.cand[,l2+1]=dbinom(z.cand[,l2+1], 1, Ez.cand[,l2], log=TRUE)
                prior.z <- prior.z + sum(ll.z[,l2+1])
                prior.z.cand <- prior.z.cand + sum(ll.z.cand[,l2+1])
              }
            }else{
              if(l<t){ ## NOTE: Don't subset with swapz
                # gamma.prime.cand[l] <- sum(zt.cand)*gammause[l] / sum(at.cand)
                gamma.primeM.cand[l]=(sum(zt.cand)*gammaMuse[l]) / sum(at.cand[sex==1])
                gamma.primeF.cand[l]=(sum(zt.cand)*gammaFuse[l]) / sum(at.cand[sex==2])
                if(gamma.primeM.cand[l] > 1 | gamma.primeF.cand[l] > 1){
                  warning("Rejected z due to low M")
                  next
                }
                # Ez.cand[,l] <- zt.cand*phiuse[l] + at.cand*gamma.prime.cand[l]
                Ez.cand[sex==1,l]=zt.cand[sex==1]*phiMuse[l] + at.cand[sex==1]*gamma.primeM.cand[l]
                Ez.cand[sex==2,l]=zt.cand[sex==2]*phiFuse[l] + at.cand[sex==2]*gamma.primeF.cand[l]
                ll.z.cand[,l+1] <- dbinom(z[,l+1], 1, Ez.cand[,l], log=TRUE)
                prior.z <- prior.z + sum(ll.z[,l+1])
                prior.z.cand <- prior.z.cand + sum(ll.z.cand[,l+1])
              }
            }
            if(reject==FALSE){
              if(runif(1) < exp((sum(ll.y.cand[i,1:J[l],l]) + prior.z.cand) - (sum(ll.y[i,1:J[l],l]) +prior.z) )) {
                if(primary[l]==1){
                  ll.y[i,1:J[l],l] <- ll.y.cand[i,1:J[l],l]
                }
                ll.z[i,l] <- ll.z.cand[i,l]
                if((fix1|fix2)&(t>3)&(l<t)){
                  z=z.cand
                  a=a.cand
                  ll.z[,l:t]=ll.z.cand[,l:t]
                  Ez[,l:(t-1)]=Ez.cand[,l:(t-1)]
                  gamma.primeM= gamma.primeM.cand
                  gamma.primeF= gamma.primeF.cand
                }else{
                  z[,l] <- zt.cand
                  a[,l] <- at.cand
                  if(l < t) {
                    gamma.primeM <- gamma.primeM.cand
                    gamma.primeF <- gamma.primeF.cand
                    ll.z[,l+1] <- ll.z.cand[,l+1]
                    Ez[,l] <- Ez.cand[,l]
                  }
                }
                N[l] <- sum(z[,l])
              }
            }
          }
        }
      }else{
        # jointZ update
        ll_z_possibleM[,1]=ll_z_possibleF[,1]=dbinom(zpossible[,1], 1, psi,log=TRUE)
        #Get likelihood for all possible z histories. Need to update when accepted
        for(l in 2:t){
          EzpossibleM[,l-1]=zpossible[,l-1]*phiMuse[l-1] + apossible[,l-1]*gamma.primeM[l-1]
          EzpossibleF[,l-1]=zpossible[,l-1]*phiFuse[l-1] + apossible[,l-1]*gamma.primeF[l-1]
          ll_z_possibleM[,l]=dbinom(zpossible[,l], 1, EzpossibleM[,l-1],log=TRUE)
          ll_z_possibleF[,l]=dbinom(zpossible[,l], 1, EzpossibleF[,l-1],log=TRUE)
        }
        for(i in upz3){
          #new z stuff
          if(sex[i]==1){
            propto1=rowSums(exp(ll_z_possibleM))*(1*(cancel[i,]==1))
          }else{
            propto1=rowSums(exp(ll_z_possibleF))*(1*(cancel[i,]==1))
          }
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
          reject=FALSE
          for(l in 2:t){
            gamma.primeM.cand[l-1]=(Ntmp[l-1]*gammaMuse[l-1]) / sum(a.cand[sex==1,l-1])
            gamma.primeF.cand[l-1]=(Ntmp[l-1]*gammaFuse[l-1]) / sum(a.cand[sex==2,l-1])
            if(gamma.primeM.cand[l-1] > 1 | gamma.primeF.cand[l-1] > 1) { # E(Recruits) must be < nAvailable
              reject=TRUE
            }else{
              Ez.cand[sex==1,l-1]=z.cand[sex==1,l-1]*phiMuse[l-1] + a.cand[sex==1,l-1]*gamma.primeM.cand[l-1]
              Ez.cand[sex==2,l-1]=z.cand[sex==2,l-1]*phiFuse[l-1] + a.cand[sex==2,l-1]*gamma.primeF.cand[l-1]
              ll.z.cand[,l]=dbinom(z.cand[,l], 1, Ez.cand[,l-1], log=TRUE)
            }
          }
          if(reject){
            warning("Rejected z due to low M")
          }else{
            #update ll.y
            if(obstype=="bernoulli"){
              for(l in 1:t){
                if(primary[l]==1){
                  ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,],pd[i,1:J[l],l]*z.cand[i,l],log=TRUE)
                }
              }
            }else{
              for(l in 1:t){
                if(primary[l]==1){
                  ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd[i,1:J[l],l]*z.cand[i,l],log=TRUE)
                }
              }
            }
            #MH step
            if(runif(1)<exp((sum(ll.y.cand[i,,])+ll.z.cand[i,1]+sum(ll.z.cand[,2:t]))-(sum(ll.y[i,,])+ll.z[i,1]+sum(ll.z[,2:t])))*(back.prob/prop.prob)){
              z[i,]=zprop
              a[i,]=aprop
              gamma.primeM=gamma.primeM.cand
              gamma.primeF=gamma.primeF.cand
              N=Ntmp
              Ez=Ez.cand
              ll.z[i,1]=ll.z.cand[i,1]
              ll.z[,2:t]=ll.z.cand[,2:t]
              ll.y[i,,]=ll.y.cand[i,,]
              #Update likelihood for all possible z histories
              for(l in 2:t){
                EzpossibleM[,l-1]=zpossible[,l-1]*phiMuse[l-1] + apossible[,l-1]*gamma.primeM[l-1]
                EzpossibleF[,l-1]=zpossible[,l-1]*phiFuse[l-1] + apossible[,l-1]*gamma.primeF[l-1]
                ll_z_possibleM[,l]=dbinom(zpossible[,l], 1, EzpossibleM[,l-1],log=TRUE)
                ll_z_possibleF[,l]=dbinom(zpossible[,l], 1, EzpossibleF[,l-1],log=TRUE)
              }
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

      #update sexes
      ll.sex.cand=ll.sex
      upsex=sample(which(known.sex==0),proppars$sex)
      for(i in upsex){
        sex.cand=sex
        # if(known.sex[i]==1)next
        if(sex[i]==1){
          sex.cand[i]=2
        }else{
          sex.cand[i]=1
        }
        for(l in 1:t){
          if(primary[l]==1){
            if(length(lam0)==2&length(sigma)==2){
              lamd.cand[i,1:J[l],l]=lam0[sex.cand[i]]*exp(-D[i,1:J[l],l]^2/(2*sigma[sex.cand[i]]*sigma[sex.cand[i]]))
            }else if(length(lam0)==2&length(sigma)==1){
              lamd.cand[i,1:J[l],l]=lam0[sex.cand[i]]*exp(-D[i,1:J[l],l]^2/(2*sigma*sigma))
            }else if(length(lam0)==1&length(sigma)==2){
              lamd.cand[i,1:J[l],l]=lam0*exp(-D[i,1:J[l],l]^2/(2*sigma[sex.cand[i]]*sigma[sex.cand[i]]))
            }else{
              lamd.cand[i,1:J[l],l]=lam0*exp(-D[i,1:J[l],l]^2/(2*sigma*sigma))
            }
            if(obstype=="bernoulli"){
              pd.cand[i,1:J[l],l]=1-exp(-lamd.cand[i,1:J[l],l])
              ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,], pd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
            }else{
              ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
            }
          }else{#might cause a problem later so fix now
            lamd.cand[i,,l]=0
          }
        }
        #ll.z. changing 1 sex changed gamma.primeM and gamma.primeF all the way through to ll.z for *everyone*
        reject=FALSE
        for(l in 2:t){
          gamma.primeM.cand[l-1]=(N[l-1]*gammaMuse[l-1]) / sum(a[sex.cand==1,l-1])
          gamma.primeF.cand[l-1]=(N[l-1]*gammaFuse[l-1]) / sum(a[sex.cand==2,l-1])
          if(gamma.primeM.cand[l-1] > 1 | gamma.primeF.cand[l-1] > 1) { # E(Recruits) must be < nAvailable
            reject=TRUE
            next
          }
          Ez.cand[sex.cand==1,l-1]=z[sex.cand==1,l-1]*phiMuse[l-1] + a[sex.cand==1,l-1]*gamma.primeM.cand[l-1]
          Ez.cand[sex.cand==2,l-1]=z[sex.cand==2,l-1]*phiFuse[l-1] + a[sex.cand==2,l-1]*gamma.primeF.cand[l-1]
          ll.z.cand[,l]=dbinom(z[,l], 1, Ez.cand[,l-1], log=TRUE)
        }
        if(reject){
          warning("Rejected z due to low M")
          next
        }
        ll.sex.cand[i]=dbinom(sex.cand[i]-1,1,psex,log=TRUE)
        if(ACtype%in%c("metamu","metamu2")){
          if(ACtype=="metamu2"|(ACtype=="metamu"&useverts)){
            ll.s2.cand[i,]=(dnorm(s2[i,,1],s1[i,1],sigma_t[sex.cand[i]],log=TRUE)+dnorm(s2[i,,2],s1[i,2],sigma_t[sex.cand[i]],log=TRUE))
          }else{
            ll.s2.cand[i,]=log(dnorm(s2[i,,1],s1[i,1],sigma_t[sex.cand[i]])/
                        (pnorm(xlim[2],s1[i,1],sigma_t[sex.cand[i]])-pnorm(xlim[1],s1[i,1],sigma_t[sex.cand[i]])))
            ll.s2.cand[i,]=ll.s2.cand[i,]+log(dnorm(s2[i,,2],s1[i,2],sigma_t[sex.cand[i]])/
                              (pnorm(ylim[2],s1[i,2],sigma_t[sex.cand[i]])-pnorm(ylim[1],s1[i,2],sigma_t[sex.cand[i]])))
          }
        }else if(ACtype%in%c("markov")){
          for(l in 2:t){
            if(!useverts){#truncate
              ll.s2.cand[i,l-1]=log(dnorm(s2[i,l,1],s2[i,l-1,1],sigma_t[sex.cand[i]])/
                                (pnorm(xlim[2],s2[i,l-1,1],sigma_t[sex.cand[i]])-pnorm(xlim[1],s2[i,l-1,1],sigma_t[sex.cand[i]])))
              ll.s2.cand[i,l-1]=ll.s2.cand[i,l-1]+log(dnorm(s2[i,l,2],s2[i,l-1,2],sigma_t[sex.cand[i]])/
                                            (pnorm(ylim[2],s2[i,l-1,2],sigma_t[sex.cand[i]])-pnorm(ylim[1],s2[i,l-1,2],sigma_t[sex.cand[i]])))
            }else{
              ll.s2.cand[i,l-1]=dnorm(s2[i,l,1],s2[i,l-1,1],sigma_t[sex.cand[i]],log=TRUE)+dnorm(s2[i,l,2],s2[i,l-1,2],sigma_t[sex.cand[i]],log=TRUE)
            }
          }
        }else if(ACtype=="markov2"){
          for(l in 2:t){
            dists=distances[s2.cell[i,l-1],]
            probs=dexp(dists,1/sigma_t[sex.cand[i]])#cells you can pick
            probs=probs/sum(probs)
            pick=rep(0,NdSS)
            pick[s2.cell[i,l]]=1#cells you did pick
            ll.s2.cand[i,l-1]=dmultinom(pick,1,probs,log=TRUE)
          }
        }
        if(ACtype%in%c("metamu","metamu2","markov","markov2")){
          if(runif(1) < exp((ll.sex.cand[i]+sum(ll.y.cand[i,,])+sum(ll.z.cand[,2:t])+sum(ll.s2.cand[i,])-(ll.sex[i]+sum(ll.y[i,,])+sum(ll.z[,2:t])+sum(ll.s2[i,]))))){
            sex[i]=sex.cand[i]
            ll.sex[i]=ll.sex.cand[i]
            lamd[i,,]=lamd.cand[i,,]
            if(obstype=="bernoulli"){
              pd[i,,]=pd.cand[i,,]
            }
            ll.y[i,,]=ll.y.cand[i,,]
            ll.z[,2:t]=ll.z.cand[,2:t]
            ll.s2[i,]=ll.s2.cand[i,]
            gamma.primeM=gamma.primeM.cand
            gamma.primeF=gamma.primeF.cand
            Ez=Ez.cand
          }
        }else{
          if(runif(1) < exp((ll.sex.cand[i]+sum(ll.y.cand[i,,])+sum(ll.z.cand[,2:t]))-(ll.sex[i]+sum(ll.y[i,,])+sum(ll.z[,2:t])))){
            sex[i]=sex.cand[i]
            ll.sex[i]=ll.sex.cand[i]
            lamd[i,,]=lamd.cand[i,,]
            if(obstype=="bernoulli"){
              pd[i,,]=pd.cand[i,,]
            }
            ll.y[i,,]=ll.y.cand[i,,]
            ll.z[,2:t]=ll.z.cand[,2:t]
            gamma.primeM=gamma.primeM.cand
            gamma.primeF=gamma.primeF.cand
            Ez=Ez.cand
          }
        }
      }

      #update
      nfemale=sum(sex==2&z[,1]==1)
      psex <- rbeta(1, 1+nfemale, 1+N[1]-nfemale)
      ll.sex=dbinom(sex-1, 1, psex, log=TRUE)
      Nm=colSums(z==1&sex==1)
      Nf=colSums(z==1&sex==2)

      #Update phi
      if(sexparms$phi=="sex"){#if sex-specific survival
        surviveM=sum(z[sex==1,-t]==1&z[sex==1,-1]==1)
        deadM=sum(z[sex==1,-t]==1&z[sex==1,-1]==0)
        phi[1]=rbeta(1, 1+surviveM, 1+deadM)
        surviveF=sum(z[sex==2,-t]==1&z[sex==2,-1]==1)
        deadF=sum(z[sex==2,-t]==1&z[sex==2,-1]==0)
        phi[2]=rbeta(1, 1+surviveF, 1+deadF)
      }else{
        survive=sum(z[,-t]==1&z[,-1]==1)
        dead=sum(z[,-t]==1&z[,-1]==0)
        phi=rbeta(1, 1+survive, 1+dead)
      }
      if(sexparms$phi=="sex"){
        phiMuse=rep(phi[1],t-1)
        phiFuse=rep(phi[2],t-1)
      }else{
        phiMuse=rep(phi,t-1)
        phiFuse=rep(phi,t-1)
      }

      # Update gamma
      # NOTE: Must update ll.z, Ez, etc...
      if(sexparms$gamma=="fixed"){
        gamma.cand <- rnorm(1, gamma, proppars$gamma)
        gamma.cand.ok <- TRUE
        for(l in 2:t) {
          gamma.primeM.cand[l-1] <- (N[l-1]*gamma.cand) / sum(a[sex==1,l-1])
          gamma.primeF.cand[l-1] <- (N[l-1]*gamma.cand) / sum(a[sex==2,l-1])
          if(gamma.primeM.cand[l-1] > 1| gamma.primeF.cand[l-1] > 1){ ## Note don't break loop b/c ll.z needs updating because phi changed
            gamma.cand.ok=FALSE
          }
          Ez[sex==1,l-1] <- z[sex==1,l-1]*phiMuse[l-1] + a[sex==1,l-1]*gamma.primeM[l-1]
          Ez[sex==2,l-1] <- z[sex==2,l-1]*phiFuse[l-1] + a[sex==2,l-1]*gamma.primeF[l-1]
          ll.z[,l] <- dbinom(z[,l], 1, Ez[,l-1], log=TRUE)
        }
        if(gamma.cand>0 & gamma.cand.ok) {
          #Only update ll.z for a=1 cases. originally. Changed from Chandler. changed back.
          for(l in 2:t) {
            Ez.cand[sex==1,l-1] <- z[sex==1,l-1]*phiMuse[l-1] + a[sex==1,l-1]*gamma.primeM.cand[l-1]
            Ez.cand[sex==2,l-1] <- z[sex==2,l-1]*phiFuse[l-1] + a[sex==2,l-1]*gamma.primeF.cand[l-1]
            ll.z.cand[,l] <- dbinom(z[,l], 1, Ez.cand[,l-1], log=TRUE)
          }
          if(runif(1) < exp(sum(ll.z.cand[,-1]) - sum(ll.z[,-1]))) {
            gamma <- gamma.cand
            gamma.primeM <- gamma.primeM.cand
            gamma.primeF <- gamma.primeF.cand
            Ez <- Ez.cand
            ll.z[,-1] <- ll.z.cand[,-1]
          }
        }
      }else{ #sex specific
        #males first
        gamma.cand <- rnorm(1, gamma[1], proppars$gamma[1])
        gamma.cand.ok <- TRUE
        for(l in 2:t) {
          gamma.primeM.cand[l-1] <- (N[l-1]*gamma.cand) / sum(a[sex==1,l-1])
          if(gamma.primeM.cand[l-1] > 1){ ## Note don't break loop b/c ll.z needs updating because phi changed
            gamma.cand.ok=FALSE
          }
          #still need to update these for both sexes
          Ez[sex==1,l-1] <- z[sex==1,l-1]*phiMuse[l-1] + a[sex==1,l-1]*gamma.primeM[l-1]
          Ez[sex==2,l-1] <- z[sex==2,l-1]*phiFuse[l-1] + a[sex==2,l-1]*gamma.primeF[l-1]
          ll.z[,l] <- dbinom(z[,l], 1, Ez[,l-1], log=TRUE)
        }
        if(gamma.cand>0 & gamma.cand.ok) {
          #Only update ll.z for a=1 cases. originally. Changed from Chandler. changed back.
          for(l in 2:t) {
            Ez.cand[sex==1,l-1] <- z[sex==1,l-1]*phiMuse[l-1] + a[sex==1,l-1]*gamma.primeM.cand[l-1]
            Ez.cand[sex==2,l-1] <- z[sex==2,l-1]*phiFuse[l-1] + a[sex==2,l-1]*gamma.primeF[l-1]#curr
            ll.z.cand[,l] <- dbinom(z[,l], 1, Ez.cand[,l-1], log=TRUE)
          }
          if(runif(1) < exp(sum(ll.z.cand[,-1]) - sum(ll.z[,-1]))) {
            gamma[1] <- gamma.cand
            gamma.primeM <- gamma.primeM.cand
            Ez <- Ez.cand
            ll.z[,-1] <- ll.z.cand[,-1]
          }
        }
        #now females
        gamma.cand <- rnorm(1, gamma[2], proppars$gamma[2])
        gamma.cand.ok <- TRUE
        for(l in 2:t) {
          gamma.primeF.cand[l-1] <- (N[l-1]*gamma.cand) / sum(a[sex==2,l-1])
          if(gamma.primeF.cand[l-1] > 1){ ## Note don't break loop b/c ll.z needs updating because phi changed
            gamma.cand.ok=FALSE
          }
        }
        if(gamma.cand>0 & gamma.cand.ok) {
          #Only update ll.z for a=1 cases. originally. Changed from Chandler. changed back.
          for(l in 2:t) {
            Ez.cand[sex==1,l-1] <- z[sex==1,l-1]*phiMuse[l-1] + a[sex==1,l-1]*gamma.primeM[l-1]#curr
            Ez.cand[sex==2,l-1] <- z[sex==2,l-1]*phiFuse[l-1] + a[sex==2,l-1]*gamma.primeF.cand[l-1]
            ll.z.cand[,l] <- dbinom(z[,l], 1, Ez.cand[,l-1], log=TRUE)
          }
          if(runif(1) < exp(sum(ll.z.cand[,-1]) - sum(ll.z[,-1]))) {
            gamma[2] <- gamma.cand
            gamma.primeF <- gamma.primeF.cand
            Ez <- Ez.cand
            ll.z[,-1] <- ll.z.cand[,-1]
          }
        }

      }
      #Update gamma use
      if(sexparms$gamma=="sex"){
        gammaMuse=rep(gamma[1],t-1)
        gammaFuse=rep(gamma[2],t-1)
      }else{
        gammaMuse=gammaFuse=rep(gamma,t-1)
      }

      ## Now we have to update the activity centers
      #Do s2s first
      if(ACtype=="fixed"){#Stationary ACs
        for (i in 1:M) {
          if(usedSS){
            #dist based proposal
            dists=distances[s1.cell[i],]
            prop.probs=exp(-dists*dists/(2*sigma[sex[i]]*sigma[sex[i]]))
            prop.probs=prop.probs/sum(prop.probs)
            s1.cell.cand=sample(1:length(dists),1,prob=prop.probs)
            Scand=dSS[s1.cell.cand,1:2]
            dists2=distances[s1.cell.cand,]
            back.probs=exp(-dists2*dists2/(2*sigma[sex[i]]*sigma[sex[i]]))
            back.probs=back.probs/sum(back.probs)
            inbox=TRUE
            MHratio=back.probs[s1.cell[i]]/prop.probs[s1.cell.cand]
          }else{
            Scand <- c(rnorm(1, s1[i, 1], proppars$s1x), rnorm(1, s1[i, 2], proppars$s1y))
            if(useverts==FALSE){
              inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
            }else{
              inbox=any(unlist(lapply(vertices,function(x){inout(Scand,x)})))
            }
            MHratio=1
          }
          if(inbox){
            dtmp=matrix(Inf,maxJ,t)
            for(l in 1:t){
              if(primary[l]==1){
                dtmp[1:J[l],l] <- sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
                if(length(lam0)==1&length(sigma)==1){
                  lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp[1:J[l],l]*dtmp[1:J[l],l]/(2*sigma*sigma))
                }else if(length(lam0)==2&length(sigma)==1){
                  lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp[1:J[l],l]*dtmp[1:J[l],l]/(2*sigma*sigma))
                }else if(length(lam0)==1&length(sigma)==2){
                  lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp[1:J[l],l]*dtmp[1:J[l],l]/(2*sigma[sex[i]]*sigma[sex[i]]))
                }else{
                  lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp[1:J[l],l]*dtmp[1:J[l],l]/(2*sigma[sex[i]]*sigma[sex[i]]))
                }
                if(obstype=="bernoulli"){
                  pd.cand[i,1:J[l],l]=1-exp(-lamd.cand[i,1:J[l],l])
                  ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,], pd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                }else{
                  ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                }
              }else{
                lamd.cand[i,,l]=0
              }
            }
            if(runif(1) < exp(sum(ll.y.cand[i,,]) -sum(ll.y[i,,]))*MHratio){
              s1[i, ]=Scand
              s2[i,,]=matrix(rep(Scand,t),byrow=TRUE,ncol=2)
              for(l in 1:t){
                if(primary[l]==1){
                  D[i,1:J[l],l] <- dtmp[1:J[l],l]
                  lamd[i,1:J[l], l] <- lamd.cand[i,1:J[l],l]
                  if(obstype=="bernoulli"){
                    pd[i,1:J[l],l]=pd.cand[i,1:J[l],l]
                  }
                  ll.y[i,1:J[l],l]=ll.y.cand[i,1:J[l],l]
                }
              }
              if(usedSS){
                s1.cell[i]=s1.cell.cand
              }
            }
          }
        }
      }else if(ACtype=="independent"){#independent ACS
        for(l in 1:t){
          for (i in 1:M) {
            if(usedSS){
              dists=distances[s2.cell[i,l],]
              prop.probs=exp(-dists*dists/(2*sigma[sex[i]]*sigma[sex[i]]))
              prop.probs=prop.probs/sum(prop.probs)
              s2.cell.cand=sample(1:length(dists),1,prob=prop.probs)
              Scand=dSS[s2.cell.cand,1:2]
              dists2=distances[s2.cell.cand,]
              back.probs=exp(-dists2*dists2/(2*sigma[sex[i]]*sigma[sex[i]]))
              back.probs=back.probs/sum(back.probs)
              inbox=TRUE
              MHratio=back.probs[s2.cell[i,l]]/prop.probs[s2.cell.cand]
            }else{
              Scand <- c(rnorm(1, s2[i,l,1], proppars$s2x), rnorm(1, s2[i,l,2], proppars$s2y))
              if(useverts==FALSE){
                inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
              }else{
                inbox=any(unlist(lapply(vertices,function(x){inout(Scand,x)})))
              }
              MHratio=1
            }
            if(inbox) {
              if(primary[l]==1){
                dtmp<- sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
                if(length(lam0)==1&length(sigma)==1){
                  lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
                }else if(length(lam0)==2&length(sigma)==1){
                  lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp*dtmp/(2*sigma*sigma))
                }else if(length(lam0)==1&length(sigma)==2){
                  lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma[sex[i]]*sigma[sex[i]]))
                }else{
                  lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp*dtmp/(2*sigma[sex[i]]*sigma[sex[i]]))
                }
                if(obstype=="bernoulli"){
                  pd.cand[i,1:J[l],l]=1-exp(-lamd.cand[i,1:J[l],l])
                  ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,], pd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                }else{
                  ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                }
              }
              if(runif(1) < exp(sum(ll.y.cand[i,,l]) -sum(ll.y[i,,l]))*MHratio){
                s2[i,l, ] <- Scand
                if(primary[l]==1){
                  D[i,1:J[l],l] <- dtmp
                  lamd[i,1:J[l], l] <- lamd.cand[i,1:J[l],l]
                  if(obstype=="bernoulli"){
                    pd[i,1:J[l],l]=pd.cand[i,1:J[l],l]
                  }
                  ll.y[i,1:J[l],l]=ll.y.cand[i,1:J[l],l]
                }
                if(usedSS){
                  s2.cell[i,l]=s2.cell.cand
                }
              }
            }
          }
        }
      }else if(ACtype%in%c("metamu","metamu2")){
        #Update within year ACs
        for (i in 1:M){
          for(l in 1:t){
            if(usedSS){
              #dist based proposal
              dists=distances[s2.cell[i,l],]
              prop.probs=exp(-dists*dists/(2*sigma[sex[i]]*sigma[sex[i]]))
              prop.probs=prop.probs/sum(prop.probs)
              s2.cell.cand=sample(1:length(dists),1,prob=prop.probs)
              Scand=dSS[s2.cell.cand,1:2]
              dists2=distances[s2.cell.cand,]
              back.probs=exp(-dists2*dists2/(2*sigma[sex[i]]*sigma[sex[i]]))
              back.probs=back.probs/sum(back.probs)
              inbox=TRUE
              MHratio=back.probs[s2.cell[i,l]]/prop.probs[s2.cell.cand]
            }else{
              Scand=c(rnorm(1, s2[i,l,1], proppars$s2x), rnorm(1, s2[i,l,2], proppars$s2y))
              if(useverts==FALSE){
                inbox=Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
              }else{
                inbox=any(unlist(lapply(vertices,function(x){inout(Scand,x)})))
              }
              if(ACtype=="metamu2"){
                inbox=TRUE
              }
              MHratio=1
            }
            if(inbox){
              if(primary[l]==1){
                dtmp=sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
                if(length(lam0)==1&length(sigma)==1){
                  lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
                }else if(length(lam0)==2&length(sigma)==1){
                  lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp*dtmp/(2*sigma*sigma))
                }else if(length(lam0)==1&length(sigma)==2){
                  lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma[sex[i]]*sigma[sex[i]]))
                }else{
                  lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp*dtmp/(2*sigma[sex[i]]*sigma[sex[i]]))
                }
                if(obstype=="bernoulli"){
                  pd.cand[i,1:J[l],l]=1-exp(-lamd.cand[i,1:J[l],l])
                  ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,], pd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                }else{
                  ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                }
              }
              if(ACtype=="metamu2"|(ACtype=="metamu"&useverts)){
                ll.s2.cand[i,l]=(dnorm(Scand[1],s1[i,1],sigma_t[sex[i]],log=TRUE)+dnorm(Scand[2],s1[i,2],sigma_t[sex[i]],log=TRUE))
              }else{
                ll.s2.cand[i,l]=log(dnorm(Scand[1],s1[i,1],sigma_t[sex[i]])/
                                     (pnorm(xlim[2],s1[i,1],sigma_t[sex[i]])-pnorm(xlim[1],s1[i,1],sigma_t[sex[i]])))
                ll.s2.cand[i,l]=ll.s2.cand[i,l]+log(dnorm(Scand[2],s1[i,2],sigma_t[sex[i]])/
                                                    (pnorm(ylim[2],s1[i,2],sigma_t[sex[i]])-pnorm(ylim[1],s1[i,2],sigma_t[sex[i]])))
              }
              if(runif(1) < exp((sum(ll.y.cand[i,,l])+ll.s2.cand[i,l]) -(sum(ll.y[i,,l])+ll.s2[i,l]))*MHratio){
                s2[i,l,] <- Scand
                if(primary[l]==1){
                  D[i,1:J[l],l] <- dtmp
                  lamd[i,1:J[l], l] <- lamd.cand[i,1:J[l],l]
                  if(obstype=="bernoulli"){
                    pd[i,1:J[l],l]=pd.cand[i,1:J[l],l]
                  }
                  ll.y[i,1:J[l],l]=ll.y.cand[i,1:J[l],l]
                }
                ll.s2[i,l]=ll.s2.cand[i,l]
                if(usedSS){
                  s2.cell[i,l]=s2.cell.cand
                }
              }
            }
          }
        }
      }else if(ACtype%in%c("markov","markov2")){#Markov ACs
        for (i in 1:M){
          for(l in 1:t){
            if(usedSS){
              #dist based proposal
              dists=distances[s2.cell[i,l],]
              prop.probs=exp(-dists*dists/(2*sigma[sex[i]]*sigma[sex[i]]))
              prop.probs=prop.probs/sum(prop.probs)
              s2.cell.cand=sample(1:length(dists),1,prob=prop.probs)
              Scand=dSS[s2.cell.cand,1:2]
              dists2=distances[s2.cell.cand,]
              back.probs=exp(-dists2*dists2/(2*sigma[sex[i]]*sigma[sex[i]]))
              back.probs=back.probs/sum(back.probs)
              inbox=TRUE
              MHratio=back.probs[s2.cell[i,l]]/prop.probs[s2.cell.cand]
            }else{
              Scand=c(rnorm(1, s2[i,l,1], proppars$s2x), rnorm(1, s2[i,l,2], proppars$s2y))
              if(useverts==FALSE){
                inbox=Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
              }else{
                inbox=any(unlist(lapply(vertices,function(x){inout(Scand,x)})))
              }
              MHratio=1
            }
            if(inbox){
              if(primary[l]==1){
                dtmp=sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
                if(length(lam0)==1&length(sigma)==1){
                  lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
                }else if(length(lam0)==2&length(sigma)==1){
                  lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp*dtmp/(2*sigma*sigma))
                }else if(length(lam0)==1&length(sigma)==2){
                  lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma[sex[i]]*sigma[sex[i]]))
                }else{
                  lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp*dtmp/(2*sigma[sex[i]]*sigma[sex[i]]))
                }
                if(obstype=="bernoulli"){
                  pd.cand[i,1:J[l],l]=1-exp(-lamd.cand[i,1:J[l],l])
                  ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,], pd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                }else{
                  ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                }
              }
              ll.s2.cand[i,]=ll.s2[i,]
              if(ACtype=="markov"){
                if(!useverts){#truncate
                  if(l==1){#only ll.s2[i,1] matters
                    #time 1 to 2
                    ll.s2.cand[i,1]=log(dnorm(s2[i,2,1],Scand[1],sigma_t[sex[i]])/(pnorm(xlim[2],Scand[1],sigma_t[sex[i]])-
                                     pnorm(xlim[1],Scand[1],sigma_t[sex[i]])))
                    ll.s2.cand[i,1]=ll.s2.cand[i,1]+log(dnorm(s2[i,2,2],Scand[2],sigma_t[sex[i]])/
                                  (pnorm(ylim[2],Scand[2],sigma_t[sex[i]])-pnorm(ylim[1],Scand[2],sigma_t[sex[i]])))
                  }else if(l>1&l<t){#ll.s2[i,l-1] and ll.s2[i,l] matter
                    #time l-1 to time l
                    ll.s2.cand[i,l-1]=log(dnorm(Scand[1],s2[i,l-1,1],sigma_t[sex[i]])/(pnorm(xlim[2],s2[i,l-1,1],sigma_t[sex[i]])-
                                       pnorm(xlim[1],s2[i,l-1,1],sigma_t[sex[i]])))
                    ll.s2.cand[i,l-1]=ll.s2.cand[i,l-1]+log(dnorm(Scand[2],s2[i,l-1,2],sigma_t[sex[i]])/
                                  (pnorm(ylim[2],s2[i,l-1,2],sigma_t[sex[i]])-pnorm(ylim[1],s2[i,l-1,2],sigma_t[sex[i]])))
                    #time l to l+1
                    ll.s2.cand[i,l]=log(dnorm(s2[i,l+1,1],Scand[1],sigma_t[sex[i]])/(pnorm(xlim[2],Scand[1],sigma_t[sex[i]])-
                                   pnorm(xlim[1],Scand[1],sigma_t[sex[i]])))
                    ll.s2.cand[i,l]=ll.s2.cand[i,l]+log(dnorm(s2[i,l+1,2],Scand[2],sigma_t[sex[i]])/
                                (pnorm(ylim[2],Scand[2],sigma_t[sex[i]])-pnorm(ylim[1],Scand[2],sigma_t[sex[i]])))
                  }else{#only ll.s2[i,t-1] matters
                    #time t-1 to t
                    ll.s2.cand[i,t-1]=log(dnorm(Scand[1],s2[i,t-1,1],sigma_t[sex[i]])/(pnorm(xlim[2],s2[i,t-1,1],sigma_t[sex[i]])-
                                      pnorm(xlim[1],s2[i,t-1,1],sigma_t[sex[i]])))
                    ll.s2.cand[i,t-1]=ll.s2.cand[i,t-1]+log(dnorm(Scand[2],s2[i,t-1,2],sigma_t[sex[i]])/
                                      (pnorm(ylim[2],s2[i,t-1,2],sigma_t[sex[i]])-pnorm(ylim[1],s2[i,t-1,2],sigma_t[sex[i]])))
                  }
                }else{
                  if(l==1){#only ll.s2[i,1] matters
                    #time 1 to 2
                    ll.s2.cand[i,1]=dnorm(s2[i,2,1],Scand[1],sigma_t[sex[i]],log=TRUE)+dnorm(s2[i,2,2],Scand[2],sigma_t[sex[i]],log=TRUE)
                  }else if(l>1&l<t){#ll.s2[i,l-1] and ll.s2[i,l] matter
                    #time l-1 to time l
                    ll.s2.cand[i,l-1]=dnorm(Scand[1],s2[i,l-1,1],sigma_t[sex[i]],log=TRUE)+dnorm(Scand[2],s2[i,l-1,2],sigma_t[sex[i]],log=TRUE)
                    #time l to l+1
                    ll.s2.cand[i,l]=dnorm(s2[i,l+1,1],Scand[1],sigma_t[sex[i]],log=TRUE)+dnorm(s2[i,l+1,2],Scand[2],sigma_t[sex[i]],log=TRUE)
                  }else{#only ll.s2[i,t-1] matters
                    #time t-1 to t
                    ll.s2.cand[i,t-1]=dnorm(Scand[1],s2[i,t-1,1],sigma_t[sex[i]],log=TRUE)+dnorm(Scand[2],s2[i,t-1,2],sigma_t[sex[i]],log=TRUE)
                  }
                }
              }else{#markov 2
                if(l==1){#only ll.s2[i,1] matters
                  #time 1 to 2
                  dists=distances[s2.cell.cand,]
                  probs=dexp(dists,1/sigma_t[sex[i]])#cells you can pick
                  probs=probs/sum(probs)
                  pick=rep(0,NdSS)
                  pick[s2.cell[i,2]]=1#you picked this one
                  ll.s2.cand[i,1]=dmultinom(pick,1,probs,log=TRUE)
                }else if(l>1&l<t){#ll.s2[i,l-1] and ll.s2[i,l] matter
                  #time l-1 to time l
                  dists=distances[s2.cell[i,l-1],]
                  probs=dexp(dists,1/sigma_t[sex[i]])#cells you can pick
                  probs=probs/sum(probs)
                  pick=rep(0,NdSS)
                  pick[s2.cell.cand]=1#propose to pick this one
                  ll.s2.cand[i,l-1]=dmultinom(pick,1,probs,log=TRUE)
                  #time l to l+1
                  dists=distances[s2.cell.cand,]
                  probs=dexp(dists,1/sigma_t[sex[i]])#cells you can pick
                  probs=probs/sum(probs)
                  pick=rep(0,NdSS)
                  pick[s2.cell[i,l+1]]=1#propose to pick this one
                  ll.s2.cand[i,l]=dmultinom(pick,1,probs,log=TRUE)
                }else{#only ll.s2[i,t-1] matters
                  #time t-1 to t
                  dists=distances[s2.cell[i,l-1],]
                  probs=dexp(dists,1/sigma_t[sex[i]])#cells you can pick
                  probs=probs/sum(probs)
                  pick=rep(0,NdSS)
                  pick[s2.cell.cand]=1#propose to pick this one
                  ll.s2.cand[i,t-1]=dmultinom(pick,1,probs,log=TRUE)
                }
              }
              if(runif(1) < exp((sum(ll.y.cand[i,,l])+sum(ll.s2.cand[i,])) -(sum(ll.y[i,,l])+sum(ll.s2[i,])))*MHratio){
                s2[i,l,] <- Scand
                if(primary[l]==1){
                  D[i,1:J[l],l] <- dtmp
                  lamd[i,1:J[l], l] <- lamd.cand[i,1:J[l],l]
                  if(obstype=="bernoulli"){
                    pd[i,1:J[l],l]=pd.cand[i,1:J[l],l]
                  }
                  ll.y[i,1:J[l],l]=ll.y.cand[i,1:J[l],l]
                }
                ll.s2[i,]=ll.s2.cand[i,]
                if(usedSS){
                  s2.cell[i,l]=s2.cell.cand
                }
              }
            }
          }
        }
      }
      #Second AC update
      if(dualACup){
        if(iter%%proppars$dualAC==0){
          if(ACtype=="fixed"){
            for(i in 1:M){
              if(sum(known.matrix[i,])==1) next
              currpatch=dSS[s1.cell[i],3]
              idx2=which((dSS[,3]!=currpatch))
              dists=prop.probs=rep(0,NdSS)
              dists[idx2]=distances[s1.cell[i],idx2]
              prop.probs[idx2]=exp(-dists[idx2]*dists[idx2]/(2*sigma[sex[i]]*sigma[sex[i]]))
              prop.probs=prop.probs/sum(prop.probs)
              s1.cell.cand=sample(1:length(dists),1,prob=prop.probs)
              Scand=dSS[s1.cell.cand,1:2]
              backpatch=dSS[s1.cell.cand,3]
              idx2=which((dSS[,3]!=backpatch))
              dists=back.probs=rep(0,NdSS)
              dists[idx2]=distances[s1.cell.cand,idx2]
              back.probs[idx2]=exp(-dists[idx2]*dists[idx2]/(2*sigma[sex[i]]*sigma[sex[i]]))
              back.probs=back.probs/sum(back.probs)
              dtmp=matrix(Inf,maxJ,t)
              for(l in 1:t){
                if(z[i,l]==0)
                  next
                if(primary[l]==1){
                  dtmp[1:J[l],l] <- sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
                  if(length(lam0)==1&length(sigma)==1){
                    lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp[1:J[l],l]*dtmp[1:J[l],l]/(2*sigma*sigma))
                  }else if(length(lam0)==2&length(sigma)==1){
                    lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp[1:J[l],l]*dtmp[1:J[l],l]/(2*sigma*sigma))
                  }else if(length(lam0)==1&length(sigma)==2){
                    lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp[1:J[l],l]*dtmp[1:J[l],l]/(2*sigma[sex[i]]*sigma[sex[i]]))
                  }else{
                    lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp[1:J[l],l]*dtmp[1:J[l],l]/(2*sigma[sex[i]]*sigma[sex[i]]))
                  }
                  if(obstype=="bernoulli"){
                    pd.cand[i,1:J[l],l]=1-exp(-lamd.cand[i,1:J[l],l])
                    ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,], pd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                  }else{
                    ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                  }
                }
              }
              if(runif(1) < exp(sum(ll.y.cand[i,,]) -sum(ll.y[i,,]))*(back.probs[s1.cell[i]]/prop.probs[s1.cell.cand])){
                s1[i, ]=Scand
                # s2[i,,]=rep(Scand,t)
                for(l in 1:t){
                  if(primary[l]==1){
                    D[i,1:J[l],l] <- dtmp[1:J[l],l]
                    lamd[i,1:J[l], l] <- lamd.cand[i,1:J[l],l]
                    if(obstype=="bernoulli"){
                      pd[i,1:J[l],]=pd.cand[i,1:J[l],]
                    }
                    ll.y[i,1:J[l],]=ll.y.cand[i,1:J[l],]
                  }
                }
                s1.cell[i]=s1.cell.cand
              }
            }
          }else if(ACtype=="independent"){
            for(i in 1:M){
              for(l in 1:t){
                if(known.matrix[i,l]==1) next
                currpatch=dSS[s2.cell[i,l],3]
                idx2=which((dSS[,3]!=currpatch))
                dists=prop.probs=rep(0,NdSS)
                dists[idx2]=distances[s2.cell[i,l],idx2]
                prop.probs[idx2]=exp(-dists[idx2]*dists[idx2]/(2*sigma[sex[i]]*sigma[sex[i]]))
                prop.probs=prop.probs/sum(prop.probs)
                s2.cell.cand=sample(1:length(dists),1,prob=prop.probs)
                Scand=dSS[s2.cell.cand,1:2]
                backpatch=dSS[s2.cell.cand,3]
                idx2=which((dSS[,3]!=backpatch))
                dists=back.probs=rep(0,NdSS)
                dists[idx2]=distances[s2.cell.cand,idx2]
                back.probs[idx2]=exp(-dists[idx2]*dists[idx2]/(2*sigma[sex[i]]*sigma[sex[i]]))
                back.probs=back.probs/sum(back.probs)
                if(primary[l]==1){
                  dtmp<- sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
                  if(length(lam0)==1&length(sigma)==1){
                    lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
                  }else if(length(lam0)==2&length(sigma)==1){
                    lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp*dtmp/(2*sigma*sigma))
                  }else if(length(lam0)==1&length(sigma)==2){
                    lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma[sex[i]]*sigma[sex[i]]))
                  }else{
                    lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp*dtmp/(2*sigma[sex[i]]*sigma[sex[i]]))
                  }
                  if(obstype=="bernoulli"){
                    pd.cand[i,,l]=1-exp(-lamd.cand[i,,l])
                    ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,], pd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                  }else{
                    ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                  }
                }
                if(runif(1) < exp(sum(ll.y.cand[i,,l]) -sum(ll.y[i,,l]))*(back.probs[s2.cell[i,l]]/prop.probs[s2.cell.cand])){
                  s2[i,l, ] <- Scand
                  if(primary[l]==1){
                    D[i,1:J[l],l] <- dtmp
                    lamd[i,, l] <- lamd.cand[i,,l]
                    if(obstype=="bernoulli"){
                      pd[i,,l]=pd.cand[i,,l]
                    }
                    ll.y[i,1:J[l],l]=ll.y.cand[i,1:J[l],l]
                  }
                  s2.cell[i,l]=s2.cell.cand
                }
              }
            }
          }else if(ACtype%in%c("metamu","metamu2")){
            for(i in 1:M){
              for(l in 1:t){
                if(known.matrix[i,l]==1) next
                currpatch=dSS[s2.cell[i,l],3]
                idx2=which((dSS[,3]!=currpatch))
                dists=prop.probs=rep(0,NdSS)
                dists[idx2]=distances[s2.cell[i,l],idx2]
                prop.probs[idx2]=exp(-dists[idx2]*dists[idx2]/(2*sigma[sex[i]]*sigma[sex[i]]))
                prop.probs=prop.probs/sum(prop.probs)
                s2.cell.cand=sample(1:length(dists),1,prob=prop.probs)
                Scand=dSS[s2.cell.cand,1:2]
                backpatch=dSS[s2.cell.cand,3]
                idx2=which((dSS[,3]!=backpatch))
                dists=back.probs=rep(0,NdSS)
                dists[idx2]=distances[s2.cell.cand,idx2]
                back.probs[idx2]=exp(-dists[idx2]*dists[idx2]/(2*sigma[sex[i]]*sigma[sex[i]]))
                back.probs=back.probs/sum(back.probs)
                if(primary[l]==1){
                  dtmp=sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
                  if(length(lam0)==1&length(sigma)==1){
                    lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
                  }else if(length(lam0)==2&length(sigma)==1){
                    lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp*dtmp/(2*sigma*sigma))
                  }else if(length(lam0)==1&length(sigma)==2){
                    lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma[sex[i]]*sigma[sex[i]]))
                  }else{
                    lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp*dtmp/(2*sigma[sex[i]]*sigma[sex[i]]))
                  }
                  if(obstype=="bernoulli"){
                    pd.cand[i,1:J[l],l]=1-exp(-lamd.cand[i,1:J[l],l])
                    ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,], pd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                  }else{
                    ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                  }
                }
                if(ACtype=="metamu2"|(ACtype=="metamu"&useverts)){
                  ll.s2.cand[i,l]=(dnorm(Scand[1],s1[i,1],sigma_t[sex[i]],log=TRUE)+dnorm(Scand[2],s1[i,2],sigma_t[sex[i]],log=TRUE))
                }else{
                  ll.s2.cand[i,l]=log(dnorm(Scand[1],s1[i,1],sigma_t[sex[i]])/
                                        (pnorm(xlim[2],s1[i,1],sigma_t[sex[i]])-pnorm(xlim[1],s1[i,1],sigma_t[sex[i]])))
                  ll.s2.cand[i,l]=ll.s2.cand[i,l]+log(dnorm(Scand[2],s1[i,2],sigma_t[sex[i]])/
                                                       (pnorm(ylim[2],s1[i,2],sigma_t[sex[i]])-pnorm(ylim[1],s1[i,2],sigma_t[sex[i]])))
                }
                if(runif(1) < exp((sum(ll.y.cand[i,,l])+ll.s2.cand[i,l]) -(sum(ll.y[i,,l])+ll.s2[i,l]))*(back.probs[s2.cell[i,l]]/prop.probs[s2.cell.cand])){
                  s2[i,l,] <- Scand
                  if(primary[l]==1){
                    D[i,1:J[l],l] <- dtmp
                    lamd[i,1:J[l], l] <- lamd.cand[i,1:J[l],l]
                    if(obstype=="bernoulli"){
                      pd[i,1:J[l],]=pd.cand[i,1:J[l],]
                    }
                    ll.y[i,1:J[l],]=ll.y.cand[i,1:J[l],]
                  }
                  ll.s2[i,l]=ll.s2.cand[i,l]
                  s2.cell[i,l]=s2.cell.cand
                }
              }
            }
          }else if(ACtype%in%c("markov","markov2")){
            for(i in 1:M){
              for(l in 1:t){
                if(known.matrix[i,l]==1) next
                currpatch=dSS[s2.cell[i,l],3]
                idx2=which((dSS[,3]!=currpatch))
                dists=prop.probs=rep(0,NdSS)
                dists[idx2]=distances[s2.cell[i,l],idx2]
                prop.probs[idx2]=exp(-dists[idx2]*dists[idx2]/(2*sigma[sex[i]]*sigma[sex[i]]))
                prop.probs=prop.probs/sum(prop.probs)
                s2.cell.cand=sample(1:length(dists),1,prob=prop.probs)
                dists=sqrt((s2[i,l,1]-dSS[,1])^2+(s2[i,l,2]-dSS[,2])^2)
                # s2.cell.cand=sample(1:length(dists),1)
                Scand=dSS[s2.cell.cand,1:2]
                backpatch=dSS[s2.cell.cand,3]
                idx2=which((dSS[,3]!=backpatch))
                dists=back.probs=rep(0,NdSS)
                dists[idx2]=distances[s2.cell.cand,idx2]
                back.probs[idx2]=exp(-dists[idx2]*dists[idx2]/(2*sigma[sex[i]]*sigma[sex[i]]))
                back.probs=back.probs/sum(back.probs)
                if(primary[l]==1){
                  dtmp=sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
                  if(length(lam0)==1&length(sigma)==1){
                    lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
                  }else if(length(lam0)==2&length(sigma)==1){
                    lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp*dtmp/(2*sigma*sigma))
                  }else if(length(lam0)==1&length(sigma)==2){
                    lamd.cand[i,1:J[l],l]<- lam0*exp(-dtmp*dtmp/(2*sigma[sex[i]]*sigma[sex[i]]))
                  }else{
                    lamd.cand[i,1:J[l],l]<- lam0[sex[i]]*exp(-dtmp*dtmp/(2*sigma[sex[i]]*sigma[sex[i]]))
                  }
                  if(obstype=="bernoulli"){
                    pd.cand[i,1:J[l],l]=1-exp(-lamd.cand[i,1:J[l],l])
                    ll.y.cand[i,1:J[l],l] <- dbinom(y[i,1:J[l],l], tf[[l]][i,], pd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                  }else{
                    ll.y.cand[i,1:J[l],l] <- dpois(y[i,1:J[l],l], tf[[l]][i,]*lamd.cand[i,1:J[l],l]*z[i,l], log=TRUE)
                  }
                }
                if(!is.finite(sum(ll.y.cand[i,,]))) next #jump too far for obs mod ll
                ll.s2.cand[i,]=ll.s2[i,]
                if(ACtype=="markov"){
                  if(!useverts){#truncate
                    if(l==1){#only ll.s2[i,1] matters
                      #time 1 to 2
                      ll.s2.cand[i,1]=log(dnorm(s2[i,2,1],Scand[1],sigma_t[sex[i]])/(pnorm(xlim[2],Scand[1],sigma_t[sex[i]])-
                                                                                       pnorm(xlim[1],Scand[1],sigma_t[sex[i]])))
                      ll.s2.cand[i,1]=ll.s2.cand[i,1]+log(dnorm(s2[i,2,2],Scand[2],sigma_t[sex[i]])/
                                                            (pnorm(ylim[2],Scand[2],sigma_t[sex[i]])-pnorm(ylim[1],Scand[2],sigma_t[sex[i]])))
                    }else if(l>1&l<t){#ll.s2[i,l-1] and ll.s2[i,l] matter
                      #time l-1 to time l
                      ll.s2.cand[i,l-1]=log(dnorm(Scand[1],s2[i,l-1,1],sigma_t[sex[i]])/(pnorm(xlim[2],s2[i,l-1,1],sigma_t[sex[i]])-
                                                                                           pnorm(xlim[1],s2[i,l-1,1],sigma_t[sex[i]])))
                      ll.s2.cand[i,l-1]=ll.s2.cand[i,l-1]+log(dnorm(Scand[2],s2[i,l-1,2],sigma_t[sex[i]])/
                                                                (pnorm(ylim[2],s2[i,l-1,2],sigma_t[sex[i]])-pnorm(ylim[1],s2[i,l-1,2],sigma_t[sex[i]])))
                      #time l to l+1
                      ll.s2.cand[i,l]=log(dnorm(s2[i,l+1,1],Scand[1],sigma_t[sex[i]])/(pnorm(xlim[2],Scand[1],sigma_t[sex[i]])-
                                                                                         pnorm(xlim[1],Scand[1],sigma_t[sex[i]])))
                      ll.s2.cand[i,l]=ll.s2.cand[i,l]+log(dnorm(s2[i,l+1,2],Scand[2],sigma_t[sex[i]])/
                                                            (pnorm(ylim[2],Scand[2],sigma_t[sex[i]])-pnorm(ylim[1],Scand[2],sigma_t[sex[i]])))
                    }else{#only ll.s2[i,t-1] matters
                      #time t-1 to t
                      ll.s2.cand[i,t-1]=log(dnorm(Scand[1],s2[i,t-1,1],sigma_t[sex[i]])/(pnorm(xlim[2],s2[i,t-1,1],sigma_t[sex[i]])-
                                                                                           pnorm(xlim[1],s2[i,t-1,1],sigma_t[sex[i]])))
                      ll.s2.cand[i,t-1]=ll.s2.cand[i,t-1]+log(dnorm(Scand[2],s2[i,t-1,2],sigma_t[sex[i]])/
                                                                (pnorm(ylim[2],s2[i,t-1,2],sigma_t[sex[i]])-pnorm(ylim[1],s2[i,t-1,2],sigma_t[sex[i]])))
                    }
                  }else{
                    if(l==1){#only ll.s2[i,1] matters
                      #time 1 to 2
                      ll.s2.cand[i,1]=dnorm(s2[i,2,1],Scand[1],sigma_t[sex[i]],log=TRUE)+dnorm(s2[i,2,2],Scand[2],sigma_t[sex[i]],log=TRUE)
                    }else if(l>1&l<t){#ll.s2[i,l-1] and ll.s2[i,l] matter
                      #time l-1 to time l
                      ll.s2.cand[i,l-1]=dnorm(Scand[1],s2[i,l-1,1],sigma_t[sex[i]],log=TRUE)+dnorm(Scand[2],s2[i,l-1,2],sigma_t[sex[i]],log=TRUE)
                      #time l to l+1
                      ll.s2.cand[i,l]=dnorm(s2[i,l+1,1],Scand[1],sigma_t[sex[i]],log=TRUE)+dnorm(s2[i,l+1,2],Scand[2],sigma_t[sex[i]],log=TRUE)
                    }else{#only ll.s2[i,t-1] matters
                      #time t-1 to t
                      ll.s2.cand[i,t-1]=dnorm(Scand[1],s2[i,t-1,1],sigma_t[sex[i]],log=TRUE)+dnorm(Scand[2],s2[i,t-1,2],sigma_t[sex[i]],log=TRUE)
                    }
                  }
                }else{#markov 2
                  if(l==1){#only ll.s2[i,1] matters
                    #time 1 to 2
                    dists=distances[s2.cell.cand,]
                    probs=dexp(dists,1/sigma_t[sex[i]])#cells you can pick
                    probs=probs/sum(probs)
                    pick=rep(0,NdSS)
                    pick[s2.cell[i,2]]=1#you picked this one
                    ll.s2.cand[i,1]=dmultinom(pick,1,probs,log=TRUE)
                  }else if(l>1&l<t){#ll.s2[i,l-1] and ll.s2[i,l] matter
                    #time l-1 to time l
                    dists=distances[s2.cell[i,l-1],]
                    probs=dexp(dists,1/sigma_t[sex[i]])#cells you can pick
                    probs=probs/sum(probs)
                    pick=rep(0,NdSS)
                    pick[s2.cell.cand]=1#propose to pick this one
                    ll.s2.cand[i,l-1]=dmultinom(pick,1,probs,log=TRUE)
                    #time l to l+1
                    dists=distances[s2.cell.cand,]
                    probs=dexp(dists,1/sigma_t[sex[i]])#cells you can pick
                    probs=probs/sum(probs)
                    pick=rep(0,NdSS)
                    pick[s2.cell[i,l+1]]=1#propose to pick this one
                    ll.s2.cand[i,l]=dmultinom(pick,1,probs,log=TRUE)
                  }else{#only ll.s2[i,t-1] matters
                    #time t-1 to t
                    dists=distances[s2.cell[i,l-1],]
                    probs=dexp(dists,1/sigma_t[sex[i]])#cells you can pick
                    probs=probs/sum(probs)
                    pick=rep(0,NdSS)
                    pick[s2.cell.cand]=1#propose to pick this one
                    ll.s2.cand[i,t-1]=dmultinom(pick,1,probs,log=TRUE)
                  }
                }
                if(runif(1) < exp((sum(ll.y.cand[i,,l])+sum(ll.s2.cand[i,])) -(sum(ll.y[i,,l])+sum(ll.s2[i,])))*(back.probs[s2.cell[i,l]]/prop.probs[s2.cell.cand])){
                  s2[i,l,] <- Scand
                  if(primary[l]==1){
                    D[i,1:J[l],l] <- dtmp
                    lamd[i,1:J[l], l] <- lamd.cand[i,1:J[l],l]
                    if(obstype=="bernoulli"){
                      pd[i,1:J[l],l]=pd.cand[i,1:J[l],l]
                    }
                    ll.y[i,1:J[l],l]=ll.y.cand[i,1:J[l],l]
                  }
                  ll.s2[i,]=ll.s2.cand[i,]
                  s2.cell[i,l]=s2.cell.cand
                }
              }
            }
          }
        }
      }
      #update metamus and sigma_t
      if(ACtype%in%c("metamu","metamu2")){
        #Update meta mus
        for (i in 1:M){
          Scand <- c(rnorm(1,s1[i,1],proppars$s1x), rnorm(1,s1[i,2],proppars$s1y))
          # if(usedSS){
          #   dists=sqrt((Scand[1]-dSS[,1])^2+(Scand[2]-dSS[,2])^2)
          #   s1.cell.cand=which(dists==min(dists))
          #   Scand=dSS[s1.cell.cand,1:2]
          # }
          if(useverts==FALSE){
            inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
          }else{
            inbox=any(unlist(lapply(vertices,function(x){inout(Scand,x)})))
          }
          if(inbox){
            ll.s2.cand[i,]<-dnorm(s2[i,,1],Scand[1],sigma_t[sex[i]],log=TRUE)+dnorm(s2[i,,2],Scand[2],sigma_t[sex[i]],log=TRUE)
            if(ACtype=="metamu2"|(ACtype=="metamu"&useverts)){
              ll.s2.cand[i,]=dnorm(s2[i,,1],Scand[1],sigma_t[sex[i]],log=TRUE)+dnorm(s2[i,,2],Scand[2],sigma_t[sex[i]],log=TRUE)
            }else{#truncate
              ll.s2.cand[i,]=log(dnorm(s2[i,,1],Scand[1],sigma_t[sex[i]])/
                                   (pnorm(xlim[2],Scand[1],sigma_t[sex[i]])-pnorm(xlim[1],Scand[1],sigma_t[sex[i]])))
              ll.s2.cand[i,]=ll.s2.cand[i,]+log(dnorm(s2[i,,2],Scand[2],sigma_t[sex[i]])/
                                                  (pnorm(ylim[2],Scand[2],sigma_t[sex[i]])-pnorm(ylim[1],Scand[2],sigma_t[sex[i]])))
            }
            if(runif(1) < exp(sum(ll.s2.cand[i,]) - sum(ll.s2[i,]))) {
              s1[i, ]=Scand
              ll.s2[i,]=ll.s2.cand[i,]
            }
          }
        }
        #Update sigma_t
        if(sexparms$sigma_t=="fixed"){
          sigma_t.cand <- rnorm(1,sigma_t[1],proppars$sigma_t[1])
          if(sigma_t.cand > 0){
            if(ACtype=="metamu2"|(ACtype=="metamu"&useverts)){
              ll.s2.cand=dnorm(s2[,,1],s1[,1],sigma_t.cand,log=TRUE)+dnorm(s2[,,2],s1[,2],sigma_t.cand,log=TRUE)
            }else{
              ll.s2.cand=log(dnorm(s2[,,1],s1[,1],sigma_t.cand)/
                               (pnorm(xlim[2],s1[,1],sigma_t.cand)-pnorm(xlim[1],s1[,1],sigma_t.cand)))
              ll.s2.cand=ll.s2.cand+log(dnorm(s2[,,2],s1[,2],sigma_t.cand)/
                                          (pnorm(ylim[2],s1[,2],sigma_t.cand)-pnorm(ylim[1],s1[,2],sigma_t.cand)))
            }
            if(runif(1) < exp(sum(ll.s2.cand) - sum(ll.s2))) {
              sigma_t[1:2]=sigma_t.cand
              ll.s2=ll.s2.cand
            }
          }
        }else{
          for(i2 in 1:2){
            sigma_t.cand <- rnorm(1,sigma_t[i2],proppars$sigma_t[i2])
            if(sigma_t.cand > 0){
              ll.s2.cand=ll.s2
              ll.s2.cand[sex==i2,]=dnorm(s2[sex==i2,,1],s1[sex==i2,1],sigma_t.cand,log=TRUE)+dnorm(s2[sex==i2,,2],s1[sex==i2,2],sigma_t.cand,log=TRUE)
              if(ACtype=="metamu2"|(ACtype=="metamu"&useverts)){
                ll.s2.cand[sex==i2,]=dnorm(s2[sex==i2,,1],s1[sex==i2,1],sigma_t.cand,log=TRUE)+dnorm(s2[sex==i2,,2],s1[sex==i2,2],sigma_t.cand,log=TRUE)
              }else{
                ll.s2.cand[sex==i2,]=log(dnorm(s2[sex==i2,,1],s1[sex==i2,1],sigma_t.cand)/
                                 (pnorm(xlim[2],s1[sex==i2,1],sigma_t.cand)-pnorm(xlim[1],s1[sex==i2,1],sigma_t.cand)))
                ll.s2.cand[sex==i2,]=ll.s2.cand[sex==i2,]+log(dnorm(s2[sex==i2,,2],s1[sex==i2,2],sigma_t.cand)/
                                            (pnorm(ylim[2],s1[sex==i2,2],sigma_t.cand)-pnorm(ylim[1],s1[sex==i2,2],sigma_t.cand)))
              }
              if (runif(1) < exp(sum(ll.s2.cand) - sum(ll.s2))) {
                sigma_t[i2]=sigma_t.cand
                ll.s2=ll.s2.cand
              }
            }
          }
        }
      }else if(ACtype%in%c("markov")){#Markov ACs
        #Update sigma_t
        if(sexparms$sigma_t=="fixed"){
          sigma_t.cand <- rnorm(1,sigma_t[1],proppars$sigma_t[1])
          if(sigma_t.cand > 0){
            for(l in 2:t){
              for(i in 1:M){
                if(!useverts){#truncate
                  ll.s2.cand[i,l-1]=log(dnorm(s2[i,l,1],s2[i,l-1,1],sigma_t.cand)/(pnorm(xlim[2],s2[i,l-1,1],sigma_t.cand)-
                                                                                     pnorm(xlim[1],s2[i,l-1,1],sigma_t.cand)))
                  ll.s2.cand[i,l-1]=ll.s2.cand[i,l-1]+log(dnorm(s2[i,l,2],s2[i,l-1,2],sigma_t.cand)/
                                                            (pnorm(ylim[2],s2[i,l-1,2],sigma_t.cand)-pnorm(ylim[1],s2[i,l-1,2],sigma_t.cand)))
                }else{
                  ll.s2.cand[i,l-1]=dnorm(s2[i,l,1],s2[i,l-1,1],sigma_t.cand,log=TRUE)+dnorm(s2[i,l,2],s2[i,l-1,2],sigma_t.cand,log=TRUE)
                }
              }
            }
            if(runif(1) < exp(sum(ll.s2.cand) - sum(ll.s2))){
              sigma_t[1:2]=sigma_t.cand
              ll.s2=ll.s2.cand
            }
          }
        }else{
          for(i2 in 1:2){
            sigma_t.cand <- rnorm(1,sigma_t[i2],proppars$sigma_t[i2])
            if(sigma_t.cand > 0){
              ll.s2.cand=ll.s2
              for(l in 2:t){
                for(i in 1:M){
                  if(sex[i]==i2){
                    if(!useverts){#truncate
                      ll.s2.cand[i,l-1]=log(dnorm(s2[i,l,1],s2[i,l-1,1],sigma_t.cand)/
                                              (pnorm(xlim[2],s2[i,l-1,1],sigma_t.cand)- pnorm(xlim[1],s2[i,l-1,1],sigma_t.cand)))
                      ll.s2.cand[i,l-1]=ll.s2.cand[i,l-1]+log(dnorm(s2[i,l,2],s2[i,l-1,2],sigma_t.cand)/
                                                                (pnorm(ylim[2],s2[i,l-1,2],sigma_t.cand)-pnorm(ylim[1],s2[i,l-1,2],sigma_t.cand)))
                    }else{
                      ll.s2.cand[i,l-1]=dnorm(s2[i,l,1],s2[i,l-1,1],sigma_t.cand,log=TRUE)+dnorm(s2[i,l,2],s2[i,l-1,2],sigma_t.cand,log=TRUE)
                    }
                  }
                }
              }
              if (runif(1) < exp(sum(ll.s2.cand) - sum(ll.s2))) {
                sigma_t[i2]=sigma_t.cand
                ll.s2=ll.s2.cand
              }
            }
          }
        }
      }else if(ACtype%in%c("markov2")){#Markov2 ACs
        #Update sigma_t
        if(sexparms$sigma_t=="fixed"){
          sigma_t.cand <- rnorm(1,sigma_t[1],proppars$sigma_t[1])
          if(sigma_t.cand > 0){
            for(l in 2:t){
              dists=distances[s2.cell[,l-1],]#cells each ind can choose
              for(i in 1:M){
                probs=dexp(dists[i,],1/sigma_t.cand)
                probs=probs/sum(probs)
                pick=rep(0,NdSS)
                pick[s2.cell[i,l]]=1#you picked this one
                ll.s2.cand[i,l-1]=dmultinom(pick,1,probs,log=TRUE)
              }
            }
            if (runif(1) < exp(sum(ll.s2.cand) - sum(ll.s2))) {
              sigma_t[1:2]=sigma_t.cand
              ll.s2=ll.s2.cand
            }
          }
        }else{
          for(i2 in 1:2){
            sigma_t.cand= rnorm(1,sigma_t[i2],proppars$sigma_t[i2])
            if(sigma_t.cand > 0){
              ll.s2.cand=ll.s2
              for(l in 2:t){
                dists=distances[s2.cell[,l-1],]#cells each ind can choose
                for(i in 1:M){
                  if(sex[i]==i2){
                    probs=dexp(dists[i,],1/sigma_t.cand)
                    probs=probs/sum(probs)
                    pick=rep(0,NdSS)
                    pick[s2.cell[i,l]]=1#you picked this one
                    ll.s2.cand[i,l-1]=dmultinom(pick,1,probs,log=TRUE)
                  }
                }
              }
              if (runif(1) < exp(sum(ll.s2.cand) - sum(ll.s2))) {
                sigma_t[i2]=sigma_t.cand
                ll.s2=ll.s2.cand
              }
            }
          }
        }
      }

      #Do we record output on this iteration?
      if(iter>(nburn)&iter%%nthin==0){
        if(storeLatent){
          zout[idx,,]<- z
          sexout[idx,]=sex
        }
        if(ACtype%in%c("markov","markov2")){
          if(sexparms$sigma_t=="sex"){
            out[idx,]<- c(lam0,sigma ,gamma,phi,N,Nm,Nf,sigma_t,psex,psi)
          }else{
            out[idx,]<- c(lam0,sigma ,gamma,phi,N,Nm,Nf,sigma_t[1],psex,psi)
          }
          if(storeLatent){
            s2xout[idx,,]<- s2[,,1]
            s2yout[idx,,]<- s2[,,2]
          }
        }else if(ACtype%in%c("metamu","metamu2")){
          if(sexparms$sigma_t=="sex"){
            out[idx,]<- c(lam0,sigma ,gamma,phi,N,Nm,Nf,sigma_t,psex,psi)
          }else{
            out[idx,]<- c(lam0,sigma ,gamma,phi,N,Nm,Nf,sigma_t[1],psex,psi)
          }
          if(storeLatent){
            s1xout[idx,]<- s1[,1]
            s1yout[idx,]<- s1[,2]
            s2xout[idx,,]<- s2[,,1]
            s2yout[idx,,]<- s2[,,2]
          }
        }else if(ACtype=="independent"){
          out[idx,]<- c(lam0,sigma ,gamma,phi,N,Nm,Nf,psex,psi)
          if(storeLatent){
            s2xout[idx,,]<- s2[,,1]
            s2yout[idx,,]<- s2[,,2]
          }
        }else{
          out[idx,]<- c(lam0,sigma ,gamma,phi,N,Nm,Nf,psex,psi)
          if(storeLatent){
            s1xout[idx,]<- s1[,1]
            s1yout[idx,]<- s1[,2]
          }
        }
        idx=idx+1
      }
    }  # end of MCMC algorithm
    if(storeLatent==TRUE){
      if(ACtype%in%c("markov","markov2","independent")){
        list(out=out,s2xout=s2xout, s2yout=s2yout, zout=zout,sexout=sexout)
      }else if(ACtype%in%c("metamu","metamu2")){
        list(out=out, s1xout=s1xout, s1yout=s1yout,s2xout=s2xout, s2yout=s2yout, zout=zout,sexout=sexout)
      }else{
        list(out=out, s1xout=s1xout, s1yout=s1yout, zout=zout,sexout=sexout)
      }
    }else{
      list(out=out)
    }
  }
