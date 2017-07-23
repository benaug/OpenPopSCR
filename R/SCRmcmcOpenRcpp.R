SCRmcmcOpenRcpp <-
  function(data,niter=2400,nburn=1200, nthin=5,M = 200, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),
           jointZ=TRUE,storeLatent=TRUE,ACtype="fixed",obstype="bernoulli",dSS=NA,dualACup=FALSE){
    library(abind)
    t=dim(data$y)[3]
    y<-data$y
    X<-data$X
    #make sure X list elements are matrices
    for(i in 1:length(X)){
      X[[i]]=as.matrix(X[[i]])
    }
    dSS=as.matrix(dSS)
    if(!is.na(dSS[1])&"vertices"%in%names(data)){
      rem=which(names(data)=="vertices")
      data[[rem]]=NULL
      warning("Discarding vertices since dSS supplied")
      if(max(dSS[,1])>xlim[2]|min(dSS[,1])<xlim[1]|max(dSS[,2])>ylim[2]|min(dSS[,2])<ylim[1]){
        stop("dSS dimensions exceed xlim or ylim. Change dSS or buff")
      }
    }
    if(length(X)!=t){
      stop("must input traps for each year")
    }
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
    if(!ACtype%in%c("fixed","independent","metamu","metamu2","markov","markov2")){
      stop("ACtype must be 'fixed','independent','metamu', 'metamu2', 'markov' or 'markov2'")
    }
    if(ACtype%in%c("metamu","metamu2","markov","markov2")){
      if(!"sigma_t"%in%names(proppars)){
        stop("must supply proppars$sigma_t if ACtype is metamu, metamu2, markov, or markov2")
      }
      if(is.null(sigma_t)){
        stop("must supply inits$sigma_t if ACtype is metamu, metamu2, markov, or markov2")
      }
    }
    if(ACtype%in%c("metamu","metamu2","fixed")){
      if(!"s1x"%in%names(proppars)|!"s1y"%in%names(proppars)){
        stop("must supply proppars$s1x and proppars$s1y if ACtype is metamu, metamu2, or fixed")
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
    #Initialize s1 and s2
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
              if(ACtype%in%c("metamu","metamu2")){
                inside=FALSE
                while(inside==FALSE){
                  s2[i,l,]=c(rnorm(1,s1[i,1],sigma_t),rnorm(1,s1[i,2],sigma_t))
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
                  s2[i,l,1]=rnorm(1,s2[i,l-1,1],sigma_t)
                  s2[i,l,2]=rnorm(1,s2[i,l-1,2],sigma_t)
                  inside=any(unlist(lapply(vertices,function(x){inout(s2[i,l,],x)})))
                }
              }
            }
          }
        }
        #calculate likelihoods
        if(ACtype%in%c("metamu","metamu2")){
          if(ACtype=="metamu2"|(ACtype=="metamu"&useverts)){
            ll.s2=(dnorm(s2[,,1],s1[,1],sigma_t,log=TRUE)+dnorm(s2[,,2],s1[,2],sigma_t,log=TRUE))
          }else{
            ll.s2=log(dnorm(s2[,,1],s1[,1],sigma_t)/
                        (pnorm(xlim[2],s1[,1],sigma_t)-pnorm(xlim[1],s1[,1],sigma_t)))
            ll.s2=ll.s2+log(dnorm(s2[,,2],s1[,2],sigma_t)/
                              (pnorm(ylim[2],s1[,2],sigma_t)-pnorm(ylim[1],s1[,2],sigma_t)))
          }
          ll.s2.cand=ll.s2
        }else if(ACtype%in%c("markov")){
          ll.s2=matrix(NA,nrow=M,ncol=t-1)
          for(l in 2:t){
            if(!useverts){#truncate
              ll.s2[,l-1]=log(dnorm(s2[,l,1],s2[,l-1,1],sigma_t)/
                                (pnorm(xlim[2],s2[,l-1,1],sigma_t)-pnorm(xlim[1],s2[,l-1,1],sigma_t)))
              ll.s2[,l-1]=ll.s2[,l-1]+log(dnorm(s2[,l,2],s2[,l-1,2],sigma_t)/
                                            (pnorm(ylim[2],s2[,l-1,2],sigma_t)-pnorm(ylim[1],s2[,l-1,2],sigma_t)))
            }else{
              ll.s2[,l-1]=(dnorm(s2[,l,1],s2[,l-1,1],sigma_t,log=TRUE)+dnorm(s2[,l,2],s2[,l-1,2],sigma_t,log=TRUE))
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
                  s2[i,l,]=c(rnorm(1,s1[i,1],sigma_t),rnorm(1,s1[i,2],sigma_t))
                  inside=any(unlist(lapply(vertices,function(x){inout(s2[i,l,],x)})))
                }
              }
            }
          }
          if(ACtype=="metamu2"|(ACtype=="metamu"&useverts)){
            ll.s2=(dnorm(s2[,,1],s1[,1],sigma_t,log=TRUE)+dnorm(s2[,,2],s1[,2],sigma_t,log=TRUE))
          }else{
            ll.s2=log(dnorm(s2[,,1],s1[,1],sigma_t)/
                        (pnorm(xlim[2],s1[,1],sigma_t)-pnorm(xlim[1],s1[,1],sigma_t)))
            ll.s2=ll.s2+log(dnorm(s2[,,2],s1[,2],sigma_t)/
                              (pnorm(ylim[2],s1[,2],sigma_t)-pnorm(ylim[1],s1[,2],sigma_t)))
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
                  s2[i,l,1]=rnorm(1,s2[i,l-1,1],sigma_t)
                  s2[i,l,2]=rnorm(1,s2[i,l-1,2],sigma_t)
                  inside=any(unlist(lapply(vertices,function(x){inout(s2[i,l,],x)})))
                }
              }
            }
          }
          for(l in 2:t){
            if(!useverts){#truncate
              ll.s2[,l-1]=log(dnorm(s2[,l,1],s2[,l-1,1],sigma_t)/
                                (pnorm(xlim[2],s2[,l-1,1],sigma_t)-pnorm(xlim[1],s2[,l-1,1],sigma_t)))
              ll.s2[,l-1]=ll.s2[,l-1]+log(dnorm(s2[,l,2],s2[,l-1,2],sigma_t)/
                                            (pnorm(ylim[2],s2[,l-1,2],sigma_t)-pnorm(ylim[1],s2[,l-1,2],sigma_t)))
            }else{
              ll.s2[,l-1]=(dnorm(s2[,l,1],s2[,l-1,1],sigma_t,log=TRUE)+dnorm(s2[,l,2],s2[,l-1,2],sigma_t,log=TRUE))
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
              probs=dexp(dists[i,],1/sigma_t)
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
    D=lamd=ll.y=ll.y.cand=array(0,dim=c(M,maxJ,t))
    D[is.na(D)]=Inf  #hack to allow years with different J and K to fit in one array
    for(l in 1:t){
      if(primary[l]==1){
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
    }else{
      #make up some fake stuff to feed to rcpp
      zpossible=apossible=cancel=matrix(0,nrow=2,ncol=2)
    }
    Xidx=unlist(lapply(X,dim))[seq(1,2*length(X)-1,2)]
    Xcpp=array(NA,dim=c(t,max(Xidx),2))
    for(l in 1:t){
      if(primary[l]){
        for(j in 1:Xidx[l]){
          Xcpp[l,j,]=as.numeric(X[[l]][j,])
        }
      }
    }
    each=unlist(lapply(inits,length))[1:4]
    npar=sum(each)+t+1
    if(ACtype%in%c("metamu","metamu2","markov","markov2")){
      npar=npar+1
    }
    if(!dualACup){#dummy for Rcpp
      proppars$dualAC=1
    }
    if(ACtype%in%c("fixed","independent")){#dummy for Rcpp
      sigma_t=1
      proppars$sigma_t=1
    }
    if(!ACtype%in%c("metamu","metamu2","fixed")){#dummy for Rcpp
      proppars$s1x=proppars$s1y=1
    }
    if(is.null(proppars$s2x)|is.null(proppars$s2x)){#dummy for Rcpp
      proppars$s2x=proppars$s2y=1
    }
    if(!dualACup){#dummy for Rcpp
      proppars$dualAC=1
    }
    if(is.null(proppars$propz)){#dummy for Rcpp
      proppars$propz=10
    }

    #So these aren't modified by Rcpp
    lam0in=lam0
    sigmain=sigma
    gammain=gamma
    phiin=phi
    if(ACtype=="fixed"){
      ACtype=1
    }else if(ACtype=="metamu"){
      ACtype=2
    }else if(ACtype=="markov"){
      ACtype=3
    }else if(ACtype=="independent"){
      ACtype=4#independent
    }else if(ACtype=="metamu2"){#metamu2
      ACtype=5
    }else{
      ACtype=6
    }
    if(obstype=="bernoulli"){
      obstype2=1
    }else{
      obstype2=2
    }
    tf2=matrix(NA,nrow=maxJ,ncol=t)
    for(l in 1:t){
      if(primary[l]){
        tf2[1:J[l],l]=tf[[l]][1,]
      }
    }
    if(!usedSS){
      dSS=matrix(0.5,nrow=2,ncol=2)#dummy to fool rcpp
    }
    if(is.null(sigma_t)){
      sigma_t=0.5 #dummy for Rcpp
    }
    if(is.null(proppars$propz)){
      proppars$propz=rep(10,t-1)
    }
    primaryin=primary==1
    if(usedSS==FALSE){
      s1.cell=rep(1,M)
      s2.cell=matrix(1,nrow=M,ncol=t)
      distances=matrix(c(0,0,0,0),nrow=2,ncol=2)
    }

    store=mcmc_Open(lam0in,  sigmain,  gammain, gamma.prime, phiin, D,lamd, y, z, a,s1,s2,
                    ACtype, useverts, vertices, xlim, ylim, known.matrix, Xidx, Xcpp, K, Ez,  psi,
                    N, proppars$lam0, proppars$sigma, proppars$propz,  proppars$gamma, proppars$s1x,  proppars$s1y,
                    proppars$s2x,proppars$s2y,proppars$sigma_t,sigma_t,niter,nburn,nthin,npar,each,jointZ,
                    zpossible,apossible,cancel,obstype2,tf2,s2.cell-1,s1.cell-1,dSS,usedSS,primaryin,
                    dualACup,proppars$dualAC,distances,storeLatent)

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
    # # Sanity check to make sure there are no shenanigans in z or a
    # if(t>3){
    #   for(i in 1:niter){
    #     sanity=any(rowSums(zout[i,,])==0&rowSums(aout[i,,-t])!=(t-1))
    #     for(l2 in 2:(t-2)){
    #       sanity=c(sanity,any(aout[i,,l2]==0&aout[i,,l2-1]==1&aout[i,,l2+1]==1))
    #     }
    #     sanity=c(sanity,any((zout[i,,-t]+aout[i,,-t])>1))
    #     if(any(sanity)){stop("insanity")}
    #   }
    # }else if(t==3){
    #   for(i in 1:niter){
    #     sanity=any(rowSums(zout[i,,])==0&rowSums(aout[i,,-t])!=(t-1))
    #     sanity=c(sanity,any((zout[i,,-t]+aout[i,,-t])>1))
    #     if(any(sanity)){stop("insanity")}
    #   }
    # }

    # name out
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
    if(ACtype%in%c(2,3,5,6)){
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,"sigma_t","psi")
    }else{
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,"psi")
    }
    if(storeLatent==TRUE){
      if(ACtype%in%c(3,4,6)){
        list(out=out,s2xout=s2xout, s2yout=s2yout, zout=zout)
      }else if(ACtype%in%c(2,5)){
        list(out=out, s1xout=s1xout, s1yout=s1yout,s2xout=s2xout, s2yout=s2yout, zout=zout)
      }else{
        list(out=out, s1xout=s1xout, s1yout=s1yout, zout=zout)
      }
    }else{
      list(out=out)
    }
  }
