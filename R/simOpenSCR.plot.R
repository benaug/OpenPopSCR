#' Plot the open population simulation. Will document later.
#' @description Plot the open population simulation. Will document later.
#' @param data an object created by simOpenSCR
#' @return a plot
#' @author Ben Augustine
#' @export
simOpenSCR.plot=function(data,animated=FALSE,delay=2){
  X=data$X
  s=data$sfull
  z=data$z
  y=data$yfull
  t=dim(y)[3]
  known.matrix=1*(apply(y,c(1,3),sum)>1)
  xlim=c(min(unlist(lapply(X,function(x){min(x[,1])})))-buff,max(unlist(lapply(X,function(x){max(x[,1])})))+buff)
  ylim=c(min(unlist(lapply(X,function(x){min(x[,2])})))-buff,max(unlist(lapply(X,function(x){max(x[,2])})))+buff)
  M=nrow(z)
  ACtype=data$ACtype
  recruit=die=rep(NA,M)
  for(i in 2:t){
    recruit[z[,i-1]==0&z[,i]==1]=i
    die[z[,i-1]==1&z[,i]==0]=i
  }
  if(ACtype=="fixed"){
    plot(X[[1]],xlim=xlim,ylim=ylim,pch=4,main=paste("Year",1),xlab="X",ylab="Y")
    points(s[z[,1]==1,1,],pch=16)
    for(i in 2:t){
      Sys.sleep(delay)
      plot(X[[1]],xlim=xlim,ylim=ylim,pch=4,main=paste("Year",i),xlab="X",ylab="Y")
      points(s[z[,i-1]==1,i-1,],pch=16)
      points(s[which(recruit==i),i,],pch=16,col="green")
      points(s[which(die==i),i-1,],pch=16,col="red")
      Sys.sleep(delay)
      plot(X[[1]],xlim=xlim,ylim=ylim,pch=4,main=paste("Year",i),xlab="X",ylab="Y")
      points(s[z[,i]==1,i,],pch=16)
    }
  }else if(ACtype=="metamu"){
    if(animated){
    plot(X[[1]],xlim=xlim,ylim=ylim,pch=4,main=paste("Year",1),xlab="X",ylab="Y")
    points(s[z[,1]==1,1,],pch=16)
    for(i in 2:t){
      Sys.sleep(delay)
      plot(X[[1]],xlim=xlim,ylim=ylim,pch=4,main=paste("Year",i),xlab="X",ylab="Y")
      points(s[z[,i-1]==1,i-1,],pch=16)
      points(s[which(recruit==i),i-1,],pch=16,col="green")
      points(s[which(die==i),i-1,],pch=16,col="red")
      Sys.sleep(delay)
      plot(X[[1]],xlim=xlim,ylim=ylim,pch=4,main=paste("Year",i),xlab="X",ylab="Y")
      points(s[z[,i]==1,i-1,],pch=16)
      Sys.sleep(delay)
      plot(X[[1]],xlim=xlim,ylim=ylim,pch=4,main=paste("Year",i),xlab="X",ylab="Y")
      points(s[z[,i-1]==1&z[,i]==1,i-1,],pch=16,col=rgb(128,128,128,100,max=255))
      Sys.sleep(delay)
      plot(X[[1]],xlim=xlim,ylim=ylim,pch=4,main=paste("Year",i),xlab="X",ylab="Y")
      points(s[z[,i-1]==1&z[,i]==1,i-1,],pch=16,col=rgb(128,128,128,100,max=255))
      points(s[z[,i]==1,i,],pch=16)
      for(j in 1:M){
        if(z[j,i-1]==1&z[j,i]==1){
          lines(x=c(s[j,i-1,1],s[j,i,1]),y=c(s[j,i-1,2],s[j,i,2]))
        }
      }
      Sys.sleep(delay)
      plot(X[[1]],xlim=xlim,ylim=ylim,pch=4,main=paste("Year",i),xlab="X",ylab="Y")
      points(s[z[,i]==1,i,],pch=16)
    }
    }else{
      mu=data$mufull
      plot(X[[1]],xlim=xlim,ylim=ylim,pch=4,xlab="X",ylab="Y")
      zAny=rowSums(z)>0
      points(mu[which(zAny),],pch=16)
      for(i in 1:t){
        points(s[z[,i]==1,i,],pch=16,col=rgb(128,128,128,100,max=255))
        for(j in 1:M){
          if(z[j,i]==1){
            lines(x=c(mu[j,1],s[j,i,1]),y=c(mu[j,2],s[j,i,2]))
          }
        }
      }
    }
  }else if(ACtype=="markov"){
    if(animated){
      stop("No animated plot for markov movement, yet")
    }else{
      plot(X[[1]],xlim=xlim,ylim=ylim,pch=4,xlab="X",ylab="Y")
      points(s[z[,1]==1,1,],pch=16,col="black")
      for(i in 2:t){
        points(s[z[,i]==1,i,],pch=16,col="black")
        for(j in 1:M){
          if(z[j,i-1]==1&z[j,i]==1){
            lines(x=c(s[j,i-1,1],s[j,i,1]),y=c(s[j,i-1,2],s[j,i,2]))
          }
        }
      }
    }


  }else{
    stop("No plot function for independent activity centers, yet")
  }
}
