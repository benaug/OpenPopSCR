library(OpenPopSCR)
library(coda)
#Rcpp markov sigma_t not updating
#same but t=3
t=3
N=50
p0=0.3
lam0=-log(1-p0)
sigma=1
sigma_t=0.5 #dispersal sigma
phi=0.7
gamma=0.3
buff=2
X=list(expand.grid(4:11,4:11),expand.grid(4:11,4:11),expand.grid(4:11,4:11))
K=c(10,10,10)
M=200
ACtype="metamu2"
data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,sigma_t=sigma_t,K=K,X=X,t=t,M=M,buff=buff,
                ACtype=ACtype,poly=poly)
inits=list(lam0=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M,sigma_t=sigma_t)
niter=2000
nburn=0
nthin=1
proppars=list(lam0=0.01,sigma=0.02,gamma=0.1,s1x=0.05,s1y=0.05,s2x=0.2,s2y=0.2,sigma_t=0.09)

a=Sys.time()
out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,ACtype=ACtype,Rcpp=FALSE)
b=Sys.time()
b-a
plot(mcmc(out$out))
summary(mcmc(out$out))

#multi poly
X=list(rbind(expand.grid(3:8,3:8),expand.grid(17:23,17:23)),rbind(expand.grid(3:8,3:8),expand.grid(17:23,17:23)),rbind(expand.grid(3:8,3:8),expand.grid(17:23,17:23)))
poly=vector("list",2)
poly[[1]]=rbind(c(1,1),c(1,10),c(10,10),c(10,1))
poly[[2]]=rbind(c(15,15),c(15,25),c(25,25),c(25,15))
data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,sigma_t=sigma_t,K=K,X=X,t=t,M=M,buff=buff,
                ACtype=ACtype,vertices=poly)

#same but t=3
t=3
N=50
p0=0.5
lam0=-log(1-p0)
sigma=0.50
sigma_t=0.4 #dispersal sigma
phi=0.95
gamma=0.1
buff=3
X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9))
K=c(10,10,10)
M=175
obstype="bernoulli"
ACtype="metamu2"
data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,sigma_t=sigma_t,K=K,X=X,t=t,M=M,buff=buff,
                ACtype=ACtype,obstype="bernoulli")
inits=list(lam01=lam0,sigma=sigma,gamma=gamma,phi=phi,psi=N[1]/M,sigma_t=sigma_t)
niter=2000
nburn=0
nthin=1
proppars=list(lam0=0.01,sigma=0.02,gamma=0.1,s1x=0.05,s1y=0.05,s2x=0.2,s2y=0.2,sigma_t=0.04)

a=Sys.time()
out=mcmc.OpenSCR(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,ACtype=ACtype,
                 Rcpp=TRUE,keepACs=TRUE)
b=Sys.time()
b-a
plot(mcmc(out$out))
summary(mcmc(out$out))













X=X[[1]]
xlim=c(min(X[,1])-buff,max(X[,1])+buff)
ylim=c(min(X[,2])-buff,max(X[,2])+buff)

k=1 #both guy
i=data$simID_L[k] #Real individual
offset=0.15
plot(X[X[,3]==1,1:2],xlim=xlim,ylim=ylim,pch=4,cex=1.5,lwd=2,xlab="X",ylab="Y")
points(X[X[,3]==2,1]-offset,X[X[,3]==2,2],xlim=xlim,ylim=ylim,pch=4,cex=1.5,lwd=2)
points(X[X[,3]==2,1]+offset,X[X[,3]==2,2],xlim=xlim,ylim=ylim,pch=4,cex=1.5,lwd=2)
who=out$ID_Lout[,k]
for(j in 1:(niter-nburn)){
  points(out$s1xout[j,who[j]],out$s1yout[j,who[j]],pch=20,col=rgb(1,0,0,0.1),cex=0.1)
}
# Nfixed=max(data$IDknown)
Nfixed=0
capleft=which(rowSums(data$left)>0&1:nrow(data$left)>Nfixed)
capright=which(rowSums(data$right)>0&1:nrow(data$right)>Nfixed)
notboth=which(rowSums(data$z)>0&1:M>Nfixed)
points(data$s[data$IDknown,1,],col="forestgreen",pch=20,cex=2)
points(data$s[notboth,1,1],data$s[notboth,1,2],col="black",pch=20,cex=1)
singles= unique(c(data$simID_L[capleft],data$simID_R[capright]))
points(data$s[singles,1,1],data$s[singles,1,2],col="gold",pch=20,cex=2)
points(data$s[i,1,1],data$s[i,1,2],col="red",pch=20,cex=2)

lefts=setdiff(data$simID_L,data$simID_R)
rights=setdiff(data$simID_R,data$simID_L)
LRs=data$simID_L[(Nfixed+1):length(data$simID_L)][data$simID_L[(Nfixed+1):length(data$simID_L)]%in%data$simID_R[(Nfixed+1):length(data$simID_R)]]
text(x=data$s[lefts,1,1],y=data$s[lefts,1,2],rep("L",length(lefts)),cex=0.75)
text(x=data$s[rights,1,1],y=data$s[rights,1,2],rep("R",length(lefts)),cex=0.75)
text(x=data$s[LRs,1,1],y=data$s[LRs,1,2],rep("LR",length(lefts)),cex=0.75)
text(x=data$s[1:Nfixed,1,1],y=data$s[1:Nfixed,1,2],rep("B",Nfixed),cex=0.75)

#plot left caps
offset2=0.25
jit=0.2
l=which(data$simID_L==i)
if(length(l)>0){
  leftcaps=X[data$ID_L[which(apply(data$left,c(1,3),sum)[l,]>0)],]
  text(x=jitter(leftcaps[,1],jit),y=jitter(leftcaps[,2]-offset2,jit),labels=rep("L",nrow(leftcaps)),cex=1.25,col="darkblue")
}
r=which(data$simID_R==i)
if(length(r)>0){
  rightcaps=X[data$ID_R[which(apply(data$right,c(1,3),sum)[r,]>0)],]
  text(x=jitter(rightcaps[,1],jit),y=jitter(rightcaps[,2]+offset2,jit),labels=rep("R",nrow(rightcaps)),cex=1.25,col="darkblue")
}
