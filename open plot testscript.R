t=4
N=25
M=450
p0=rep(0.5,1)
lam0=-log(1-p0)
sigma=rep(0.750,1)
sigma_t=NULL #dispersal sigma
gamma=0.3
phi=rep(0.7,1)
buff=2
X=list(expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9),expand.grid(4:9,4:9))
K=c(8,8,8,8)


data=simOpenSCR(N=N,gamma=gamma,phi=phi,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,obstype="bernoulli")
simOpenSCR.plot(data,delay=2)
#individual or time covariates for per cap recruitment or survival
#Density independent or density dependent growth
#N_t-1 all or female only
#habitat covariates for recruitment or survival

#movement
sigma_t=0.5
buff=4
data=simOpenSCR(N=N,gamma=gamma,phi=phi,lam0=lam0,sigma=sigma,sigma_t=sigma_t,K=K,X=X,t=t,M=M,buff=buff,obstype="bernoulli")
simOpenSCR.plot(data,animated=FALSE,delay=2)

sigma_t=1
buff=4
data=simOpenSCR(N=N,gamma=gamma,phi=phi,lam0=lam0,sigma=sigma,sigma_t=sigma_t,K=K,X=X,t=t,M=M,buff=buff,obstype="bernoulli")
simOpenSCR.plot(data,animated=FALSE,delay=2)


sigma_t=1
buff=2
data=simOpenSCR(N=N,gamma=gamma,phi=phi,lam0=lam0,sigma=sigma,sigma_t=sigma_t,K=K,X=X,t=t,M=M,buff=buff,obstype="bernoulli")
simOpenSCR.plot(data,animated=TRUE,delay=2)

#apparent vs. true survival
#meta mu vs. markov movement
#recruitment dispersal kernel
#habitat movement covariates
