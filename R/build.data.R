#' Format "secr-like" data for use in mcmc.OpenSCR or mcmc.OpenSCR.sex
#' @param dataIn a matrix with 4 or 5 columns.  If observations are already summed over occasions (K index), the
#' columns indicate the individual number, trap number, primary period number, and number of observations for
#' that individual at that trap on that primary period.  If the observations are not summed over occasion, the
#' colums indicate the individual number, trap number, secondary occasion number, primary period number and number of
#' observation for that individual at that trap on that secondary occasion on that primary period.  All individuals, traps,
#' secondary occasions, and primary periods should be numbered sequentially starting at 1.
#' @param X a list with elements that consists of J[l] x 2 matrices housing the X and Y trap locations for each primary period, l.
#' If the traps do not vary across primary periods, just repeat the same traps in each list element. If population was not
#' observed in primary period l, enter NA instead of a matrix.  The row numbers of J[l] correspond to the trap
#' numbers in dataIn
#' @param K a vector of integers indicating the number of secondary occasions within each primary period
#' @param obstype a character indicating the observation model "bernoulli" or "poisson"
#' @param buff an optional numeric specifying the buffer for the traps to produce the state space.  It is applied to the minimum and maximum
#' X and Y trap locations across primary periods, producing a square or rectangular state space.
#' @param vertices an optional *list* of matrices with the X and Y coordinates
#' of a polygonal state space with one polygon in each list element.  If there is just one polygon, the list is of length 1.
#'  @param primary an optional vector of length t with entries 1 if data was recorded in that primary period and 0 if not. This allows
#' population dynamics to occur at equal interval primary periods even if data was not recorded at each primary period.
#' If not entered, the population is assumed to have been sampled in all primary periods.
#' @param sex an vector of sexes for use in the sex-specific sampler.  Each of the n individuals should have an
#' entry with 1 indicating male, 2 female, and NA unknown.
#' @param tf an optional list of vectors of length t containing the trap operation information.  Each vector has one element for each trap
#' and indicates how many occasions each trap was operational in that primary period.
#' @author Ben Augustine
#' @description
#' @examples
#' \dontrun{
#' #simulate some data
#' t=3
#' N=50
#' p0=0.3
#' lam0=-log(1-p0)
#' sigma=1
#' phi=0.7
#' gamma=0.3
#' buff=2
#' X=list(expand.grid(4:11,4:11),expand.grid(4:11,4:11),expand.grid(4:11,4:11))
#' K=c(10,10,10)
#' M=125
#' obstype="bernoulli"
#' data=simOpenSCR(N=N,phi=phi,gamma=gamma,lam0=lam0,sigma=sigma,K=K,X=X,t=t,M=M,buff=buff,
#'               obstype=obstype,)
#' #Create "secr-like" input matrix
#' idx=which(data$y>0,arr.ind=TRUE)
#' dataIn=matrix(NA,nrow=nrow(idx),ncol=4)
#' for(i in 1:nrow(idx)){
#'   dataIn[i,]=c(idx[i,],data$y[idx[i,1],idx[i,2],idx[i,3]])
#' }
#' #Your input data should have this format if it is summed across secondary occasions
#' str(dataIn)
#' data.formatted=build.data(dataIn,X, K, obstype="bernoulli",buff=3)
#' str(data.formatted)
#' all(data.formatted$y==data$y)
#'}
#'@export
build.data=function(dataIn=NA,X=NA, K=NA, obstype=NA,buff=NA,vertices=NA,primary=NA,sex=NA,tf=NA){
  #entry check
  if(any(is.na(X))){
    stop("Must enter X")
  }
  if(any(e2dist(X[[1]],X[[1]])>1000,na.rm=TRUE)){
    warning("Are trap units in meters? Consider switching to km.")
  }
  if(!is.list(X)){
    stop("X must be a list with one element per primary period")
  }
  if(is.na(K[1])){
    stop("Must enter K")
  }
  if(is.na(obstype)){
    stop("Must enter obstype")
  }
  if(!obstype%in%c("bernoulli","poisson")){
    stop("obstype must be bernoulli or poisson")
  }
  if(any(unlist(lapply(X,function(x){any(is.na(x))})))){
    stop("no missing values allowed in the trap locations")
  }
  J=unlist(lapply(X,nrow))
  maxJ=max(J)
  n=max(dataIn[,1])
  if(ncol(dataIn)==4){
    t=max(dataIn[,3])
  }else{
    t=max(dataIn[,4])
  }
  nobs=nrow(dataIn)
  if(length(J)!=t){
    stop("J must be of length t")
  }
  if(length(K)!=t){
    stop("K must be of length t")
  }
  #data check
  if(!all(round(dataIn)==dataIn)){
    stop("All entries of dataIn must be integers")
  }
  if(!is.numeric(dataIn)){
    stop("dataIn must all be sequentially numbered. e.g ID number, trap number, year number")
  }
  if(!length(unique(sort(dataIn[,1])))==(max(dataIn[,1])-min(dataIn[,1])+1)){
    stop("Must number individuals consecutively from 1 to n")
  }
  if(ncol(dataIn)==5){
    for(l in 1:t){
      if(!all(dataIn[dataIn[,4]==l,3]<=K[l])){
        stop(paste("some occasion numbers in primary period",l,"are larger than K[l]"))
      }
    }
    for(l in 1:t){
      if(any(dataIn[,4]==l&dataIn[,2]>J[l])){
        stop(paste("Trap number in dataIn on primary period",l,"is larger than J[l]"))
      }
    }
  }else{
    for(l in 1:t){
      if(any(dataIn[,3]==l&dataIn[,2]>J[l])){
        stop(paste("Trap number in dataIn on primary period",l,"is larger than J[l]"))
      }
    }
  }
  #sex check
  if(!is.na(sex[1])){
    if(length(sex)!=n)
      stop("if using sex, must enter a sex for all n captured individuals. Enter NA if unknown")
  }
  #primary check
  if(!is.na(primary[1])){
    if(primary[1]==0){
      stop("first primary period must be a 1 (data recorded)")
    }
    if(any(!primary%in%c(0,1))){
      stop("primary must have entries 0 or 1")
    }
    if(any(primary[unique(dataIn[,3])]==0)){
      stop("observations recorded in year with primary=0")
    }
    if(any(is.na(X)&primary==0)){
      stop("Traps entered for primary=0 period. Switch to NA.")
    }
  }else{
    primary=rep(1,t)
  }
  #vertices check
  if(!is.na(vertices[1])){
    if(!is.list(vertices)){
      stop("vertices must be a list")
    }
  }
  #buffer check
  if(is.na(buff)&is.na(vertices[1])){
    stop("Must enter buff or vertices")
  }
  #trap file check
  if(!is.na(tf[1])){
    if(length(tf)!=length(X)){
      stop("If using a trap operation file, must input one for each year")
    }
    for(i in 1:t){
      if(primary[l]==1){
        if(is.matrix(tf[[l]])){
          stop("If using trap operation file, must enter vector of number of occasions operational, not matrix of trap by occasion operation")
        }
        if(nrow(X[[l]])!=length(tf[[l]])){
          stop("If using trap operation file, must enter operation for every trap in every primary=1 year")
        }
      }
      if(any(tf[[l]]>K[l])){
        stop(paste("Some trap file operation occasions are larger than K[l] on occasion",l))
      }
    }
  }
  #build y array
  if(ncol(dataIn)==4){
    dataOut=array(0,dim=c(n,maxJ,t))
    for(i in 1:nobs){
      dataOut[dataIn[i,1],dataIn[i,2],dataIn[i,3]]=dataOut[dataIn[i,1],dataIn[i,2],dataIn[i,3]]+dataIn[i,4]
    }
  }else{#4D data
    warning("collapsing n x j x k x l data to n x j x l")
    maxK=max(K)
    dataOut=array(0,dim=c(n,maxJ,maxK,t))
    for(i in 1:nobs){
      dataOut[dataIn[i,1],dataIn[i,2],dataIn[i,3],dataIn[i,4]]=dataOut[dataIn[i,1],dataIn[i,2],dataIn[i,3],dataIn[i,4]]+dataIn[i,5]
    }
    dataOut=apply(dataOut,c(1,2,4),sum)
  }
  data=list(y=dataOut,X=X,K=K,n=n,J=J,obstype=obstype,primary=primary,sex=sex,tf=tf,vertices=vertices,buff=buff)
  if(is.na(vertices[1])){
    data$vertices=NULL
  }
  if(is.na(tf[1])){
    data$tf=NULL
  }
  if(length(sex)==1){
    data$sex=NULL
  }
  if(is.na(buff)){
    data$buff=NULL
  }
  return(data)
}
