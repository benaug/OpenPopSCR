////////////////////////////////////in polygon functions/////////////////////
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
bool intersectCppOpen(NumericVector sx, NumericVector sy,NumericVector vertex1, NumericVector vertex2) {
  NumericVector swap;
  double m_blue;
  double m_red;
  bool out;
  if((sy(0)==vertex1[1])|(sy(0)==vertex2[1])){
    sy(0)=sy(0)+0.000001;
  }
  if(vertex1[1]>vertex2[1]){
    swap=vertex1;
    vertex1=vertex2;
    vertex2=swap;
  }
  if((sy(0)<vertex1[1])|(sy(0)>vertex2[1])){
    out=FALSE;
  }else if((sx(0) > vertex1[0]) & (sx(0)> vertex2[0])){
    out=FALSE;
  }else{
    if((sx(0) < vertex1[0]) & (sx(0) < vertex2[0])){
      out=TRUE;
    }else{
      if(vertex1[0]!=vertex2[0]){
        m_red=(vertex2(1)-vertex1(1))/(vertex2(0)-vertex1(0));
      }else{
        m_red=1000000000000000000;
      }
      if(vertex1[0]!=sx(0)){
        m_blue=(sy(0)-vertex1(1))/(sx(0)-vertex1(0));
      }else{
        m_blue=1000000000000000000;
      }
      if(m_blue>=m_red){
        out=TRUE;
      }else{
        out=FALSE;
      }
    }
  }
  return out;
}


#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
bool inoutCppOpen(NumericVector sx,NumericVector sy,NumericMatrix vertices) {
  bool to_return=false;
  int count=0;
  int I=vertices.nrow();
  for(int i=0; i<(I-1); i++) {
    if(intersectCppOpen(sx,sy,vertices(i,_),vertices(i+1,_))){
      count=count+1;
    }
  }
  if(count % 2 != 0){
    to_return = true;
  }
  return to_return;
}


//Full open pop MCMC sampler
//[[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::export]]
List mcmc_Open(NumericVector lam0, NumericVector sigma, NumericVector gamma,NumericVector gammaprime, NumericVector phi,
               arma::cube D,arma::cube lamd, arma::cube y,IntegerMatrix z,IntegerMatrix a, NumericMatrix s1,arma::cube s2,
               int ACtype, bool useverts,List vertices,NumericVector xlim,NumericVector ylim,
               IntegerMatrix knownmatrix,IntegerVector Xidx, arma::cube Xcpp,IntegerVector K,NumericMatrix Ez, double psi,
               IntegerVector N,NumericVector proplam0, NumericVector propsig,NumericVector propz, NumericVector propgamma,double props1x,
               double props1y,double props2x,double props2y, double propsigma_t,NumericVector sigma_t,
               int niter, int nburn, int nthin,int npar,IntegerVector each,bool jointZ,IntegerMatrix zpossible,
               IntegerMatrix apossible,IntegerMatrix cancel,int obstype,IntegerVector tf,NumericMatrix dSS,bool usedSS) {
  RNGScope scope;
  int M = size(lamd)[0];
  int J = size(lamd)[1];
  int t = size(lamd)[2];
  //Preallocate detection function
  NumericVector lam0cand(t);
  NumericVector sigmacand(t);
  NumericVector rand;
  NumericVector rand2;
  arma::cube pd=zeros<cube>(M,J,t);
  arma::cube ll_y_curr=zeros<cube>(M,J,t);
  arma::cube lamdcand=zeros<cube>(M,J,t);
  arma::cube pdcand=zeros<cube>(M,J,t);
  arma::cube ll_y_cand=zeros<cube>(M,J,t);
  //Structures for dealing with year specific parms
  NumericVector likcurr2D(t);
  NumericVector likcand2D(t);
  NumericVector sigmause(t);
  NumericVector lam0use(t);
  for(int l=0; l<t; l++) {
    if(sigma.size()==1){
      sigmause(l)=sigma(0);
    }else{
      sigmause(l)=sigma(l);
    }
  }
  for(int l=0; l<t; l++) {
    if(lam0.size()==1){
      lam0use(l)=lam0(0);
    }else{
      lam0use(l)=lam0(l);
    }
  }
  //Preallocate Z update
  IntegerVector z1cand(M);
  IntegerVector a1cand(M);
  NumericMatrix ll_z(M,t);
  NumericMatrix Ezcand(M,t);
  NumericMatrix zcand(M,t);
  NumericMatrix acand(M,t);

  NumericMatrix ll_z_cand(M,t);
  double llysum=0;
  double llycandsum=0;
  double llzsum=0;
  double llzcandsum=0;
  bool fix1=FALSE;
  bool fix2=FALSE;
  int sumz=0;
  int sumz1tmp=0;
  int suma1tmp=0;
  int suma=0;
  LogicalVector warn(M,FALSE);
  int warncount=0;
  LogicalVector upz(M);
  LogicalVector upz3(M);
  IntegerVector upz2(M);
  IntegerVector swapz(M);
  IntegerVector latecaps(M,0);
  IntegerVector Ntmp(t);
  //Preallocate z[,2+]
  IntegerVector at_cand(M);
  IntegerVector zt_cand(M);
  int navail=0;
  int idx=0;
  int propzuse=0;
  //Preallocate jointZ
  int nzpossible=zpossible.nrow();
  NumericMatrix Ezpossible(nzpossible,t-1);
  NumericMatrix llzpossible(nzpossible,t);
  LogicalVector fixed(t);
  NumericVector propto(nzpossible);
  NumericVector propto1(nzpossible);
  NumericVector propto2(nzpossible);
  IntegerVector zchoose(1);
  IntegerVector choose=Rcpp::seq(0,(nzpossible-1));
  IntegerVector zprop(t);
  IntegerVector aprop(t);
  double propprob;
  double backprob;
  int currz=0;
  double sumpropto=0;

  //Preallocate phi and gamma
  int survive=0;
  int dead=0;
  bool gamma_cand_ok=TRUE;
  double gamma_cand=0;
  NumericVector gammaprimecand(t-1);
  //Fixed or year-specific stuff
  NumericVector gammacand(t-1);
  NumericVector phicand(t-1);
  NumericVector gammause(t-1);
  NumericVector phiuse(t-1);
  for(int l=0; l<(t-1); l++) {
    if(gamma.size()==1){
      gammause(l)=gamma(0);
    }else{
      gammause(l)=gamma(l);
    }
    if(phi.size()==1){
      phiuse(l)=phi(0);
    }else{
      phiuse(l)=phi(l);
    }
  }
  //Preallocate for updating activity centers
  LogicalVector inbox(1);
  NumericVector inbox2(1);
  NumericMatrix dtmp(J,t);
  NumericVector ScandX(1);
  NumericVector ScandY(1);
  //metamu stuff
  // NumericMatrix ll_s2(M,t);
  // NumericMatrix ll_s2_cand(M,t);
  double lls2sum=0;
  double lls2candsum=0;
  NumericVector sigma_t_cand(1);
  //Structures to record output
  int nstore=(niter-nburn)/nthin;
  if(nburn % nthin!=0){
    nstore=nstore+1;
  }
  NumericMatrix out(nstore,npar);
  NumericMatrix s1xout(nstore,M);
  NumericMatrix s1yout(nstore,M);
  arma::cube s2xout(nstore,M,t);
  arma::cube s2yout(nstore,M,t);
  arma::cube zout(nstore,M,t);

  //Can update anyone who wasn't captured on every occasion
  for(int i=0; i<M; i++){
    sumz=0;//reusing sum z to sum known.matrix
    for(int l=0; l<t; l++){
      sumz+=knownmatrix(i,l);
    }
    if(sumz<t){
      upz3(i)=TRUE;
    }else{
      upz3(i)=FALSE;
    }
  }
  int polys=vertices.size();//number of polygons for SS


  int iteridx=0;
  //////Calculate starting log likelihoods///////
  //ll.s2
  NumericMatrix ll_s2(M,t);
  NumericMatrix ll_s2_cand(M,t);
  if((ACtype==2)|(ACtype==5)){
    for(int l=0; l<t; l++){
      for(int i=0; i<M; i++) { //X and Y normal log-likelihood simplified
        ll_s2(i,l)=-log(pow(sigma_t(0),2.0))-(1/(2*pow(sigma_t(0),2.0)))*(pow(s2(i,l,0)-s1(i,0),2.0)+pow(s2(i,l,1)-s1(i,1),2.0));
      }
    }
  }else if(ACtype==3){
    for(int l=1; l<t; l++){
      for(int i=0; i<M; i++) { //X and Y normal log-likelihood simplified
        ll_s2(i,l-1)=-log(pow(sigma_t(0),2.0))-(1/(2*pow(sigma_t(0),2.0)))*(pow(s2(i,l,0)-s2(i,l-1,0),2.0)+pow(s2(i,l,1)-s2(i,l-1,1),2.0));
        ll_s2_cand(i,l-1)=ll_s2(i,l-1);
      }
    }
  }
  //  Detection function
  if(obstype==1){
    for(int l=0; l<t; l++){
      likcurr2D(l)=0;
      for(int i=0; i<M; i++) {
        for(int j=0; j<Xidx(l); j++){
          pd(i,j,l)=1-exp(-lamd(i,j,l));
          ll_y_curr(i,j,l)=z(i,l)*(y(i,j,l)*log(pd(i,j,l))+(tf(j,l)-y(i,j,l))*log(1-pd(i,j,l)));
          if(ll_y_curr(i,j,l)==ll_y_curr(i,j,l)){
            likcurr2D(l)+=ll_y_curr(i,j,l); //year-specific components
          }
        }
      }
      llysum+=likcurr2D(l); //full likelihood sum
    }
  }else{
    for(int l=0; l<t; l++){
      likcurr2D(l)=0;
      for(int i=0; i<M; i++) {
        for(int j=0; j<Xidx(l); j++){
          ll_y_curr(i,j,l)=z(i,l)*(y(i,j,l)*log(tf(j,l)*lamd(i,j,l))-tf(j,l)*lamd(i,j,l));
          if(ll_y_curr(i,j,l)==ll_y_curr(i,j,l)){
            likcurr2D(l)+=ll_y_curr(i,j,l); //year-specific components
          }
        }
      }
      llysum+=likcurr2D(l); //full likelihood sum
    }
  }
  //ll.z. Some Ez are so small we're close to log(0) which is NaN in Rcpp
  for(int i=0; i<M; i++){//z1
    ll_z(i,0)= z(i,0)*log(psi)+(1-z(i,0))*log(1-psi);
  }
  for(int l=1; l<t; l++){//z2+
    for(int i=0; i<M; i++){
      ll_z(i,l)= z(i,l)*log(Ez(i,l-1))+(1-z(i,l))*log(1-Ez(i,l-1));
      if(ll_z(i,l)!=ll_z(i,l)){//Turn NaNs to 0s
        ll_z(i,l)=0;
      }
    }
  }
  //Set up discrete state space
  int NdSS=dSS.nrow();
  NumericVector dists(NdSS);

  //Here we go!
  for(int iter=0; iter<niter; iter++){
    /////////////////Detection function update///////////////////////
    // Need to resum the ll_y on each iter
    llysum=0;
    for(int l=0; l<t; l++){
      likcurr2D(l)=0;
      for(int i=0; i<M; i++) {
        for(int j=0; j<Xidx(l); j++){
          if(ll_y_curr(i,j,l)==ll_y_curr(i,j,l)){
            likcurr2D(l)+=ll_y_curr(i,j,l);
          }
        }
      }
      llysum+=likcurr2D(l); //full likelihood sum
    }
    //Update lam0
    if(lam0.size()==1){//fixed lambda
      rand=Rcpp::rnorm(1,lam0(0),proplam0(0));
      if(rand(0) > 0){
        llycandsum=0;
        lam0cand(0)=rand(0);
        //  Update lamd and calculate cand likelihood
        for(int l=0; l<t; l++){
          likcand2D(l)=0;
          for(int i=0; i<M; i++) {
            for(int j=0; j<Xidx(l); j++){
              lamdcand(i,j,l)=lam0cand(0)*exp(-D(i,j,l)*D(i,j,l)/(2*sigmause(l)*sigmause(l)));
              if(obstype==1){
                pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(tf(j,l)-y(i,j,l))*log(1-pdcand(i,j,l)));
              }else{
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(tf(j,l)*lamdcand(i,j,l))-tf(j,l)*lamdcand(i,j,l));
              }
              if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                likcand2D(l)+=ll_y_cand(i,j,l);
              }
            }
          }
          llycandsum+=likcand2D(l);
        }
        rand2=Rcpp::runif(1);
        if(rand2(0)<exp(llycandsum-llysum)){
          lam0(0)=lam0cand(0);
          for(int l=0; l<t; l++){
            for(int i=0; i<M; i++) {
              for(int j=0; j<Xidx(l); j++){
                lamd(i,j,l)=lamdcand(i,j,l);
                if(obstype==1){
                  pd(i,j,l)=pdcand(i,j,l);
                }
                ll_y_curr(i,j,l)=ll_y_cand(i,j,l);
              }
            }
            likcurr2D(l)=likcand2D(l);
          }
          llysum=llycandsum;
        }
      }
    }else{//year-specific lambdas
      llysum=0;
      for(int l=0; l<t; l++){
        likcand2D(l)=0;
        rand=Rcpp::rnorm(1,lam0(l),proplam0(l));
        if(rand(0) > 0){
          lam0cand(l)=rand(0);
          //  Calculate likelihood
          for(int i=0; i<M; i++) {
            for(int j=0; j<Xidx(l); j++){
              lamdcand(i,j,l)=lam0cand(l)*exp(-D(i,j,l)*D(i,j,l)/(2*sigmause(l)*sigmause(l)));
              if(obstype==1){
                pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(tf(j,l)-y(i,j,l))*log(1-pdcand(i,j,l)));
              }else{
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(tf(j,l)*lamdcand(i,j,l))-tf(j,l)*lamdcand(i,j,l));
              }
              if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                likcand2D(l)+=ll_y_cand(i,j,l);
              }
            }
          }
          rand2=Rcpp::runif(1);
          if(rand2(0)<exp(likcand2D(l)-likcurr2D(l))){
            lam0(l)=lam0cand(l);
            for(int i=0; i<M; i++) {
              for(int j=0; j<Xidx(l); j++){
                ll_y_curr(i,j,l)=ll_y_cand(i,j,l);
                if(obstype==1){
                  pd(i,j,l)=pdcand(i,j,l);
                }
                lamd(i,j,l)=lamdcand(i,j,l);
              }
            }
            likcurr2D(l)=likcand2D(l);
          }
        }
        llysum+=likcurr2D(l);
      }
    }
    //fill lam0use
    for(int l=0; l<t; l++) {
      if(lam0.size()==1){
        lam0use(l)=lam0(0);
      }else{
        lam0use(l)=lam0(l);
      }
    }
    // Update sigma
    if(sigma.size()==1){//fixed sigma
      rand=Rcpp::rnorm(1,sigma(0),propsig(0));
      if(rand(0) > 0){
        sigmacand(0)=rand(0);
        llycandsum=0;
        //  Update lamd and calculate cand likelihood
        for(int l=0; l<t; l++){
          likcand2D(l)=0;
          for(int i=0; i<M; i++) {
            for(int j=0; j<Xidx(l); j++){
              lamdcand(i,j,l)=lam0use(l)*exp(-D(i,j,l)*D(i,j,l)/(2*sigmacand(0)*sigmacand(0)));
              if(obstype==1){
                pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(tf(j,l)-y(i,j,l))*log(1-pdcand(i,j,l)));
              }else{
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(tf(j,l)*lamdcand(i,j,l))-tf(j,l)*lamdcand(i,j,l));
              }
              if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                likcand2D(l)+=ll_y_cand(i,j,l);
              }
            }
          }
          llycandsum+=likcand2D(l);
        }
        rand=Rcpp::runif(1);
        if(rand(0)<exp(llycandsum-llysum)){
          sigma(0)=sigmacand(0);
          for(int l=0; l<t; l++){
            for(int i=0; i<M; i++) {
              for(int j=0; j<Xidx(l); j++){
                lamd(i,j,l)=lamdcand(i,j,l);
                if(obstype==1){
                  pd(i,j,l)=pdcand(i,j,l);
                }
                ll_y_curr(i,j,l)=ll_y_cand(i,j,l);
              }
            }
            likcurr2D(l)=likcand2D(l);
          }
          llysum=llycandsum;
        }
      }
    }else{//year-specific sigmas
      llysum=0;
      for(int l=0; l<t; l++){
        rand=Rcpp::rnorm(1,sigma(l),propsig(l));
        likcand2D(l)=0;
        if(rand(0) > 0){
          sigmacand(l)=rand(0);
          //  Calculate likelihood
          for(int i=0; i<M; i++) {
            for(int j=0; j<Xidx(l); j++){
              lamdcand(i,j,l)=lam0use(l)*exp(-D(i,j,l)*D(i,j,l)/(2*sigmacand(l)*sigmacand(l)));
              if(obstype==1){
                pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(tf(j,l)-y(i,j,l))*log(1-pdcand(i,j,l)));
              }else{
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(tf(j,l)*lamdcand(i,j,l))-tf(j,l)*lamdcand(i,j,l));

              }
              if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                likcand2D(l)+=ll_y_cand(i,j,l);
              }
            }
          }
          rand2=Rcpp::runif(1);
          if(rand2(0)<exp(likcand2D(l)-likcurr2D(l))){
            sigma(l)=sigmacand(l);
            likcurr2D(l)=likcand2D(l);
            for(int i=0; i<M; i++) {
              for(int j=0; j<Xidx(l); j++){
                ll_y_curr(i,j,l)=ll_y_cand(i,j,l);
                if(obstype==1){
                  pd(i,j,l)=pdcand(i,j,l);
                }
                lamd(i,j,l)=lamdcand(i,j,l);
              }
            }
          }
        }
        llysum+=likcurr2D(l);
      }
    }
    //fill sigmause
    for(int l=0; l<t; l++) {
      if(sigma.size()==1){
        sigmause(l)=sigma(0);
      }else{
        sigmause(l)=sigma(l);
      }
    }
    if(jointZ==FALSE){
      ////////////////Z1 stuff////////////////////
      // Figure out who can be updated
      if(t==2){
        for(int i=0; i<M; i++){
          upz(i)=(knownmatrix(i,0)==0);
        }
      }else{
        for(int i=0; i<M; i++){
          latecaps(i)=0;
          for(int l=2; l<t; l++){
            latecaps(i)+=z(i,l);
          }
          upz(i)=(!((z(i,0)==0)&(z(i,1)==0)&(latecaps(i)>0)))&(knownmatrix(i,0)==0);//Don't turn on a guy that is turned on later, but not the next occasion
        }
      }
      //update z[,1]
      N(0)=0;
      for(int i=0; i<M; i++){
        if(upz(i)){
          for(int l=1; l<t; l++) {
            gammaprimecand(l-1)=gammaprime(l-1);
          }
          for(int i2=0; i2<M; i2++){
            if(z(i2,0)==1){
              z1cand(i2) = 1;
            }else{
              z1cand(i2) = 0;
            }
            if(a(i2,0)==1){
              a1cand(i2) = 1;
            }else{
              a1cand(i2) = 0;
            }
          }
          if(z1cand(i)==1){
            z1cand(i)=0;
            a1cand(i)=1;
          }else{
            z1cand(i)=1;
            a1cand(i)=0;
          }
          sumz1tmp=0;
          for(int i2=0; i2<M; i2++){
            sumz1tmp+=z1cand(i2);
          }
          //sum y ll across j dimension for each i and l=0
          llysum=0;
          llycandsum=0;
          for(int j=0; j<Xidx(0); j++){
            if(obstype==1){
              ll_y_cand(i,j,0) =z1cand(i)*(y(i,j,0)*log(pd(i,j,0))+(tf(j,0)-y(i,j,0))*log(1-pd(i,j,0)));
            }else{
              ll_y_cand(i,j,0)=z1cand(i)*(y(i,j,0)*log(tf(j,0)*lamd(i,j,0))-tf(j,0)*lamd(i,j,0));
            }
            if(ll_y_cand(i,j,0)==ll_y_cand(i,j,0)){
              llycandsum+=ll_y_cand(i,j,0);
            }
            if(ll_y_curr(i,j,0)==ll_y_curr(i,j,0)){
              llysum+=ll_y_curr(i,j,0);
            }
          }
          sumz=0;
          for(int l=0; l<t; l++){
            sumz+= z(i,l);
          }
          if((((z1cand(i)==1)&(sumz==0))|((z1cand(i)==0)&(z(i,0)==1)&(sumz==1)))&(t>2)){//Are we turning on a guy that was never on before? or turning off a guy that was only on on z1?
            for(int l=0; l<t; l++){
              for(int i2=0; i2<M; i2++){
                if(z(i2,l)==1){
                  zcand(i2,l) = 1;
                }else{
                  zcand(i2,l) = 0;
                }
                if(a(i2,l)==1){
                  acand(i2,l) = 1;
                }else{
                  acand(i2,l) = 0;
                }
              }
            }
            if(zcand(i,0)==1){
              zcand(i,0)=0;
            }else{
              zcand(i,0)=1;
            }
            if((z1cand(i)==1)&(sumz==0)){//if caught on 1st occasion, turn availability all off
              for(int l=0; l<t; l++){
                acand(i,l)=0;
              }
            }else{  //if never caught, turn availability all on
              for(int l=0; l<t; l++){
                acand(i,l)=1;
              }
            }
            for(int l=1; l<t; l++){
              Ntmp(l)=N(l);
            }
            Ntmp(0)=sumz1tmp;
            for(int l=1; l<t; l++){
              suma=0;
              for(int i2=0; i2<M; i2++){
                suma+=acand(i2,l-1);
              }
              gammaprimecand(l-1)=(Ntmp(l-1)*gammause(l-1))/suma;
              warn(i)=FALSE;
              if(gammaprimecand(l-1) > 1) { // E(Recruits) must be < nAvailable
                warn(i)=TRUE;
              }
            }
            if(warn(i)==FALSE){
              ll_z_cand(i,0) = zcand(i,0)*log(psi)+(1-zcand(i,0))*log(1-psi);
              llzcandsum=0;
              llzsum=0;
              llzcandsum+=ll_z_cand(i,0);//add ll.z[i,1] for focal guy
              llzsum+=ll_z(i,0);
              for(int l=1; l<t; l++){
                for(int i2=0; i2<M; i2++){
                  Ezcand(i2,l-1)=zcand(i2,l-1)*phiuse(l-1) + acand(i2,l-1)*gammaprimecand(l-1);
                  ll_z_cand(i2,l) = zcand(i2,l)*log(Ezcand(i2,l-1))+(1-zcand(i2,l))*log(1-Ezcand(i2,l-1));
                  if(ll_z_cand(i2,l)!=ll_z_cand(i2,l)){
                    ll_z_cand(i2,l)=0;
                  }
                  llzcandsum+=ll_z_cand(i2,l);
                  llzsum+=ll_z(i2,l);
                }
              }
              rand=Rcpp::runif(1);
              if(rand(0) < exp((llycandsum+ llzcandsum)-(llysum+llzsum ))) {
                for(int j=0; j<Xidx(0); j++){
                  ll_y_curr(i,j,0) = ll_y_cand(i,j,0);
                }
                if(zcand(i,0)==1){
                  z(i,0)=1;
                }else{
                  z(i,0)=0;
                }
                for(int l=0; l<t; l++){//Only changed focal individual
                  if(acand(i,l)==1){
                    a(i,l)=1;
                  }else{
                    a(i,l)=0;
                  }
                }
                ll_z(i,0)=ll_z_cand(i,0);
                for(int l=1; l<t; l++){
                  for(int i2=0; i2<M; i2++){
                    Ez(i2,l-1) =Ezcand(i2,l-1);
                    ll_z(i2,l) = ll_z_cand(i2,l);
                  }
                  gammaprime(l-1)=gammaprimecand(l-1);
                }
              }
            }
          }else{//Don't need to modify more than 1 year
            suma1tmp=0;
            for(int i2=0; i2<M; i2++){
              suma1tmp+=a1cand(i2);
            }
            gammaprimecand(0)=sumz1tmp*gammause(0)/suma1tmp;
            warn(i)=FALSE;
            if(gammaprimecand(0) > 1) { // E(Recruits) must be < nAvailable
              warn(i)=TRUE;
            }
            if(warn(i)==FALSE){
              ll_z_cand(i,0) = z1cand(i)*log(psi)+(1-z1cand(i))*log(1-psi);
              //sum z ll
              llzcandsum=0;
              llzsum=0;
              llzcandsum+=ll_z_cand(i,0);//add ll.z[i,1] for focal guy
              llzsum+=ll_z(i,0);
              for(int i2=0; i2<M; i2++){
                Ezcand(i2,0)=z1cand(i2)*phiuse(0) + a1cand(i2)*gammaprimecand(0);
                ll_z_cand(i2,1) = z(i2,1)*log(Ezcand(i2,0))+(1-z(i2,1))*log(1-Ezcand(i2,0));
                if(ll_z_cand(i2,1)!=ll_z_cand(i2,1)){
                  ll_z_cand(i2,1)=0;
                }
                llzcandsum+=ll_z_cand(i2,1);
                llzsum+=ll_z(i2,1);
              }
              rand=Rcpp::runif(1);
              if(rand(0) < exp((llycandsum+ llzcandsum)-(llysum+llzsum ))) {
                for(int j=0; j<Xidx(0); j++){
                  ll_y_curr(i,j,0) = ll_y_cand(i,j,0);
                }
                ll_z(i,0)=ll_z_cand(i,0);
                for(int i2=0; i2<M; i2++){
                  Ez(i2,0) = Ezcand(i2,0);
                  ll_z(i2,1) = ll_z_cand(i2,1);
                }
                if(z1cand(i)==1){
                  z(i,0)=1;
                }else{
                  z(i,0)=0;
                }
                if(a1cand(i)==1){
                  a(i,0)=1;
                }else{
                  a(i,0)=0;
                }
                gammaprime(0) = gammaprimecand(0);
              }
            }
          }
        }
        N(0)+=z(i,0);
      }
      // update z[,2+]
      for(int l=1; l<t; l++){
        //figure out who can be updated with upz
        //always remove dead guys except t=2
        if(t==2){
          for(int i=0; i<M; i++){
            upz(i)=(knownmatrix(i,l)==0);
          }
        }else if(t==3){
          if(l==1){
            for(int i=0; i<M; i++){
              upz(i)=(!((z(i,l-1)==1)&(z(i,l+1)==1)))&(knownmatrix(i,l)==0);  //remove guys that are on before and after l. can't be dead guys
            }
          }else{
            for(int i=0; i<M; i++){
              upz(i)=(knownmatrix(i,l)==0)&(!( (z(i,l-1)==0) & (a(i,l-1)==0)));  //remove dead guys
            }
          }
        }else{//t>3
          if(l==(t-1)){ //can update anyone on last occasion unless they're dead
            for(int i=0; i<M; i++){
              upz(i)=(knownmatrix(i,l)==0)&(!( (z(i,l-1)==0) & (a(i,l-1)==0)));
            }
          }else if(l==(t-2)){ //second to last occasion
            for(int i=0; i<M; i++){
              upz(i)=(!((z(i,l-1)==1)&(z(i,l+1)==1)))&(knownmatrix(i,l)==0)&(!( (z(i,l-1)==0) & (a(i,l-1)==0)));  //remove guys that are on before and after l and dead guys
            }
          }else{//l between 2 and t-1
            for(int i=0; i<M; i++){
              latecaps(i)=0;  //used to identify guys you can't turn on because they're on later, but not in next occasion.
              for(int l2=(l+2); l2<t; l2++){
                latecaps(i)+=z(i,l2);
              }
              upz(i)=(!((z(i,l-1)==1)&(z(i,l+1)==1)))&(knownmatrix(i,l)==0)&(!((z(i,l)==0)&(z(i,l+1)==0)&(latecaps(i)>0)))&(!( (z(i,l-1)==0) & (a(i,l-1)==0))); //remove guys that are on before and after l or 0,0,1 and dead guys
            }
          }
        }

        //Can we swap?
        navail=0;
        for(int i=0; i<M; i++){
          if(upz(i)){
            navail+=1;
          }
        }
        // propzuse=0;
        //How many to swap?
        if(navail < propz(l-1)) {
          propzuse=navail;
        }else{
          propzuse=propz(l-1);
        }
        //get upz2
        // IntegerVector upz2(navail,0);
        idx=0;
        for(int i=0; i<M; i++){
          if(upz(i)){
            upz2(idx)=i;
            idx+=1;
          }
        }
        if(navail>0){//if so, proceed with swapping
          //Who to swap?
          IntegerVector swapzidx(propzuse,0);
          IntegerVector choose=Rcpp::seq(0,(navail-1));
          swapzidx=Rcpp::RcppArmadillo::sample(choose,propzuse,FALSE);
          // IntegerVector swapz(propzuse,0);
          for(int i=0; i<propzuse; i++){
            swapz(i)=upz2(swapzidx(i));
            // storeswapz(iter,i,l-1)=swapz(i);
          }
          //Update swapz 1 at a time
          for(int i=0; i<propzuse; i++){
            for(int i2=0; i2<M; i2++){
              if(z(i2,l)==1){
                zt_cand(i2)=1;
              }else{
                zt_cand(i2)=0;
              }
              if(a(i2,l)==1){
                at_cand(i2)=1;
              }else{
                at_cand(i2)=0;
              }
            }
            if(zt_cand(swapz(i))==1){
              zt_cand(swapz(i))=0;
            }else{
              zt_cand(swapz(i))=1;
            }
            if((a(swapz(i),l-1)==1)&(zt_cand(swapz(i))==0)){//who was available on last occasion and not proposed to be captured?
              at_cand(swapz(i))=1;
            }else{
              at_cand(swapz(i))=0;
            }
            //sum ll across j for each l and chosen i
            llycandsum=0;
            llysum=0;
            for(int j=0; j<Xidx(l); j++){
              if(obstype==1){
                ll_y_cand(swapz(i),j,l)=zt_cand(swapz(i))*(y(swapz(i),j,l)*log(pd(swapz(i),j,l))+(tf(j,l)-y(swapz(i),j,l))*log(1-pd(swapz(i),j,l)));
              }else{
                ll_y_cand(swapz(i),j,l)=zt_cand(swapz(i))*(y(swapz(i),j,l)*log(tf(j,l)*lamd(swapz(i),j,l))-tf(j,l)*lamd(swapz(i),j,l));
              }
              if(ll_y_cand(swapz(i),j,l)==ll_y_cand(swapz(i),j,l)){
                llycandsum+=ll_y_cand(swapz(i),j,l);
              }
              if(ll_y_curr(swapz(i),j,l)==ll_y_curr(swapz(i),j,l)){
                llysum+=ll_y_curr(swapz(i),j,l);
              }
            }
            //Add up the z likelihood contributions curr and cand for swapped guys
            ll_z_cand(swapz(i),l) = zt_cand(swapz(i))*log(Ez(swapz(i),l-1))+(1-zt_cand(swapz(i)))*log(1-Ez(swapz(i),l-1));
            if(ll_z_cand(swapz(i),l)!=ll_z_cand(swapz(i),l)){//Turn NaNs to 0s
              ll_z_cand(swapz(i),l)=0;
            }
            llzcandsum=0;
            llzsum=0;
            llzcandsum+=ll_z_cand(swapz(i),l);//prior.z.cand
            llzsum+=ll_z(swapz(i),l);//prior.z
            //If we make this change to z[,l] and a[,l], how does it change ll.z[,l+1]?
            if(t>3){
              fix1=FALSE;
              fix2=FALSE;
              sumz=0;
              for(int l2=0; l2<t; l2++){
                sumz+=z(swapz(i),l2);
              }
              if((zt_cand(swapz(i))==1)&(sumz==0)){//are we turning a guy on that was previously never on?
                fix1=TRUE;
              }
              if((sumz==1)&(zt_cand(swapz(i))==0)&(z(swapz(i),l)==1)){
                fix2=TRUE;
              }
            }
            if((fix1|fix2)&(l<(t-1))){//fix1 and 2 only true if t>3
              for(int l2=0; l2<t; l2++){
                for(int i2=0; i2<M; i2++){
                  if(z(i2,l2)==1){
                    zcand(i2,l2)=1;
                  }else{
                    zcand(i2,l2)=0;
                  }
                  if(a(i2,l2)==1){
                    acand(i2,l2)=1;
                  }else{
                    acand(i2,l2)=0;
                  }
                }
              }
              if(zt_cand(swapz(i))==1){
                zcand(swapz(i),l)=1;
              }else{
                zcand(swapz(i),l)=0;
              }
              if(fix1){
                for(int l2=l; l2<t; l2++){
                  acand(swapz(i),l2)=0; //all l:t a off
                }
              }
              if(fix2){
                for(int l2=l; l2<t; l2++){
                  acand(swapz(i),l2)=1; //all l:t a on
                }
              }
              sumz=0;
              for(int i2=0; i2<M; i2++){
                sumz+=zt_cand(i2);//same as zcand[,1]
              }
              for(int l2=0; l2<t; l2++){
                Ntmp(l2)=N(l2);
              }
              Ntmp(l)=sumz;
              for(int l2=l; l2<(t-1); l2++){
                suma=0;
                for(int i2=0; i2<M; i2++){
                  suma+=acand(i2,l2);
                }
                gammaprimecand(l2)=(Ntmp(l2)*gammause(l2)) / suma;
                if(gammaprimecand(l2) < 1){ //Is this a valid probability?
                  for(int i2=0; i2<M; i2++){
                    //Add on contributions to z ll from l+1 for cand and curr
                    Ezcand(i2,l2) = zcand(i2,l2)*phiuse(l2) + acand(i2,l2)*gammaprimecand(l2);
                    ll_z_cand(i2,l2+1) = z(i2,l2+1)*log(Ezcand(i2,l2))+(1-z(i2,l2+1))*log(1-Ezcand(i2,l2));
                    llzsum+=ll_z(i2,l2+1);//prior.z
                    if(ll_z_cand(i2,l2+1)==ll_z_cand(i2,l2+1)){
                      llzcandsum+=ll_z_cand(i2,l2+1);//prior.z.cand
                    }else{
                      ll_z_cand(i2,l2+1)=0; //Turn NaN to 0
                    }
                  }
                }
              }
            }else{//future years not affected
              if(l<(t-1)){
                sumz=0;
                suma=0;
                for(int i2=0; i2<M; i2++){
                  sumz+=zt_cand(i2);
                  suma+=at_cand(i2);
                }
                gammaprimecand(l)=sumz*gammause(l)/suma;
                if(gammaprimecand(l) < 1){
                  for(int i2=0; i2<M; i2++){
                    Ezcand(i2,l)=zt_cand(i2)*phiuse(l) + at_cand(i2)*gammaprimecand(l);
                    ll_z_cand(i2,l+1) = z(i2,l+1)*log(Ezcand(i2,l))+(1-z(i2,l+1))*log(1-Ezcand(i2,l));
                    llzsum+=ll_z(i2,l+1);//prior.z
                    if(ll_z_cand(i2,l+1)==ll_z_cand(i2,l+1)){
                      llzcandsum+=ll_z_cand(i2,l+1);//prior.z.cand
                    }else{
                      ll_z_cand(i2,l+1)=0; //Turn NaN to 0
                    }
                  }
                }
              }
            }
            rand=Rcpp::runif(1);
            if(rand(0) < exp(( llycandsum+llzcandsum)-( llysum+llzsum) )) {
              for(int j=0; j<Xidx(l); j++){
                ll_y_curr(swapz(i),j,l) = ll_y_cand(swapz(i),j,l);
              }
              ll_z(swapz(i),l)=ll_z_cand(swapz(i),l);
              if(zt_cand(swapz(i))==1){//only changed z for this l
                z(swapz(i),l)=1;
              }else{
                z(swapz(i),l)=0;
              }
              if((fix1|fix2)&(l<(t-1))){
                for(int l2=l; l2<t; l2++){ //need to fill in l:t a for this swapz
                  if(acand(swapz(i),l2)==1){
                    a(swapz(i),l2)=1;
                  }else{
                    a(swapz(i),l2)=0;
                  }
                }
                for(int l2=l; l2<(t-1); l2++){
                  for(int i2=0; i2<M; i2++){
                    Ez(i2,l2)=Ezcand(i2,l2);
                    ll_z(i2,l2+1)=ll_z_cand(i2,l2+1);
                  }
                  gammaprime(l2)=gammaprimecand(l2);
                }
              }else{
                if(at_cand(swapz(i))==1){//only need to fill in a for this l
                  a(swapz(i),l)=1;
                }else{
                  a(swapz(i),l)=0;
                }
                if(l<(t-1)){
                  for(int i2=0; i2<M; i2++){
                    Ez(i2,l)=Ezcand(i2,l);
                    ll_z(i2,l+1)=ll_z_cand(i2,l+1);
                  }
                  gammaprime(l)=gammaprimecand(l);
                }
              }
            }
          }
        }
        // Update N[,2+]
        N(l)=0;
        for(int i=0; i<M; i++){
          N(l)+=z(i,l);
        }
      }
    }else{
      //jointZ update
      //ll.z[,1] won't change across i
      for(int i=0; i<nzpossible; i++){
        llzpossible(i,0) = zpossible(i,0)*log(psi)+(1-zpossible(i,0))*log(1-psi);
      }
      //Get likelihood for all possible z histories. must update if accepted
      for(int l=1; l<t; l++){
        for(int i2=0; i2<nzpossible; i2++){
          Ezpossible(i2,l-1)=zpossible(i2,l-1)*phiuse(l-1) + apossible(i2,l-1)*gammaprime(l-1);
          llzpossible(i2,l)=zpossible(i2,l)*log(Ezpossible(i2,l-1))+(1-zpossible(i2,l))*log(1-Ezpossible(i2,l-1));
          if(llzpossible(i2,l)!=llzpossible(i2,l)){//fix NaNs
            llzpossible(i2,l)=0;
          }
        }
      }
      for(int i=0; i<M; i++){
        if(upz3){
          //new z stuff
          sumpropto=0;
          for(int i2=0; i2<nzpossible; i2++){
            propto1(i2)=0;
            if(cancel(i,i2)==1){
              for(int l=0; l<t; l++){
                propto1(i2)+=exp(llzpossible(i2,l));
              }
              sumpropto+=propto1(i2);
            }
          }
          for(int i2=0; i2<nzpossible; i2++){
            propto(i2)=propto1(i2)/sumpropto;
            propto2(i2)=propto1(i2)/sumpropto;
          }
          //sample screws up ordering so feeding in second copy not used later
          zchoose=Rcpp::RcppArmadillo::sample(choose,1,FALSE,propto2);
          sumz=0; //using here to see if prop z is same as curr z
          for(int l=0; l<t; l++){
            zprop(l)=zpossible(zchoose(0),l);
            if(zprop(l)==z(i,l)){
              sumz+=1;
            }
          }
          if(sumz!=t){
            propprob=propto(zchoose(0));
            //old z stuff
            for(int i2=0; i2<nzpossible; i2++){//find which zpossible matches current z
              sumz=0;
              for(int l=0; l<t; l++){
                if(zpossible(i2,l)==z(i,l)){
                  sumz+=1;
                }
              }
              if(sumz==t){
                currz=i2;
              }
            }
            backprob=propto(currz);
            //Because a and z changes, must update gamma.prime and Ez
            //Don't need to update all years every time, but not figuring that out for now
            for(int l=0; l<t; l++){
              aprop(l)=apossible(zchoose(0),l);
            }

            //Update z and a cand, Ntmp
            for(int l=0; l<t; l++){
              for(int i2=0; i2<M; i2++){
                if(z(i2,l)==1){
                  zcand(i2,l)=1;
                }else{
                  zcand(i2,l)=0;
                }
                if(a(i2,l)==1){
                  acand(i2,l)=1;
                }else{
                  acand(i2,l)=0;
                }
              }
            }
            for(int l=0; l<t; l++){
              zcand(i,l)=zprop(l);
              acand(i,l)=aprop(l);
            }
            for(int l=0; l<t; l++){
              Ntmp(l)=0;
              for(int i2=0; i2<M; i2++){
                Ntmp(l)+=zcand(i2,l);
              }
            }
            //ll.z[i,1]
            llzcandsum=0;
            llzsum=0;
            ll_z_cand(i,0) = zcand(i,0)*log(psi)+(1-zcand(i,0))*log(1-psi);
            llzcandsum+=ll_z_cand(i,0);
            llzsum+=ll_z(i,0);
            //Calculate gamma.prime, Ez, and ll.z[,2+] candidates
            for(int l=0; l<(t-1); l++){
              suma=0;
              for(int i2=0; i2<M; i2++){
                suma+=acand(i2,l);
              }
              gammaprimecand(l)=(Ntmp(l)*gammause(l)) / suma;
              if(gammaprimecand(l) < 1){ //Is this a valid probability?
                for(int i2=0; i2<M; i2++){
                  //Add on contributions to z ll from l+1 for cand and curr
                  Ezcand(i2,l) = zcand(i2,l)*phiuse(l) + acand(i2,l)*gammaprimecand(l);
                  ll_z_cand(i2,l+1) = zcand(i2,l+1)*log(Ezcand(i2,l))+(1-zcand(i2,l+1))*log(1-Ezcand(i2,l));
                  llzsum+=ll_z(i2,l+1);//prior.z
                  if(ll_z_cand(i2,l+1)==ll_z_cand(i2,l+1)){
                    llzcandsum+=ll_z_cand(i2,l+1);//prior.z.cand
                  }else{
                    ll_z_cand(i2,l+1)=0; //Turn NaN to 0
                  }
                }
              }
            }
            //update ll.y
            llycandsum=0;
            llysum=0;
            for(int l=0; l<t; l++){
              for(int j=0; j<Xidx(l); j++){
                if(obstype==1){
                  ll_y_cand(i,j,l)=zcand(i,l)*(y(i,j,l)*log(pd(i,j,l))+(tf(j,l)-y(i,j,l))*log(1-pd(i,j,l)));
                }else{
                  ll_y_cand(i,j,l)=zcand(i,l)*(y(i,j,l)*log(tf(j,l)*lamd(i,j,l))-tf(j,l)*lamd(i,j,l));
                }
                if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                  llycandsum+=ll_y_cand(i,j,l);
                }
                if(ll_y_curr(i,j,l)==ll_y_curr(i,j,l)){
                  llysum+=ll_y_curr(i,j,l);
                }
              }
            }
            //MH step
            rand=Rcpp::runif(1);
            if(rand(0)<exp((llycandsum+llzcandsum)-(llysum+llzsum))*(backprob/propprob)){
              for(int l=0; l<t; l++){
                z(i,l)=zprop(l);
                a(i,l)=aprop(l);
                N(l)=Ntmp(l);
                ll_z(i,0)=ll_z_cand(i,0);
                for(int j=0; j<Xidx(l); j++){
                  ll_y_curr(i,j,l)=ll_y_cand(i,j,l);
                }
              }
              for(int l=0; l<(t-1); l++){
                gammaprime(l)=gammaprimecand(l);
                for(int i2=0; i2<M; i2++){
                  Ez(i2,l)=Ezcand(i2,l);
                  ll_z(i2,l+1)=ll_z_cand(i2,l+1);
                }
              }
              //Update likelihood for all possible z histories
              for(int l=1; l<t; l++){
                for(int i2=0; i2<nzpossible; i2++){
                  Ezpossible(i2,l-1)=zpossible(i2,l-1)*phiuse(l-1) + apossible(i2,l-1)*gammaprime(l-1);
                  llzpossible(i2,l)=zpossible(i2,l)*log(Ezpossible(i2,l-1))+(1-zpossible(i2,l))*log(1-Ezpossible(i2,l-1));
                  if(llzpossible(i2,l)!=llzpossible(i2,l)){//fix NaNs
                    llzpossible(i2,l)=0;
                  }
                }
              }
            }
          }
        }
      }
    }
    //Update psi
    rand=Rcpp::rbeta(1, 1+N(0), 1+M-N(0));
    psi=rand(0);
    for(int i=0; i<M; i++){
      ll_z(i,0)= z(i,0)*log(psi)+(1-z(i,0))*log(1-psi);
    }
    //Update phi
    if(phi.size()==(t-1)){//if time-specific survival
      for(int l=1; l<t; l++){
        survive=0;
        dead=0;
        for(int i=0; i<M; i++){
          survive+=(z(i,l-1)==1)&(z(i,l)==1);
          dead+=(z(i,l-1)==1)&(z(i,l)==0);
        }
        rand=Rcpp::rbeta(1, 1+survive, 1+dead);
        phi(l-1)=rand(0);
        phiuse(l-1)=rand(0);
      }
    }else{
      survive=0;
      dead=0;
      for(int l=1; l<t; l++){
        for(int i=0; i<M; i++){
          survive+=(z(i,l-1)==1)&(z(i,l)==1);
          dead+=(z(i,l-1)==1)&(z(i,l)==0);
        }
      }
      rand=Rcpp::rbeta(1, 1+survive, 1+dead);
      phi(0)=rand(0);
      for(int l=1; l<t; l++){
        phiuse(l-1)=rand(0);
      }
    }
    //  Update gamma
    // NOTE: Must update ll.z, Ez, etc...
    if(gamma.size()==1){
      rand=Rcpp::rnorm(1, gamma(0), propgamma(0));
      gamma_cand=rand(0);
      gamma_cand_ok=TRUE;
      llzsum=0;
      for(int l=1; l<t; l++){
        suma=0;
        for(int i=0; i<M; i++){
          suma+=a(i,l-1);
        }
        gammaprimecand(l-1)=(N(l-1)*gamma_cand) / suma;
        if(gammaprimecand(l-1) > 1){  //Note don't break loop b/c ll.z needs updating because phi changed
          gamma_cand_ok=FALSE;
        }
        for(int i=0; i<M; i++){
          Ez(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gammaprime(l-1);
          ll_z(i,l)= z(i,l)*log(Ez(i,l-1))+(1-z(i,l))*log(1-Ez(i,l-1));
          if(ll_z(i,l)!=ll_z(i,l)){//Turn NaNs to 0s
            ll_z(i,l)=0;
          }
          llzsum+=ll_z(i,l);
        }
      }
      if((gamma_cand>0)&(gamma_cand_ok)){
        llzcandsum=0;
        for(int l=1; l<t; l++){
          for(int i=0; i<M; i++){
            Ezcand(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gammaprimecand(l-1);
            ll_z_cand(i,l)= z(i,l)*log(Ezcand(i,l-1))+(1-z(i,l))*log(1-Ezcand(i,l-1));
            if(ll_z_cand(i,l)!=ll_z_cand(i,l)){//Turn NaNs to 0s
              ll_z_cand(i,l)=0;
            }
            llzcandsum+=ll_z_cand(i,l);
          }
        }
        rand=Rcpp::runif(1);
        if(rand(0) < exp(llzcandsum - llzsum)){
          gamma(0)=gamma_cand;
          for(int l=1; l<t; l++){
            gammaprime(l-1)=gammaprimecand(l-1);
            gammause(l-1)=gamma_cand;
            for(int i=0; i<M; i++){
              Ez(i,l-1)=Ezcand(i,l-1);
              ll_z(i,l)=ll_z_cand(i,l);
            }
          }
        }
      }
    }else{
      for(int l=1; l<t; l++){
        rand=Rcpp::rnorm(1, gamma(l-1), propgamma(l-1));
        gamma_cand=rand(0);
        gamma_cand_ok=TRUE;
        suma=0;
        for(int i=0; i<M; i++){
          suma+=a(i,l-1);
        }
        gammaprimecand(l-1)=(N(l-1)*gamma_cand) / suma;
        if(gammaprimecand(l-1) > 1){  //Note don't break loop b/c ll.z needs updating because phi changed
          gamma_cand_ok=FALSE;
        }
        llzsum=0;
        for(int i=0; i<M; i++){
          Ez(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gammaprime(l-1);
          ll_z(i,l)= z(i,l)*log(Ez(i,l-1))+(1-z(i,l))*log(1-Ez(i,l-1));
          if(ll_z(i,l)!=ll_z(i,l)){//Turn NaNs to 0s
            ll_z(i,l)=0;
          }
          llzsum+=ll_z(i,l);
        }
        if((gamma_cand>0)&(gamma_cand_ok)){
          llzcandsum=0;
          for(int i=0; i<M; i++){
            Ezcand(i,l-1)=z(i,l-1)*phiuse(l-1) + a(i,l-1)*gammaprimecand(l-1);
            ll_z_cand(i,l)= z(i,l)*log(Ezcand(i,l-1))+(1-z(i,l))*log(1-Ezcand(i,l-1));
            if(ll_z_cand(i,l)!=ll_z_cand(i,l)){//Turn NaNs to 0s
              ll_z_cand(i,l)=0;
            }
            llzcandsum+=ll_z_cand(i,l);
          }
          rand=Rcpp::runif(1);
          if(rand(0) < exp(llzcandsum - llzsum)){
            gamma(l-1)=gamma_cand;
            gammaprime(l-1)=gammaprimecand(l-1);
            gammause(l-1)=gamma_cand;
            for(int i=0; i<M; i++){
              Ez(i,l-1)=Ezcand(i,l-1);
              ll_z(i,l)=ll_z_cand(i,l);
            }
          }
        }
      }
    }
    for(int i=0; i<M; i++) {
      warncount+=warn(i);
    }
    //// Now we have to update the activity centers//////////////////
    if((ACtype==2)|(ACtype==5)){//metamu
      // Update within year ACs
      for(int i=0; i<M; i++) {
        for(int l=0; l<t; l++) {
          ScandX=Rcpp::rnorm(1,s2(i,l,0),props2x);
          ScandY=Rcpp::rnorm(1,s2(i,l,1),props2y);
          if(usedSS){
            double mindist=100000;
            int idxdist=0;
            for(int i2=0; i2<NdSS;i2++){
              dists(i2)=pow(pow(ScandX(0)-dSS(i2,0),2)+pow(ScandY(0)-dSS(i2,1),2),0.5);
              if(dists(i2)<mindist){
                mindist=dists(i2);
                idxdist=i2;
              }
            }
            ScandX(0)=dSS(idxdist,0);
            ScandY(0)=dSS(idxdist,1);
          }
          if(ACtype==2){
            if(useverts==FALSE){
              inbox=(ScandX<xlim(1)) & (ScandX>xlim(0)) & (ScandY<ylim(1)) & (ScandY>ylim(0));
            }else{
              // inbox=inoutCppOpen(ScandX,ScandY,vertices);
              inbox(0)=FALSE;
              for(int p=0; p<polys; p++){
                inbox2(0)=inoutCppOpen(ScandX,ScandY,vertices[p]);
                if(inbox2(0)){
                  inbox(0)=TRUE;
                }
              }
            }
          }else{
            inbox(0)=TRUE;
          }
          if(inbox(0)){
            //sum ll across j for each i and l
            llysum=0;
            llycandsum=0;
            lls2sum=0;
            lls2candsum=0;
            for(int j=0; j<Xidx(l); j++){
              dtmp(j,l)=pow( pow(ScandX(0) - Xcpp(l,j,0), 2.0) + pow(ScandY(0)-Xcpp(l,j,1), 2.0), 0.5 );
              lamdcand(i,j,l)=lam0use(l)*exp(-dtmp(j,l)*dtmp(j,l)/(2*sigmause(l)*sigmause(l)));
              pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
              if(obstype==1){
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(tf(j,l)-y(i,j,l))*log(1-pdcand(i,j,l)));
              }else{
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(tf(j,l)*lamdcand(i,j,l))-tf(j,l)*lamdcand(i,j,l));
              }
              if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                llycandsum+=ll_y_cand(i,j,l);
              }
              if(ll_y_curr(i,j,l)==ll_y_curr(i,j,l)){
                llysum+=ll_y_curr(i,j,l);
              }
            }
            ll_s2_cand(i,l)=-log(pow(sigma_t(0),2.0))-(1/(2*pow(sigma_t(0),2.0)))*(pow(ScandX(0)-s1(i,0),2.0)+pow(ScandY(0)-s1(i,1),2.0));
            rand=Rcpp::runif(1);
            if(rand(0)<exp((llycandsum+ll_s2_cand(i,l))-(llysum+ll_s2(i,l)))){
              s2(i,l,0)=ScandX(0);
              s2(i,l,1)=ScandY(0);
              ll_s2(i,l)=ll_s2_cand(i,l);
              for(int j=0; j<Xidx(l); j++){
                D(i,j,l) = dtmp(j,l);
                lamd(i,j,l) = lamdcand(i,j,l);
                if(obstype==1){
                  pd(i,j,l) = pdcand(i,j,l);
                }
                ll_y_curr(i,j,l) = ll_y_cand(i,j,l);
              }
            }
          }

        }
      }
      // Update meta mus
      for(int i=0; i<M; i++) {
        ScandX=Rcpp::rnorm(1,s1(i,0),props1x);
        ScandY=Rcpp::rnorm(1,s1(i,1),props1y);
        if(usedSS){
          double mindist=100000;
          int idxdist=0;
          for(int i2=0; i2<NdSS;i2++){
            dists(i2)=pow(pow(ScandX(0)-dSS(i2,0),2)+pow(ScandY(0)-dSS(i2,1),2),0.5);
            if(dists(i2)<mindist){
              mindist=dists(i2);
              idxdist=i2;
            }
          }
          ScandX(0)=dSS(idxdist,0);
          ScandY(0)=dSS(idxdist,1);
        }
        if(useverts==FALSE){
          inbox=(ScandX<xlim(1)) & (ScandX>xlim(0)) & (ScandY<ylim(1)) & (ScandY>ylim(0));
        }else{
          // inbox=inoutCppOpen(ScandX,ScandY,vertices);
          inbox(0)=FALSE;
          for(int p=0; p<polys; p++){
            inbox2(0)=inoutCppOpen(ScandX,ScandY,vertices[p]);
            if(inbox2(0)){
              inbox(0)=TRUE;
            }
          }
        }
        if(inbox(0)){
          lls2sum=0;
          lls2candsum=0;
          for(int l=0; l<t; l++) {
            ll_s2_cand(i,l)=-log(pow(sigma_t(0),2.0))-(1/(2*pow(sigma_t(0),2.0)))*(pow(s2(i,l,0)-ScandX(0),2.0)+pow(s2(i,l,1)-ScandY(0),2.0));
            lls2sum+=ll_s2(i,l);
            lls2candsum+=ll_s2_cand(i,l);
          }
          rand=Rcpp::runif(1);
          if (rand(0) < exp(lls2candsum - lls2sum)) {
            s1(i,0)=ScandX(0);
            s1(i,1)=ScandY(0);
            for(int l=0; l<t; l++) {
              ll_s2(i,l)=ll_s2_cand(i,l);
            }
          }
        }
      }
      // Update sigma_t
      sigma_t_cand=rnorm(1,sigma_t(0),propsigma_t);
      if(sigma_t_cand(0) > 0){
        lls2sum=0;
        lls2candsum=0;
        for(int i=0; i<M; i++) {
          for(int l=0; l<t; l++) {
            ll_s2_cand(i,l)=-log(pow(sigma_t_cand(0),2.0))-(1/(2*pow(sigma_t_cand(0),2.0)))*(pow(s2(i,l,0)-s1(i,0),2.0)+pow(s2(i,l,1)-s1(i,1),2.0));
            lls2sum+=ll_s2(i,l);
            lls2candsum+=ll_s2_cand(i,l);
          }
        }
        rand=Rcpp::runif(1);
        if (rand(0) < exp(lls2candsum - lls2sum)) {
          sigma_t(0)=sigma_t_cand(0);
          for(int i=0; i<M; i++) {
            for(int l=0; l<t; l++) {
              ll_s2(i,l)=ll_s2_cand(i,l);
            }
          }
        }
      }
    }else if(ACtype==1){//Stationary ACs
      //Update Activity Centers
      for(int i=0; i<M; i++) {
        ScandX=Rcpp::rnorm(1,s1(i,0),props2x);
        ScandY=Rcpp::rnorm(1,s1(i,1),props2y);
        if(usedSS){
          double mindist=100000;
          int idxdist=0;
          for(int i2=0; i2<NdSS;i2++){
            dists(i2)=pow(pow(ScandX(0)-dSS(i2,0),2)+pow(ScandY(0)-dSS(i2,1),2),0.5);
            if(dists(i2)<mindist){
              mindist=dists(i2);
              idxdist=i2;
            }
          }
          ScandX(0)=dSS(idxdist,0);
          ScandY(0)=dSS(idxdist,1);
        }
        if(useverts==FALSE){
          inbox=(ScandX<xlim(1)) & (ScandX>xlim(0)) & (ScandY<ylim(1)) & (ScandY>ylim(0));
        }else{
          // inbox=inoutCppOpen(ScandX,ScandY,vertices);
          inbox(0)=FALSE;
          for(int p=0; p<polys; p++){
            inbox2(0)=inoutCppOpen(ScandX,ScandY,vertices[p]);
            if(inbox2(0)){
              inbox(0)=TRUE;
            }
          }
        }
        if(inbox(0)){
          //sum ll across j for each i and l
          llysum=0;
          llycandsum=0;
          for(int l=0; l<t; l++){
            for(int j=0; j<Xidx(l); j++){
              dtmp(j,l)=pow( pow(ScandX(0) - Xcpp(l,j,0), 2.0) + pow(ScandY(0)-Xcpp(l,j,1), 2.0), 0.5 );
              lamdcand(i,j,l)=lam0use(l)*exp(-dtmp(j,l)*dtmp(j,l)/(2*sigmause(l)*sigmause(l)));
              if(obstype==1){
                pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(tf(j,l)-y(i,j,l))*log(1-pdcand(i,j,l)));
              }else{
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(tf(j,l)*lamdcand(i,j,l))-tf(j,l)*lamdcand(i,j,l));
              }
              if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                llycandsum+=ll_y_cand(i,j,l);
              }
              if(ll_y_curr(i,j,l)==ll_y_curr(i,j,l)){
                llysum+=ll_y_curr(i,j,l);
              }
            }
          }
          rand=Rcpp::runif(1);
          if((rand(0)<exp(llycandsum-llysum))){
            s1(i,0)=ScandX(0);
            s1(i,1)=ScandY(0);
            for(int l=0; l<t; l++){
              s2(i,l,0)=ScandX(0);
              s2(i,l,1)=ScandY(0);
              for(int j=0; j<J; j++){
                D(i,j,l) = dtmp(j,l);
                lamd(i,j,l) = lamdcand(i,j,l);
                if(obstype==1){
                  pd(i,j,l) = pdcand(i,j,l);
                }
                ll_y_curr(i,j,l) = ll_y_cand(i,j,l);
              }
            }
          }
        }
      }
    }else if(ACtype==3){//markov ACs
      // Update within year ACs
      for(int i=0; i<M; i++) {
        for(int l=0; l<t; l++) {
          ScandX=Rcpp::rnorm(1,s2(i,l,0),props2x);
          ScandY=Rcpp::rnorm(1,s2(i,l,1),props2y);
          if(usedSS){
            double mindist=100000;
            int idxdist=0;
            for(int i2=0; i2<NdSS;i2++){
              dists(i2)=pow(pow(ScandX(0)-dSS(i2,0),2)+pow(ScandY(0)-dSS(i2,1),2),0.5);
              if(dists(i2)<mindist){
                mindist=dists(i2);
                idxdist=i2;
              }
            }
            ScandX(0)=dSS(idxdist,0);
            ScandY(0)=dSS(idxdist,1);
          }
          if(useverts==FALSE){
            inbox=(ScandX<xlim(1)) & (ScandX>xlim(0)) & (ScandY<ylim(1)) & (ScandY>ylim(0));
          }else{
            // inbox=inoutCppOpen(ScandX,ScandY,vertices);
            inbox(0)=FALSE;
            for(int p=0; p<polys; p++){
              inbox2(0)=inoutCppOpen(ScandX,ScandY,vertices[p]);
              if(inbox2(0)){
                inbox(0)=TRUE;
              }
            }
          }
          if(inbox(0)){
            //sum ll across j for each i and l
            llysum=0;
            llycandsum=0;
            for(int j=0; j<Xidx(l); j++){
              dtmp(j,l)=pow( pow(ScandX(0) - Xcpp(l,j,0), 2.0) + pow(ScandY(0)-Xcpp(l,j,1), 2.0), 0.5 );
              lamdcand(i,j,l)=lam0use(l)*exp(-dtmp(j,l)*dtmp(j,l)/(2*sigmause(l)*sigmause(l)));
              if(obstype==1){
                pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(tf(j,l)-y(i,j,l))*log(1-pdcand(i,j,l)));
              }else{
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(tf(j,l)*lamdcand(i,j,l))-tf(j,l)*lamdcand(i,j,l));
              }
              if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                llycandsum+=ll_y_cand(i,j,l);
              }
              if(ll_y_curr(i,j,l)==ll_y_curr(i,j,l)){
                llysum+=ll_y_curr(i,j,l);
              }
            }
            lls2sum=0;
            lls2candsum=0;
            if(l==0){ //first occasion
              //step from 1 to 2
              ll_s2_cand(i,0)=-log(pow(sigma_t(0),2.0))-(1/(2*pow(sigma_t(0),2.0)))*(pow(s2(i,1,0)-ScandX(0),2.0)+pow(s2(i,1,1)-ScandY(0),2.0));
              lls2sum+=ll_s2(i,0);
              lls2candsum+=ll_s2_cand(i,0);
            }else if((l>0)&(l<(t-1))){//middle occasions
              //step from l-1 to l
              ll_s2_cand(i,l-1)=-log(pow(sigma_t(0),2.0))-(1/(2*pow(sigma_t(0),2.0)))*(pow(ScandX(0)-s2(i,l-1,0),2.0)+pow(ScandY(0)-s2(i,l-1,1),2.0));
              //step from l to l+1
              ll_s2_cand(i,l)=-log(pow(sigma_t(0),2.0))-(1/(2*pow(sigma_t(0),2.0)))*(pow(s2(i,l+1,0)-ScandX(0),2.0)+pow(s2(i,l+1,1)-ScandY(0),2.0));
              lls2sum+=ll_s2(i,l-1);
              lls2sum+=ll_s2(i,l);
              lls2candsum+=ll_s2_cand(i,l-1);
              lls2candsum+=ll_s2_cand(i,l);
            }else{//final occasion
              //step from t-1 to t
              ll_s2_cand(i,l-1)=-log(pow(sigma_t(0),2.0))-(1/(2*pow(sigma_t(0),2.0)))*(pow(ScandX(0)-s2(i,l-1,0),2.0)+pow(ScandY(0)-s2(i,l-1,1),2.0));
              lls2sum+=ll_s2(i,l-1);
              lls2candsum+=ll_s2_cand(i,l-1);
            }
            rand=Rcpp::runif(1);
            if(rand(0)<exp((llycandsum+lls2candsum)-(llysum+lls2sum))){
              s2(i,l,0)=ScandX(0);
              s2(i,l,1)=ScandY(0);
              if(l==0){
                ll_s2(i,0)=ll_s2_cand(i,0);
              }else if((l>0)&(l<(t-1))){
                ll_s2(i,l)=ll_s2_cand(i,l);
                ll_s2(i,l-1)=ll_s2_cand(i,l-1);
              }else{
                ll_s2(i,l-1)=ll_s2_cand(i,l-1);
              }
              for(int j=0; j<Xidx(l); j++){
                D(i,j,l) = dtmp(j,l);
                lamd(i,j,l) = lamdcand(i,j,l);
                if(obstype==1){
                  pd(i,j,l) = pdcand(i,j,l);
                }
                ll_y_curr(i,j,l) = ll_y_cand(i,j,l);
              }
            }
          }
        }
      }
      // Update sigma_t
      sigma_t_cand=rnorm(1,sigma_t(0),propsigma_t);
      if(sigma_t_cand(0) > 0){
        lls2sum=0;
        lls2candsum=0;
        for(int i=0; i<M; i++) {
          for(int l=1; l<t; l++) {
            ll_s2_cand(i,l-1)=-log(pow(sigma_t_cand(0),2.0))-(1/(2*pow(sigma_t_cand(0),2.0)))*(pow(s2(i,l,0)-s2(i,l-1,0),2.0)+pow(s2(i,l,1)-s2(i,l-1,1),2.0));
            lls2sum+=ll_s2(i,l-1);
            lls2candsum+=ll_s2_cand(i,l-1);
          }
        }
        rand=Rcpp::runif(1);
        if (rand(0) < exp(lls2candsum - lls2sum)) {
          sigma_t(0)=sigma_t_cand(0);
          for(int i=0; i<M; i++) {
            for(int l=1; l<t; l++) {
              ll_s2(i,l-1)=ll_s2_cand(i,l-1);
            }
          }
        }
      }
    }else{//independent ACs
      for(int l=0; l<t; l++){
        for(int i=0; i<M; i++) {
          ScandX=Rcpp::rnorm(1,s2(i,l,0),props2x);
          ScandY=Rcpp::rnorm(1,s2(i,l,1),props2y);
          if(usedSS){
            double mindist=100000;
            int idxdist=0;
            for(int i2=0; i2<NdSS;i2++){
              dists(i2)=pow(pow(ScandX(0)-dSS(i2,0),2)+pow(ScandY(0)-dSS(i2,1),2),0.5);
              if(dists(i2)<mindist){
                mindist=dists(i2);
                idxdist=i2;
              }
            }
            ScandX(0)=dSS(idxdist,0);
            ScandY(0)=dSS(idxdist,1);
          }
          if(useverts==FALSE){
            inbox=(ScandX<xlim(1)) & (ScandX>xlim(0)) & (ScandY<ylim(1)) & (ScandY>ylim(0));
          }else{
            // inbox=inoutCppOpen(ScandX,ScandY,vertices);
            inbox(0)=FALSE;
            for(int p=0; p<polys; p++){
              inbox2(0)=inoutCppOpen(ScandX,ScandY,vertices[p]);
              if(inbox2(0)){
                inbox(0)=TRUE;
              }
            }
          }
          if(inbox(0)){
            //sum ll across j for each i and l
            llysum=0;
            llycandsum=0;
            for(int j=0; j<Xidx(l); j++){
              dtmp(j,l)=pow( pow(ScandX(0) - Xcpp(l,j,0), 2.0) + pow(ScandY(0)-Xcpp(l,j,1), 2.0), 0.5 );
              lamdcand(i,j,l)=lam0use(l)*exp(-dtmp(j,l)*dtmp(j,l)/(2*sigmause(l)*sigmause(l)));
              if(obstype==1){
                pdcand(i,j,l)=1-exp(-lamdcand(i,j,l));
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(pdcand(i,j,l))+(tf(j,l)-y(i,j,l))*log(1-pdcand(i,j,l)));
              }else{
                ll_y_cand(i,j,l)=z(i,l)*(y(i,j,l)*log(tf(j,l)*lamdcand(i,j,l))-tf(j,l)*lamdcand(i,j,l));
              }
              if(ll_y_cand(i,j,l)==ll_y_cand(i,j,l)){
                llycandsum+=ll_y_cand(i,j,l);
              }
              if(ll_y_curr(i,j,l)==ll_y_curr(i,j,l)){
                llysum+=ll_y_curr(i,j,l);
              }
            }
            rand=Rcpp::runif(1);
            if((rand(0)<exp(llycandsum-llysum))){
              s2(i,l,0)=ScandX(0);
              s2(i,l,1)=ScandY(0);
              for(int j=0; j<J; j++){
                D(i,j,l) = dtmp(j,l);
                lamd(i,j,l) = lamdcand(i,j,l);
                if(obstype==1){
                  pd(i,j,l) = pdcand(i,j,l);
                }
                ll_y_curr(i,j,l) = ll_y_cand(i,j,l);
              }
            }
          }
        }
      }
    }



    //Record output ll_y_curr.subcube(i,0,0,i,maxJ-1,0) s2yout(nstore,M,t)
    if(((iter+1)>nburn)&((iter+1) % nthin==0)){
      for(int i=0; i<M; i++){
        s1xout(iteridx,i)= s1(i,0);
        s1yout(iteridx,i)= s1(i,1);
        for(int l=0; l<t; l++){
          s2xout(iteridx,i,l)=s2(i,l,0);
          s2yout(iteridx,i,l)=s2(i,l,1);
          zout(iteridx,i,l)= z(i,l);
        }
      }
      idx=0;
      //fill in lam0
      for(int l=0; l<each(0); l++){
        out(iteridx,idx)=lam0(l);
        idx=idx+1;
      }
      //fill in sigma
      for(int l=0; l<each(1); l++){
        out(iteridx,idx)=sigma(l);
        idx=idx+1;
      }
      //fill in gamma
      for(int l=0; l<each(2); l++){
        out(iteridx,idx)=gamma(l);
        idx=idx+1;
      }
      //fill in pki
      for(int l=0; l<each(3); l++){
        out(iteridx,idx)=phi(l);
        idx=idx+1;
      }
      //fill in N
      for(int l=0; l<t; l++){
        out(iteridx,idx)=N(l);
        idx=idx+1;
      }
      if((ACtype==2)|(ACtype==3)|(ACtype==5)){
        out(iteridx,idx)=sigma_t(0);
      }
      iteridx=iteridx+1;
    }
  }
  List to_return(10);
  to_return[0] = out;
  to_return[1] = s1xout;
  to_return[2] = s1yout;
  to_return[3] = s2xout;
  to_return[4] = s2yout;
  to_return[5] = zout;
  to_return[6] = warncount;
  to_return[7] = a;
  to_return[8] = ll_s2;
  to_return[9] = ll_s2_cand;
  return to_return;
}
