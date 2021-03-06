// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// inoutCppOpen
bool inoutCppOpen(NumericVector sx, NumericVector sy, NumericMatrix vertices);
RcppExport SEXP _OpenPopSCR_inoutCppOpen(SEXP sxSEXP, SEXP sySEXP, SEXP verticesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sx(sxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sy(sySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type vertices(verticesSEXP);
    rcpp_result_gen = Rcpp::wrap(inoutCppOpen(sx, sy, vertices));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_Open
List mcmc_Open(NumericVector lam0in, NumericVector sigmain, NumericVector gammain, NumericVector gammaprime, NumericVector phiin, arma::cube D, arma::cube lamd, arma::cube y, IntegerMatrix z, IntegerMatrix a, NumericMatrix s1, arma::cube s2, int ACtype, bool useverts, List vertices, NumericVector xlim, NumericVector ylim, IntegerMatrix knownmatrix, IntegerVector Xidx, arma::cube Xcpp, IntegerVector K, NumericMatrix Ez, double psi, IntegerVector N, NumericVector proplam0, NumericVector propsig, NumericVector propz, NumericVector propgamma, double props1x, double props1y, double props2x, double props2y, double propsigma_t, NumericVector sigma_tin, int niter, int nburn, int nthin, int npar, IntegerVector each, bool jointZ, IntegerMatrix zpossible, IntegerMatrix apossible, IntegerMatrix cancel, int obstype, IntegerMatrix tf, IntegerMatrix s2cell, IntegerVector s1cell, NumericMatrix dSS, bool usedSS, LogicalVector primary, bool dualACup, int propdualAC, NumericMatrix distances, bool storeLatent);
RcppExport SEXP _OpenPopSCR_mcmc_Open(SEXP lam0inSEXP, SEXP sigmainSEXP, SEXP gammainSEXP, SEXP gammaprimeSEXP, SEXP phiinSEXP, SEXP DSEXP, SEXP lamdSEXP, SEXP ySEXP, SEXP zSEXP, SEXP aSEXP, SEXP s1SEXP, SEXP s2SEXP, SEXP ACtypeSEXP, SEXP usevertsSEXP, SEXP verticesSEXP, SEXP xlimSEXP, SEXP ylimSEXP, SEXP knownmatrixSEXP, SEXP XidxSEXP, SEXP XcppSEXP, SEXP KSEXP, SEXP EzSEXP, SEXP psiSEXP, SEXP NSEXP, SEXP proplam0SEXP, SEXP propsigSEXP, SEXP propzSEXP, SEXP propgammaSEXP, SEXP props1xSEXP, SEXP props1ySEXP, SEXP props2xSEXP, SEXP props2ySEXP, SEXP propsigma_tSEXP, SEXP sigma_tinSEXP, SEXP niterSEXP, SEXP nburnSEXP, SEXP nthinSEXP, SEXP nparSEXP, SEXP eachSEXP, SEXP jointZSEXP, SEXP zpossibleSEXP, SEXP apossibleSEXP, SEXP cancelSEXP, SEXP obstypeSEXP, SEXP tfSEXP, SEXP s2cellSEXP, SEXP s1cellSEXP, SEXP dSSSEXP, SEXP usedSSSEXP, SEXP primarySEXP, SEXP dualACupSEXP, SEXP propdualACSEXP, SEXP distancesSEXP, SEXP storeLatentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lam0in(lam0inSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmain(sigmainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gammain(gammainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gammaprime(gammaprimeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phiin(phiinSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type lamd(lamdSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< int >::type ACtype(ACtypeSEXP);
    Rcpp::traits::input_parameter< bool >::type useverts(usevertsSEXP);
    Rcpp::traits::input_parameter< List >::type vertices(verticesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xlim(xlimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ylim(ylimSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type knownmatrix(knownmatrixSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Xidx(XidxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Xcpp(XcppSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type K(KSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ez(EzSEXP);
    Rcpp::traits::input_parameter< double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type proplam0(proplam0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propsig(propsigSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propz(propzSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propgamma(propgammaSEXP);
    Rcpp::traits::input_parameter< double >::type props1x(props1xSEXP);
    Rcpp::traits::input_parameter< double >::type props1y(props1ySEXP);
    Rcpp::traits::input_parameter< double >::type props2x(props2xSEXP);
    Rcpp::traits::input_parameter< double >::type props2y(props2ySEXP);
    Rcpp::traits::input_parameter< double >::type propsigma_t(propsigma_tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma_tin(sigma_tinSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< int >::type nthin(nthinSEXP);
    Rcpp::traits::input_parameter< int >::type npar(nparSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type each(eachSEXP);
    Rcpp::traits::input_parameter< bool >::type jointZ(jointZSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type zpossible(zpossibleSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type apossible(apossibleSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type cancel(cancelSEXP);
    Rcpp::traits::input_parameter< int >::type obstype(obstypeSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type tf(tfSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type s2cell(s2cellSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type s1cell(s1cellSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dSS(dSSSEXP);
    Rcpp::traits::input_parameter< bool >::type usedSS(usedSSSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type primary(primarySEXP);
    Rcpp::traits::input_parameter< bool >::type dualACup(dualACupSEXP);
    Rcpp::traits::input_parameter< int >::type propdualAC(propdualACSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type distances(distancesSEXP);
    Rcpp::traits::input_parameter< bool >::type storeLatent(storeLatentSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_Open(lam0in, sigmain, gammain, gammaprime, phiin, D, lamd, y, z, a, s1, s2, ACtype, useverts, vertices, xlim, ylim, knownmatrix, Xidx, Xcpp, K, Ez, psi, N, proplam0, propsig, propz, propgamma, props1x, props1y, props2x, props2y, propsigma_t, sigma_tin, niter, nburn, nthin, npar, each, jointZ, zpossible, apossible, cancel, obstype, tf, s2cell, s1cell, dSS, usedSS, primary, dualACup, propdualAC, distances, storeLatent));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_Open_sex
List mcmc_Open_sex(NumericVector lam0in, NumericVector sigmain, NumericVector gammain, NumericVector gammaprimeM, NumericVector gammaprimeF, NumericVector phiin, double psex, arma::cube D, arma::cube lamd, arma::cube y, IntegerMatrix z, IntegerMatrix a, NumericMatrix s1, arma::cube s2, int ACtype, bool useverts, List vertices, NumericVector xlim, NumericVector ylim, IntegerVector sex, IntegerMatrix knownmatrix, IntegerVector Xidx, arma::cube Xcpp, IntegerVector K, NumericMatrix Ez, double psi, IntegerVector N, NumericVector proplam0, NumericVector propsig, NumericVector propz, NumericVector propgamma, double props1x, double props1y, double props2x, double props2y, NumericVector propsigma_t, int propsex, NumericVector sigma_tin, int niter, int nburn, int nthin, int npar, IntegerVector each, bool jointZ, IntegerMatrix zpossible, IntegerMatrix apossible, IntegerMatrix cancel, int obstype, IntegerMatrix tf, NumericMatrix dSS, bool usedSS, LogicalVector sexparms, IntegerVector choosesex, LogicalVector primary, IntegerMatrix s2cell, IntegerVector s1cell, bool dualACup, int propdualAC, NumericMatrix distances, bool storeLatent);
RcppExport SEXP _OpenPopSCR_mcmc_Open_sex(SEXP lam0inSEXP, SEXP sigmainSEXP, SEXP gammainSEXP, SEXP gammaprimeMSEXP, SEXP gammaprimeFSEXP, SEXP phiinSEXP, SEXP psexSEXP, SEXP DSEXP, SEXP lamdSEXP, SEXP ySEXP, SEXP zSEXP, SEXP aSEXP, SEXP s1SEXP, SEXP s2SEXP, SEXP ACtypeSEXP, SEXP usevertsSEXP, SEXP verticesSEXP, SEXP xlimSEXP, SEXP ylimSEXP, SEXP sexSEXP, SEXP knownmatrixSEXP, SEXP XidxSEXP, SEXP XcppSEXP, SEXP KSEXP, SEXP EzSEXP, SEXP psiSEXP, SEXP NSEXP, SEXP proplam0SEXP, SEXP propsigSEXP, SEXP propzSEXP, SEXP propgammaSEXP, SEXP props1xSEXP, SEXP props1ySEXP, SEXP props2xSEXP, SEXP props2ySEXP, SEXP propsigma_tSEXP, SEXP propsexSEXP, SEXP sigma_tinSEXP, SEXP niterSEXP, SEXP nburnSEXP, SEXP nthinSEXP, SEXP nparSEXP, SEXP eachSEXP, SEXP jointZSEXP, SEXP zpossibleSEXP, SEXP apossibleSEXP, SEXP cancelSEXP, SEXP obstypeSEXP, SEXP tfSEXP, SEXP dSSSEXP, SEXP usedSSSEXP, SEXP sexparmsSEXP, SEXP choosesexSEXP, SEXP primarySEXP, SEXP s2cellSEXP, SEXP s1cellSEXP, SEXP dualACupSEXP, SEXP propdualACSEXP, SEXP distancesSEXP, SEXP storeLatentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lam0in(lam0inSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigmain(sigmainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gammain(gammainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gammaprimeM(gammaprimeMSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gammaprimeF(gammaprimeFSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phiin(phiinSEXP);
    Rcpp::traits::input_parameter< double >::type psex(psexSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type lamd(lamdSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< int >::type ACtype(ACtypeSEXP);
    Rcpp::traits::input_parameter< bool >::type useverts(usevertsSEXP);
    Rcpp::traits::input_parameter< List >::type vertices(verticesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xlim(xlimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ylim(ylimSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sex(sexSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type knownmatrix(knownmatrixSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Xidx(XidxSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Xcpp(XcppSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type K(KSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ez(EzSEXP);
    Rcpp::traits::input_parameter< double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type proplam0(proplam0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propsig(propsigSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propz(propzSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propgamma(propgammaSEXP);
    Rcpp::traits::input_parameter< double >::type props1x(props1xSEXP);
    Rcpp::traits::input_parameter< double >::type props1y(props1ySEXP);
    Rcpp::traits::input_parameter< double >::type props2x(props2xSEXP);
    Rcpp::traits::input_parameter< double >::type props2y(props2ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propsigma_t(propsigma_tSEXP);
    Rcpp::traits::input_parameter< int >::type propsex(propsexSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma_tin(sigma_tinSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type nburn(nburnSEXP);
    Rcpp::traits::input_parameter< int >::type nthin(nthinSEXP);
    Rcpp::traits::input_parameter< int >::type npar(nparSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type each(eachSEXP);
    Rcpp::traits::input_parameter< bool >::type jointZ(jointZSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type zpossible(zpossibleSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type apossible(apossibleSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type cancel(cancelSEXP);
    Rcpp::traits::input_parameter< int >::type obstype(obstypeSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type tf(tfSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dSS(dSSSEXP);
    Rcpp::traits::input_parameter< bool >::type usedSS(usedSSSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type sexparms(sexparmsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type choosesex(choosesexSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type primary(primarySEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type s2cell(s2cellSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type s1cell(s1cellSEXP);
    Rcpp::traits::input_parameter< bool >::type dualACup(dualACupSEXP);
    Rcpp::traits::input_parameter< int >::type propdualAC(propdualACSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type distances(distancesSEXP);
    Rcpp::traits::input_parameter< bool >::type storeLatent(storeLatentSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_Open_sex(lam0in, sigmain, gammain, gammaprimeM, gammaprimeF, phiin, psex, D, lamd, y, z, a, s1, s2, ACtype, useverts, vertices, xlim, ylim, sex, knownmatrix, Xidx, Xcpp, K, Ez, psi, N, proplam0, propsig, propz, propgamma, props1x, props1y, props2x, props2y, propsigma_t, propsex, sigma_tin, niter, nburn, nthin, npar, each, jointZ, zpossible, apossible, cancel, obstype, tf, dSS, usedSS, sexparms, choosesex, primary, s2cell, s1cell, dualACup, propdualAC, distances, storeLatent));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_OpenPopSCR_inoutCppOpen", (DL_FUNC) &_OpenPopSCR_inoutCppOpen, 3},
    {"_OpenPopSCR_mcmc_Open", (DL_FUNC) &_OpenPopSCR_mcmc_Open, 54},
    {"_OpenPopSCR_mcmc_Open_sex", (DL_FUNC) &_OpenPopSCR_mcmc_Open_sex, 60},
    {NULL, NULL, 0}
};

RcppExport void R_init_OpenPopSCR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
