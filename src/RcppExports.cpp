// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// multiSigma
NumericVector multiSigma(NumericMatrix signal, int deg = 3);
RcppExport SEXP mwaved_multiSigma(SEXP signalSEXP, SEXP degSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type signal(signalSEXP );
        Rcpp::traits::input_parameter< int >::type deg(degSEXP );
        NumericVector __result = multiSigma(signal, deg);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// waveletThresh
List waveletThresh(List betaList, NumericVector thr, String shrinkType = "hard");
RcppExport SEXP mwaved_waveletThresh(SEXP betaListSEXP, SEXP thrSEXP, SEXP shrinkTypeSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type betaList(betaListSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type thr(thrSEXP );
        Rcpp::traits::input_parameter< String >::type shrinkType(shrinkTypeSEXP );
        List __result = waveletThresh(betaList, thr, shrinkType);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// multiProj
NumericVector multiProj(NumericVector beta, int j0 = 3, int j1 = NA_INTEGER, int deg = 3);
RcppExport SEXP mwaved_multiProj(SEXP betaSEXP, SEXP j0SEXP, SEXP j1SEXP, SEXP degSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< int >::type j0(j0SEXP );
        Rcpp::traits::input_parameter< int >::type j1(j1SEXP );
        Rcpp::traits::input_parameter< int >::type deg(degSEXP );
        NumericVector __result = multiProj(beta, j0, j1, deg);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// multiThresh
NumericVector multiThresh(NumericMatrix signal, NumericMatrix G, NumericVector alpha = NumericVector::create(), String blur = "direct", int j0 = 3, int j1 = NA_INTEGER, double eta = NA_REAL, int deg = 3);
RcppExport SEXP mwaved_multiThresh(SEXP signalSEXP, SEXP GSEXP, SEXP alphaSEXP, SEXP blurSEXP, SEXP j0SEXP, SEXP j1SEXP, SEXP etaSEXP, SEXP degSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type signal(signalSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type G(GSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< String >::type blur(blurSEXP );
        Rcpp::traits::input_parameter< int >::type j0(j0SEXP );
        Rcpp::traits::input_parameter< int >::type j1(j1SEXP );
        Rcpp::traits::input_parameter< double >::type eta(etaSEXP );
        Rcpp::traits::input_parameter< int >::type deg(degSEXP );
        NumericVector __result = multiThresh(signal, G, alpha, blur, j0, j1, eta, deg);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// multiEstimate
NumericVector multiEstimate(NumericMatrix signal, NumericMatrix G, NumericVector alpha = NumericVector::create(), String blur = "direct", NumericVector sigma = NumericVector::create(), int j0 = 3, int j1 = NA_INTEGER, double eta = NA_REAL, NumericVector thresh = NumericVector::create(), String shrinkage = "Hard", int deg = 3);
RcppExport SEXP mwaved_multiEstimate(SEXP signalSEXP, SEXP GSEXP, SEXP alphaSEXP, SEXP blurSEXP, SEXP sigmaSEXP, SEXP j0SEXP, SEXP j1SEXP, SEXP etaSEXP, SEXP threshSEXP, SEXP shrinkageSEXP, SEXP degSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type signal(signalSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type G(GSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< String >::type blur(blurSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< int >::type j0(j0SEXP );
        Rcpp::traits::input_parameter< int >::type j1(j1SEXP );
        Rcpp::traits::input_parameter< double >::type eta(etaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type thresh(threshSEXP );
        Rcpp::traits::input_parameter< String >::type shrinkage(shrinkageSEXP );
        Rcpp::traits::input_parameter< int >::type deg(degSEXP );
        NumericVector __result = multiEstimate(signal, G, alpha, blur, sigma, j0, j1, eta, thresh, shrinkage, deg);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// multiCoef
List multiCoef(NumericMatrix signal, NumericMatrix G, NumericVector alpha = NumericVector::create(), String blur = "direct", int j0 = 3, int j1 = NA_INTEGER, NumericVector thresh = NumericVector::create(), double eta = NA_REAL, int deg = 3);
RcppExport SEXP mwaved_multiCoef(SEXP signalSEXP, SEXP GSEXP, SEXP alphaSEXP, SEXP blurSEXP, SEXP j0SEXP, SEXP j1SEXP, SEXP threshSEXP, SEXP etaSEXP, SEXP degSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type signal(signalSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type G(GSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< String >::type blur(blurSEXP );
        Rcpp::traits::input_parameter< int >::type j0(j0SEXP );
        Rcpp::traits::input_parameter< int >::type j1(j1SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type thresh(threshSEXP );
        Rcpp::traits::input_parameter< double >::type eta(etaSEXP );
        Rcpp::traits::input_parameter< int >::type deg(degSEXP );
        List __result = multiCoef(signal, G, alpha, blur, j0, j1, thresh, eta, deg);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// mWaveD
List mWaveD(NumericMatrix signal, NumericMatrix G, NumericVector alpha = NumericVector::create(), int j0 = 3, int j1 = NA_INTEGER, String blur = "direct", NumericVector thresh = NumericVector::create(), double eta = NA_REAL, String shrinkage = "Hard", int deg = 3);
RcppExport SEXP mwaved_mWaveD(SEXP signalSEXP, SEXP GSEXP, SEXP alphaSEXP, SEXP j0SEXP, SEXP j1SEXP, SEXP blurSEXP, SEXP threshSEXP, SEXP etaSEXP, SEXP shrinkageSEXP, SEXP degSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type signal(signalSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type G(GSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< int >::type j0(j0SEXP );
        Rcpp::traits::input_parameter< int >::type j1(j1SEXP );
        Rcpp::traits::input_parameter< String >::type blur(blurSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type thresh(threshSEXP );
        Rcpp::traits::input_parameter< double >::type eta(etaSEXP );
        Rcpp::traits::input_parameter< String >::type shrinkage(shrinkageSEXP );
        Rcpp::traits::input_parameter< int >::type deg(degSEXP );
        List __result = mWaveD(signal, G, alpha, j0, j1, blur, thresh, eta, shrinkage, deg);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
