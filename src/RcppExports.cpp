// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// RCPPSampleFromRIM
arma::mat RCPPSampleFromRIM(NumericVector numSamplesVec, arma::mat rimNodesMat);
RcppExport SEXP _RIM_RCPPSampleFromRIM(SEXP numSamplesVecSEXP, SEXP rimNodesMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type numSamplesVec(numSamplesVecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type rimNodesMat(rimNodesMatSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPSampleFromRIM(numSamplesVec, rimNodesMat));
    return rcpp_result_gen;
END_RCPP
}
// RCPPAverageDiscMatrix
arma::mat RCPPAverageDiscMatrix(arma::imat samples);
RcppExport SEXP _RIM_RCPPAverageDiscMatrix(SEXP samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat >::type samples(samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPAverageDiscMatrix(samples));
    return rcpp_result_gen;
END_RCPP
}
// RCPPStructByDP
arma::mat RCPPStructByDP(arma::mat aveDiscMatrix, arma::ivec refRanking, bool makeCanonical);
RcppExport SEXP _RIM_RCPPStructByDP(SEXP aveDiscMatrixSEXP, SEXP refRankingSEXP, SEXP makeCanonicalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type aveDiscMatrix(aveDiscMatrixSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type refRanking(refRankingSEXP);
    Rcpp::traits::input_parameter< bool >::type makeCanonical(makeCanonicalSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPStructByDP(aveDiscMatrix, refRanking, makeCanonical));
    return rcpp_result_gen;
END_RCPP
}
// RCPPSASearch
arma::mat RCPPSASearch(arma::mat aveDiscMatrix, arma::ivec refRanking, double inverseTemp, int maxIter, bool makeCanonical, bool verbose);
RcppExport SEXP _RIM_RCPPSASearch(SEXP aveDiscMatrixSEXP, SEXP refRankingSEXP, SEXP inverseTempSEXP, SEXP maxIterSEXP, SEXP makeCanonicalSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type aveDiscMatrix(aveDiscMatrixSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type refRanking(refRankingSEXP);
    Rcpp::traits::input_parameter< double >::type inverseTemp(inverseTempSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< bool >::type makeCanonical(makeCanonicalSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPSASearch(aveDiscMatrix, refRanking, inverseTemp, maxIter, makeCanonical, verbose));
    return rcpp_result_gen;
END_RCPP
}
// RCPPLogProbRIM
double RCPPLogProbRIM(arma::mat rimNodesMat, arma::mat aveDiscMatrix);
RcppExport SEXP _RIM_RCPPLogProbRIM(SEXP rimNodesMatSEXP, SEXP aveDiscMatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type rimNodesMat(rimNodesMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aveDiscMatrix(aveDiscMatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPLogProbRIM(rimNodesMat, aveDiscMatrix));
    return rcpp_result_gen;
END_RCPP
}
// RCPPthetaMLERIM
arma::mat RCPPthetaMLERIM(arma::mat rimNodesMat, arma::mat aveDiscMatrix);
RcppExport SEXP _RIM_RCPPthetaMLERIM(SEXP rimNodesMatSEXP, SEXP aveDiscMatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type rimNodesMat(rimNodesMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aveDiscMatrix(aveDiscMatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(RCPPthetaMLERIM(rimNodesMat, aveDiscMatrix));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RIM_RCPPSampleFromRIM", (DL_FUNC) &_RIM_RCPPSampleFromRIM, 2},
    {"_RIM_RCPPAverageDiscMatrix", (DL_FUNC) &_RIM_RCPPAverageDiscMatrix, 1},
    {"_RIM_RCPPStructByDP", (DL_FUNC) &_RIM_RCPPStructByDP, 3},
    {"_RIM_RCPPSASearch", (DL_FUNC) &_RIM_RCPPSASearch, 6},
    {"_RIM_RCPPLogProbRIM", (DL_FUNC) &_RIM_RCPPLogProbRIM, 2},
    {"_RIM_RCPPthetaMLERIM", (DL_FUNC) &_RIM_RCPPthetaMLERIM, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_RIM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
