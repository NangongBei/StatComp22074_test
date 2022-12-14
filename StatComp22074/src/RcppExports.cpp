// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Gauss_Seidel
NumericVector Gauss_Seidel(NumericMatrix A, NumericVector b, NumericVector x0, int N);
RcppExport SEXP _StatComp22074_Gauss_Seidel(SEXP ASEXP, SEXP bSEXP, SEXP x0SEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(Gauss_Seidel(A, b, x0, N));
    return rcpp_result_gen;
END_RCPP
}
// RandomC
NumericMatrix RandomC(double x0, double y0, int N);
RcppExport SEXP _StatComp22074_RandomC(SEXP x0SEXP, SEXP y0SEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(RandomC(x0, y0, N));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StatComp22074_Gauss_Seidel", (DL_FUNC) &_StatComp22074_Gauss_Seidel, 4},
    {"_StatComp22074_RandomC", (DL_FUNC) &_StatComp22074_RandomC, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_StatComp22074(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
