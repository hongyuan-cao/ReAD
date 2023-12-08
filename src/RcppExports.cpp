// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// em_hmm
SEXP em_hmm(SEXP pa_in, SEXP pb_in, SEXP pi0a_in, SEXP pi0b_in);
RcppExport SEXP _ReAD_em_hmm(SEXP pa_inSEXP, SEXP pb_inSEXP, SEXP pi0a_inSEXP, SEXP pi0b_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pa_in(pa_inSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pb_in(pb_inSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pi0a_in(pi0a_inSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pi0b_in(pi0b_inSEXP);
    rcpp_result_gen = Rcpp::wrap(em_hmm(pa_in, pb_in, pi0a_in, pi0b_in));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ReAD_em_hmm", (DL_FUNC) &_ReAD_em_hmm, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_ReAD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
