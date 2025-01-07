// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/EpiSky.h"
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// systematic_sample_cpp
NumericVector systematic_sample_cpp(int n_particles, NumericVector norm_weights);
static SEXP _EpiSky_systematic_sample_cpp_try(SEXP n_particlesSEXP, SEXP norm_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type n_particles(n_particlesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type norm_weights(norm_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(systematic_sample_cpp(n_particles, norm_weights));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _EpiSky_systematic_sample_cpp(SEXP n_particlesSEXP, SEXP norm_weightsSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_EpiSky_systematic_sample_cpp_try(n_particlesSEXP, norm_weightsSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _EpiSky_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("NumericVector(*systematic_sample_cpp)(int,NumericVector)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _EpiSky_RcppExport_registerCCallable() { 
    R_RegisterCCallable("EpiSky", "_EpiSky_systematic_sample_cpp", (DL_FUNC)_EpiSky_systematic_sample_cpp_try);
    R_RegisterCCallable("EpiSky", "_EpiSky_RcppExport_validate", (DL_FUNC)_EpiSky_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_EpiSky_systematic_sample_cpp", (DL_FUNC) &_EpiSky_systematic_sample_cpp, 2},
    {"_EpiSky_RcppExport_registerCCallable", (DL_FUNC) &_EpiSky_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_EpiSky(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
