// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;


RcppExport SEXP _rcpp_module_boot_stan_fit4FA_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4corr_est_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4gFA_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4gNBFA_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4gPFA_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_stan_fit4FA_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4FA_mod, 0},
    {"_rcpp_module_boot_stan_fit4corr_est_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4corr_est_mod, 0},
    {"_rcpp_module_boot_stan_fit4gFA_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4gFA_mod, 0},
    {"_rcpp_module_boot_stan_fit4gNBFA_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4gNBFA_mod, 0},
    {"_rcpp_module_boot_stan_fit4gPFA_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4gPFA_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_gfac(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
