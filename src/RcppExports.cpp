// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/fido.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// conjugateLinearModel
List conjugateLinearModel(const Eigen::Map<Eigen::MatrixXd> Y, const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> Theta, const Eigen::Map<Eigen::MatrixXd> Gamma, const Eigen::Map<Eigen::MatrixXd> Xi, const double upsilon, int n_samples);
RcppExport SEXP _fido_conjugateLinearModel(SEXP YSEXP, SEXP XSEXP, SEXP ThetaSEXP, SEXP GammaSEXP, SEXP XiSEXP, SEXP upsilonSEXP, SEXP n_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Xi(XiSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(conjugateLinearModel(Y, X, Theta, Gamma, Xi, upsilon, n_samples));
    return rcpp_result_gen;
END_RCPP
}
// loglikMaltipooCollapsed
double loglikMaltipooCollapsed(const Eigen::ArrayXXd Y, const double upsilon, const Eigen::MatrixXd Theta, const Eigen::MatrixXd X, const Eigen::MatrixXd KInv, const Eigen::MatrixXd U, Eigen::MatrixXd eta, Eigen::VectorXd ell, bool sylv);
RcppExport SEXP _fido_loglikMaltipooCollapsed(SEXP YSEXP, SEXP upsilonSEXP, SEXP ThetaSEXP, SEXP XSEXP, SEXP KInvSEXP, SEXP USEXP, SEXP etaSEXP, SEXP ellSEXP, SEXP sylvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type KInv(KInvSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type U(USEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type ell(ellSEXP);
    Rcpp::traits::input_parameter< bool >::type sylv(sylvSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikMaltipooCollapsed(Y, upsilon, Theta, X, KInv, U, eta, ell, sylv));
    return rcpp_result_gen;
END_RCPP
}
// gradMaltipooCollapsed
Eigen::VectorXd gradMaltipooCollapsed(const Eigen::ArrayXXd Y, const double upsilon, const Eigen::MatrixXd Theta, const Eigen::MatrixXd X, const Eigen::MatrixXd KInv, const Eigen::MatrixXd U, Eigen::MatrixXd eta, Eigen::VectorXd ell, bool sylv);
RcppExport SEXP _fido_gradMaltipooCollapsed(SEXP YSEXP, SEXP upsilonSEXP, SEXP ThetaSEXP, SEXP XSEXP, SEXP KInvSEXP, SEXP USEXP, SEXP etaSEXP, SEXP ellSEXP, SEXP sylvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type KInv(KInvSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type U(USEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type ell(ellSEXP);
    Rcpp::traits::input_parameter< bool >::type sylv(sylvSEXP);
    rcpp_result_gen = Rcpp::wrap(gradMaltipooCollapsed(Y, upsilon, Theta, X, KInv, U, eta, ell, sylv));
    return rcpp_result_gen;
END_RCPP
}
// hessMaltipooCollapsed
Eigen::MatrixXd hessMaltipooCollapsed(const Eigen::ArrayXXd Y, const double upsilon, const Eigen::MatrixXd Theta, const Eigen::MatrixXd X, const Eigen::MatrixXd KInv, const Eigen::MatrixXd U, Eigen::MatrixXd eta, Eigen::VectorXd ell, bool sylv);
RcppExport SEXP _fido_hessMaltipooCollapsed(SEXP YSEXP, SEXP upsilonSEXP, SEXP ThetaSEXP, SEXP XSEXP, SEXP KInvSEXP, SEXP USEXP, SEXP etaSEXP, SEXP ellSEXP, SEXP sylvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type KInv(KInvSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type U(USEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type ell(ellSEXP);
    Rcpp::traits::input_parameter< bool >::type sylv(sylvSEXP);
    rcpp_result_gen = Rcpp::wrap(hessMaltipooCollapsed(Y, upsilon, Theta, X, KInv, U, eta, ell, sylv));
    return rcpp_result_gen;
END_RCPP
}
// optimMaltipooCollapsed
List optimMaltipooCollapsed(const Eigen::ArrayXXd Y, const double upsilon, const Eigen::MatrixXd Theta, const Eigen::MatrixXd X, const Eigen::MatrixXd KInv, const Eigen::MatrixXd U, Eigen::MatrixXd init, Eigen::VectorXd ellinit, int n_samples, bool calcGradHess, double b1, double b2, double step_size, double epsilon, double eps_f, double eps_g, int max_iter, bool verbose, int verbose_rate, String decomp_method, double eigvalthresh, double jitter);
RcppExport SEXP _fido_optimMaltipooCollapsed(SEXP YSEXP, SEXP upsilonSEXP, SEXP ThetaSEXP, SEXP XSEXP, SEXP KInvSEXP, SEXP USEXP, SEXP initSEXP, SEXP ellinitSEXP, SEXP n_samplesSEXP, SEXP calcGradHessSEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP step_sizeSEXP, SEXP epsilonSEXP, SEXP eps_fSEXP, SEXP eps_gSEXP, SEXP max_iterSEXP, SEXP verboseSEXP, SEXP verbose_rateSEXP, SEXP decomp_methodSEXP, SEXP eigvalthreshSEXP, SEXP jitterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type KInv(KInvSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type U(USEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type init(initSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type ellinit(ellinitSEXP);
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< bool >::type calcGradHess(calcGradHessSEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type step_size(step_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type eps_f(eps_fSEXP);
    Rcpp::traits::input_parameter< double >::type eps_g(eps_gSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type verbose_rate(verbose_rateSEXP);
    Rcpp::traits::input_parameter< String >::type decomp_method(decomp_methodSEXP);
    Rcpp::traits::input_parameter< double >::type eigvalthresh(eigvalthreshSEXP);
    Rcpp::traits::input_parameter< double >::type jitter(jitterSEXP);
    rcpp_result_gen = Rcpp::wrap(optimMaltipooCollapsed(Y, upsilon, Theta, X, KInv, U, init, ellinit, n_samples, calcGradHess, b1, b2, step_size, epsilon, eps_f, eps_g, max_iter, verbose, verbose_rate, decomp_method, eigvalthresh, jitter));
    return rcpp_result_gen;
END_RCPP
}
// loglikPibbleCollapsed
double loglikPibbleCollapsed(const Eigen::ArrayXXd Y, const double upsilon, const Eigen::MatrixXd ThetaX, const Eigen::MatrixXd KInv, const Eigen::MatrixXd AInv, Eigen::MatrixXd eta, bool sylv);
RcppExport SEXP _fido_loglikPibbleCollapsed(SEXP YSEXP, SEXP upsilonSEXP, SEXP ThetaXSEXP, SEXP KInvSEXP, SEXP AInvSEXP, SEXP etaSEXP, SEXP sylvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type ThetaX(ThetaXSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type KInv(KInvSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type AInv(AInvSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< bool >::type sylv(sylvSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikPibbleCollapsed(Y, upsilon, ThetaX, KInv, AInv, eta, sylv));
    return rcpp_result_gen;
END_RCPP
}
// gradPibbleCollapsed
Eigen::VectorXd gradPibbleCollapsed(const Eigen::ArrayXXd Y, const double upsilon, const Eigen::MatrixXd ThetaX, const Eigen::MatrixXd KInv, const Eigen::MatrixXd AInv, Eigen::MatrixXd eta, bool sylv);
RcppExport SEXP _fido_gradPibbleCollapsed(SEXP YSEXP, SEXP upsilonSEXP, SEXP ThetaXSEXP, SEXP KInvSEXP, SEXP AInvSEXP, SEXP etaSEXP, SEXP sylvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type ThetaX(ThetaXSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type KInv(KInvSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type AInv(AInvSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< bool >::type sylv(sylvSEXP);
    rcpp_result_gen = Rcpp::wrap(gradPibbleCollapsed(Y, upsilon, ThetaX, KInv, AInv, eta, sylv));
    return rcpp_result_gen;
END_RCPP
}
// hessPibbleCollapsed
Eigen::MatrixXd hessPibbleCollapsed(const Eigen::ArrayXXd Y, const double upsilon, const Eigen::MatrixXd ThetaX, const Eigen::MatrixXd KInv, const Eigen::MatrixXd AInv, Eigen::MatrixXd eta, bool sylv);
RcppExport SEXP _fido_hessPibbleCollapsed(SEXP YSEXP, SEXP upsilonSEXP, SEXP ThetaXSEXP, SEXP KInvSEXP, SEXP AInvSEXP, SEXP etaSEXP, SEXP sylvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type ThetaX(ThetaXSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type KInv(KInvSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type AInv(AInvSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< bool >::type sylv(sylvSEXP);
    rcpp_result_gen = Rcpp::wrap(hessPibbleCollapsed(Y, upsilon, ThetaX, KInv, AInv, eta, sylv));
    return rcpp_result_gen;
END_RCPP
}
// optimPibbleCollapsed
List optimPibbleCollapsed(const Eigen::ArrayXXd Y, const double upsilon, const Eigen::MatrixXd ThetaX, const Eigen::MatrixXd KInv, const Eigen::MatrixXd AInv, Eigen::MatrixXd init, int n_samples, bool calcGradHess, double b1, double b2, double step_size, double epsilon, double eps_f, double eps_g, int max_iter, bool verbose, int verbose_rate, String decomp_method, String optim_method, double eigvalthresh, double jitter, double multDirichletBoot, bool useSylv, int ncores, long seed);
RcppExport SEXP _fido_optimPibbleCollapsed(SEXP YSEXP, SEXP upsilonSEXP, SEXP ThetaXSEXP, SEXP KInvSEXP, SEXP AInvSEXP, SEXP initSEXP, SEXP n_samplesSEXP, SEXP calcGradHessSEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP step_sizeSEXP, SEXP epsilonSEXP, SEXP eps_fSEXP, SEXP eps_gSEXP, SEXP max_iterSEXP, SEXP verboseSEXP, SEXP verbose_rateSEXP, SEXP decomp_methodSEXP, SEXP optim_methodSEXP, SEXP eigvalthreshSEXP, SEXP jitterSEXP, SEXP multDirichletBootSEXP, SEXP useSylvSEXP, SEXP ncoresSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type ThetaX(ThetaXSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type KInv(KInvSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type AInv(AInvSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type init(initSEXP);
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< bool >::type calcGradHess(calcGradHessSEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type step_size(step_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type eps_f(eps_fSEXP);
    Rcpp::traits::input_parameter< double >::type eps_g(eps_gSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type verbose_rate(verbose_rateSEXP);
    Rcpp::traits::input_parameter< String >::type decomp_method(decomp_methodSEXP);
    Rcpp::traits::input_parameter< String >::type optim_method(optim_methodSEXP);
    Rcpp::traits::input_parameter< double >::type eigvalthresh(eigvalthreshSEXP);
    Rcpp::traits::input_parameter< double >::type jitter(jitterSEXP);
    Rcpp::traits::input_parameter< double >::type multDirichletBoot(multDirichletBootSEXP);
    Rcpp::traits::input_parameter< bool >::type useSylv(useSylvSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< long >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(optimPibbleCollapsed(Y, upsilon, ThetaX, KInv, AInv, init, n_samples, calcGradHess, b1, b2, step_size, epsilon, eps_f, eps_g, max_iter, verbose, verbose_rate, decomp_method, optim_method, eigvalthresh, jitter, multDirichletBoot, useSylv, ncores, seed));
    return rcpp_result_gen;
END_RCPP
}
// uncollapsePibble
List uncollapsePibble(const Eigen::Map<Eigen::VectorXd> eta, const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> Theta, const Eigen::Map<Eigen::MatrixXd> Gamma, const Eigen::Map<Eigen::MatrixXd> Xi, const double upsilon, long seed, bool ret_mean, int ncores);
RcppExport SEXP _fido_uncollapsePibble(SEXP etaSEXP, SEXP XSEXP, SEXP ThetaSEXP, SEXP GammaSEXP, SEXP XiSEXP, SEXP upsilonSEXP, SEXP seedSEXP, SEXP ret_meanSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Xi(XiSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< long >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type ret_mean(ret_meanSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(uncollapsePibble(eta, X, Theta, Gamma, Xi, upsilon, seed, ret_mean, ncores));
    return rcpp_result_gen;
END_RCPP
}
// rMatNormalCholesky_test
Eigen::MatrixXd rMatNormalCholesky_test(Eigen::MatrixXd M, Eigen::MatrixXd LU, Eigen::MatrixXd LV, int discard);
RcppExport SEXP _fido_rMatNormalCholesky_test(SEXP MSEXP, SEXP LUSEXP, SEXP LVSEXP, SEXP discardSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type M(MSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type LU(LUSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type LV(LVSEXP);
    Rcpp::traits::input_parameter< int >::type discard(discardSEXP);
    rcpp_result_gen = Rcpp::wrap(rMatNormalCholesky_test(M, LU, LV, discard));
    return rcpp_result_gen;
END_RCPP
}
// rInvWishRevCholesky_test
Eigen::MatrixXd rInvWishRevCholesky_test(int v, Eigen::MatrixXd Psi);
RcppExport SEXP _fido_rInvWishRevCholesky_test(SEXP vSEXP, SEXP PsiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Psi(PsiSEXP);
    rcpp_result_gen = Rcpp::wrap(rInvWishRevCholesky_test(v, Psi));
    return rcpp_result_gen;
END_RCPP
}
// rInvWishRevCholesky_thread_test
Eigen::MatrixXd rInvWishRevCholesky_thread_test(int v, Eigen::MatrixXd Psi, int discard);
RcppExport SEXP _fido_rInvWishRevCholesky_thread_test(SEXP vSEXP, SEXP PsiSEXP, SEXP discardSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< int >::type discard(discardSEXP);
    rcpp_result_gen = Rcpp::wrap(rInvWishRevCholesky_thread_test(v, Psi, discard));
    return rcpp_result_gen;
END_RCPP
}
// rInvWishRevCholesky_thread_inplace_test
Eigen::MatrixXd rInvWishRevCholesky_thread_inplace_test(int v, Eigen::MatrixXd Psi, int discard);
RcppExport SEXP _fido_rInvWishRevCholesky_thread_inplace_test(SEXP vSEXP, SEXP PsiSEXP, SEXP discardSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< int >::type discard(discardSEXP);
    rcpp_result_gen = Rcpp::wrap(rInvWishRevCholesky_thread_inplace_test(v, Psi, discard));
    return rcpp_result_gen;
END_RCPP
}
// rMatUnitNormal_test1
Eigen::MatrixXd rMatUnitNormal_test1(int n, int m);
RcppExport SEXP _fido_rMatUnitNormal_test1(SEXP nSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(rMatUnitNormal_test1(n, m));
    return rcpp_result_gen;
END_RCPP
}
// rMatUnitNormal_test2
Eigen::MatrixXd rMatUnitNormal_test2(int n);
RcppExport SEXP _fido_rMatUnitNormal_test2(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(rMatUnitNormal_test2(n));
    return rcpp_result_gen;
END_RCPP
}
// uncollapsePibble_sigmaKnown
List uncollapsePibble_sigmaKnown(const Eigen::Map<Eigen::VectorXd> eta, const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> Theta, const Eigen::Map<Eigen::MatrixXd> Gamma, const Eigen::Map<Eigen::MatrixXd> GammaComb, const Eigen::Map<Eigen::MatrixXd> Xi, const Eigen::Map<Eigen::VectorXd> sigma, const double upsilon, long seed, bool ret_mean, bool linear, int ncores);
RcppExport SEXP _fido_uncollapsePibble_sigmaKnown(SEXP etaSEXP, SEXP XSEXP, SEXP ThetaSEXP, SEXP GammaSEXP, SEXP GammaCombSEXP, SEXP XiSEXP, SEXP sigmaSEXP, SEXP upsilonSEXP, SEXP seedSEXP, SEXP ret_meanSEXP, SEXP linearSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type GammaComb(GammaCombSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Xi(XiSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< long >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type ret_mean(ret_meanSEXP);
    Rcpp::traits::input_parameter< bool >::type linear(linearSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(uncollapsePibble_sigmaKnown(eta, X, Theta, Gamma, GammaComb, Xi, sigma, upsilon, seed, ret_mean, linear, ncores));
    return rcpp_result_gen;
END_RCPP
}
// lmvgamma
double lmvgamma(double a, int p);
RcppExport SEXP _fido_lmvgamma(SEXP aSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(lmvgamma(a, p));
    return rcpp_result_gen;
END_RCPP
}
// lmvgamma_deriv
double lmvgamma_deriv(double a, int p);
RcppExport SEXP _fido_lmvgamma_deriv(SEXP aSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(lmvgamma_deriv(a, p));
    return rcpp_result_gen;
END_RCPP
}
// eigen_lap_test
Eigen::MatrixXd eigen_lap_test(int n_samples, Eigen::VectorXd m, Eigen::MatrixXd S, double eigvalthresh);
RcppExport SEXP _fido_eigen_lap_test(SEXP n_samplesSEXP, SEXP mSEXP, SEXP SSEXP, SEXP eigvalthreshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type m(mSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type eigvalthresh(eigvalthreshSEXP);
    rcpp_result_gen = Rcpp::wrap(eigen_lap_test(n_samples, m, S, eigvalthresh));
    return rcpp_result_gen;
END_RCPP
}
// cholesky_lap_test
Eigen::MatrixXd cholesky_lap_test(int n_samples, Eigen::VectorXd m, Eigen::MatrixXd S, double eigvalthresh);
RcppExport SEXP _fido_cholesky_lap_test(SEXP n_samplesSEXP, SEXP mSEXP, SEXP SSEXP, SEXP eigvalthreshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type m(mSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type eigvalthresh(eigvalthreshSEXP);
    rcpp_result_gen = Rcpp::wrap(cholesky_lap_test(n_samples, m, S, eigvalthresh));
    return rcpp_result_gen;
END_RCPP
}
// LaplaceApproximation_test
Eigen::MatrixXd LaplaceApproximation_test(int n_samples, Eigen::VectorXd m, Eigen::MatrixXd S, String decomp_method, double eigvalthresh);
RcppExport SEXP _fido_LaplaceApproximation_test(SEXP n_samplesSEXP, SEXP mSEXP, SEXP SSEXP, SEXP decomp_methodSEXP, SEXP eigvalthreshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type m(mSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type S(SSEXP);
    Rcpp::traits::input_parameter< String >::type decomp_method(decomp_methodSEXP);
    Rcpp::traits::input_parameter< double >::type eigvalthresh(eigvalthreshSEXP);
    rcpp_result_gen = Rcpp::wrap(LaplaceApproximation_test(n_samples, m, S, decomp_method, eigvalthresh));
    return rcpp_result_gen;
END_RCPP
}
// alrInv_default_test
Eigen::MatrixXd alrInv_default_test(Eigen::MatrixXd eta);
RcppExport SEXP _fido_alrInv_default_test(SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(alrInv_default_test(eta));
    return rcpp_result_gen;
END_RCPP
}
// alr_default_test
Eigen::MatrixXd alr_default_test(Eigen::MatrixXd pi);
RcppExport SEXP _fido_alr_default_test(SEXP piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type pi(piSEXP);
    rcpp_result_gen = Rcpp::wrap(alr_default_test(pi));
    return rcpp_result_gen;
END_RCPP
}
// rDirichlet_test
Eigen::MatrixXd rDirichlet_test(int n_samples, Eigen::VectorXd alpha);
RcppExport SEXP _fido_rDirichlet_test(SEXP n_samplesSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(rDirichlet_test(n_samples, alpha));
    return rcpp_result_gen;
END_RCPP
}
// MultDirichletBoot_test
Eigen::MatrixXd MultDirichletBoot_test(int n_samples, Eigen::MatrixXd eta, Eigen::ArrayXXd Y, double pseudocount);
RcppExport SEXP _fido_MultDirichletBoot_test(SEXP n_samplesSEXP, SEXP etaSEXP, SEXP YSEXP, SEXP pseudocountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Eigen::ArrayXXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type pseudocount(pseudocountSEXP);
    rcpp_result_gen = Rcpp::wrap(MultDirichletBoot_test(n_samples, eta, Y, pseudocount));
    return rcpp_result_gen;
END_RCPP
}
// fillUnitNormal_test
void fillUnitNormal_test(Eigen::Map<Eigen::MatrixXd>& Z);
RcppExport SEXP _fido_fillUnitNormal_test(SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type Z(ZSEXP);
    fillUnitNormal_test(Z);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fido_conjugateLinearModel", (DL_FUNC) &_fido_conjugateLinearModel, 7},
    {"_fido_loglikMaltipooCollapsed", (DL_FUNC) &_fido_loglikMaltipooCollapsed, 9},
    {"_fido_gradMaltipooCollapsed", (DL_FUNC) &_fido_gradMaltipooCollapsed, 9},
    {"_fido_hessMaltipooCollapsed", (DL_FUNC) &_fido_hessMaltipooCollapsed, 9},
    {"_fido_optimMaltipooCollapsed", (DL_FUNC) &_fido_optimMaltipooCollapsed, 22},
    {"_fido_loglikPibbleCollapsed", (DL_FUNC) &_fido_loglikPibbleCollapsed, 7},
    {"_fido_gradPibbleCollapsed", (DL_FUNC) &_fido_gradPibbleCollapsed, 7},
    {"_fido_hessPibbleCollapsed", (DL_FUNC) &_fido_hessPibbleCollapsed, 7},
    {"_fido_optimPibbleCollapsed", (DL_FUNC) &_fido_optimPibbleCollapsed, 25},
    {"_fido_uncollapsePibble", (DL_FUNC) &_fido_uncollapsePibble, 9},
    {"_fido_rMatNormalCholesky_test", (DL_FUNC) &_fido_rMatNormalCholesky_test, 4},
    {"_fido_rInvWishRevCholesky_test", (DL_FUNC) &_fido_rInvWishRevCholesky_test, 2},
    {"_fido_rInvWishRevCholesky_thread_test", (DL_FUNC) &_fido_rInvWishRevCholesky_thread_test, 3},
    {"_fido_rInvWishRevCholesky_thread_inplace_test", (DL_FUNC) &_fido_rInvWishRevCholesky_thread_inplace_test, 3},
    {"_fido_rMatUnitNormal_test1", (DL_FUNC) &_fido_rMatUnitNormal_test1, 2},
    {"_fido_rMatUnitNormal_test2", (DL_FUNC) &_fido_rMatUnitNormal_test2, 1},
    {"_fido_uncollapsePibble_sigmaKnown", (DL_FUNC) &_fido_uncollapsePibble_sigmaKnown, 12},
    {"_fido_lmvgamma", (DL_FUNC) &_fido_lmvgamma, 2},
    {"_fido_lmvgamma_deriv", (DL_FUNC) &_fido_lmvgamma_deriv, 2},
    {"_fido_eigen_lap_test", (DL_FUNC) &_fido_eigen_lap_test, 4},
    {"_fido_cholesky_lap_test", (DL_FUNC) &_fido_cholesky_lap_test, 4},
    {"_fido_LaplaceApproximation_test", (DL_FUNC) &_fido_LaplaceApproximation_test, 5},
    {"_fido_alrInv_default_test", (DL_FUNC) &_fido_alrInv_default_test, 1},
    {"_fido_alr_default_test", (DL_FUNC) &_fido_alr_default_test, 1},
    {"_fido_rDirichlet_test", (DL_FUNC) &_fido_rDirichlet_test, 2},
    {"_fido_MultDirichletBoot_test", (DL_FUNC) &_fido_MultDirichletBoot_test, 4},
    {"_fido_fillUnitNormal_test", (DL_FUNC) &_fido_fillUnitNormal_test, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_fido(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
