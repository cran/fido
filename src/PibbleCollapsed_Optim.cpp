#include <fido.h>
#include <Rcpp/Benchmark/Timer.h>

// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppEigen)]]



using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using Eigen::VectorXd;

//' Function to Optimize the Collapsed Pibble Model
//' 
//' See details for model. Should likely be followed by function 
//' \code{\link{uncollapsePibble}}. Notation: \code{N} is number of samples,
//' \code{D} is number of multinomial categories, and \code{Q} is number
//' of covariates. 
//' 
//' @param Y D x N matrix of counts
//' @param upsilon (must be > D)
//' @param ThetaX D-1 x N matrix formed by Theta*X (Theta is Prior mean 
//'    for regression coefficients) 
//' @param KInv D-1 x D-1 precision matrix (inverse of Xi)
//' @param AInv N x N precision matrix given by \eqn{(I_N + X'*Gamma*X)^{-1}}
//' @param init D-1 x N matrix of initial guess for eta used for optimization
//' @param n_samples number of samples for Laplace Approximation (=0 very fast
//'    as no inversion or decomposition of Hessian is required)
//' @param calcGradHess if n_samples=0 should Gradient and Hessian 
//'   still be calculated using closed form solutions?
//' @param b1 (ADAM) 1st moment decay parameter (recommend 0.9) "aka momentum"
//' @param b2 (ADAM) 2nd moment decay parameter (recommend 0.99 or 0.999)
//' @param step_size (ADAM) step size for descent (recommend 0.001-0.003)
//' @param epsilon (ADAM) parameter to avoid divide by zero
//' @param eps_f (ADAM) normalized function improvement stopping criteria 
//' @param eps_g (ADAM) normalized gradient magnitude stopping criteria
//' @param max_iter (ADAM) maximum number of iterations before stopping
//' @param verbose (ADAM) if true will print stats for stopping criteria and 
//'   iteration number
//' @param verbose_rate (ADAM) rate to print verbose stats to screen
//' @param decomp_method decomposition of hessian for Laplace approximation
//'   'eigen' (more stable-slightly, slower) or 'cholesky' (less stable, faster, default)
//' @param optim_method (default:"lbfgs") or "adam"
//' @param eigvalthresh threshold for negative eigenvalues in 
//'   decomposition of negative inverse hessian (should be <=0)
//' @param jitter (default: 0) if >=0 then adds that factor to diagonal of Hessian 
//' before decomposition (to improve matrix conditioning)
//' @param multDirichletBoot if >0 then it overrides laplace approximation and samples
//'  eta efficiently at MAP estimate from pseudo Multinomial-Dirichlet posterior.
//' @param useSylv (default: true) if N<D-1 uses Sylvester Determinant Identity
//'   to speed up calculation of log-likelihood and gradients. 
//' @param ncores (default:-1) number of cores to use, if ncores==-1 then 
//' uses default from OpenMP typically to use all available cores. 
//' @param seed (random seed for Laplace approximation -- integer)
//'  
//' @details Notation: Let \eqn{Z_j} denote the J-th row of a matrix Z.
//' Model:
//'    \deqn{Y_j \sim Multinomial(\pi_j)}{Y_j \sim Multinomial(Pi_j)}
//'    \deqn{\pi_j = \Phi^{-1}(\eta_j)}{Pi_j = Phi^(-1)(Eta_j)}
//'    \deqn{\eta \sim T_{D-1, N}(\upsilon, \Theta X, K, A)}{Eta \sim T_{D-1, N}(upsilon, Theta*X, K, A)}
//' Where \eqn{A = I_N + X  \Gamma X'}{A = I_N + X  Gamma  X'}, K is a (D-1)x(D-1) covariance 
//' matrix, \eqn{\Gamma}{Gamma} is a Q x Q covariance matrix, and \eqn{\Phi^{-1}}{Phi^(-1)} is ALRInv_D 
//' transform. 
//' 
//' Gradient and Hessian calculations are fast as they are computed using closed
//' form solutions. That said, the Hessian matrix can be quite large 
//' \[N*(D-1) x N*(D-1)\] and storage may be an issue. 
//' 
//' Note: Warnings about large negative eigenvalues can either signal 
//' that the optimizer did not reach an optima or (more commonly in my experience)
//' that the prior / degrees of freedom for the covariance (given by parameters
//' \code{upsilon} and \code{KInv}) were too specific and at odds with the observed data.
//' If you get this warning try the following. 
//' 1. Try restarting the optimization using a different initial guess for eta
//' 2. Try decreasing (or even increasing )\code{step_size} (by increments of 0.001 or 0.002) 
//'   and increasing \code{max_iter} parameters in optimizer. Also can try 
//'   increasing \code{b1} to 0.99 and decreasing \code{eps_f} by a few orders
//'   of magnitude
//' 3. Try relaxing prior assumptions regarding covariance matrix. (e.g., may want
//' to consider decreasing parameter \code{upsilon} closer to a minimum value of 
//' D)
//' 4. Try adding small amount of jitter (e.g., set \code{jitter=1e-5}) to address
//'   potential floating point errors. 
//' @return List containing (all with respect to found optima)
//' 1. LogLik - Log Likelihood of collapsed model (up to proportionality constant)
//' 2. Gradient - (if \code{calcGradHess}=true)
//' 3. Hessian - (if \code{calcGradHess}=true) of the POSITIVE LOG POSTERIOR
//' 4. Pars - Parameter value of eta at optima
//' 5. Samples - (D-1) x N x n_samples array containing posterior samples of eta 
//'    based on Laplace approximation (if n_samples>0)
//' 6. Timer - Vector of Execution Times
//' 7. logInvNegHessDet - the log determinant of the covariacne of the Laplace 
//'    approximation, useful for calculating marginal likelihood 
//' 8. logMarginalLikelihood - A calculation of the log marginal likelihood based on
//'    the laplace approximation
//' @md 
//' @export
//' @name optimPibbleCollapsed
//' @references S. Ruder (2016) \emph{An overview of gradient descent 
//' optimization algorithms}. arXiv 1609.04747
//' 
//' JD Silverman K Roche, ZC Holmes, LA David, S Mukherjee. 
//'   \emph{Bayesian Multinomial Logistic Normal Models through Marginally Latent Matrix-T Processes}. 
//'   2022, Journal of Machine Learning
//' @seealso \code{\link{uncollapsePibble}}
//' @examples
//' sim <- pibble_sim()
//' 
//' # Fit model for eta
//' fit <- optimPibbleCollapsed(sim$Y, sim$upsilon, sim$Theta%*%sim$X, sim$KInv, 
//'                              sim$AInv, random_pibble_init(sim$Y))  
// [[Rcpp::export]]
List optimPibbleCollapsed(const Eigen::ArrayXXd Y, 
               const double upsilon, 
               const Eigen::MatrixXd ThetaX, 
               const Eigen::MatrixXd KInv, 
               const Eigen::MatrixXd AInv, 
               Eigen::MatrixXd init, 
               int n_samples=2000, 
               bool calcGradHess = true,
               double b1 = 0.9,         
               double b2 = 0.99,        
               double step_size = 0.003, // was called eta in ADAM code
               double epsilon = 10e-7, 
               double eps_f=1e-10,       
               double eps_g=1e-4,       
               int max_iter=10000,      
               bool verbose=false,      
               int verbose_rate=10,
               String decomp_method="cholesky",
               String optim_method="lbfgs",
               double eigvalthresh=0, 
               double jitter=0,
               double multDirichletBoot = -1.0, 
               bool useSylv = true, 
               int ncores=-1, 
               long seed=-1){  
  #ifdef FIDO_USE_PARALLEL 
    Eigen::initParallel();
    if (ncores > 0) Eigen::setNbThreads(ncores);
  #endif 
  Timer timer;
  timer.step("Overall_start");
  int N = Y.cols();
  int D = Y.rows();
  PibbleCollapsed cm(Y, upsilon, ThetaX, KInv, AInv, useSylv);
  Map<VectorXd> eta(init.data(), init.size()); // will rewrite by optim
  double nllopt; // NEGATIVE LogLik at optim
  List out(8);
  out.names() = CharacterVector::create("LogLik", "Gradient", "Hessian",
					"Pars", "Samples", "Timer", "logInvNegHessDet",
					"logMarginalLikelihood");
  
  // Pick optimizer (ADAM - without perturbation appears to be best)
  //   ADAM with perturbations not fully implemented
  timer.step("Optimization_start");
  int status;
  if (optim_method=="lbfgs"){
    status = Numer::optim_lbfgs(cm, eta, nllopt, max_iter, eps_f, eps_g);
  } else if (optim_method=="adam"){
    status = adam::optim_adam(cm, eta, nllopt, b1, b2, step_size, epsilon, 
                                  eps_f, eps_g, max_iter, verbose, verbose_rate);  
  } else {
    Rcpp::stop("unrecognized optimization method");
  }
   
  timer.step("Optimization_stop");

  if (status<0)
    Rcpp::warning("Max Iterations Hit, May not be at optima");
  Map<MatrixXd> etamat(eta.data(), D-1, N);
  out[0] = -nllopt; // Return (positive) LogLik
  out[3] = etamat;
  
  if (n_samples > 0 || calcGradHess){
    if (verbose) Rcout << "Allocating for Gradient" << std::endl;
    VectorXd grad(N*(D-1));
    MatrixXd hess; // don't preallocate this thing could be unneeded
    if (verbose) Rcout << "Calculating Gradient" << std::endl;
    grad = cm.calcGrad(); // should have eta at optima already
    
    // "Multinomial-Dirichlet" option
    if (multDirichletBoot>=0.0){
      timer.step("MultDirichletBoot_start");
      if (verbose) Rcout << "Performing Multinomial Dirichlet Bootstrap" << std::endl;
      MatrixXd samp = MultDirichletBoot::MultDirichletBoot(n_samples, etamat, Y, 
                                                           multDirichletBoot);
      timer.step("MultDirichletBoot_stop");
      out[1] = R_NilValue;
      out[2] = R_NilValue;
      IntegerVector d = IntegerVector::create(D-1, N, n_samples);
      NumericVector samples = wrap(samp);
      samples.attr("dim") = d; // convert to 3d array for return to R
      out[4] = samples;
      timer.step("Overall_stop");
      NumericVector t(timer);
      out[5] = t;
      out[6] = R_NilValue;
      out[7] = R_NilValue;
      return out;
    }
    // "Multinomial-Dirchlet" option 
    if (verbose) Rcout << "Calculating Hessian" << std::endl;
    timer.step("HessianCalculation_start");
    hess = -cm.calcHess(); // should have eta at optima already
    timer.step("HessianCalculation_Stop");
    out[1] = grad;
    if ((N * (D-1)) > 44750){
      Rcpp::warning("Hessian is to large to return to R");
    } else {
      if (calcGradHess)
        out[2] = hess;    
    }

    if (n_samples>0){
      // Laplace Approximation
      int status;
      timer.step("LaplaceApproximation_start");
      MatrixXd samp = MatrixXd::Zero(N*(D-1), n_samples);
      double logInvNegHessDet;
      status = lapap::LaplaceApproximation(samp, eta, hess, 
                                           decomp_method, eigvalthresh, 
                                           jitter, 
                                           logInvNegHessDet, 
                                           seed);
      timer.step("LaplaceApproximation_stop");
      if (status != 0){
        Rcpp::warning("Decomposition of Hessian Failed, returning MAP Estimate only");
        return out;
      }
      out[6] = logInvNegHessDet;

      // (log)marginallikelihood calculation
      double normConstT; 
      double normConstMult;
      Eigen::LLT<MatrixXd> KInvSqrt; 
      Eigen::LLT<MatrixXd> AInvSqrt; 
      KInvSqrt.compute(KInv);
      AInvSqrt.compute(AInv);
      // // log |K|^(-(D-1)/2) MPN: Comments are wrong, code is correct.
      normConstT = N*0.5*KInvSqrt.matrixLLT().diagonal().array().log().sum();
      // // log |A|^(-N/2) MPN: Comments are wrong, code is correct.
      normConstT = normConstT + (D-1)*0.5*AInvSqrt.matrixLLT().diagonal().array().log().sum();
      // // TODO normConstT needs to account for the multivariate gamma terms...
      ArrayXd numer = Y.colwise().sum();
      numer += ArrayXd::Ones(N);
      Eigen::ArrayXXd denom = Y;
      denom += ArrayXXd::Ones(D, N);
      for (int n=0; n<N; n++){
	numer(n) = lmvgamma(numer(n), 1);
	for (int d=0; d<D; d++){
	  denom(d,n) = lmvgamma(denom(d,n), 1); 
	}
      }
      normConstMult = numer.sum()-denom.sum(); 
      out[7] = -nllopt + 0.5*logInvNegHessDet  + (D-1)*N/2*log(2*3.1415)  + normConstT + normConstMult - 0.5*log((D-1)*N);
      IntegerVector d = IntegerVector::create(D-1, N, n_samples);
      NumericVector samples = wrap(samp);
      samples.attr("dim") = d; // convert to 3d array for return to R
      out[4] = samples;
    } // endif n_samples || calcGradHess
  } // endif n_samples || calcGradHess
  timer.step("Overall_stop");
  NumericVector t(timer);
  out[5] = t;
  return out;
}
