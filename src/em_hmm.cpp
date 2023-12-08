#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

vec my_pava(vec &values, vec &weight, int decreasing);
vec replis(vec &pi, mat &A, vec &f1, vec &f2);
vec fdr(vec &rLIS);

//' @title EM algorithm in combination with a non-parametric algorithm for estimation of the rLIS statistic.
//' @description Estimate the rLIS values accounting for the linkage disequilibrium across two genome-wide association studies via the four-state hidden Markov model. Apply a step-up procedure to control the FDR of replicability null.
//' @param pa_in A numeric vector of p-values from study 1.
//' @param pb_in A numeric vector of p-values from study 2.
//' @param pi0a_in An initial estimate of the null probability in study 1.
//' @param pi0b_in An initial estimate of the null probability in study 2.
//'
//' @return
//' \item{rLIS}{The estimated rLIS for replicability null.}
//' \item{fdr}{The adjusted values based on rLIS for FDR control.}
//' \item{loglik}{The log-likelihood value with converged estimates of the unknowns.}
//' \item{pi}{An estimate of the stationary probabilities of four states {(0,0), (0,1), (1,0), (1,1)}.}
//' \item{A}{An estimate of the 4-by-4 transition matrix.}
//' \item{f1}{A non-parametric estimate for the non-null probability density function in study 1.}
//' \item{f2}{A non-parametric estimate for the non-null probability density function in study 2.}
//'
//' @export
// [[Rcpp::export]]
SEXP em_hmm(SEXP pa_in, SEXP pb_in, SEXP pi0a_in, SEXP pi0b_in) {
  try{
    const int maxIter = 200;
    const double tol = 1e-3;
    const double pvalue_cutoff = 1e-15;
    const double f_cutoff = 1e-15;

    vec pa = as<arma::vec>(pa_in), pb = as<arma::vec>(pb_in);
    const double pi0_pa = Rcpp::as<double>(pi0a_in), pi0_pb = Rcpp::as<double>(pi0b_in);

    uword J = pa.size();
    double min_a = pa.elem(find(pa>0)).min(), min_b = pb.elem(find(pb>0)).min();
    pa.elem(find(pa<=0)).fill(pvalue_cutoff < min_a ? pvalue_cutoff : min_a);
    pb.elem(find(pb<=0)).fill(pvalue_cutoff < min_b ? pvalue_cutoff : min_b);

    vec p1 = sort(pa), p2 = sort(pb);
    uvec ix1 = sort_index(pa), ix2 = sort_index(pb);

    vec p1_diff(J), p2_diff(J);
    p1_diff(0) = p1(0);
    p2_diff(0) = p2(0);
    for(uword i = 1; i<J; ++i){
      p1_diff(i) = p1(i) - p1(i-1);
      p2_diff(i) = p2(i) - p2(i-1);
    }

    // Initialization
    vec pi(4);
    pi(0) = pi0_pa * pi0_pb;
    pi(1) = pi0_pa * (1 - pi0_pb);
    pi(2) = (1 - pi0_pa) * pi0_pb;
    pi(3) = (1 - pi0_pa) * (1 - pi0_pb);

    mat A = {{0.9, 0.025, 0.025, 0.05},
    {0.025, 0.9, 0.025, 0.05},
    {0.025, 0.025, 0.9, 0.05},
    {0.1, 0.1, 0.1, 0.7}};

    vec f0 = ones(J,1), f1(J), f2(J);
    f1 = 1 - pa;
    f2 = 1 - pb;

    mat f = join_rows(f0, f2, f1, f1%f2);

    vec loglik(maxIter);
    loglik(0) = -datum::inf;

    mat alpha = zeros(J, 4), beta = zeros(J, 4), gamma = zeros(J, 4);
    cube xi = zeros(4, 4, J-1);

    vec xi_kl(J-1), Q1(J), Q2(J), y1(J), y2(J), res1(J), res2(J);
    mat xi_k(4,J-1);
    mat f_jp1(4,4);
    double sub_loglik, loglik_delta;

    // std::cout << "EM begins:" << std::endl;
    Rcout << "EM begins:" << "\n";

    for (int i = 1; i < maxIter; i++){
      // E-step
      // calculate the forward, backward and posterior probabilities based on current HMM
      alpha.row(0) = pi.t() % f.row(0);
      alpha.row(0) = alpha.row(0)/sum(alpha.row(0));
      
      for (uword j = 1; j < J; j++){
        alpha(j,0) = sum(alpha.row(j-1).t() % A.col(0)) * f(j,0);
        alpha(j,1) = sum(alpha.row(j-1).t() % A.col(1)) * f(j,1);
        alpha(j,2) = sum(alpha.row(j-1).t() % A.col(2)) * f(j,2);
        alpha(j,3) = sum(alpha.row(j-1).t() % A.col(3)) * f(j,3);

        alpha.row(j) = alpha.row(j)/sum(alpha.row(j));
      }

      beta.row(J-1).fill(0.25);
      for(uword j = J-2; j > 0; j--){
        beta(j,0) = sum(beta.row(j+1) % A.row(0) % f.row(j+1));
        beta(j,1) = sum(beta.row(j+1) % A.row(1) % f.row(j+1));
        beta(j,2) = sum(beta.row(j+1) % A.row(2) % f.row(j+1));
        beta(j,3) = sum(beta.row(j+1) % A.row(3) % f.row(j+1));
        beta.row(j) = beta.row(j)/sum(beta.row(j));
      }
      beta(0,0) = sum(beta.row(1) % A.row(0) % f.row(1));
      beta(0,1) = sum(beta.row(1) % A.row(1) % f.row(1));
      beta(0,2) = sum(beta.row(1) % A.row(2) % f.row(1));
      beta(0,3) = sum(beta.row(1) % A.row(3) % f.row(1));
      beta.row(0) = beta.row(0)/sum(beta.row(0));

      for(uword j = 0; j < J; j++){
        gamma.row(j) = alpha.row(j) % beta.row(j) / sum(alpha.row(j) % beta.row(j));
      }

      for(uword j = 0; j < J-1; j++){
        f_jp1 = repmat(f.row(j+1), 4, 1);
        xi.slice(j) = alpha.row(j).t() * beta.row(j+1) % A % f_jp1/
          accu(alpha.row(j).t() * beta.row(j+1) % A % f_jp1);
      }

      // M-step
      // update the parameters pi and A
      pi = gamma.row(0).t();
      for(int k = 0; k < 4; k++){
        for(int l = 0; l < 4; l++){
          xi_kl = xi.tube(k,l);
          xi_k = xi.row(k);
          A(k,l) = sum(xi_kl)/accu(xi_k);
        }
      }

      // update f1 and f2
      Q1 = gamma.col(2) + gamma.col(3);
      Q2 = gamma.col(1) + gamma.col(3);
      Q1 = Q1(ix1);
      Q2 = Q2(ix2);

      y1 = - p1_diff * sum(Q1) / Q1;
      y2 = - p2_diff * sum(Q2) / Q2;

      y1.elem(find_nonfinite(y1)).fill(y1.elem(find_finite(y1)).min());
      y2.elem(find_nonfinite(y2)).fill(y2.elem(find_finite(y2)).min());

      res1 = my_pava(y1, Q1, 1);
      res2 = my_pava(y2, Q2, 1);

      f1 = -1 / res1;
      f1 = f1 / sum(f1 % p1_diff);
      f1(ix1) = f1;
      f1.elem(find_nan(f1)).fill(f1.min());

      f2 = -1 / res2;
      f2 = f2 / sum(f2 % p2_diff);
      f2(ix2) = f2;
      f2.elem(find_nan(f2)).fill(f2.min());

      double min_f1 = f1.elem(find(f1>0)).min(), min_f2 = f2.elem(find(f2>0)).min();
      f1.elem(find(f1<=0)).fill(f_cutoff < min_f1 ? f_cutoff : min_f1);
      f2.elem(find(f2<=0)).fill(f_cutoff < min_f2 ? f_cutoff : min_f2);

      f = join_rows(f0, f2, f1, f1%f2);

      // calculate the updated log-likelihood
      sub_loglik = 0;
      for(uword j = 0; j < J - 1; j++){
        sub_loglik = sub_loglik + accu(log(A) % xi.slice(j));
      }

      loglik(i) = sum(log(pi).t() % gamma.row(1)) + sub_loglik + accu(gamma % log(f));
      loglik_delta = std::fabs((loglik(i) / loglik(i-1)) - 1);

      // std::cout<<i<<". "<< loglik(i) << ", delta = " << loglik_delta << std::endl;
      Rcout << i << ". " << loglik(i) << ", delta = " << loglik_delta << "\n";

      if(loglik_delta < tol){
        break;
      }
    }

    vec rLIS = replis(pi, A, f1, f2);
    vec rlis_adj = fdr(rLIS);

    return Rcpp::List::create(Rcpp::Named("rLIS") = rLIS,
                              Rcpp::Named("fdr") = rlis_adj,
                              Rcpp::Named("loglik") = loglik,
                              Rcpp::Named("pi") = pi,
                              Rcpp::Named("A") = A,
                              Rcpp::Named("f1") = f1,
                              Rcpp::Named("f2") = f2);
  } catch( std::exception &ex ) {
    forward_exception_to_r(ex);
    return Rcpp::List::create();
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)..." );
    return Rcpp::List::create();
  }
}

double na_rm(vec &p){
  double _min = std::numeric_limits<double>::max();

  p.for_each([&_min](double &val) { if(val > 0 && val < _min) _min = val; });

  return _min;
}

vec replis(vec &pi, mat &A, vec &f1, vec &f2){
  uword J = f1.size();
  vec f0 = ones(J,1);
  mat f = join_rows(f0, f2, f1, f1%f2);

  mat alpha = zeros(J, 4), beta = zeros(J, 4), f_jp1(4,4);

  alpha.row(0) = pi.t() % f.row(0);
  alpha.row(0) = alpha.row(0)/sum(alpha.row(0));

  for (uword j = 1; j < J; j++){
    alpha(j,0) = sum(alpha.row(j-1).t() % A.col(0)) * f(j,0);
    alpha(j,1) = sum(alpha.row(j-1).t() % A.col(1)) * f(j,1);
    alpha(j,2) = sum(alpha.row(j-1).t() % A.col(2)) * f(j,2);
    alpha(j,3) = sum(alpha.row(j-1).t() % A.col(3)) * f(j,3);

    alpha.row(j) = alpha.row(j)/sum(alpha.row(j));
  }

  beta.row(J-1).fill(0.25);
  for(uword j = J-2; j > 0; j--){
    beta(j,0) = sum(beta.row(j+1) % A.row(0) % f.row(j+1));
    beta(j,1) = sum(beta.row(j+1) % A.row(1) % f.row(j+1));
    beta(j,2) = sum(beta.row(j+1) % A.row(2) % f.row(j+1));
    beta(j,3) = sum(beta.row(j+1) % A.row(3) % f.row(j+1));
    beta.row(j) = beta.row(j)/sum(beta.row(j));
  }
  beta(0,0) = sum(beta.row(1) % A.row(0) % f.row(1));
  beta(0,1) = sum(beta.row(1) % A.row(1) % f.row(1));
  beta(0,2) = sum(beta.row(1) % A.row(2) % f.row(1));
  beta(0,3) = sum(beta.row(1) % A.row(3) % f.row(1));
  beta.row(0) = beta.row(0)/sum(beta.row(0));

  vec rLIS(J);
  for(uword j = 0; j < J; j++){
    rLIS(j) = sum(alpha.row(j).head(3) % beta.row(j).head(3))/
      sum(alpha.row(j) % beta.row(j));
  }

  return rLIS;
}

vec fdr(vec &rLIS){
  uword J = rLIS.size();

  vec ordered_lis = sort(rLIS), s = linspace(1,J,J);
  uvec ix_lis = sort_index(rLIS);

  vec radj = cumsum(ordered_lis)/s;
  radj(ix_lis) = radj;

  return radj;
}

vec my_pava(vec &values, vec &weight, int decreasing)
{
  if(decreasing == 1){
    values = reverse(values);
    weight = reverse(weight);
  }
  vec w(values.size(), fill::zeros);
  vec x(values.size(), fill::zeros);
  x(0) = values(0);
  w(0) = weight(0);
  uword j = 0;
  vec s(values.size(), fill::zeros);

  for (uword i = 1; i < values.size(); i++) {
    j += 1;
    x(j) = values(i);
    w(j) = weight(i);
    while (j > 0 && x(j - 1)>x(j)) {
      x(j - 1) = (w(j) * x(j) + w(j - 1) * x(j - 1)) / (w(j) + w(j - 1));
      w(j - 1) += w(j);
      j -= 1;
    }
    s(j + 1) = i + 1;
  }

  vec ww(values.size(), fill::zeros);
  vec xx(values.size(), fill::zeros);
  for (uword k = 0; k < j + 1; k++) {
    for (uword i = s(k); i < s(k + 1); i++) {
      ww(i) = w(k);
      xx(i) = x(k);
    }
  }

  if(decreasing == 1){
    xx = reverse(xx);
  }

  return xx;
}
