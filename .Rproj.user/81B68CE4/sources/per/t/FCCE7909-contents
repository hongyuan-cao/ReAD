############################################################################################
# Package: ReAD
# Version: 1.0.0
# Data: 2023-06-12
# Authors: Y. Li, H. Lei, X. Wen and H. Cao.
############################################################################################
#' @title Replicability analysis across two genome-wide association studies accounting for the linkage disequilibrium structure.
#'
#' @param pa A numeric vector of p-values from study 1.
#' @param pb A numeric vector of p-values from study 2.
#'
#' @return A list:
#' \item{rLIS}{The estimated rLIS for replicability null.}
#' \item{fdr}{The adjusted values based on rLIS for FDR control.}
#' \item{loglik}{The log-likelihood value with converged estimates of the unknowns.}
#' \item{pi}{An estimate of the stationary probabilities of four states {(0,0), (0,1), (1,0), (1,1)}.}
#' \item{A}{An estimate of the 4-by-4 transition matrix.}
#' \item{f1}{A non-parametric estimate for the non-null probability density function in study 1.}
#' \item{f2}{A non-parametric estimate for the non-null probability density function in study 2.}
#'
#' @importFrom Rcpp, RcppArmadillo, qvalue
#'
#' @export
#'
#' @examples
#' # Simulate p-values in two studies locally dependent via a four-state hidden Markov model
#' data <- SimuData(m = 10000)
#' p1 = data$pa; p2 = data$pb; theta1 = data$theta1; theta2 = data$theta2
#' # Run ReAD to identify replicable signals
#' res.read = ReAD(p1, p2)
#' sig.idx = which(res.read$fdr <= 0.05)
#'
ReAD <- function(pa, pb){
  pvals.cutoff = 1e-15
  pa[pa == 0] <- min(min(pa[pa != 0]), pvals.cutoff)
  pb[pb == 0] <- min(min(pb[pb != 0]), pvals.cutoff)

  pi0_pa <- min(pi0est(pa)$pi0, 0.999)
  pi0_pb <- min(pi0est(pb)$pi0, 0.999)

  res <- em_hmm(pa, pb, pi0_pa, pi0_pb)

  return(res)
}
