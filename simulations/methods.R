
# ----------- Functions for performing the MaRR procedure ------------------------

# get.SS: a function to produce a vector of sums of squared differences between observed
#         and actual survival functions for k-hat=0,...,n-1
# INPUTS: max.rank = a vector of maximum rank statistics
# OUTPUTS: mySS = a vector of SS(i/k) values as described in manuscript
get.SS = function(max.rank)
{
  n = length(max.rank)
  x = 1:n / n
  mySS = rep(NA, n)
  my.W = sapply(1:n, function(i) {
    sum(max.rank == i)
  })
  surv.function = (n - cumsum(my.W)) / n
  for (k in 0:(n - 1))
  {
    i = k + 1
    pi1 = k / n
    pi0 = 1 - pi1
    temp.W = surv.function[i:n]
    Sn = pi0 * (1 - (x - pi1) ^ 2 / pi0 ^ 2)[i:n]
    temp.diff = temp.W - Sn
    sq.diff = (temp.diff * temp.diff) / n
    mySS[i] = sum(sq.diff) / pi0
  }
  return(mySS)
}

# est.fdr: a function to produce a vector of estimated false discovery rats based on
#          a vector of maximum rank statistics and a khat.
# INPUTS: max.rank = a vector of maximum rank statistics
#         khat = a value of khat calculated from the argmin of mySS
# OUTPUTS: a vector of estimated FDR values for each potential threshold Nhat=1,...,n
est.fdr = function(khat, max.rank)
{
  n = length(max.rank)
  Nhat = (khat + 1):n
  Q.khat = sum(max.rank <= khat)
  R.Nhat = sapply(Nhat, function(i) {
    sum(max.rank <= i)
  })
  indicate.R = as.numeric(R.Nhat > 0)
  R.Nhat[indicate.R == 0] = 1
  temp = Nhat - khat
  numer = temp * temp
  denom = (n - khat) * R.Nhat
  FDR.Nhatkhat = (numer / denom) * indicate.R
  return(c(rep(0, khat), FDR.Nhatkhat))
}

# MaRR: a function that performs the MaRR procedure based on a vector of maximum rank statistics.
# INPUTS: max.rank = a vector of maximum rank statistics
#         cutoff = a value between 0 and 1 that provides the maximum allowed value for pi-hat.
#         alpha = desired level of FDR control
#         khat.to.zero = T/F, whether or not to set k-hat to zero for mFDR calculation (recommended for very small pi1)
# OUTPUTS: khat = n*pi-hat, discrete estimate of where irreproducible signals begin
#          Nhat = estimated cut-off for maximum ranks that will control FDR at level alpha
#          est.fdr = the estimated fdr value for each of potential N-hat
#          SS = vector of values for SS loss function evaluated at i/n =0, 1/n, 2/n,..., 1
#          which.sig = vector of indices for maximum ranks declared to be reproducible
MaRR = function(max.rank,
                cutoff = .9,
                alpha = .05,
                khat.to.zero = F) {
  maxx = floor(cutoff * length(max.rank))
  mySS = get.SS(max.rank)
  khat = which(mySS[1:maxx] == min(mySS[1:maxx])) - 1
  if (khat.to.zero == T) {
    khat = 0
  }
  temp.fdr = est.fdr(khat, max.rank)
  Nhat = max(which(temp.fdr <= alpha))
  which.sig = which(max.rank <= Nhat)
  return(list(
    Nhat = Nhat,
    khat = khat,
    est.fdr = temp.fdr,
    SS = mySS,
    which.sig = which.sig
  ))
}


#QL.randomstarts - a function that performs the IDR-procedure (copula mixture model) a set number of times and
#                  reports the results with the highest likelihood
QL.randomstarts = function(signals,
                           starts = 25,
                           omegas = c(.16, .6),
                           means = c(1.5, 3.5),
                           ps = c(0.01, 1))
{
  r.omegas = runif(starts, min = omegas[1], max = omegas[2])
  r.means = runif(starts, min = means[1], max = means[2])
  r.ps = runif(starts, min = ps[1], max = ps[2])
  likes = rep(NA, starts)
  for (i in 1:starts)
  {
    likes[i] = -est.IDR(
      signals,
      mu = r.means[i],
      sigma = 1,
      rho = r.omegas[i],
      p = r.ps[i],
      eps = 0.001,
      max.ite = 30
    )$loglik
  }
  chosen = which(likes == min(likes))
  idr.out = est.IDR(
    signals,
    mu = r.means[chosen],
    sigma = 1,
    rho = r.omegas[chosen],
    p = r.ps[chosen],
    eps = .001,
    max.ite = 30
  )
  # idr.sig = which(idr.out$IDR<=alpha)
  return(idr.out$IDR)
}

## CHMM (repLIS)
bwfw.hmm.Cartesian <- function(x1, x2, pii, A, f0, f1, f2)
{
  ###################################################################
  ## Initialize
  
  NUM <- length(x1)
  
  ## Densities
  
  f0x1 <- dnorm(x1, f0[1], f0[2])
  f0x2 <- dnorm(x2, f0[1], f0[2])
  f1x1 <- dnorm(x1, f1[1], f1[2])
  f1x2 <- dnorm(x2, f2[1], f2[2])
  
  ## the backward-forward procedure
  
  # a. the backward variables
  # --rescaled
  
  alpha <- array(rep(0, NUM * 4), c(NUM, 2, 2))
  # scaling variable c_0
  c0 <- rep(0, NUM)
  
  alpha[1, 1, 1] <- pii[1] * f0x1[1] * f0x2[1]
  alpha[1, 2, 1] <- pii[2] * f1x1[1] * f0x2[1]
  alpha[1, 1, 2] <- pii[3] * f0x1[1] * f1x2[1]
  alpha[1, 2, 2] <- pii[4] * f1x1[1] * f1x2[1]
  
  # rescaling alpha
  c0[1] <-
    1 / (alpha[1, 1, 1] + alpha[1, 2, 1] + alpha[1, 1, 2] + alpha[1, 2, 2])
  for (i in 1:2)
    for (j in 1:2)
      alpha[1, i, j] <- c0[1] * alpha[1, i, j]
  
  for (k in 1:(NUM - 1))
  {
    alpha[k + 1, 1, 1] <- (alpha[k, 1, 1] * A[1, 1] + alpha[k, 2, 1] * A[2, 1] +
                             alpha[k, 1, 2] * A[3, 1] + alpha[k, 2, 2] * A[4, 1]) *
      f0x1[k + 1] * f0x2[k + 1]
    alpha[k + 1, 2, 1] <-
      (alpha[k, 1, 1] * A[1, 2] + alpha[k, 2, 1] * A[2, 2] +
         alpha[k, 1, 2] * A[3, 2] + alpha[k, 2, 2] * A[4, 2]) *
      f1x1[k + 1] * f0x2[k + 1]
    alpha[k + 1, 1, 2] <-
      (alpha[k, 1, 1] * A[1, 3] + alpha[k, 2, 1] * A[2, 3] +
         alpha[k, 1, 2] * A[3, 3] + alpha[k, 2, 2] * A[4, 3]) *
      f0x1[k + 1] * f1x2[k + 1]
    alpha[k + 1, 2, 2] <-
      (alpha[k, 1, 1] * A[1, 4] + alpha[k, 2, 1] * A[2, 4] +
         alpha[k, 1, 2] * A[3, 4] + alpha[k, 2, 2] * A[4, 4]) *
      f1x1[k + 1] * f1x2[k + 1]
    # rescaling alpha
    c0[k + 1] <-
      1 / (alpha[k + 1, 1, 1] + alpha[k + 1, 2, 1] + alpha[k + 1, 1, 2] + alpha[k +
                                                                                  1, 2, 2])
    for (i in 1:2)
      for (j in 1:2)
        alpha[k + 1, i, j] <- c0[k + 1] * alpha[k + 1, i, j]
  }
  
  # b. the forward variables
  # --rescaled
  
  beta <- array(rep(0, NUM * 4), c(NUM, 2, 2))
  
  beta[NUM, 1, 1] <- 1 / 4
  beta[NUM, 2, 1] <- 1 / 4
  beta[NUM, 1, 2] <- 1 / 4
  beta[NUM, 2, 2] <- 1 / 4
  b0 <- rep(0, NUM)
  
  for (k in (NUM - 1):1)
  {
    beta[k, 1, 1] <-
      (
        beta[k + 1, 1, 1] * f0x1[k + 1] * f0x2[k + 1] * A[1, 1] + beta[k + 1, 2, 1] *
          f1x1[k + 1] * f0x2[k + 1] * A[1, 2] +
          beta[k + 1, 1, 2] * f0x1[k + 1] * f1x2[k + 1] * A[1, 3] +
          beta[k + 1, 2, 2] * f1x1[k + 1] * f1x2[k + 1] * A[1, 4]
      )
    beta[k, 2, 1] <-
      (
        beta[k + 1, 1, 1] * f0x1[k + 1] * f0x2[k + 1] * A[2, 1] + beta[k + 1, 2, 1] *
          f1x1[k + 1] * f0x2[k + 1] * A[2, 2] +
          beta[k + 1, 1, 2] * f0x1[k + 1] * f1x2[k + 1] * A[2, 3] +
          beta[k + 1, 2, 2] * f1x1[k + 1] * f1x2[k + 1] * A[2, 4]
      )
    beta[k, 1, 2] <-
      (
        beta[k + 1, 1, 1] * f0x1[k + 1] * f0x2[k + 1] * A[3, 1] + beta[k + 1, 2, 1] *
          f1x1[k + 1] * f0x2[k + 1] * A[3, 2] +
          beta[k + 1, 1, 2] * f0x1[k + 1] * f1x2[k + 1] * A[3, 3] +
          beta[k + 1, 2, 2] * f1x1[k + 1] * f1x2[k + 1] * A[3, 4]
      )
    beta[k, 2, 2] <-
      (
        beta[k + 1, 1, 1] * f0x1[k + 1] * f0x2[k + 1] * A[4, 1] + beta[k + 1, 2, 1] *
          f1x1[k + 1] * f0x2[k + 1] * A[4, 2] +
          beta[k + 1, 1, 2] * f0x1[k + 1] * f1x2[k + 1] * A[4, 3] +
          beta[k + 1, 2, 2] * f1x1[k + 1] * f1x2[k + 1] * A[4, 4]
      )
    # rescaling beta
    # using the same scaling factors as alpha
    b0[k] <-
      1 / (beta[k, 1, 1] + beta[k, 2, 1] + beta[k, 1, 2] + beta[k, 2, 2])
    for (i in 1:2)
      for (j in 1:2)
        beta[k, i, j] <- b0[k] * beta[k, i, j]
  }
  
  # c. lfdr variables
  # --original
  # --the same formula hold for the rescaled alpha and beta
  
  lfdr <- rep(0, NUM)
  
  for (k in 1:NUM)
  {
    q1 <-
      alpha[k, 1, 1] * beta[k, 1, 1] + alpha[k, 2, 1] * beta[k, 2, 1] + alpha[k, 1, 2] *
      beta[k, 1, 2]
    q2 <- alpha[k, 2, 2] * beta[k, 2, 2]
    lfdr[k] <- q1 / (q1 + q2)
  }
  
  lfdr1 <- rep(0, NUM)
  for (k in 1:NUM)
  {
    q1 <- alpha[k, 1, 1] * beta[k, 1, 1] + alpha[k, 1, 2] * beta[k, 1, 2]
    q2 <- alpha[k, 2, 2] * beta[k, 2, 2] + alpha[k, 2, 1] * beta[k, 2, 1]
    lfdr1[k] <- q1 / (q1 + q2)
  }
  
  lfdr2 <- rep(0, NUM)
  for (k in 1:NUM)
  {
    q1 <- alpha[k, 2, 2] * beta[k, 2, 2] + alpha[k, 1, 2] * beta[k, 1, 2]
    q2 <- alpha[k, 1, 1] * beta[k, 1, 1] + alpha[k, 2, 1] * beta[k, 2, 1]
    lfdr2[k] <- q1 / (q1 + q2)
  }
  
  # d. probabilities of hidden states
  # -- and transition variables
  # -- both are rescaled
  
  pdtheta <- array(rep(0, 16 * (NUM - 1)), c((NUM - 1), 2, 2, 2, 2))
  
  b1 <- rep(0, NUM - 1)
  for (k in 1:(NUM - 1))
  {
    pdtheta[k, 1, 1, 1, 1] <-
      alpha[k, 1, 1] * beta[k + 1, 1, 1] * f0x1[k + 1] * f0x2[k + 1] * A[1, 1]
    pdtheta[k, 2, 1, 1, 1] <-
      alpha[k, 2, 1] * beta[k + 1, 1, 1] * f0x1[k + 1] * f0x2[k + 1] * A[2, 1]
    pdtheta[k, 1, 2, 1, 1] <-
      alpha[k, 1, 2] * beta[k + 1, 1, 1] * f0x1[k + 1] * f0x2[k + 1] * A[3, 1]
    pdtheta[k, 2, 2, 1, 1] <-
      alpha[k, 2, 2] * beta[k + 1, 1, 1] * f0x1[k + 1] * f0x2[k + 1] * A[4, 1]
    
    pdtheta[k, 1, 1, 2, 1] <-
      alpha[k, 1, 1] * beta[k + 1, 2, 1] * f1x1[k + 1] * f0x2[k + 1] * A[1, 2]
    pdtheta[k, 2, 1, 2, 1] <-
      alpha[k, 2, 1] * beta[k + 1, 2, 1] * f1x1[k + 1] * f0x2[k + 1] * A[2, 2]
    pdtheta[k, 1, 2, 2, 1] <-
      alpha[k, 1, 2] * beta[k + 1, 2, 1] * f1x1[k + 1] * f0x2[k + 1] * A[3, 2]
    pdtheta[k, 2, 2, 2, 1] <-
      alpha[k, 2, 2] * beta[k + 1, 2, 1] * f1x1[k + 1] * f0x2[k + 1] * A[4, 2]
    
    pdtheta[k, 1, 1, 1, 2] <-
      alpha[k, 1, 1] * beta[k + 1, 1, 2] * f0x1[k + 1] * f1x2[k + 1] * A[1, 3]
    pdtheta[k, 2, 1, 1, 2] <-
      alpha[k, 2, 1] * beta[k + 1, 1, 2] * f0x1[k + 1] * f1x2[k + 1] * A[2, 3]
    pdtheta[k, 1, 2, 1, 2] <-
      alpha[k, 1, 2] * beta[k + 1, 1, 2] * f0x1[k + 1] * f1x2[k + 1] * A[3, 3]
    pdtheta[k, 2, 2, 1, 2] <-
      alpha[k, 2, 2] * beta[k + 1, 1, 2] * f0x1[k + 1] * f1x2[k + 1] * A[4, 3]
    
    pdtheta[k, 1, 1, 2, 2] <-
      alpha[k, 1, 1] * beta[k + 1, 2, 2] * f1x1[k + 1] * f1x2[k + 1] * A[1, 4]
    pdtheta[k, 2, 1, 2, 2] <-
      alpha[k, 2, 1] * beta[k + 1, 2, 2] * f1x1[k + 1] * f1x2[k + 1] * A[2, 4]
    pdtheta[k, 1, 2, 2, 2] <-
      alpha[k, 1, 2] * beta[k + 1, 2, 2] * f1x1[k + 1] * f1x2[k + 1] * A[3, 4]
    pdtheta[k, 2, 2, 2, 2] <-
      alpha[k, 2, 2] * beta[k + 1, 2, 2] * f1x1[k + 1] * f1x2[k + 1] * A[4, 4]
    
    b1[k] <-
      1 / (
        pdtheta[k, 1, 1, 1, 1] + pdtheta[k, 1, 2, 1, 1] + pdtheta[k, 2, 1, 1, 1] +
          pdtheta[k, 2, 2, 1, 1] +
          pdtheta[k, 1, 1, 2, 1] + pdtheta[k, 1, 2, 2, 1] + pdtheta[k, 2, 1, 2, 1] +
          pdtheta[k, 2, 2, 2, 1] +
          pdtheta[k, 1, 1, 1, 2] + pdtheta[k, 1, 2, 1, 2] + pdtheta[k, 2, 1, 1, 2] +
          pdtheta[k, 2, 2, 1, 2] +
          pdtheta[k, 1, 1, 2, 2] + pdtheta[k, 1, 2, 2, 2] + pdtheta[k, 2, 1, 2, 2] +
          pdtheta[k, 2, 2, 2, 2]
      )
    
    for (i in 1:2)
      for (j in 1:2)
        for (p in 1:2)
          for (q in 1:2)
            pdtheta[k, i, j, p, q] <- b1[k] * pdtheta[k, i, j, p, q]
  }
  
  # f. return the results of the bwfw proc.
  bwfw.var <-
    list(
      bw = alpha,
      fw = beta,
      lis = lfdr,
      lfdr1 = lfdr1,
      lfdr2 = lfdr2,
      pdthe = pdtheta,
      c0 = c0
    )
  
  return(bwfw.var)
}

em.hmm.Cartesian <- function(x1, x2, maxiter = 200)
{
  #####################################################################################
  NUM <- length(x1)
  # precision tolerance level
  ptol <- 1e-4
  niter <- 0
  
  ### initializing model parameters
  
  pii.new <- c(1, 0, 0, 0)
  
  A.new <- matrix(c(
    c(0.25, 0.25, 0.25, 0.25),
    c(0.25, 0.25, 0.25, 0.25),
    c(0.25, 0.25, 0.25, 0.25),
    c(0.25, 0.25, 0.25, 0.25)
  ), 4, 4, byrow = TRUE)
  f0 <- c(0, 1)
  f1.new <- c(2, 1)
  f2.new <- c(2, 1)
  
  diff <- 10
  Loglikelihood.new <- -10000
  
  ### The E-M Algorithm
  
  while (diff > ptol && niter < maxiter)
  {
    niter <- niter + 1
    
    pii.old <- pii.new
    A.old <- A.new
    f1.old <- f1.new
    f2.old <- f2.new
    Loglikelihood.old <- Loglikelihood.new
    ## updating the weights and probabilities of hidden states
    
    bwfw.res <-
      bwfw.hmm.Cartesian(x1, x2, pii.old, A.old, f0, f1.old, f2.old)
    
    # the backward-forward variable
    
    alpha <- bwfw.res$bw
    beta <- bwfw.res$fw
    
    # the hidden states probabilities
    gamma <- bwfw.res$pr
    # the transition variables
    dgamma <- bwfw.res$ts
    
    ## updating the parameter estimates
    # a. initial state distribution
    # transform variable of theta
    pdtheta <- bwfw.res$pdthe
    
    # b. transition matrix of Main effect
    
    for (i in 1:2)
    {
      for (j in 1:2)
      {
        for (p in 1:2)
        {
          for (q in 1:2)
          {
            q1 <- sum(pdtheta[, i, j, p, q])
            q2 <-
              sum(pdtheta[, i, j, 1, 1]) + sum(pdtheta[, i, j, 2, 1]) + sum(pdtheta[, i, j, 1, 2]) +
              sum(pdtheta[, i, j, 2, 2])
            if (j == 1)
              index_j <- j - 1
            else
              index_j <- j
            if (q == 1)
              index_q <- q - 1
            else
              index_q <- q
            A.new[i + index_j, p + index_q] <- q1 / q2
          }
        }
      }
    }
    # c. non-null distribution
    
    ptheta1 <- array(rep(0, 2 * NUM), c(NUM, 2))
    ptheta1[, 1] <-
      alpha[, 1, 1] * beta[, 1, 1] + alpha[, 1, 2] * beta[, 1, 2]
    ptheta1[, 2] <-
      alpha[, 2, 2] * beta[, 2, 2] + alpha[, 2, 1] * beta[, 2, 1]
    lf.s1 <- ptheta1[, 2] / (ptheta1[, 1] + ptheta1[, 2])
    q1 <- sum(lf.s1)
    q2 <- sum(lf.s1 * x1)
    mu1 <- q2 / q1
    q3 <- sum(lf.s1 * (x1 - mu1) * (x1 - mu1))
    sd1 <- sqrt(q3 / q1)
    f1.new <- c(mu1, sd1)
    
    ptheta2 <- array(rep(0, 2 * NUM), c(NUM, 2))
    ptheta2[, 1] <-
      alpha[, 1, 1] * beta[, 1, 1] + alpha[, 2, 1] * beta[, 2, 1]
    ptheta2[, 2] <-
      alpha[, 2, 2] * beta[, 2, 2] + alpha[, 1, 2] * beta[, 1, 2]
    lf.s2 <- ptheta2[, 2] / (ptheta2[, 1] + ptheta2[, 2])
    q1 <- sum(lf.s2)
    q2 <- sum(lf.s2 * x2)
    mu2 <- q2 / q1
    q3 <- sum(lf.s2 * (x2 - mu2) * (x2 - mu2))
    sd2 <- sqrt(q3 / q1)
    f2.new <- c(mu2, sd2)
    
    pii.new <- c(1, 0, 0, 0)
    
    c0 <- bwfw.res$c0
    Loglikelihood.new <- -sum(log(c0))
    df1 <- abs(Loglikelihood.old - Loglikelihood.new)
    diff <- df1
  }
  # g. return the results of the E-M algorithm
  em.var <- list(
    pii = pii.new,
    A = A.new,
    f1 = f1.new,
    f2 = f2.new,
    ni = niter
  )
  return (em.var)
}

mt.hmm <- function(lfdr, q)
{
  ## USAGE
  # mt.hmm(lfdr, q)
  
  ## ARGUMENTS
  # lfdr: local false discovery rate sequence
  # q: the FDR level
  
  ## DETAILS
  # mt1.hmm gives a multiple testing rule in a hidden markov model
  # --that controls the FDR at level q, based on sequence of lfdr
  
  ## VALUES
  # nr: the number of rejected hypotheses
  # th: the threshold
  # re: the rejected hypotheses
  # ac: the accepted hypotheses
  # de: the decision rule
  
  m = length(lfdr)
  st.lfdr <- sort(lfdr)
  hps <- rep(0, m)
  if (min(lfdr) > q)
  {
    k <- 0
    threshold <- 1
    reject <- NULL
    accept <- 1:m
  }
  else
  {
    k = 1
    while (k < m && (1 / k) * sum(st.lfdr[1:k]) < q) {
      k = k + 1
    }
    k <- k - 1
    threshold <- st.lfdr[k]
    reject <- which(lfdr <= threshold)
    accept <- which(lfdr > threshold)
    hps[reject] <- 1
  }
  y <- list(
    nr = k,
    th = threshold,
    re = reject,
    ac = accept,
    de = hps
  )
  return (y)
}
