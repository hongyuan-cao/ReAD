# ReAD
 Powerful replicability analysis accounting for local dependence structure

## Overview

A robust and powerful approach is developed for replicability analysis of two Genome-wide association studies (GWASs) accounting for the linkage disequilibrium (LD) among genetic variants. The LD structure in two GWASs is captured by a four-state hidden Markov model (HMM). The unknowns involved in the HMM are estimated by an efficient expectation-maximization (EM) algorithm in combination with a non-parametric estimation of functions. By incorporating information from adjacent locations via the HMM, this approach identifies the entire clusters of genotype-phenotype associated signals, improving the power of replicability analysis while effectively controlling the false discovery rate.

## Installation

```R
## Install dependency packages if necessary
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")


## Install ReAD
install.packages("devtools")
devtools::install_github("YanLi15/ReAD")
# An alternative installation method from CRAN
install.packages("ReAD")

## Load ReAD
library(ReAD)
```

## An numerical example

We illustrate with simulated data the application of ReAD to the replicability analysis of two large-scale multiple testing problems, both of which consider local dependency structures based on Hidden Markov Models.

```R
## Pre-specify the number of hypotheses, m, the prior probabilities of the joint hidden states, pi, the transition matrix, A, and the alternative settings
m = 10000; pi = c(0.9, 0.025, 0.025, 0.05); A = 0.6 * diag(4) + 0.1;
mu1 = 2; mu2 = 2.5; sigma1 = 1; sigma2 = 1

## Generate the hidden states and corresponding p-values in two studies based on a four-state hidden Markov model
s <- c()
s[1] <- sample(0:3, 1, prob = pi)
for (j in 2:J){
  s[j] <- sample(0:3, 1, prob = A[s[j-1]+1,])
}

states1 = rep(0, J)
states1[c(which(s == 2), which(s == 3))] = 1
states2 = rep(0, J)
states2[c(which(s == 1), which(s == 3))] = 1

xa <- rnorm(J, mean = muA * states1, sd = sdA)
xb <- rnorm(J, mean = muB * states2, sd = sdB)

pa <- 1 - pnorm(xa)
pb <- 1 - pnorm(xb)

## Replicability analysis
library(ReAD)
alpha = 0.05
rep.obj <- ReAD(pa, pb)
rep.snps <- which(rep.obj$fdr <= alpha)
```

