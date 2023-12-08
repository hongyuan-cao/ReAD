library(radjust)
library(STAREG)
library(ReAD)
source("./R/SimuFunc.R")
source("./R/methods.R")

nrep = 1
J  = 10000
pi = c(0.9,0.025,0.025,0.05)
A00 = 0.8
A11 = A00
A22 = A00
A33 = 0.4
A = matrix(c(A00, (1-A00)/3, (1-A00)/3, (1-A00)/3,
             (1-A11)/3, A11, (1-A11)/3, (1-A11)/3,
             (1-A22)/3, (1-A22)/3, A22, (1-A22)/3,
             (1-A33)/3, (1-A33)/3, (1-A33)/3, A33), 4, 4, byrow=TRUE)
q = 0.05

methods <- c("BH", "MaxP", "MaRR", "radjust", "STAREG", "repLIS", "ReAD")
res <- list()

for (i in 1: length(methods)){
  res[[methods[i]]]$fdr <- c()
  res[[methods[i]]]$power <- c()
}

for (i in 1:nrep){
  cat(paste0("Replication ", i, "\n"))
  
  data.obj <- SimuData_norm_mix(J  = J,
                                pi = pi,
                                A  = A,
                                muA = c(1, 6),
                                prA = c(0.5, 0.5),
                                muB = c(1, 6),
                                prB = c(0.4, 0.6))
  
  pa = data.obj$pa
  pb = data.obj$pb
  theta1 = data.obj$theta1
  theta2 = data.obj$theta2
  truth <- theta1*theta2
  
  # ad hoc BH
  tic <- Sys.time()
  padj1.bh <- p.adjust(pa, method = "BH")
  padj2.bh <- p.adjust(pb, method = "BH")
  toc <- Sys.time()
  toc - tic
  res$BH$fdr[i] <- sum(padj1.bh <= q & padj2.bh <= q & !truth)/max(sum(padj1.bh <= q & padj2.bh <= q), 1)
  res$BH$power[i] <- sum(padj1.bh <= q & padj2.bh <= q & truth) / sum(truth)
  
  # MaxP
  tic <- Sys.time()
  maxp <- pmax(pa, pb)
  padj.maxp <- p.adjust(maxp, method = "BH")
  toc <- Sys.time()
  toc - tic
  res$MaxP$fdr[i] <- sum(padj.maxp <= q & !truth)/ max(sum(padj.maxp <= q), 1)
  res$MaxP$power[i] <- sum(padj.maxp <= q & truth) / sum(truth)
  
  # MaRR
  x <- rank(pa)
  y <- rank(pb)
  max.rank = apply(cbind(x,y),1,max)
  res.marr <- MaRR(max.rank, alpha = q)
  marr.rej <- rep(0, J)
  marr.rej[res.marr$which.sig] = 1
  res$MaRR$fdr[i] <- sum(marr.rej & !truth)/ max(sum(marr.rej), 1)
  res$MaRR$power[i] <- sum(marr.rej & truth) / sum(truth)
  
  # radjust Bogomolov & Heller 2018 (Biometrika)
  library(radjust)
  p1 = pa
  p2 = pb
  p1[which(p1==0)] <- 1e-15
  p2[which(p2==0)] <- 1e-15
  res.rv18 <- radjust_sym(p1, p2, input_type = "all", directional_rep_claim = FALSE,
                          variant = "adaptive", alpha = q)
  rv18 <- rep(1, J)
  rv18[as.numeric(res.rv18$results_table$name)] <- res.rv18$results_table$r_value
  res$radjust$fdr[i] <- sum(rv18 <= q & !truth)/max(sum(rv18 <= q), 1)
  res$radjust$power[i] <- sum(rv18 <= q & truth) / sum(truth)
  
  # STAREG
  res.eb <- stareg(pa, pb)
  padj.eb <- res.eb$fdr
  res$STAREG$fdr[i] <- sum(padj.eb <= q & !truth)/max(sum(padj.eb <= q), 1)
  res$STAREG$power[i] <- sum(padj.eb <= q & truth) / sum(truth)
  
  # repLIS (Wang and Zhu, 2019)
  x1 <- qnorm(p1,lower.tail = FALSE)
  x2 <- qnorm(p2,lower.tail = FALSE)
  f0 <- c(0, 1)
  em.var <- em.hmm.Cartesian(x1, x2)
  bwfw.var <- bwfw.hmm.Cartesian(x1, x2, em.var$pii, em.var$A, f0, em.var$f1, em.var$f2)
  lis <- bwfw.var$lis
  res.replis <- mt.hmm(lis, q)
  res$repLIS$fdr[i] <- sum(res.replis$de & !truth)/max(sum(res.replis$de), 1)
  res$repLIS$power[i] <- sum(res.replis$de & truth) / sum(truth)
  
  # ReAD
  library(ReAD)
  res.hmm <- ReAD(pa, pb)
  padj.hmm <- res.hmm$fdr
  res$ReAD$fdr[i] <- sum(padj.hmm <= q & !truth)/max(sum(padj.hmm <= q), 1)
  res$ReAD$power[i] <- sum(padj.hmm <= q & truth) / sum(truth)
}

for (k in 1:length(methods)){
  res[[k]]$fdr <- mean(res[[k]]$fdr)
  res[[k]]$power <- mean(res[[k]]$power)
}