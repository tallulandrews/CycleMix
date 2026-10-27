test_that("Discrete All Phases Regression Works", {
  set.seed(3892)
  classification <- list(phase = rep(c("G1S", "G2M", "None"),each=100))
  # discrete gaussian expression
  expr1 <- rnorm(100,mean=5)
  expr2 <- rnorm(100, mean=7)
  expr3 <- rnorm(100, mean=10)

  expr_mat <- rbind(c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2))
  colnames(expr_mat) <- paste0("cell", 1:300)

  out <- regressCyclePartial(expr_mat, classification, type="norm", method="phase", phases=c("None", "G1S","G2M"), subsample_cells=1000, allow_negative=TRUE)

  expect_equal(round(apply(out, 1, mean), digits=0), c(10,5,7))
})

test_that("Discrete Partial Phases Regression Works", {
  set.seed(3892)
  classification <- list(phase = rep(c("G1S", "G2M", "None"),each=100))
  # discrete gaussian expression
  expr1 <- rnorm(100,mean=5)
  expr2 <- rnorm(100, mean=7)
  expr3 <- rnorm(100, mean=10)

  expr_mat <- rbind(c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2))
  colnames(expr_mat) <- paste0("cell", 1:300)

  out <- regressCyclePartial(expr_mat, classification, type="norm", method="phase", phases=c( "G1S","G2M"), subsample_cells=1000, allow_negative=TRUE)

  expect_equal(round(apply(out[,1:200], 1, mean), digits=0), c(5,7,10))
})

test_that("Continuous Full Phases Regression Works", {
  set.seed(3892)
  G1S_scores <- 1:300 * (1+runif(300)/50)
  G2M_scores <- c(1:150,1:150) * (1+runif(300)/50)
  classification <- list(scores=cbind(G1S_scores,G2M_scores))
  colnames(classification$scores) <- c("G1S", "G2M")
  # continuous gaussian expression
  expr1 <- G1S_scores*0.3 + rnorm(300, mean=5, sd=1)
  expr2 <- G2M_scores*1.5 - rnorm(300, mean=6, sd=1)
  expr3 <- G2M_scores*0.7 - rnorm(300, mean=1, sd=1)

  expr_mat <- rbind(expr1, expr2, expr3)
  colnames(expr_mat) <- paste0("cell", 1:300)

  out <- regressCyclePartial(expr_mat, classification, type="norm", method="scores", phases=c( "G1S","G2M"), subsample_cells=1000, allow_negative=TRUE)

  stat <- round(cor(out[1,], G1S_scores)+cor(out[2,], G2M_scores)+cor(out[3,], G2M_scores), digits=5)
  expect_equal(stat, 0)
})

test_that("Continuous Partial Phases Regression Works", {
  set.seed(3892)
  G1S_scores <- 1:300 * (1+runif(300)/50)
  G2M_scores <- c(1:150,1:150) * (1+runif(300)/50)
  classification <- list(scores=cbind(G1S_scores,G2M_scores))
  colnames(classification$scores) <- c("G1S", "G2M")
  # continuous gaussian expression
  expr1 <- G1S_scores*0.3 + rnorm(300, mean=5, sd=1)
  expr2 <- G2M_scores*1.5 - rnorm(300, mean=6, sd=1)
  expr3 <- G2M_scores*0.7 - rnorm(300, mean=1, sd=1)

  expr_mat <- rbind(expr1, expr2, expr3)
  colnames(expr_mat) <- paste0("cell", 1:300)

  out <- regressCyclePartial(expr_mat, classification, type="norm", method="scores", phases=c( "G2M"), subsample_cells=1000, allow_negative=TRUE)

  stat <- round(cor(out[1,], G1S_scores)+cor(out[2,], G2M_scores)+cor(out[3,], G2M_scores), digits=0)
  expect_equal(stat, 1)
})

## Counts - Fixed

test_that("Discrete All Phases Counts Regression Works", {
  set.seed(2891)
  classification <- list(phase = factor(rep(c("G1S", "G2M", "None"),each=100), levels=c("None", "G1S", "G2M")))
  # discrete negative binomial expression
  expr1 <- rnbinom(100, mu=1, size=100)
  expr2 <- rnbinom(100, mu=5, size=100)
  expr3 <- rnbinom(100, mu=50, size=100)

  expr_mat <- rbind(c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2))
  colnames(expr_mat) <- paste0("cell", 1:300)

  out <- regressCyclePartial(expr_mat, classification, type="counts", method="phase", phases=c("None", "G1S","G2M"), subsample_cells=1000, allow_negative=TRUE)

  expect_equal(round(apply(out, 1, mean), digits=0), c(50,1,5))
})

test_that("Discrete Partial Phases Counts Regression Works", {
  set.seed(2891)
  classification <- list(phase = factor(rep(c("G1S", "G2M", "None"),each=100), levels=c("None", "G1S", "G2M")))
  # discrete negative binomial expression
  expr1 <- rnbinom(100, mu=1, size=100)
  expr2 <- rnbinom(100, mu=5, size=100)
  expr3 <- rnbinom(100, mu=50, size=100)

  expr_mat <- rbind(c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2))
  colnames(expr_mat) <- paste0("cell", 1:300)

  out <- regressCyclePartial(expr_mat, classification, type="norm", method="phase", phases=c( "G1S","G2M"), subsample_cells=1000, allow_negative=TRUE)

  expect_equal(round(apply(out[,1:200], 1, mean), digits=0), c(1,5,50))
})

test_that("Continuous Full Phases Counts Regression Works", {
  set.seed(2891)
  G1S_scores <- 1:300 * (1+runif(300)/50)
  G2M_scores <- c(1:150,1:150) * (1+runif(300)/50)
  classification <- list(scores=cbind(G1S_scores,G2M_scores))
  colnames(classification$scores) <- c("G1S", "G2M")
  # continuous negative binomial expression
  expr1 <- rnbinom(300, mu=G1S_scores*0.3, size=100)
  expr2 <- rnbinom(300, mu=G2M_scores*2.5, size=100)
  expr3 <- rnbinom(300, mu=G2M_scores*1.2, size=100)

  expr_mat <- rbind(expr1, expr2, expr3)
  colnames(expr_mat) <- paste0("cell", 1:300)

  out <- regressCyclePartial(expr_mat, classification, type="norm", method="scores", phases=c( "G1S","G2M"), subsample_cells=1000, allow_negative=TRUE)

  stat <- round(cor(out[1,], G1S_scores)+cor(out[2,], G2M_scores)+cor(out[3,], G2M_scores), digits=5)
  expect_equal(stat, 0)
})

test_that("Continuous Partial Phases Counts Regression Works", {
  set.seed(2891)
  G1S_scores <- 1:300 * (1+runif(300)/50)
  G2M_scores <- c(1:150,1:150) * (1+runif(300)/50)
  classification <- list(scores=cbind(G1S_scores,G2M_scores))
  colnames(classification$scores) <- c("G1S", "G2M")
  # continuous negative binomial expression
  expr1 <- rnbinom(300, mu=G1S_scores*0.3, size=100)
  expr2 <- rnbinom(300, mu=G2M_scores*2.5, size=100)
  expr3 <- rnbinom(300, mu=G2M_scores*1.2, size=100)

  expr_mat <- rbind(expr1, expr2, expr3)
  colnames(expr_mat) <- paste0("cell", 1:300)

  out <- regressCyclePartial(expr_mat, classification, type="norm", method="scores", phases=c( "G2M"), subsample_cells=1000, allow_negative=TRUE)

  stat <- round(cor(out[1,], G1S_scores)+cor(out[2,], G2M_scores)+cor(out[3,], G2M_scores), digits=0)
  expect_equal(stat, 1)
})



