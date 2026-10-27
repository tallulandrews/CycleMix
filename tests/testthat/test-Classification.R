test_that("assignPhase works", {
  set.seed(182)
  truth <- rep(c("None", "G1S", "G2M"), each=300)
  expr1 <- rnorm(300, mean=1)
  expr2 <- rnorm(300, mean=4)
  expr3 <- rnorm(300, mean=10)
  expr_mat <- rbind(c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2),
  			c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2),
  			c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2),
  			c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2),
  			c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2),
  			c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2),
  			c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2))
  colnames(expr_mat) <- paste0("cell",1:ncol(expr_mat))
  rownames(expr_mat) <- paste0("gene",1:nrow(expr_mat))
 
  CC_tab <- data.frame(Gene=rownames(expr_mat)[-seq(from=3, to=nrow(expr_mat), by=3)], 
			Stage=rep(c("G2M", "G1S"), nrow(expr_mat)/3), Dir=1) 
  out <- assignPhase(expr_mat, CC_tab, phase="G2M")
  truth[truth != "G2M"] <- ""
  expect_identical(out$phase, truth)
})

test_that("classifyCells works", {
  set.seed(182)
  truth <- rep(c("None", "G1S", "G2M"), each=300)
  expr1 <- rnorm(300, mean=1)
  expr2 <- rnorm(300, mean=4)
  expr3 <- rnorm(300, mean=10)
  expr_mat <- rbind(c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2),
  			c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2),
  			c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2),
  			c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2),
  			c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2),
  			c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2),
  			c(expr1,expr2,expr3),c(expr2,expr3,expr1), c(expr3,expr1,expr2))
  colnames(expr_mat) <- paste0("cell",1:ncol(expr_mat))
  rownames(expr_mat) <- paste0("gene",1:nrow(expr_mat))
 
  CC_tab <- data.frame(Gene=rownames(expr_mat)[-seq(from=3, to=nrow(expr_mat), by=3)], 
			Stage=rep(c("G2M", "G1S"), nrow(expr_mat)/3), Dir=1) 
  out <- classifyCells(expr_mat, CC_tab)
  expect_identical(out$phase, truth)
})

