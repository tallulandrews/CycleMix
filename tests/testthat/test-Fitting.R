test_that("checkFit matching distribution", {
  set.seed(3892)
  dat1 <- rnorm(1000, mean=5, sd=sqrt(4))
  dat2 <- rnorm(1000, mean=10, sd=sqrt(4))
  classification <- list(fits=list("G2M"=list(parameters=list(mean=c(5,10), variance=list(sigmasq=c(4,4))),
				   classification=c(rep("A",length(dat1)), rep("B",length(dat2))), 
			data=matrix(c(dat1,dat2), ncol=1))))
  expect_true(checkFit(classification)$good.fit)
})

test_that("checkFit non-matching distribution", {
  set.seed(3892)
  dat1 <- rnorm(1000, mean=3, sd=sqrt(6))
  dat2 <- rnorm(1000, mean=5, sd=sqrt(1))
  classification <- list(fits=list("G2M"=list(parameters=list(mean=c(5,10), variance=list(sigmasq=c(4,4))),
				   classification=c(rep("A",length(dat1)), rep("B",length(dat2))), 
			data=matrix(c(dat1,dat2), ncol=1))))
  expect_false(checkFit(classification)$good.fit)
})
