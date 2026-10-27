test_that("knnSmooth all proliferating", {
  set.seed(191)
  coords <- cbind(runif(300), runif(300))
  classification <- list(phase = rep(c("G1S", "G2M", "None"),100))
  smoothed <- knnSmooth(classification, dims=coords, k=20, min.threshold=0.1)			
  expect_equal(sum(smoothed$phase == "None"), 0)
})

test_that("knnSmooth no change", {
  set.seed(191)
  coords <- cbind(runif(300), runif(300))
  classification <- list(phase = rep(c("G1S", "G2M", "None"),100))
  smoothed <- knnSmooth(classification, dims=coords, clusters=1:300, min.threshold=0.9)			
  expect_true(identical(smoothed$phase, classification$phase))
})
