test_that("Reads SCE", {
  set.seed(3892)
  counts <- matrix(rpois(2000, lambda=2),ncol=20);
  obj <- SingleCellExperiment::SingleCellExperiment(counts)
  SummarizedExperiment::assay(obj, "logcounts") <- counts
  expect_identical(prepData(obj), counts)
})

test_that("Reads Seurat", {
  set.seed(3892)
  counts <- matrix(rpois(2000, lambda=2),ncol=20);
  rownames(counts) <- paste0("Gene", 1:nrow(counts))
  obj <- Seurat::CreateSeuratObject(counts=counts)
  obj <- Seurat::NormalizeData(obj)
  expect_identical(prepData(obj, expr_name="RNA"), Seurat::GetAssayData(obj, assay="RNA",layer="data"))
  TRUE
})

test_that("Reads Matrix", {
  set.seed(3892)
  counts <- matrix(rpois(2000, lambda=2),ncol=20);
  counts <- Matrix::Matrix(counts)
  expect_identical(prepData(counts), counts)
})

test_that("Reads matrix", {
  set.seed(3892)
  counts <- matrix(rpois(2000, lambda=2),ncol=20);
  expect_identical(prepData(counts), counts)
})
