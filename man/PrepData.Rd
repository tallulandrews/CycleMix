\name{prepData}
\alias{prepData}
\title{Prepare Data}
\description{
	Extract the normalized expression matrix from a variety of objects used in scRNAseq.
}
\usage{
  prepData(obj, expr_name="logcounts", do.scale=FALSE, symbol_column="feature_symbol")
}
\arguments{
  \item{obj}{a SingleCellExperiment object or a Seurat object or a matrix/Matrix.}
  \item{expr_name}{name of the expression matrix to use, must be an element of assays(obj), should be at minimum library-size normalized and log-transformed. If obj is a Seurat object the "assay" from which to obtain the "data" matrix from.}
  \item{do.scale}{whether to mean-center and scale variance of each gene.}
  \item{symbol_column}{name of the column of rowData(obj) which contains HGNC gene symbols, if set to NULL uses rownames.}
}
\details{ 
Extracts the log normalized expression matrix from the provided object - currently supports Seurat version 5, and SingleCellExperiment objects. \code{expr_name} should be the assay name for the object, e.g. "logcounts" for SingleCellExperiment or "RNA" for Seurat. 

symbol_column is either the column of the rowData (SingleCellExperiment) containing the appropriate gene symbols, or should be set to NULL to use rownames, or can be a vector of gene symbols. Duplicates are removed and these are added as the rownames of the normalized expression matrix.
}
\value{ 
	a "matrix" or "Matrix" of the specified log-normalized expression values with the specified gene symbols as the rownames.
}
\examples{
	expr_mat <- prepData(Ex) # Ex = SingleCellExperiment Object
}
