\name{regressCyclePartial}
\alias{regressCyclePartial}
\title{Regress out cell-cycle}
\description{
Uses a GLM to remove differences between specified cell-cycle phases.
}
\usage{
  knnSmooth(classification, dims=NULL, clusters=NULL, k=20)
}
\arguments{
  \item{classification}{output from \code{classifyCells}.}
  \item{dims}{a dataframe or matrix of the lower dimensional space to be used for calculating nearest neighbours. Rows are cells in the same order as in classification, columns are different dimension (should be PCs).}
  \item{clusters}{a vector of cell-type IDs or other clusters.}
  \item{k}{number of nearest neighbours to consider, should be no higher than 1 percent of the total dataset}
}
\details{
	Uses RcppAnnoy to identify approximate nearest neighbours based on the euclidean distance in the provided lower dimensional space. Then calculates the consensus (most common) cell-cycle phase among each cell's neighbours and reassigns cells on this basis. Preferentially reassigns cells initially assigned to "None" in order to reduce conservativeness. 

	If clusters are provided then smoothing will instead be done at the cluster level and each cluster (and all cells within) will be assigned to the most frequent cell-cycle phase within it.

}
\value{ 
	the input classification results updated with the smoothed classifications.
}
\examples{
	require(SingleCellExperiment)
	classification <- classifyCells(Ex,  MGeneSets$Cyclone, do.scale=TRUE, allow.multi=FALSE, symbol_column="feature_symbol");
	pca <- matrix(rnorm(ncol(Ex)*10), ncol=10)
	classification <- knnSmooth(classification, dims=pca, k=20)
	classification <- knnSmooth(classification, clusters=Ex$cell_type1) 
}
\keyword{cell cycle, smoothing}
\seealso{classifyCells}
