\name{regressCyclePartial}
\alias{regressCyclePartial}
\title{Regress out cell-cycle}
\description{
Uses a GLM to remove differences between specified cell-cycle phases.
}
\usage{
  regressCyclePartial(expr_mat, classification, type=c("counts","norm"), method=c("scores", "phase"), phases=c("G2M", "G1S"), allow_negative=FALSE, subsample_cells=ncol(expr_mat))
}
\arguments{
  \item{expr_mat}{a raw or normalized expression matrix to regress}
  \item{classification}{output from \code{classifyCells}.}
  \item{type}{whether the expr_mat is raw counts (uses NB-GLM) or normalized expression (uses Gaussian-GLM).}
  \item{method}{which factor to regress on.}
  \item{phases}{which phases to regress out.}
  \item{allow_negative}{whether or not to allow corrected expression values to be negative.}
}
\details{
	\code{regressCyclePartial} Creates a design matrix based on either the phase-signatures ("scores") or the
        phase assignments ("phase"). Uses a GLM to remove differences between the provided set of phases. 

}
\value{ 
	\code{regressCyclePartial} Uses a glm to regress out differences among a select set of cell cycle phases
}
\examples{
	require(SingleCellExperiment)
	classification <- classifyCells(Ex,  MGeneSets$Cyclone, do.scale=TRUE, allow.multi=FALSE);
	Ex <- regressCyclePartial(counts(Ex), classification, type="counts", method="scores", phases=c("S", "G2M"), allow_negative=TRUE, subsample_cells=200)
	Ex <- regressCyclePartial(counts(Ex), classification, type="counts", method="phase", phases=c("S", "G2M"), allow_negative=TRUE, subsample_cells=200)
	regressed_mat1 <- regressCyclePartial(logcounts(Ex), classification, type="norm", method="scores", phases=c("G1", "S", "G2M"), allow_negative=FALSE, subsample_cells=200)
	regressed_mat2 <- regressCyclePartial(logcounts(Ex), classification, type="norm", method="phase", phases=c("G1", "S", "G2M"), allow_negative=TRUE, subsample_cells=200)
}
\keyword{cell cycle, normalization}
\seealso{classifyCells}
