\name{regressCycle}
\alias{regressCycleScater}
\alias{regressCyclePartial}
\title{Regress out cell-cycle}
\description{
A wrapper for scater's normalizeExprs function to regress out the cell-cycle from a dataset. Or for glm to regress out differences between specific cell-cycle phases. 
}
\usage{
  regressCyclePartial( expr_mat, classification, expr_name="logcounts", method=c("scores", "phase"), phases=c("G2M", "G1S"), allow_negative=FALSE)
}
\arguments{
  \item{expr_mat}{a raw or normalized expression matrix to regress}
  \item{classification}{output from \code{classifyCells}.}
  \item{expr_name}{name of the expression matrix to use, must be an element of assays(SCE).}
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
	Ex <- regressCyclePartial(Ex, classification, method="scores", phases=c("G1", "S", "G2M"))
	Ex <- regressCyclePartial(Ex, classification, method="phase", phases=c("G1", "S", "G2M"))
	regressed_mat1 <- regressCyclePartial(logcounts(Ex), classification, method="scores", phases=c("G1", "S", "G2M"))
	regressed_mat2 <- regressCyclePartial(logcounts(Ex), classification, method="phase", phases=c("G1", "S", "G2M"))
}
\keyword{cell cycle, normalization}
\seealso{classifyCells}
