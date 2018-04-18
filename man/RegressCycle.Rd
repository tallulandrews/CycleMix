\name{regressCycle}
\alias{regressCycle}
\title{Regress out cell-cycle}
\description{
A wrapper for scater's normalizeExprs function to regress out the cell-cycle from a dataset.
}
\usage{
  regressCycle(SCE, classification, expr_name="logcounts", method=c("scores", "phase"))
}
\arguments{
  \item{SCE}{a SingleCellExperiment object.}
  \item{classification}{output from \code{classifyCells}.}
  \item{expr_name}{name of the expression matrix to use, must be an element of assays(SCE).}
  \item{method}{which factor to regress on.}
}
\details{
	Creates a design matrix based on either the phase-signatures ("scores") or the 
	phase assignments ("phase"). Then runs \code{normalizeExprs} from the scater package to 
	regress out the cell-cycle.

}
\value{ 
	The input SingleCellExperiment object with the regressed expression saved to the 'norm_exprs' slot.
}
\examples{
	classification <- classifyCells(Ex,  MGeneSets$Cyclone, do.scale=TRUE, allow.multi=FALSE);
	Ex <- regressCycle(Ex, classification, method="scores")
	Ex <- regressCycle(Ex, classification, method="phase")
}
\keyword{cell cycle, normalization}
\seealso{normalizeExprs, classifyCells}
