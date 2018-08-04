\name{regressCycle}
\alias{regressCycle_scater}
\alias{regressCycle_partial}
\title{Regress out cell-cycle}
\description{
A wrapper for scater's normalizeExprs function to regress out the cell-cycle from a dataset. Or for glm to regress out differences between specific cell-cycle phases. 
}
\usage{
  regressCycle_partial(SCE, classification, expr_name="logcounts", method=c("scores", "phase"), phases=c("G2M", "G1S"))
  regressCycle_scater(SCE, classification, expr_name="logcounts", method=c("scores", "phase"))
}
\arguments{
  \item{SCE}{a SingleCellExperiment object.}
  \item{classification}{output from \code{classifyCells}.}
  \item{expr_name}{name of the expression matrix to use, must be an element of assays(SCE).}
  \item{method}{which factor to regress on.}
  \item{phases}{which phases to regress out.}
}
\details{
	\code{regressCycle_scater} Creates a design matrix based on either the phase-signatures ("scores") or the 
	phase assignments ("phase"). Then runs \code{normalizeExprs} from the scater package to 
	regress out the cell-cycle.

	\code{regressCycle_partial} Creates a design matrix based on either the phase-signatures ("scores") or the
        phase assignments ("phase"). Uses a GLM to remove differences between the provided set of phases. 

}
\value{ 
	\code{regressCycle_scater} The input SingleCellExperiment object with the regressed expression saved to the 'norm_exprs' slot.
	\code{regressCycle_partial} Uses a glm to regress out differences among a select set of cell cycle phases
}
\examples{
	classification <- classifyCells(Ex,  MGeneSets$Cyclone, do.scale=TRUE, allow.multi=FALSE);
	Ex <- regressCycle_partial(Ex, classification, method="scores", phases=c("G1", "S", "G2M"))
	Ex <- regressCycle_partial(Ex, classification, method="phase", phases=c("G1", "S", "G2M"))
	Ex <- regressCycle_scater(Ex, classification, method="phase")
}
\keyword{cell cycle, normalization}
\seealso{normalizeExprs, classifyCells}
