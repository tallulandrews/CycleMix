\name{plotMixture}
\alias{plotMixture}
\title{Plot Mixture Model}
\description{
	Plots output from a single mixture model.
}
\usage{
  plotMixture(fit, BIC=FALSE)
}
\arguments{
  \item{fit}{output from \code{assignPhase}}
  \item{BIC}{whether to also plot the BICs for all the models considered.}
}
\details{
	Plots a histogram of the gene expression signature of a single cell-cycle phase across cells. Superimposes the fitted gaussian curves for the optimal model according to BIC. Optionally plots the BIC for all the models considered in the fitting.

}
\value{ 
	The input SingleCellExperiment object with the regressed expression saved to the 'norm_exprs' slot.
}
\examples{
	fit <- assignPhase(Ex, MGeneSets$Cyclone, "G2M")
	plotMixture(fit, BIC=TRUE)

	classification <- classifyCells(Ex,  MGeneSets$Cyclone, do.scale=TRUE, allow.multi=FALSE);
	plotMixture(classification$fits[[1]], BIC=FALSE)
}
\keyword{cell cycle, mixture model, plotting}
\seealso{classifyCells, assignPhase}
