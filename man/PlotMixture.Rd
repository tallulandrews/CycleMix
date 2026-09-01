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
  \item{fit}{output from \code{assignPhase} or the "fits" slot from the output from classifyCells}
  \item{BIC}{whether to also plot the BICs for all the models considered.}
}
\details{
	Plots a histogram of the gene expression signature of a single cell-cycle phase across cells. Superimposes the fitted gaussian curves for the optimal model according to BIC. Optionally plots the BIC for all the models considered in the fitting.

}
\value{ 
	No output.
}
\examples{
	data(ileum_10X_example)
	data(esc_Smartseq_example)
	data(MGeneSets)
	data(HGeneSets)
	fit <- assignPhase(ileum_10X_example$counts, HGeneSets$Tirosh, "G2M")
	plotMixture(fit, BIC=TRUE)

	classification <- classifyCells(esc_Smartseq_example,  MGeneSets$Cyclone, do.scale=TRUE, allow.multi=FALSE, symbol_column="feature_symbol")

	plotMixture(classification$fits[[1]], BIC=FALSE)
}
\keyword{cell cycle, mixture model, plotting}
\seealso{classifyCells, assignPhase}
