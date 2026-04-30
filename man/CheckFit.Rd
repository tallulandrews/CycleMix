\name{checkFit}
\alias{checkFit}
\title{Check Goodness of Fit of the GMM}
\description{
 Checks the goodness of fit of the gaussian mixture model using a bootstrap appraoch.
}
\usage{
	checkFit(classification, nbootstrap=100, summarize=FALSE, D.threshold=0.1) 
}
\arguments{
  \item{classification}{output from \code{classifyCells}.}
  \item{nbootstrap}{number of random samples from the GMM to perform}
  \item{summarize}{whether to return all results of the random samples or just the average.}
  \item{D.threshold}{the threshold for average D statistic to consider the GMM a good fit to the data.}
}
\details{
	Randomly samples from the fitted GMM then performs a KS test to assess the discrepancy between the fit distribution and the observed distribution of phase scores for each phase.

}
\value{ 
	if summarize = TRUE, returns a list of the average D, pval, and whether this is considered a good fit for each phase.
	if summarize = FALSE, returns a list of lists of output from the ks.test for every boostrap sampling for every phase.
}
\examples{
	require(SingleCellExperiment)
	classification <- classifyCells(Ex,  MGeneSets$Cyclone, do.scale=TRUE, allow.multi=FALSE, symbol_column="feature_symbol");
	good.fit.allres <- checkFit(classification, 100, summarize=FALSE)
	good.fit <- checkFit(classification, 100, summarize=TRUE)
}
\keyword{goodness of fit}
\seealso{classifyCells}
