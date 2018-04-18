\name{assignPhase}
\alias{assignPhase}
\alias{classifyCells}
\title{Classify Cells into Cell Cycle Stages}
\description{
	Uses mixed-gaussian distributions to assign cells to different cell-cycle stages according to provided marker/DE genes.
}
\usage{
  assignPhase(SCE, CC_table, phase="G2M", expr_name="logcounts", do.scale=FALSE, symbol_column="feature_symbol")
  classifyCells(SCE, CC_table, expr_name="logcounts", do.scale=FALSE, symbol_column="feature_symbol", allow.multi=FALSE)
}
\arguments{
  \item{SCE}{a SingleCellExperiment object.}
  \item{CC_table}{a dataframe with three columns: Gene, Stage, Dir}
  \item{phase}{which cell-cycle stage to assign cells to, must be a value found in "Stage" of CC_table}
  \item{expr_name}{name of the expression matrix to use, must be an element of assays(SCE), should be at minimum library-size normalized and log-transformed.}
  \item(do.scale}{mean-center and scale variance of each gene.}
  \item(symbol_column){name of the column of rowData(SCE) which contains HGNC gene symbols, if null uses rownames.}
  \item{allow.multi}{allow cells to be assigned to multiple cell-cycle phases.}
}
\details{ Classifies cells as being in a particular phase of the cell-cycle by modelling the 
expression signature of each stage as a mixture of guassian distributions. The gene expression 
signature of a particular stage as the weighted mean expression of the stage-specific genes 
provided in the \code{CC_table} Weights are provided in the "Dir" column of the \code{CC_table}.

Each stage specific signature is fit as a mixture of 1-3 guassian distributions and BIC is used to select the best model using \code{mclust}. If the best model is a mixture of >1 distributions then cells assigned to the distribution with the highest mean are assigned to the respective 
cell-cycle phase. 

If a cell is assigned to multiple cell-cycle phases and \code{allow.multi} is FALSE then 
signatures for each phase are centered and scaled, cells are assigned to the phase for which it
has the highest signature. 

\code{assignPhase} performs the fitting and assigns cells to a single cell-cycle phase
\code{classifyCells} runs assignPhase for all phases in the CC_table and combines the output into a single assignment vector.
}
\value{ 
  \code{assignPhase} A list of output from \code{mclust} with a extra element, "phase", which is a vector of phase assignments for each cell.
  \code{classifyCells} A list of three elements:

	"phase" : a vector of phase assignments for each cell.

	"scores" : a table of expression signatures for each cell to each phase.

	"fits" : a list of output from assignPhase.
}
\examples{
	fit <- assignPhase(Ex, MGeneSets$Cyclone, phase="G2M")
	classification <- classifyCells(Ex,  MGeneSets$Cyclone, do.scale=TRUE, allow.multi=FALSE);
	classification <- classifyCells(Ex,  MGeneSets$Cyclone, do.scale=FALSE, allow.multi=TRUE);
}
\keyword{cell cycle, mixture models}
