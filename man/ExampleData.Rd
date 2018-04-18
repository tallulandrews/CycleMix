\docType{data}
\name{Ex}
\title{Example scRNASeq Data}
\alias{Ex}
\format{ A SingleCellExperiment object containing 2002 genes and 273 cells
}
\usage{Ex}
\description{
	Filtered gene expression matrix of cell-cycle staged mESCs from Buettner et al. (2015). 
Cells were filtered to remove those with fewer than 6000 detected genes. Genes where filtered to
contain only the significant features (5\% FDR) using M3Drop and all the cell-cycle-related genes. 

	Known cell-cycle stage labels are contained in the "cell_type1" column of the colData. Gene symbols are contained in the "feature_symbols" column of the rowData. 

	Contains both counts and log2-CPM ("logcounts") matrices.
}
\value{A SingleCellExperiment object}
\keyword{datasets}
\references{
Buettner et al. (2015). Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequening data reveals hidden subpopulation of cells. Nature Biotechnology. 33:155-160.
}

