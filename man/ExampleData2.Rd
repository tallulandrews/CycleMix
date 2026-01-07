\docType{data}
\name{Ex2.2}
\title{Example scRNASeq Data}
\alias{Ex2.2}
\format{ A list containing "counts" the Matrix of gene expression for 6625 genes and 1500 cells, and data.frame of metadata as "meta".
}
\usage{Ex2.2}
\description{
	Filtered gene expression matrix of 10X Chromium data from Wang et al. (2020) J Exp Med
. 2020 Feb 3;217(2):e20191130. doi: 10.1084/jem.20191130.

Cells were randomly downsampled to 1500 total and genes where filtered to remove those detected in fewer than 100 cells. Cell-type and tissue annotations are provided from the original publication.

}
\value{A list containing "counts" and "meta" of gene expression counts and metadata}
\keyword{datasets}
\references{
 Wang et al. (2020) J Exp Med. 2020 Feb 3;217(2):e20191130. doi: 10.1084/jem.20191130
}

