\docType{data}
\name{GeneSets}
\title{Cell Cycle Gene Sets}
\alias{HGeneSets}
\alias{MGeneSets}
\format{A list of dataframes each with three columns: Gene, Stage, Dir
Genes are provided as HGNC gene symbols
Stage is which cell-cycle stage the gene is differentially expressed in
Dir is the direction of differential expression : 1=upregulated, -1=downregulated
Note: Genes may be differentially expressed in more than one stage and may appear in multiple dataframes
}
\usage{
HGeneSets
MGeneSets
}
\description{
HGeneSets includes:

	Whitfield original cell-cycle genes from Whitfield et al. 2002

	Tirosh modification of Whitfield lists used in Tirosh et al. 2016

	Macosko modification of Whitfield lists used in Macosko et al. 2015

	Quiesc combination of G0 genes from Reactome, and quiescence markers from Cheung et al. 2013

MGeneSets includes:

	Cyclone all genes used by the pre-computed mouse classifier in cyclone (Scialdone et al. 2015). Each gene was 
assigned to the stage and direction for which it exhibited the largest fold-change relative to
the overall mean in our example data.
}
\value{A list of dataframes}
\keyword{datasets}
\references{
	Cheung et al. (2013). Molecular regulation of stem cell quiescence (Table 1). Nat Rev Mol Cell Biol. 14:329-340.
	Macosko et al. (2015) Highly Parallel Genome-wide Expression Profiling of Individual Cells Using nanoliter Droplets (Table S2). Cell. 161(5):1202-1214.
	Scialdone et al. (2015) Computational assignment of cell-cycle stage from single-cell transcriptome data. Methods. 85:54-61
	Tirosh et al. (2016) Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq (Table S5). Science. 352(6282):189-96.
	Whitfield et al. (2002) Identification of genes periodically expressed in the human cell cycle and their expression in tumors. Mol Biol Cell. 13(6):1977-2000.
}

