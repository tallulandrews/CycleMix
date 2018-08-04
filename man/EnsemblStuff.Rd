\name{mapGeneNames}
\alias{mapGeneNames}
\alias{downloadEnsemblData}
\title{Map gene names}
\description{
Mapping human and mouse genes to HGNC or MGI gene symbols or across species.
}
\usage{
  downloadEnsemblData(host="www.ensembl.org")
  mapGeneNames(maps, genes, in.org=c("Hsap","Mmus"), in.name=c("symbol","ensg"), out.org=c("Hsap","Mmus"), out.name=c("symbol","ensg"))
}
\arguments{
  \item{host}{URL to download ensembl data from (e.g. archive vs current).}
  \item{maps}{output from downloadEnsemblData.}
  \item{genes}{vector of gene IDs.}
  \item{in.org}{organism of input genes.}
  \item{in.name}{input gene ID type.}
  \item{out.org}{desired organism for output genes.}
  \item{out.name}{desired ID type for output.}
}
\details{
	\code{downloadEnsemblData}: Downloads gene ids and orthologs for all protein-coding genes from ensembl biomart.
	\code{mapGeneNames}: Converts genes from one type to another or to orthologs in the other species.
}
\value{ 
	The input SingleCellExperiment object with the regressed expression saved to the 'norm_exprs' slot.
}
\examples{
	maps <- downloadEnsemblData();

	# Convert Whitfield CC table for use in mouse
	MGeneSets$Whitfield <- HGeneSets$Whitfield
	MGeneSets$Whitfield$Gene <- mapGeneNames(maps, as.character(MGeneSets$Whitfield$Gene), in.org="Hsap", in.name="symbol", out.org="Mmus", out.name="symbol")
	MGeneSets$Whitfield <- MGeneSets$Whitfield[MGeneSets$Whitfield$Gene != "",] #remove genes without orthologs.
}
