#set.seed(1973)
#
#require("SingleCellExperiment")
#require("scater")
#expr_type <- "lognorm"
#
## Load scRNASeq data
#Ex <- readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/buettner.rds");
#Ex <- scater::toSingleCellExperiment(Ex)
#Ex <- Ex[,colSums(counts(Ex)>0) > 6000]
#
## Cell-Cycle gene lists
#Whitfield <- read.table("~/Data/Whitfield_CC.txt")
#Tirosh <- read.table("~/Data/cellcycle3.txt", header=FALSE)
#Macosko <- read.table("~/Data/cellcyclegenes.csv", header=FALSE, stringsAsFactors=FALSE)
#Quiesc <- read.table("~/Data/Quiescence.txt", stringsAsFactors=FALSE) # Cheung
#G0 <- read.table("~/Data/Reactome_G0.txt", stringsAsFactors=FALSE)
#
## Build structure:
#HGeneSets <- list()
#
#HGeneSets$Whitfield <- unique(data.frame(Gene=Whitfield[,2], Stage=Whitfield[,1], Dir=rep(1, times=nrow(Whitfield))))
#
#Tirosh[,1] <- as.character(Tirosh[,1]);
#Tirosh[,1] <- sub("/", "", Tirosh[,1])
#Tirosh[,1] <- factor(Tirosh[,1]);
#HGeneSets$Tirosh <- unique(data.frame(Gene=Tirosh[,2], Stage=Tirosh[,1], Dir=rep(1, times=nrow(Tirosh))))
#
#Macosko_gene <- c(Macosko[-1,1], Macosko[-1,2], Macosko[-1,3], Macosko[-1,4], Macosko[-1,5])
#Macosko_stage <- rep(c("G1S", "S", "G2M", "M", "MG1"), each=nrow(Macosko)-1)
#junk <- is.na(Macosko_gene);
#Macosko_gene <- Macosko_gene[!junk]
#Macosko_stage <- Macosko_stage[!junk]
#HGeneSets$Macosko <- unique(data.frame(Gene=Macosko_gene, Stage=Macosko_stage, Dir=rep(1, times=length(Macosko_gene))))
#
#G0 <- cbind(G0, rep(1, times=length(G0)))
#HGeneSets$Quiesc <- unique(data.frame(Gene=c(Quiesc[,1], G0[,1]), Stage=rep("G0", times=nrow(Quiesc)+nrow(G0)), Dir=c(Quiesc[,2], G0[,2])))
#
#save(HGeneSets, file="HGenes.rda", compress="xz")
#
## mouse:
#require("scran")
#require("CellTypeProfiles")
#mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
#
## Build structure (uses buettner data to learn stage)
#MGeneSets <-list()
#
#cyclone_unique <- unique(c(mm.pairs$G1[,1], mm.pairs$G1[,2], mm.pairs$S[,1], mm.pairs$S[,2], mm.pairs$G2[,1], mm.pairs$G2[,2]))
#out <- my_row_mean_aggregate(assays(Ex)[["logcounts"]], colData(Ex)$cell_type1)
#out <- out - rowMeans(out)
#out <- out[rownames(out) %in% cyclone_unique,]
#stage <- apply(out, 1, function(x) {colnames(out)[which(abs(x) == max(abs(x)))]})
#dir <- apply(out, 1, function(x) {sign(x[which(abs(x) == max(abs(x)))])})
#gene <- names(stage)
#source("~/R-Scripts/Ensembl_Stuff.R")
#gene <- General_Map(gene, in.org="Mmus", out.org="Mmus", in.name="ensg", out.name="symbol")
#
#keep <- gene != "";
#MGeneSets$Cyclone <- data.frame(Gene=gene[keep], Stage=stage[keep], Dir=dir[keep])
#
#save(MGeneSets, file="MGenes.rda", compress="xz")
#
#
## Filter data
#
#require("M3Drop")
#keep1 <- rowData(Ex)$feature_symbol %in% MGeneSets$Cyclone$Gene
#fs <- M3DropFeatureSelection(expr_mat=assays(Ex)[["counts"]])
#keep2 <- rownames(Ex) %in% fs$Gene
#Ex <- Ex[keep1 | keep2,]
#save(Ex, file="Ex.rda", compress="xz")
#rm(list=ls())
