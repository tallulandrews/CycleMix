set.seed(1973)

require("SingleCellExperiment")
require("scater")
expr_type <- "lognorm"
Ex <- readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/buettner.rds");
Ex <- scater::toSingleCellExperiment(Ex)
Ex <- Ex[,colSums(counts(Ex)>0) > 6000]



# Cell-Cycle gene lists
Whitfield <- read.table("~/Data/Whitfield_CC.txt")
Tirosh <- read.table("~/Data/cellcycle3.txt", header=FALSE)
Macosko <- read.table("~/Data/cellcyclegenes.csv", header=FALSE, stringsAsFactors=FALSE)
Quiesc <- read.table("~/Data/Quiescence.txt", stringsAsFactors=FALSE) # Cheung
G0 <- read.table("~/Data/Reactome_G0.txt", stringsAsFactors=FALSE)

# Build structure:
HGeneSets <- list()

HGeneSets$Whitfield <- unique(data.frame(Gene=Whitfield[,2], Stage=Whitfield[,1], Dir=rep(1, times=nrow(Whitfield))))

Tirosh[,1] <- as.character(Tirosh[,1]);
Tirosh[,1] <- sub("/", "", Tirosh[,1])
Tirosh[,1] <- factor(Tirosh[,1]);
HGeneSets$Tirosh <- unique(data.frame(Gene=Tirosh[,2], Stage=Tirosh[,1], Dir=rep(1, times=nrow(Tirosh))))

Macosko_gene <- c(Macosko[-1,1], Macosko[-1,2], Macosko[-1,3], Macosko[-1,4], Macosko[-1,5])
Macosko_stage <- rep(c("G1S", "S", "G2M", "M", "MG1"), each=nrow(Macosko)-1)
junk <- is.na(Macosko_gene);
Macosko_gene <- Macosko_gene[!junk]
Macosko_stage <- Macosko_stage[!junk]
HGeneSets$Macosko <- unique(data.frame(Gene=Macosko_gene, Stage=Macosko_stage, Dir=rep(1, times=length(Macosko_gene))))

G0 <- cbind(G0, rep(1, times=length(G0)))
HGeneSets$Quiesc <- unique(data.frame(Gene=c(Quiesc[,1], G0[,1]), Stage=rep("G0", times=nrow(Quiesc)+nrow(G0)), Dir=c(Quiesc[,2], G0[,2])))

save(HGeneSets, file="data/HGenes.rda", compress="xz")

# mouse:
require("scran")
require("CellTypeProfiles")
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))

# Build structure (uses buettner data to learn stage)
MGeneSets <-list()

cyclone_unique <- unique(c(mm.pairs$G1[,1], mm.pairs$G1[,2], mm.pairs$S[,1], mm.pairs$S[,2], mm.pairs$G2[,1], mm.pairs$G2[,2]))
out <- my_row_mean_aggregate(assays(Ex)[["logcounts"]], colData(Ex)$cell_type1)
out <- out - rowMeans(out)
out <- out[rownames(out) %in% cyclone_unique,]
stage <- apply(out, 1, function(x) {colnames(out)[which(abs(x) == max(abs(x)))]})
dir <- apply(out, 1, function(x) {sign(x[which(abs(x) == max(abs(x)))])})
gene <- names(stage)
source("~/R-Scripts/Ensembl_Stuff.R")
gene <- General_Map(gene, in.org="Mmus", out.org="Mmus", in.name="ensg", out.name="symbol")

keep <- gene != "";
MGeneSets$Cyclone <- data.frame(Gene=gene[keep], Stage=stage[keep], Dir=dir[keep])

save(MGeneSets, file="data/MGenes.rda", compress="xz")


# Filter data

keep1 <- rowData(Ex)$feature_symbol %in% MGeneSets$Cyclone$Gene
fs <- M3DropFeatureSelection(expr_mat=assays(Ex)[["counts"]])
keep2 <- rownames(Ex) %in% fs$Gene
Ex <- Ex[keep1 | keep2,]
save(Ex, file="data/Ex.rda", compress="xz")

require("slalom")
gene_lists <- load_CC("all")
cellcycle <- gene_lists$Whitfield
cellcycle_simple <- gene_lists$Simple
GO_genes <- gene_lists$G0
Quiescence <- gene_lists$Quiescence
new_cellcycle <- gene_lists$Tirosh
#cellcycle <- read.table("~/Data/Whitfield_CC.txt")
#cellcycle_simple <- as.matrix(cellcycle[cellcycle[,1] != "CC",])
#cellcycle_simple[cellcycle_simple[,1] == "G2",1] = "G2M";
#cellcycle_simple[cellcycle_simple[,1] == "S",1] = "G1S";
#cellcycle_simple = cellcycle_simple[cellcycle_simple[,1] != "MG1",];
#G0_genes <- read.table("~/Data/Reactome_G0.txt", header=F) # Update this with genes from Laura/MiSigDB?
#Quiescence <- read.table("~/Data/Quiescence.txt")
#new_cellcycle <- read.table("~/Collaborations/LiverOrganoids/New_CC_171117.txt", header=FALSE) #Tirsoh2016Nature

get_prolif <- function(CC, Ex, suppress.plot=TRUE) {
	require(mixtools)
	cc_g <- fData(Ex)$Symbol %in% CC
	score <- colSums(get_exprs(Ex, expr_type)[cc_g,])
	mix <- normalmixEM(score)
	if (!suppress.plot) {
		plot(mix, which=2)
	}
	p1 <- dnorm(score, mean = mix$mu[1], sd=mix$sigma[1]) 
	p2 <- dnorm(score, mean = mix$mu[2], sd=mix$sigma[2])
	if (mix$mu[1] < mix$mu[2]) {
		assign <- p2 > p1
	} else {
		assign <- p1 > p2
	}
	return(list(score=score, assign=assign))
}

get_prolif2 <- function(CC, Ex, name="cycle genes", suppress.plot=TRUE, min_detected = 0.05) {
	require(mixtools)
        require(matrixStats)
        loggednorm <- get_exprs(Ex, expr_type)
        cc_g <- fData(Ex)$Symbol %in% CC
	loggednorm <- loggednorm[cc_g,]
	loggednorm <- loggednorm[rowVars(loggednorm) > 0,]
	loggednorm <- loggednorm[rowMeans(loggednorm > 0) >= min_detected,]
        loggednorm <- (loggednorm-rowMeans(loggednorm))/rowVars(loggednorm)
        score <- colMeans(loggednorm)
        mix <- normalmixEM(score)
        p1 <- dnorm(score, mean = mix$mu[1], sd=mix$sigma[1])
        p2 <- dnorm(score, mean = mix$mu[2], sd=mix$sigma[2])
	
	
	m1 = mix$mu[1]
	m2 = mix$mu[2]
	means_test <- abs(m1-m2) > min(mix$sigma)
	

        if (m1 < m2 & means_test) {
                assign <- p2 > p1
        } else if (m2 < m1 & means_test) {
                assign <- p1 > p2
        } else {
		assign <- rep(FALSE, times=length(p1))
		name=paste(name, "- No Difference")
	}
        if (!suppress.plot) {
                plot(mix, which=2, main2=name, xlab2="Score")
        }
        return(list(score=score, assign=assign))

}

get_prolif3 <- function(CC, Ex, name="cycle genes", suppress.plot=TRUE, min_detected = 0.05) {
        require(mixtools)
        require(matrixStats)
        loggednorm <- get_exprs(Ex, expr_type)
        cc_g <- fData(Ex)$Symbol %in% CC[,1]
        loggednorm <- loggednorm[cc_g,]
	rownames(loggednorm) <-fData(Ex)$Symbol[cc_g]; 
        loggednorm <- loggednorm[rowVars(loggednorm) > 0,]
        loggednorm <- loggednorm[rowMeans(loggednorm > 0) >= min_detected,]
        loggednorm <- (loggednorm-rowMeans(loggednorm))/rowVars(loggednorm)
	
	rownames(CC) <- CC[,1]
	weights <- CC[match(rownames(loggednorm), CC[,1]),]
	loggednorm <- loggednorm*weights[,2]
	
        score <- colMeans(loggednorm)
        mix <- normalmixEM(score)
        p1 <- dnorm(score, mean = mix$mu[1], sd=mix$sigma[1])
        p2 <- dnorm(score, mean = mix$mu[2], sd=mix$sigma[2])
        
        
        m1 = mix$mu[1]
        m2 = mix$mu[2]
        means_test <- abs(m1-m2) > min(mix$sigma)
        

        if (m1 < m2 & means_test) {
                assign <- p2 > p1
        } else if (m2 < m1 & means_test) {
                assign <- p1 > p2
        } else {
                assign <- rep(FALSE, times=length(p1))
                name=paste(name, "- No Difference")
        }
        if (!suppress.plot) {
                plot(mix, which=2, main2=name, xlab2="Score")
        }
        return(list(score=score, assign=assign))

}


prolif_analysis<- function(Ex, name="Data") {
	
	# Get Prolif
	png(paste(name, "Proliferation_MixtureModels.png", sep="_"), width=8, height=8*2/3, units="in", res=300)
	par(mfrow=c(2,3))
	sets <- levels(factor(cellcycle[,1]))
	out <- matrix(0, nrow=length(pData(Ex)[,1]), ncol=length(sets))
	for( i in 1:length(sets) ) {
		s = sets[i]
		assign<-get_prolif2(cellcycle[cellcycle[,1]==s,2], Ex, name=s, suppress.plot=F)$assign
		out[,i] = assign
	}
	dev.off()

	summary(factor(rowSums(out)))
	prolif = rowSums(out) >= 3

	# Output
	return(prolif)
}

prolif_analysis2 <- function(Ex, name="Data") {
	
#	g0 <- get_prolif2(G0_genes[,1], Ex, name="G0", suppress.plot=F)$assign
	png(paste(name, "MixtureModels2_Figure.png", sep="_"), width=7*3/2, height=7, units="in", res=300)
	par(mfrow=c(2,3))
	par(mar=c(4,4,1,1))
	g1 <- get_prolif2(cellcycle_simple[cellcycle_simple[,1] == "G1S",2], Ex, name="G1S", suppress.plot=F)
	g1_new <- get_prolif2(new_cellcycle[new_cellcycle[,2] == "G1S",1], Ex, name="(new) G1S", suppress.plot=F)
	g2 <- get_prolif2(cellcycle_simple[cellcycle_simple[,1] == "G2M",2], Ex, name="G2M", suppress.plot=F)
	g2_new <- get_prolif2(new_cellcycle[new_cellcycle[,2] == "G2M",1], Ex, name="(new) G2M", suppress.plot=F)
	g0 <- get_prolif3(Quiescence, Ex, name="G0", suppress.plot=F)
	dev.off()

	Assign <- cbind(g1$assign, g2$assign, g0$assign)
	colnames(Assign) = c("G1S", "G2M", "G0")
	Score <- cbind(g1$score, g2$score, g0$score)
	colnames(Score) = c("G1S", "G2M", "G0")
	
	state <- apply(Score, 1, function(x) {colnames(Score)[which(x == max(x))]})
	state[state=="G0" & !Assign[,3]] <- "None"
	state[state=="G1S" & !Assign[,1]] <- "None"
	state[state=="G2M" & !Assign[,2]] <- "None"

	state1 <- factor(state, levels=c("None", "G0", "G1S", "G2M"))

	Assign <- cbind(g1_new$assign, g2_new$assign, g0$assign)
	colnames(Assign) = c("G1S", "G2M", "G0")
	Score <- cbind(g1_new$score, g2_new$score, g0$score)
	colnames(Score) = c("G1S", "G2M", "G0")
	
	state <- apply(Score, 1, function(x) {colnames(Score)[which(x == max(x))]})
	state[state=="G0" & !Assign[,3]] <- "None"
	state[state=="G1S" & !Assign[,1]] <- "None"
	state[state=="G2M" & !Assign[,2]] <- "None"

	state2 <- factor(state, levels=c("None", "G0", "G1S", "G2M"))

	# Output
	return(list(state1, state2))
}

prolif_out <- prolif_analysis(Ex, args[2])
pData(Ex)$Proliferating <- prolif_out
prolif2_out <- prolif_analysis2(Ex, args[2])
pData(Ex)$CC_state <- prolif2_out[[1]]
pData(Ex)$CC_state_new <- prolif2_out[[2]]
print(table(pData(Ex)$CC_state, pData(Ex)$Proliferating))
print(table(pData(Ex)$CC_state, pData(Ex)$CC_state_new))
saveRDS(Ex, file=paste(args[2], "Prolif.rds", sep="_"))
