assignPhase <- function(SCE, CC_table, phase="G2M", expr_name="logcounts", do.scale=FALSE, symbol_column="feature_symbol") {
	if (class(SCE)[1] != "SingleCellExperiment") {
		stop("Error: Requires SingleCellExperiment object as input")
	}

	# Expression Signature
	signature <- CC_table[,"Stage"] == phase
	if (sum(signature) <= 1) {
		stop( paste("Error: Insufficient genes associated with phase:", phase) )
	}
	signature <- CC_table[signature, ]
	
	# Scale?
	exprmat <- SummarizedExperiment::assays( SCE )[[expr_name]]
	if (do.scale) {
		exprmat <- t( apply(exprmat, 1, scale ) )
	}
	if (is.null(symbol_column)) {
		rowData(SCE)$CycleMixSym <- rownames(SCE)
		symbol_column <- "CycleMixSym";
	}
	
	# Match up
	gene_names <- SummarizedExperiment::rowData(SCE)[ , symbol_column]
	matches <- base::match(signature[,"Gene"], gene_names)
	exprmat <- exprmat[matches, ]
	gene_names <- gene_names[matches]
	keep <- !is.na(gene_names)
	exprmat <- exprmat[keep,]
	gene_names <- gene_names[keep]
	signature <- signature[keep,]

	# Cell-Scores
	score <- colSums(exprmat*signature[,"Dir"])/sum(abs(signature[,"Dir"]))
	#pos_dir <- signature[,"Dir"]
	#pos_dir[pos_dir < 0] <- 0
	#pos_score <- colSums(exprmat*pos_dir)/sum(pos_dir)

	# Fit Mixture Model
	#require("mclust")
	fit <- mclust::Mclust(score, G=1:3)
	fit$phase <- as.character(fit$classification)
	
	if (fit$G > 1) {
		tag <- which(fit$parameters$mean == max(fit$parameters$mean))
		fit$phase[fit$phase == tag] <- phase
	}
	fit$phase[fit$phase != phase] <- ""
	#fit$pos_score <- pos_score
	return(fit)
}


classifyCells <- function(SCE, CC_table, expr_name="logcounts", do.scale=FALSE, symbol_column="feature_symbol", allow.multi=FALSE) {
	if (class(SCE)[1] != "SingleCellExperiment") {
		stop("Error: Requires SingleCellExperiment object as input")
	}

	out_list <- list()
	stages <- as.character(unique(CC_table[,"Stage"]))
	phases <- matrix(nrow=ncol(SCE), ncol=length(stages));
	scores <- matrix(nrow=ncol(SCE), ncol=length(stages));
	for (i in 1:length(stages)) {
		assignment <- assignPhase(SCE, CC_table, phase=stages[i], do.scale=do.scale, expr_name=expr_name, symbol_column = symbol_column)
		out_list[[stages[i]]] <- assignment
		phases[,i] <- assignment$phase
		scores[,i] <- assignment$data
	}
	if (allow.multi) {
		best <- apply(phases, 1, paste, collapse="")
	} else {
		if (!do.scale) {
			scores <- apply(scores, 2, scale)
		} else {
			scores <- apply(scores, 2, scale, center=FALSE)
		}
		best <- sapply(1:nrow(scores), function(i) {
			phases[i,which(scores[i,] == max(scores[i,]))]
			})
	}
	best[best==""] <- "None"
	scores <- as.matrix(scores)
	colnames(scores) <- stages
	return(list(phase=best, scores=scores, fits=out_list))
}

regressCycle_scater <- function(SCE, classification, expr_name="logcounts", method=c("scores", "phase")){
	if (class(SCE)[1] != "SingleCellExperiment") {
		stop("Error: Requires SingleCellExperiment object as input")
	}

	if (method[1] == "phase") {
		design <- model.matrix(~classification$phase)
	} else if (method[1] == "scores") {
		design <- model.matrix(~classification$scores)
	} else {
		stop("Error: unrecognized method.")
	}
	SCE <- scater::normalizeExprs(SCE, design=design, return_norm_as_exprs=FALSE, exprs_values=expr_name)
	return(SCE);
}

regressCycle <- function(SCE, classification, expr_name="logcounts", method=c("scores", "phase"), phases=c("G2M", "G1S")) {
	glm_fun <- function(x) {
		if (var(x) == 0) {return(x)}
		model <- stats::model.matrix(~classification[[method]])
		res <- glm(x~model[,-1])
		eff <- vector();
		mod <- vector();
		for(p in phases) {
			selected <- grep(p, names(res$coef))
			eff <- c(eff, res$coef[selected])
			mod <- cbind(mod, model[,selected])
		}
		eff <- eff*1/mean(x > 0); # adjust for not shifting zeros
		norm <- mean(eff)
		norm_factor <- -1*rowSums( t(t(mod)*eff) ) + rowSums( mod*norm );
		zeros <- which (x == 0);
		x <- x+norm_factor;
		x[zeros] <- 0;
		x[x < 0] <- 0;
		
		return(x)
	}

	corrected <- apply(assays(SCE)[[expr_name]], 1, glm_fun)
	assays(SCE)[["CC_cor"]] <- corrected;
	return(SCE);
}
