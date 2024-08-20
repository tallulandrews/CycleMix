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
	signature <- signature[signature[,"Gene"] %in% gene_names,]
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

regressCycleScater <- function(SCE, classification, expr_name="logcounts", method=c("scores", "phase")){
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

# regresses out the differences between only the specified cell cycle phases
# first phase will be used as the reference phase for discrete regression
regressCyclePartial <- function(expr_mat, classification, expr_name="logcounts", method=c("scores", "phase"), phases=c("G2M", "G1S"), allow_negative=FALSE) {
	if (method == "phase") {
		existing_phases <- unique(classification[[method]])
		if (sum(phases %in% existing_phases) < 2) {stop("Error: Insufficient stages to regress out")}
	} else if (method == "scores") {
		existing_phases <- colnames(classification[[method]])
		if (sum(phases %in% existing_phases) < 1) {stop("Error: Insufficient stages to regress out")}
	} else {
		stop("Error: not a valid method")
	}

	if (sum(phases %in% existing_phases) < length(existing_phases)){
		missing <- phases[!phases %in% existing_phases]
		phases <- phases[phases %in% existing_phases]
		warning(paste(paste(missing, collapse=","), "not found and will not be regressed."))
	}

	if (method == "phase") {
		new_phase_labels <- as.character(classification[[method]])
		new_phase_labels[!new_phase_labels %in% phases] <- "Other"
		ref_label <- phases[1]
		new_phase_labels <- factor(new_phase_labels, levels=c(ref_label, phases[2:length(phases)], "Other"))
		
		model <- stats::model.matrix(~new_phase_labels)
		model <- model[,colSums(model)>0] # remove any columns with no cells present.
		corrected <- apply(expr_mat, 1, glm_discrete, phases, model, allow_negative)
	}

	if (method == "scores") {
		model <- stats::model.matrix(~classification[[method]])
		corrected <- apply(expr_mat, 1, glm_continuous, phases, model, allow_negative)
	}
	return(corrected);
}

glm_discrete <- function(x, phases, model, allow_negative=FALSE) {
	if (var(x) == 0) {return(x)}
	if (min(x) >=0) {
		print("Assuming expression values should be >=0")
	}
	res <- MASS::glm.nb(x~model)
	# Coeffs are average difference 
	change <- rep(0, nrow(model))
	to_change <- colnames(model)
	to_change <- to_change[!(grepl("Intercept", to_change) | grepl("Other", to_change))]
	for (coeff in to_change) {
		this_eff <- res$coefficients[coeff]
		if (!allow_negative) { #adjust for not changing 0s
			adj <- mean(x[model[,coeff]==1] > 0) # proportion of non-0s - these won't be changed
			this_eff <- this_eff*1/adj; # adjust for not shifting zeros
		}
		change[model[,coeff] == 1] <- this_eff
	}
	out <- x-change
	if (min(x) >=0) {
		out[out < 0] <- 0
		out[x==0] <- 0
	}
	return(out)
}
glm_continuous <- function(x, model, allow_negative=FALSE) {
	if (var(x) == 0) {return(x)}
	res <- MASS::glm.nb(x~model)
	out <- res$residuals
	if (!allow_negative) {
		adj <- sum(out < 0)/length(out);	
		diff <- x-out
		diff <- diff/adj
		out <- x-diff
		out[out<0] <- 0
		return(out)
	} else{
		return(out)
	}
}	
