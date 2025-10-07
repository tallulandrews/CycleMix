prepData <- function(obj, expr_name="logcounts", do.scale=FALSE, symbol_column="feature_symbol") {
	gene_names <- rownames(obj)

	expr_mat <- c()
	if (class(obj)[1] == "SingleCellExperiment") {
		expr_mat <- SummarizedExperiment::assays( obj )[[expr_name]]
		gene_names <- rownames(expr_mat)
		if (!is.null(symbol_column)) {
			gene_names <- SummarizedExperiment::rowData(obj)[ , symbol_column]
		}
	}
	if (class(obj)[1] == "Seurat") {
		expr_mat <- Seurat::GetAssayData( obj, assay=expr_name, layer="data" )
		gene_names <- rownames(expr_mat)
		if (!is.null(symbol_column)) {
			gene_names <- symbol_column
		}
	}
	if (class(obj)[1] == "dgCMatrix" | class(obj)[1] == "dgTMatrix" | class(obj) == "matrix") {
		expr_mat <- obj
		gene_names <- rownames(expr_mat)
		if (!is.null(symbol_column)) {
			gene_names <- symbol_column
		}
	}
	if (do.scale) {
		expr_mat <- t( apply(expr_mat, 1, scale ) )
	}
	if (sum(duplicated(gene_names)) > 0) {
		exclude <- duplicated(gene_names)
		expr_mat <- expr_mat[!exclude,]
		gene_names <- gene_names[!exclude]
		warning(paste("Excluding", sum(exclude), "duplicated genes."))
	}
	rownames(expr_mat) <- gene_names
	return(expr_mat)
}

assignPhase <- function(expr_mat, CC_table, phase="G2M") {
	# Expression Signature
	signature <- CC_table[,"Stage"] == phase
	if (sum(signature) <= 1) {
		stop( paste("Error: Insufficient genes associated with phase:", phase) )
	}
	signature <- CC_table[signature, ]
	
	# Match up
	gene_names <- rownames(expr_mat)
	signature <- signature[signature[,"Gene"] %in% gene_names,]
	if (nrow(signature) < 5) {
		stop("Error: fewer than 5 phase genes detected in obj. Check gene names of the CC table match the gene names in the obj.")
	}
	matches <- base::match(signature[,"Gene"], gene_names)
	expr_mat <- expr_mat[matches, ]
	gene_names <- gene_names[matches]
	keep <- !is.na(gene_names)
	expr_mat <- expr_mat[keep,]
	gene_names <- gene_names[keep]
	signature <- signature[keep,]

	# Cell-Scores
	score <- Matrix::colSums(expr_mat*signature[,"Dir"])/sum(abs(signature[,"Dir"]))

	# Fit Mixture Model
	fit <- mclust::Mclust(score, G=1:3)
	fit$phase <- as.character(fit$classification)
	
	if (fit$G > 1) {
		tag <- which(fit$parameters$mean == max(fit$parameters$mean))
		fit$phase[fit$phase == tag] <- phase
	}
	fit$phase[fit$phase != phase] <- ""
	return(fit)
}

# Seurat: obj, gene_table, expr_name="RNA", symbol_column=NULL
classifyCells <- function(obj, CC_table, expr_name="logcounts", do.scale=FALSE, symbol_column="feature_symbol", allow.multi=FALSE) {

	out_list <- list()
	stages <- as.character(unique(CC_table[,"Stage"]))
	phases <- matrix(nrow=ncol(obj), ncol=length(stages));
	scores <- matrix(nrow=ncol(obj), ncol=length(stages));
	expr_mat <- prepData(obj, do.scale=do.scale, expr_name=expr_name, symbol_column=symbol_column)

	for (i in 1:length(stages)) {
		assignment <- assignPhase(expr_mat, CC_table, phase=stages[i])
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

#regressCycleScater <- function(obj, classification, expr_name="logcounts", method=c("scores", "phase")){
#	if (class(obj)[1] != "SingleCellExperiment") {
#		stop("Error: Requires SingleCellExperiment object as input")
#	}
#
#	if (method[1] == "phase") {
#		design <- model.matrix(~classification$phase)
#	} else if (method[1] == "scores") {
#		design <- model.matrix(~classification$scores)
#	} else {
#		stop("Error: unrecognized method.")
#	}
#	obj <- scater::normalizeExprs(obj, design=design, return_norm_as_exprs=FALSE, exprs_values=expr_name)
#	return(obj);
#}

# regresses out the differences between only the specified cell cycle phases
# first phase will be used as the reference phase for discrete regression
regressCyclePartial <- function(expr_mat, classification, type=c("counts","norm"), method=c("scores", "phase"), phases=c("G2M", "G1S"), allow_negative=FALSE, subsample_cells=ncol(expr_mat)) {

	if (class(classification) != "list") {
		existing_phases <- unique(classification)
		method = "phase"
	} else {
		classification <- classification[[method]]
	}
		
	if (method == "phase") {
		existing_phases <- unique(classification)
		if (sum(phases %in% existing_phases) < 2) {stop("Error: Insufficient stages to regress out")}
	} else if (method == "scores") {
		existing_phases <- colnames(classification)
		if (sum(phases %in% existing_phases) < 1) {stop("Error: Insufficient stages to regress out")}
	} else {
		stop("Error: not a valid method")
	}


	if (sum(phases %in% existing_phases) < length(phases)){
		missing <- phases[!phases %in% existing_phases]
		phases <- phases[phases %in% existing_phases]
		warning(paste(paste(missing, collapse=","), "not found and will not be regressed."))
	}


	if (method == "phase") {
		new_phase_labels <- as.character(classification)
		new_phase_labels[!new_phase_labels %in% phases] <- "Other"
		phase_counts <- table(new_phase_labels);
		# Subsample cells
		cells <- c()
		if (subsample_cells < ncol(expr_mat)) {
			if (min(phase_counts) < ceiling(subsample_cells/length(phase_counts))) {
				to_sample <- min(phase_counts)
			} else {
				to_sample <- ceiling(subsample_cells/length(phase_counts))
			}
			for (phase in names(phase_counts)) {
				cells <- c(cells, sample(colnames(expr_mat)[new_phase_labels==phase], size=to_sample))
			}
			if (length(cells) < subsample_cells) {
				cells <- c(cells, sample(colnames(expr_mat)[! colnames(expr_mat) %in% cells], size=subsample_cells-length(cells)))
			}
		} else {
			cells = colnames(expr_mat)
		}


		ref_label <- phases[1]
		new_phase_labels <- factor(new_phase_labels, levels=c(ref_label, phases[2:length(phases)], "Other"))
		
		model <- stats::model.matrix(~new_phase_labels)
		model <- model[,Matrix::colSums(model)>0] # remove any columns with no cells present.
		corrected <- apply(expr_mat, 1, glm_discrete, phases, model, type, allow_negative, colnames(expr_mat) %in% cells)
	}

	if (method == "scores") {
		to_regress <- classification[,colnames(classification) %in% phases]
		model <- stats::model.matrix(~to_regress)
		corrected <- apply(expr_mat, 1, glm_continuous,  model=model, type=type, allow_negative=allow_negative)
	}
	return(t(corrected));
}

glm_discrete <- function(x, phases, model, type=c("counts", "norm"), allow_negative=FALSE, tofit.cells=colnames(x)) {
	if (sum(tofit.cells) < length(x)) {
		x2 <- x[tofit.cells]
		model2 <- model[tofit.cells,]
	} else {
		x2 <- x
		model2 <- model
	}
	if (var(x2) == 0) {return(x)}


	if (type[1] == "counts") {
		res <- MASS::glm.nb(x2~model2)
	} else if (type[1] == "norm") {
		res <- glm(x2~model2)
	}
	# Coeffs are average difference 
	change <- rep(0, length(x))
	to_change <- names(res$coefficients)
	to_change <- to_change[!(grepl("Intercept", to_change) | grepl("Other", to_change))]
	for (coeff in to_change) {
		this_eff <- res$coefficients[coeff]
		if (!allow_negative) { #adjust for not changing 0s
			adj <- mean(x[model[,sub("model2", "",coeff)]==1] > 0) # proportion of non-0s - these won't be changed
			this_eff <- this_eff*1/adj; # adjust for not shifting zeros
		}
		change[model[,sub("model2", "",coeff)] == 1] <- this_eff
	}
	out <- x-change
	if (!allow_negative) {
		out[out < 0] <- 0
		out[x==0] <- 0
	}
	return(out)
}

glm_continuous <- function(x, model,  type=c("counts", "norm"), allow_negative=FALSE) {
	if (var(x) == 0) {return(x)}
	        if (type[1] == "counts") {
                res <- MASS::glm.nb(x~model)
        } else if (type[1] == "norm") {
                res <- glm(x~model)
        }
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
