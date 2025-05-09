plotMixture <- function(fit, BIC=FALSE) {
	mar0 <- par()$mar
	name <- unique(fit$phase[fit$phase!=""])
	if (BIC) {
		par(mfrow=c(1,2))
		par(mar=c(4,4,3,1))
		xes <- as.vector(fit$BIC)
		yes1 <- seq(from=1, to=9, by=3)
		yes2 <- seq(from=2, to=9, by=3)
		yes <- c(yes1, yes2)

		xmin <- min(xes); xmin <- signif(xmin-xmin/20, digits=2)
		plot(xes, yes, xlab="BIC", yaxt="n", bty="n", pch=21, col="black", bg=rep(c("cornflowerblue", "orchid"), each=3), ylab="", cex=2, xlim=c(xmin, max(xes)), main=name)
		arrows( rep(xmin, times=length(yes)), yes, xes, yes, col=rep(c("cornflowerblue", "orchid"), each=3), lty=2, length=0, lwd=1.5)
		axis(2, at=yes, labels=c("E-1", "E-2", "E-3", "V-1", "V-2", "V-3"), tick=FALSE, line=NA, las=2)
	}

	
	# Ylims aren't perfect but fixing them seems hard.
	par(mar=mar0)
	h <- hist(fit$data, prob=TRUE, col="grey85", xlab="Score", main="", breaks=20)
	cols <- c("forestgreen", "darkorange", "blue")
	for (i in 1:fit$G) {
		if (is.na(fit$parameters$variance$sigmasq[i])) {
			fit$parameters$variance$sigmasq[i] <- fit$parameters$variance$sigmasq[1]
		}
		distr_fxn <- function(x) { dnorm(x, mean=fit$parameters$mean[i], sd=sqrt(fit$parameters$variance$sigmasq[i]))*fit$parameters$pro[i] }
		curve(distr_fxn(x), col=cols[i], add=TRUE, lwd=1.5)
	}
	
	if (length(unique(fit$phase)) > 1) {
		thresh <- min(fit$data[fit$phase!=""])
		abline(v=thresh, col="red", lwd=3, lty=3)
		y <- max(h$density)
		x <- mean(h$breaks[h$breaks>thresh])
		text(x, y, name, pos=3, offset=0.01, font=2, cex=1.25, col="black", xpd=TRUE)
	} else {
		x <- max(h$breaks)
		y <- max(h$density)
		name <- "None"
		text(x, y, name, pos=2, offset=0.01, font=2, cex=1, col="black", xpd=TRUE)
	}
}

plot_Phases <- function(classification_out, phases, cell_colours=NULL){
	if (sum(phases %in% colnames(classification_out$scores)) < 2) {
		stop("Insufficient phases specified, or cannot find scores for specified phases.")
	}
	phases <- phases[phases %in% colnames(classification_out$scores)]
	x_dim <- ceiling(sqrt(length(phases)))
	y_dim <- ceiling(length(phases)/x_dim)

	par(mfrow=c(x_dim, y_dim))
	for (i in phases[1:(length(phases)-1)]) {
	for (j in phases[2:(length(phases))]) {
		par(mar=c(4,4,0,0))
		plot(classification_out$scores[,phases[i]], classification_out$scores[,phases[j]],
			xlab=paste(i, "scores"), ylab=paste(j, "scores"), 
			col=cell_colours, pch=18, cex=0.3)
	}}
	
}
