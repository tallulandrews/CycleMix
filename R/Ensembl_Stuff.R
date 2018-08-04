getEnsemblArchive <- function() {
	tab <- biomaRt::listEnsemblArchives()
	date <- matrix(unlist(strsplit(tab[,"date"], " ")), ncol=2, byrow=T)
	date <- as.Date(paste(1, toupper(date[,1]), date[,2], sep="-"), format="%d-%b-%Y")
	tab <- tab[order(date, decreasing=T),]
	tab <- tab[ grepl("Ensembl [1-9]", tab[,1], ignore.case=F) , ]
	return(tab[,"url"])	
}

downloadEnsemblData <- function(host="www.ensembl.org") {
	requireNamespace("biomaRt")
	archive <- getEnsemblArchive();
	biomart_mart <- "ENSEMBL_MART_ENSEMBL"
	hensembl <- tryCatch({
		message(paste("Trying :", host[1]))
		biomaRt::useMart(biomart=biomart_mart, dataset="hsapiens_gene_ensembl", host=host[1])
	}, error=function(cond){
		message("Error connecting to biomart, using archive")
		message(paste("Using :", archive[1]))
		return(biomaRt::useMart(biomart=biomart_mart, dataset="hsapiens_gene_ensembl", host=archive[1]))
	}, warning=function(cond){
		message(cond)
		return(biomaRt::useMart(biomart=biomart_mart, dataset="hsapiens_gene_ensembl", host=host[1]))
	})
	
		
	mensembl <- tryCatch({
		message(paste("Trying :", host[1]))
		biomaRt::useMart(biomart=biomart_mart, dataset="mmusculus_gene_ensembl", host=host[1])
	}, error=function(cond){
		message("Error connecting to biomart, using archive")
		message(paste("Using :", archive[1]))
		return(biomaRt::useMart(biomart=biomart_mart, dataset="mmusculus_gene_ensembl", host=archive[1]))
	}, warning=function(cond){
		message(cond)
		return(biomaRt::useMart(biomart=biomart_mart, dataset="mmusculus_gene_ensembl", host=host[1]))
	})

	ensg_name_map <- biomaRt::getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), filters = "biotype", values="protein_coding", mart=hensembl)
	ensg2musg <- biomaRt::getBM(attributes=c("ensembl_gene_id","mmusculus_homolog_ensembl_gene"),filters = "biotype", values="protein_coding", mart=hensembl)
	musg_name_map <- biomaRt::getBM(attributes=c("ensembl_gene_id","mgi_symbol"), filters = "biotype", values="protein_coding", mart=mensembl)
	return(list(Hname=ensg_name_map, Orth=ensg2musg, Mname=musg_name_map))
}

map_symbol_ensg <- function(maps, genes, is.org=c("Hsap","Mmus"), is.name=c("symbol","ensg")) {
	if (is.org[1] == "Hsap") {
		if(is.name[1]=="symbol") {
			new = as.character(maps$Hname[match(genes, maps$Hname[,2]),1])
		} else if (is.name[1] =="ensg") {
			new = as.character(maps$Hname[match(genes, maps$Hname[,1]),2])
		} else {
			stop("Unrecognized name type")
		}
	} else if (is.org[1] == "Mmus") {
		if(is.name[1]=="symbol") {
			new = as.character(maps$Mname[match(genes, maps$Mname[,2]),1])
		} else if (is.name[1] =="ensg") {
			new = as.character(maps$Mname[match(genes, maps$Mname[,1]),2])
		} else {
			stop("Unrecognized name type")
		}
	} else {
		stop("Unrecognized organism");
	}
	new[is.na(new)] = ""
	return(new);
}

map_Hsap_Mmus_one2one <- function(maps, genes, is.org=c("Hsap","Mmus")) {
	if (is.org[1] == "Hsap") {
		new = as.character(maps$Orth[match(genes, maps$Orth[,1]),2])
	} else if (is.org[1] == "Mmus") {
		new = as.character(maps$Orth[match(genes, maps$Orth[,2]),1])
	} else {
		stop("Unrecognized organism");
	}
	new[is.na(new)] = ""
	return(new);
}

mapGeneNames<- function(maps, genes, in.org=c("Hsap","Mmus"), in.name=c("symbol","ensg"), out.org=c("Hsap","Mmus"), out.name=c("symbol","ensg")) {
	genes <- as.character(genes)
	if (in.org == out.org & in.name == out.name) {
		# No change
		return(genes)
	}
	if (in.org == out.org) {
		# change names not spp
		return(map_symbol_ensg(maps, genes, is.org=in.org, is.name=in.name))

	} else {
		# to ensg names
		if (in.name == "symbol") {
			tmp <- map_symbol_ensg(maps, genes, is.org=in.org, is.name=in.name)
		} else {
			tmp <- genes
		}
		# change organism
		tmp <- map_Hsap_Mmus_one2one(maps, tmp, is.org=in.org)

		# to symbol if necessary
		if (out.name =="symbol") {
			out <- map_symbol_ensg(maps, tmp, is.org=out.org, is.name="ensg")
		} else {
			out <- tmp
		}
		return(out)	
	}
}
