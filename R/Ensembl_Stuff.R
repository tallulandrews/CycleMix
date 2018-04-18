
downloadEnsemblData <- function() {
	#require("biomaRt")
	hensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
	mensembl = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")

	ensg_name_map <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), filters = "biotype", values="protein_coding", mart=hensembl)
	ensg2musg <- getBM(attributes=c("ensembl_gene_id","mmusculus_homolog_ensembl_gene"),filters = "biotype", values="protein_coding", mart=hensembl)
	musg_name_map <- getBM(attributes=c("ensembl_gene_id","mgi_symbol"), filters = "biotype", values="protein_coding", mart=mensembl)
	return(list(Hname=ensg_name_map, Orth=ensg2musg, Mname=musg_name_map))
}

map_symbol_ensg <- function(maps, genes, is.org=c("Hsap","Mmus"), is.name=c("symbol","ensg")) {
	if (is.org[1] == "Hsap") {
		if(is.name[1]=="symbol") {
			new = as.character(maps$ensg_name_map[match(genes, ensg_name_map[,2]),1])
		} else if (is.name[1] =="ensg") {
			new = as.character(maps$ensg_name_map[match(genes, ensg_name_map[,1]),2])
		} else {
			stop("Unrecognized name type")
		}
	} else if (is.org[1] == "Mmus") {
		if(is.name[1]=="symbol") {
			new = as.character(maps$musg_name_map[match(genes, musg_name_map[,2]),1])
		} else if (is.name[1] =="ensg") {
			new = as.character(maps$musg_name_map[match(genes, musg_name_map[,1]),2])
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
		new = as.character(maps$ensg2musg[match(genes, ensg2musg[,1]),2])
	} else if (is.org[1] == "Mmus") {
		new = as.character(maps$ensg2musg[match(genes, ensg2musg[,2]),1])
	} else {
		stop("Unrecognized organism");
	}
	new[is.na(new)] = ""
	return(new);
}

mapGeneNames<- function(maps, genes, in.org=c("Hsap","Mmus"), in.name=c("symbol","ensg"), out.org=c("Hsap","Mmus"), out.name=c("symbol","ensg")) {
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
