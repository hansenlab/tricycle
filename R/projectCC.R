#' Project new data into the cell cycle pattern space
#'
#' Perform xxxx
#'
#' @param x A numeric matrix of **log-expression** values where rows are features and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param exprs_values Integer scalar or string indicating which assay of \code{x} contains the **log-expression** values. Default: "logcounts"
#' @param ref.m A custom reference projection matrix to project the new data, where rows are features and columns are dimensions.
#' Users need to use the same type of \code{gname}(or rownames of \code{x}) as for the \code{ref.m}.
#' If no custom ref.m is given, the internal reference \code{neuroRef} will be used.
#' @param gname An alternative rownames of \code{x}. If provided, this will be used to map genes within \code{x} with genes in \code{ref.m}.
#' If not provided, the rownames of \code{x} will be used instead.
#' @param gname.type The type of gene names as in \code{gname} or rownames of \code{x}. It can be either "ensembl" or "symbol". If the user uses
#' custom \code{ref.m}, this value will have no effect. Default: "ensembl"
#' @param species The type of species in \code{x}. It can be either "mouse" or "human". If the user uses
#' custom \code{ref.m}, this value will have no effect. Default: "mouse"
#' @param AnnotationDb An AnnotationDb objects. If the user uses the internal reference to project human data,
#'  and provide rownames in the format of Ensembl IDs, this object will be used to map Ensembl IDs to gene symbols.
#'  If no AnnotationDb object being given, the function will use \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}.
#' @param altexp String or integer scalar specifying an alternative experiment containing the input data.
#' @param name String specifying the name to be used to store the result in the \code{\link{reducedDims}} of the output. Default: "ccProjection"
#'





.projectCC <- function(data.m, ref.m = NULL, gname = NULL, gname.type = c("ensembl", "symbol"), species = c("mouse", "human"), AnnotationDb = NULL){
	species <- match.arg(species)
	gname.type <- match.arg(gname.type)

	if (!is.null) {
		rownames(data.m) <- gname
	}
	if (is.null(ref.m)) {
		message("No custom reference projection matrix provided. The ref learned from mouse Neuroshpere data will be used.")
		ref.m <- .getRotation(neuroRef, gname.type = gname.type, species = species)

		if (species == "human" & gname.type == "ensembl") {
			message("As the reference data was learned from mouse, we will map the human ensembl id to gene symbol.")
			rownames(data.m) <- .humanSymbol(f.id = rownames(data.m), AnnotationDb = AnnotationDb)
		}
	}
	.calProjection(data.m, rotation.m)
}



.calProjection <- function(data.m, rotation.m) {
	genes <- intersect(rownames(data.m), rownames(rotation.m))

	if (length(genes) == 0) stop("None genes found in new data. This could be caused by wrong input of rownames type.")
	message(paste0("The number of projection genes found in the new data is ", length(genes), "."))

	rotation.m <- rotation.m[genes, ]
	data.m <- data.m[genes, ]
	projection.m <- scale(t(data.m), center = T, scale = F) %*% rotation.m
	rownames(projection.m) <- colnames(data.m)
	colnames(projection.m) <- colnames(rotation.m)
	attr(projection.m, "rotation") <- rotation.m

	return(projection.m)
}

.getRotation <- function(neuroRef, gname.type, species) {
	rotation.m <- as.matrix(neuroRef[, 1:2])
	if (species == "human") {
		rownames(rotation.m) <- neuroRef$SYMBOL
	} else {
		rownames(rotation.m) <- neuroRef[, gname.type]
	}
	colnames(rotation.m) <- c("PC1", "PC2")
	return(rotation.m)
}

#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi select
.humanSymbol <- function(f.id, AnnotationDb = NULL) {
	if (is.null(AnnotationDb)) {
		AnnotationDb <- org.Hs.eg.db
		message("No AnnotationDb desginated. org.Hs.eg.db will be used to map Human ENSEMBL id to gene symbol.")
	}
	SYMBOL <- AnnotationDbi::select(AnnotationDb, keys = f.id, columns ="SYMBOL", keytype = "ENSEMBL")[["SYMBOL"]]
	return(SYMBOL)
}




#' @export
#' @rdname projectCC
setMethod("projectCC", "ANY", .projectCC)

#' @export
#' @rdname projectCC
#' @importFrom SummarizedExperiment assay
setMethod("projectCC", "SummarizedExperiment", function(x, ..., exprs_values = "logcounts") {
	.projectCC(assay(x, exprs_values), ...)
})


#' @export
#' @rdname projectCC
#' @importFrom SingleCellExperiment reducedDim<- altExp
#' @importFrom SummarizedExperiment assay
setMethod("projectCC", "SingleCellExperiment", function(x, ..., exprs_values="logcounts", altexp=NULL, name="ccProjection")
{
	if (!is.null(altexp)) {
		y <- altExp(x, altexp)
	} else {
		y <- x
	}
	reducedDim(x, name) <- .projectCC(assay(y, exprs_values), ...)
	x
})









