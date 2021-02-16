#' Project data into the cell cycle pattern space
#'
#' Project mouse and human single cell RNAseq data into a cell cycle embedding by a pre-learned reference projection matrix.
#'
#' @param x A numeric matrix of **log-expression** values where rows are features and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param exprs_values Integer scalar or string indicating which assay of \code{x} contains the **log-expression** values. Default: "logcounts"
#' @param ref.m A custom reference projection matrix to project the new data, where rows are features and columns are dimensions.
#' Users need to use the same type of \code{gname}(or rownames of \code{x}) as for the \code{ref.m}.
#' If no custom ref.m is given, the internal reference \code{neuroRef} will be used.
#' @param gname Alternative rownames of \code{x}. If provided, this will be used to map genes within \code{x} with genes in \code{ref.m}.
#' If not provided, the rownames of \code{x} will be used instead.
#' @param gname.type The type of gene names as in \code{gname} or rownames of \code{x}. It can be either "ENSEMBL" or "SYMBOL". If the user uses
#' custom \code{ref.m}, this value will have no effect. Default: "ENSEMBL"
#' @param species The type of species in \code{x}. It can be either "mouse" or "human". If the user uses
#' custom \code{ref.m}, this value will have no effect. Default: "mouse"
#' @param AnnotationDb An AnnotationDb objects. If the user uses the internal reference to project human data,
#'  and provide rownames in the format of Ensembl IDs, this object will be used to map Ensembl IDs to gene SYMBOLs.
#'  If no AnnotationDb object being given, the function will use \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}.
#' @param altexp String or integer scalar specifying an alternative experiment containing the input data.
#' @param name String specifying the name to be used to store the result in the \code{\link[SingleCellExperiment]{reducedDims}} of the output. Default: "ccProjection"
#'
#' @details
#' The function will use pre-learned cell cycle pattern to project new data to show the cell cycle progression. If the user uses internal Neuropshere reference,
#' the expression values must be **log-transformed**. Besides, we would assume the input data has been already preprocessed, library size normalized at least.
#' The projection process is to take sum of weighted mean-centered expression of chosen genes, so the mean expression of a given gene could be affected without library size normalization.
#'
#' @return
#' If the input is a numeric matrix or a \linkS4class{SummarizedExperiment}, a projection matrix with rows cells and column dimensions will be returned.
#' The actual rotation matrix used to project the data is included in the attributes with name "rotation".
#'
#' For \linkS4class{SingleCellExperiment}, an updated \linkS4class{SingleCellExperiment} is returned containing projection matrix in \code{\link[SingleCellExperiment]{reducedDims}(..., name)}.
#'
#' @name projectCC
#' @seealso
#' \code{\link{inferCCTime}}, for inferring cell cycle time.
#'
#' @author Shijie C. Zheng
#'
#' @examples
#' example_sce <- projectCC(example_sce)
#' reducedDimNames(example_sce)
#' head(reducedDim(example_sce, "ccProjection"))
#' plot(reducedDim(example_sce, "ccProjection"))
#' names(attributes(reducedDim(example_sce, "ccProjection")))
NULL




.projectCC <- function(data.m, ref.m = NULL, gname = NULL, gname.type = c("ENSEMBL", "SYMBOL"), species = c("mouse", "human"), AnnotationDb = NULL){
	species <- match.arg(species)
	gname.type <- match.arg(gname.type)

	if (!is.null(gname)) {
		rownames(data.m) <- gname
	}
	if (is.null(ref.m)) {
		message("No custom reference projection matrix provided. The ref learned from mouse Neuroshpere data will be used.")
		ref.m <- .getRotation(neuroRef, gname.type = gname.type, species = species)

		if (species == "human" & gname.type == "ENSEMBL") {
			message("As the reference data was learned from mouse, we will map the human ENSEMBL id to gene SYMBOL.")
			rownames(data.m) <- .humanSymbol(f.id = rownames(data.m), AnnotationDb = AnnotationDb)
		}
	}
	.calProjection(data.m, ref.m)
}


.calProjection <- function(data.m, rotation.m) {
	genes <- intersect(rownames(data.m), rownames(rotation.m))

	if (length(genes) == 0) stop("None genes found in new data. This could be caused by wrong input of rownames type.")
	message(paste0("The number of projection genes found in the new data is ", length(genes), "."))

	rotation.m <- rotation.m[genes, ]
	data.m <- data.m[genes, ]
	projection.m <- scale(t(as.matrix(data.m)), center = T, scale = F) %*% rotation.m
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
		rownames(rotation.m) <- neuroRef[, tolower(gname.type)]
	}
	colnames(rotation.m) <- c("PC1", "PC2")
	return(rotation.m)
}

#' @importFrom  org.Hs.eg.db org.Hs.eg.db
#' @importFrom AnnotationDbi select
.humanSymbol <- function(f.id, AnnotationDb = NULL) {
	if (is.null(AnnotationDb)) {
		AnnotationDb <- org.Hs.eg.db::org.Hs.eg.db
		message("No AnnotationDb desginated. org.Hs.eg.db will be used to map Human ENSEMBL id to gene SYMBOL.")
	}
	SYMBOL <- AnnotationDbi::select(AnnotationDb, keys = f.id, columns ="SYMBOL", keytype = "ENSEMBL")[["SYMBOL"]]
	return(SYMBOL)
}


#' @export
#' @rdname projectCC
setMethod("projectCC", "ANY", function(x, ...) {
	.projectCC(x, ...)
})

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









