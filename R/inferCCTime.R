#' Assign cell cycle time
#'
#' Assign cell cycle time by the angle formed by PC1 and PC2 in the cell cycle space.
#' If the cell cycle projection does not exist, the function will project the data.
#'
#' @param x A numeric matrix of **log-expression** values where rows are features and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param exprs_values Integer scalar or string indicating which assay of \code{x} contains the **log-expression** values, which will be used for projection.
#' If the projection already exists, you can ignore this value. Default: 'logcounts'
#' @param ... Arguments to be used by \code{\link{projectCC}}. If \code{x} is a \linkS4class{SingleCellExperiment}, and the projection is
#' already in the reducedDim with name \code{dimred}. The \code{dimred} will be directly used to assign cell cycle time withou new projection.
#' @param dimred The name of reducedDims in  \linkS4class{SingleCellExperiment} (\code{\link[SingleCellExperiment]{reducedDims}}). If the \code{dimred} already exists,
#' it will be used to assign cell cycle time. If \code{dimred} does not exist, the projection will be calculated first by \code{\link{projectCC}}
#' and stored with name \code{dimred} in {x}. Default: 'ccProjection'
#' @param center.pc1 The center of PC1 when defining the angle. Default: 0
#' @param center.pc2 The center of PC2 when defining the angle. Default: 0
#' @param altexp String or integer scalar specifying an alternative experiment containing the **log-expression** data, which will be used for projection.
#' If the projection is already calculated and stored in the \linkS4class{SingleCellExperiment} as a dimred, leave this value to default NULL.
#'
#' @details
#' The function will use assign the cell cycle time by the angle formed by the PC1 and PC2 of cell cycle projections.
#' If the input is a numeric matrix or a \linkS4class{SummarizedExperiment}, the projection will be calculated with the input **log-expression** values.
#' For \linkS4class{SingleCellExperiment}, the projection will also be calculated if the designated \code{dimred} does not exist.
#' Ohterwise, the \code{dimred} will directly be used to assign cell cycle time.
#' Therefore, this function is a wrapper function if the input is a \linkS4class{SingleCellExperiment}.
#' Refer to \code{\link{projectCC}} to all arguments during the projection.
#'
#' @return
#' If the input is a numeric matrix, the cell cycle time - a numeric vector bound between \eqn{0 \sim 2 \pi} with the same length as the number of input coumlum will be returned.
#'
#' If the input is \linkS4class{SummarizedExperiment}, the original \linkS4class{SummarizedExperiment} with cell cycle time stored in colData with name 'CCTime' will be returned.
#'
#' If the input is \linkS4class{SingleCellExperiment}, the original \linkS4class{SingleCellExperiment} with cell cycle time stored in colData with name 'CCTime' will be returned
#' and the projection will be stored in \code{\link[SingleCellExperiment]{reducedDims}(..., dimred)} if it does not exist before.
#'
#' @name inferCCTime
#' @seealso
#' \code{\link{projectCC}}, for projecting new data with a pre-learned reference
#'
#' @author Shijie C. Zheng
#'
#' @examples
#' example_sce <- inferCCTime(example_sce)
#' reducedDimNames(example_sce)
#' plot(reducedDim(example_sce, 'ccProjection'))
#' plot(example_sce$CCTime, reducedDim(example_sce, 'ccProjection')[, 1])
#' plot(example_sce$CCTime, reducedDim(example_sce, 'ccProjection')[, 2])
NULL




#' @importFrom circular coord2rad
.getTheta <- function(pc1pc2.m, center.pc1 = 0, center.pc2 = 0) {
    pc1pc2.m[, 1] <- pc1pc2.m[, 1] - center.pc1
    pc1pc2.m[, 2] <- pc1pc2.m[, 2] - center.pc2
    as.numeric(coord2rad(pc1pc2.m[, 1:2]))
}

#' @export
#' @rdname inferCCTime
setMethod("inferCCTime", "ANY", function(x, ..., center.pc1 = 0, center.pc2 = 0) {
    projection.m <- .projectCC(x, ...)
    .getTheta(projection.m, center.pc1 = center.pc1, center.pc2 = center.pc2)
})

#' @export
#' @rdname inferCCTime
#' @importFrom SummarizedExperiment assay
setMethod("inferCCTime", "SummarizedExperiment", function(x, ..., exprs_values = "logcounts", center.pc1 = 0, center.pc2 = 0) {
    projection.m <- .projectCC(assay(x, exprs_values), ...)
    x$CCTime <- .getTheta(projection.m, center.pc1 = center.pc1, center.pc2 = center.pc2)
    x
})


#' @export
#' @rdname inferCCTime
#' @importFrom SingleCellExperiment reducedDim<- altExp reducedDimNames
setMethod("inferCCTime", "SingleCellExperiment", function(x, ..., dimred = "ccProjection", center.pc1 = 0, center.pc2 = 0, altexp = NULL) {
    if (!is.null(altexp)) {
        y <- altExp(x, altexp)
        if ((dimred %in% reducedDimNames(x)) & (!(dimred %in% reducedDimNames(y)))) {
            warning(paste0(dimred, "exists in x, but not in ", altexp, ". \nThe projecion will be recalculated using ", altexp, ".\n If you intend to use pre-calculated projection, please don't set ", 
                altexp))
        }
        
    } else {
        y <- x
    }
    if (!(dimred %in% reducedDimNames(y))) {
        message(paste0("The designated dimred do not exist in the SingleCellExperiment or in altexp. projectCC will be run to calculate embedding ", dimred))
        y <- projectCC(y, name = dimred, ...)
        reducedDim(x, dimred) <- reducedDim(y, dimred)
    }
    x$CCTime <- .getTheta(reducedDim(x, dimred), center.pc1 = center.pc1, center.pc2 = center.pc2)
    x
})


