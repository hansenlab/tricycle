#' Assign cell cycle position
#'
#' Assign cell cycle position by the angle formed by PC1 and PC2 in the cell cycle space.
#' If the cell cycle projection does not exist, the function will project the data.
#'
#' @param x A numeric matrix of **log-expression** values where rows are features and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param exprs_values Integer scalar or string indicating which assay of \code{x} contains the **log-expression** values, which will be used for projection.
#' If the projection already exists, you can ignore this value. Default: 'logcounts'
#' @param ... Arguments to be used by \code{\link{project_cycle_space}}. If \code{x} is a \linkS4class{SingleCellExperiment}, and the projection is
#' already in the reducedDim with name \code{dimred}. The \code{dimred} will be directly used to assign cell cycle position withou new projection.
#' @param dimred The name of reducedDims in  \linkS4class{SingleCellExperiment} (\code{\link[SingleCellExperiment]{reducedDims}}). If the \code{dimred} already exists,
#' it will be used to assign cell cycle position. If \code{dimred} does not exist, the projection will be calculated first by \code{\link{project_cycle_space}}
#' and stored with name \code{dimred} in {x}. Default: 'tricycleEmbedding'
#' @param center.pc1 The center of PC1 when defining the angle. Default: 0
#' @param center.pc2 The center of PC2 when defining the angle. Default: 0
#' @param altexp String or integer scalar specifying an alternative experiment containing the **log-expression** data, which will be used for projection.
#' If the projection is already calculated and stored in the \linkS4class{SingleCellExperiment} as a dimred, leave this value to default NULL.
#'
#' @details
#' The function will use assign the cell cycle position by the angle formed by the PC1 and PC2 of cell cycle projections.
#' If the input is a numeric matrix or a \linkS4class{SummarizedExperiment}, the projection will be calculated with the input **log-expression** values.
#' For \linkS4class{SingleCellExperiment}, the projection will also be calculated if the designated \code{dimred} does not exist.
#' Ohterwise, the \code{dimred} will directly be used to assign cell cycle position.
#' Therefore, this function is a wrapper function if the input is a \linkS4class{SingleCellExperiment}.
#' Refer to \code{\link{project_cycle_space}} to all arguments during the projection.
#' 
#' The estimated cell cycle position is bound between 0 and 2pi.
#' Note that we strive to get high resolution of cell cycle state, and we think
#' the continuous position is more appropriate when describing the cell cycle.
#' However, to help users understand the position variable, we also note that 
#' users can approximately relate 0.5pi to be the start of S stage, pi to be the
#' start of G2M stage, 1.5pi to be the middle of M stage, and 1.75pi-0.25pi to 
#' be G1/G0 stage.
#'
#' @return
#' If the input is a numeric matrix, the cell cycle position - a numeric vector bound between \eqn{0 \sim 2 \pi} with the same length as the number of input coumlum will be returned.
#'
#' If the input is \linkS4class{SummarizedExperiment}, the original \linkS4class{SummarizedExperiment} with cell cycle position stored in colData with name 'tricyclePosition' will be returned.
#'
#' If the input is \linkS4class{SingleCellExperiment}, the original \linkS4class{SingleCellExperiment} with cell cycle position stored in colData with name 'tricyclePosition' will be returned
#' and the projection will be stored in \code{\link[SingleCellExperiment]{reducedDims}(..., dimred)} if it does not exist before.
#'
#' @name estimate_cycle_position
#' @seealso
#' \code{\link{project_cycle_space}}, for projecting new data with a pre-learned reference
#'
#' @references
#' Zheng SC, et al.
#' \emph{Universal prediction of cell cycle position using transfer learning.}
#'
#' @author Shijie C. Zheng
#'
#' @examples
#' data(neurosphere_example, package = "tricycle")
#' neurosphere_example <- estimate_cycle_position(neurosphere_example)
#' reducedDimNames(neurosphere_example)
#' plot(reducedDim(neurosphere_example, "tricycleEmbedding"))
#' plot(neurosphere_example$tricyclePosition,
#'  reducedDim(neurosphere_example, "tricycleEmbedding")[, 1])
#' plot(neurosphere_example$tricyclePosition,
#'  reducedDim(neurosphere_example, "tricycleEmbedding")[, 2])
NULL

#' @importFrom circular coord2rad
.getTheta <- function(pc1pc2.m, center.pc1 = 0, center.pc2 = 0) {
    pc1pc2.m[, 1] <- pc1pc2.m[, 1] - center.pc1
    pc1pc2.m[, 2] <- pc1pc2.m[, 2] - center.pc2
    as.numeric(coord2rad(pc1pc2.m[, seq_len(2)]))
}


#' @export
#' @rdname estimate_cycle_position
#' @importFrom methods is
#' @importFrom SingleCellExperiment reducedDim<- altExp reducedDimNames reducedDim
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' 
estimate_cycle_position <- function(x, exprs_values = "logcounts", dimred = "tricycleEmbedding", center.pc1 = 0, center.pc2 = 0, altexp = NULL, ...) {
    if (is(x, "SingleCellExperiment")) {
        if (!is.null(altexp)) {
            y <- altExp(x, altexp)
            if ((dimred %in% reducedDimNames(x)) & (!(dimred %in% reducedDimNames(y)))) {
                warning(paste0(
                    dimred, "exists in x, but not in ", altexp, ". \nThe projecion will be recalculated using ", altexp, ".\n If you intend to use pre-calculated projection, please don't set ",
                    altexp
                ))
            }
        } else {
            y <- x
        }
        if (!(dimred %in% reducedDimNames(y))) {
            message(paste0("The designated dimred do not exist in the SingleCellExperiment or in altexp. project_cycle_space will be run to calculate embedding ", dimred))
            y <- project_cycle_space(y, name = dimred, exprs_values = exprs_values, ...)
            reducedDim(x, dimred) <- reducedDim(y, dimred)
        }
        x$tricyclePosition <- .getTheta(reducedDim(x, dimred), center.pc1 = center.pc1, center.pc2 = center.pc2)
        out <- x
    } else if (is(x, "SummarizedExperiment")) {
        projection.m <- .project_cycle_space(assay(x, exprs_values), ...)
        x$tricyclePosition <- .getTheta(projection.m, center.pc1 = center.pc1, center.pc2 = center.pc2)
        out <- x
    } else {
        projection.m <- .project_cycle_space(x, ...)
        out <- .getTheta(projection.m, center.pc1 = center.pc1, center.pc2 = center.pc2)
    }
    return(out)
}


