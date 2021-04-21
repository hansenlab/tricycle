#' Pre-learned reference projection matrix from the Neurosphere dataset
#'
#' Default reference projection matrix learned from the Neurosphere dataset.
#'
#' @name neuroRef
#' @docType data
#'
#' @usage data(neuroRef)
#'
#' @format An object of class `'data.frame'`, with 5 variables. Normally, user won't call this
#' data directly as it will be automatically used in \code{\link{project_cycle_space}} if no custom reference
#' projection matrix is provided. Each row is gene, and rotation scores for PC1 and PC2, mouse ENSEMBL IDs,
#' and mouse gene SYMBOLs are included. The `SYMBOL` values are just the upper-case `symbol` values.
#'
#' @keywords datasets
#'
#' @references
#' Zheng SC, et al.
#' \emph{Universal prediction of cell cycle position using transfer learning.}
#'
#' @examples
#' data(neuroRef)
NULL
