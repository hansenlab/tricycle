#' 5 stage cell cycle gene marker list from Revelio
#'
#' This 5 stage cell cycle gene marker list is directly from \code{Revelio} package.
#' Within the list, 5 vectors corresponds to highly expressed genes at the cell cycle stage.
#' The genes are given as the human gene SYMBOLS.
#' This gene list is originally from the Whitfield et al.(2020).
#'
#' @name RevelioGeneList
#' @docType data
#'
#' @usage data(RevelioGeneList)
#'
#' @format An list of 5 string vector. The names of the elements are the names of cell cycle stages with the order of:
#' G1S, S, G2, G2M, MG1.
#'
#' @keywords datasets
#'
#' @references
#' Whitfield ML, et al.
#' \emph{Identification of Genes Periodically Expressed in the Human Cell Cycle and Their Expression in Tumors.}
#' Molecular Biology of the Cell (2002) 13: 1977-2000
#' doi:\href{https://doi.org/10.1091/mbc.02-02-0030}{10.1091/mbc.02-02-0030}.
#'
#' @importFrom utils data
#' @examples
#' data(RevelioGeneList)
NULL
