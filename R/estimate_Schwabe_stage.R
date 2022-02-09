#' Assign cell cycle stages using Schwabe method
#'
#' The function is a re-implementation of cell cycle stage assignment method
#' proposed in Schwabe et al.(2020), with a little modification. The core
#'  assignment method is not designed by the authors of this package!
#'
#' @param x A numeric matrix of **log-expression** values where rows are features and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such a matrix.
#' @param exprs_values Integer scalar or string indicating which assay of \code{x} contains the **log-expression** values, which will be used for projection.
#' If the projection already exists, you can ignore this value. Default: 'logcounts'
#' @param batch.v A string specifies which column in colData of \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} to use as the batch variable.
#' Or it can be a vector, of which the number of elements equals to the number of columns of \code{x}. The 5 stage cell cycle assignments are preformed for each batch separately.
#' No \code{NA} is permitted. Default: NULL
#' @param cycleGene.l A list contains the marker genes for each stage. The stage names should be included as names of the elements. If user feed custom list,
#' they should make sure that the same gene id type for \code{x} and \code{cycleGene.l}. If not custom list is given, \code{\link{RevelioGeneList}} will be used.
#' Default: NULL
#' @param gname Alternative rownames of \code{x}. If provided, this will be used to map genes within \code{x} with genes in \code{ref.m}.
#' If not provided, the rownames of \code{x} will be used instead. Default: NULL
#' @param gname.type The type of gene names as in \code{gname} or rownames of \code{x}. It can be either 'ENSEMBL' or 'SYMBOL'. If the user uses
#' custom \code{ref.m}, this value will have no effect. Default: 'ENSEMBL'
#' @param species The type of species in \code{x}. It can be either 'mouse' or 'human'. If the user uses
#' custom \code{cycleGene.l}, this value will have no effect. Default: 'mouse'
#' @param corThres For each batch and each stage, correlations between expression of each gene and the mean of all genes belonging to that stage
#' will be calculated to filter the final gene list used for inference. The genes with a correlation between \code{corThres} will not be used for calculating \emph{z}-scores.
#' Default: 0.2
#' @param tolerance For each cell, the function will compare the largest two \emph{z}-scores. If the difference between those two \emph{z}-scores is less than \code{tolerance},
#' the cell will be treated un-assignable with \code{NA} value returned for that cell. Default: 0.3
#' @param AnnotationDb An AnnotationDb objects. It is used to map ENSEMBL IDs to gene SYMBOLs.
#'  If no AnnotationDb object being given, the function will use \code{\link[org.Hs.eg.db]{org.Hs.eg.db}} or \code{\link[org.Mm.eg.db]{org.Mm.eg.db}} for human and mouse respectively.
#' @param altexp String or integer scalar specifying an alternative experiment containing the **log-expression** data, which will be used for projection.
#' If the projection is already calculated and stored in the \linkS4class{SingleCellExperiment} as a dimred, leave this value to default NULL.
#'
#' @details
#' The function is a re-implementation of cell cycle stage assignment method
#' proposed in Schwabe et al.(2020), with a little modification. We include 
#' this function only for the purpose of convenience. The core assignment method
#' is not designed by the authors of this package!
#' Breiefly, the function assigns cells to discretized cell cycle stages by
#'  comparing the \emph{z}-scores calculated for each stage markers.
#' Without cycleGene.l input, \code{\link{RevelioGeneList}} will be used.
#' If you use this function, you should cite Schwabe et al.(2020).
#'
#' @return
#' If the input is a numeric matrix, the discretized cell cycle stages - a factor vector corresponding to each cell will be returned.
#'
#' If the input is \linkS4class{SummarizedExperiment}, the original \linkS4class{SummarizedExperiment} with the discretized cell cycle stages stored in colData with name 'CCStage' will be returned.
#'
#' If the input is \linkS4class{SingleCellExperiment}, the original \linkS4class{SingleCellExperiment} with the discretized cell cycle stages stored in colData with name 'CCStage' will be returned.
#'
#' @name estimate_Schwabe_stage
#' @aliases estimate_Schwabe_stage
#'
#'
#' @author Shijie C. Zheng
#'
#' @references
#' Schwabe D, et al.
#' \emph{The transcriptome dynamics of single cells during the cell cycle.}
#' Molecular Systems Biology (2020) 16: e9946
#' doi:\href{https://doi.org/10.15252/msb.20209946}{10.15252/msb.20209946}.
#'
#' Zheng SC, et al.
#' \emph{Universal prediction of cell cycle position using transfer learning.}
#'
#' @examples
#' data(neurosphere_example, package = "tricycle")
#' neurosphere_example <- estimate_Schwabe_stage(neurosphere_example,
#'  gname.type = "ENSEMBL", species = "mouse")
#' neurosphere_example2 <- estimate_Schwabe_stage(neurosphere_example, batch.v = "sample")
#' neurosphere_example3 <- estimate_Schwabe_stage(neurosphere_example,
#'  batch.v = neurosphere_example$sample)
#' neurosphere_example <- project_cycle_space(neurosphere_example)
#' plot(reducedDim(neurosphere_example, "tricycleEmbedding"),
#'  col = neurosphere_example$CCStage)
NULL


#' @importFrom stats cor
#' @importMethodsFrom S4Vectors  "%in%"  do.call  levels  t
.CCStage <- function(data.m, batch.v = NULL, cycleGene.l = NULL, corThres = 0.2, tolerance = 0.3) {
    if (is.null(batch.v)) {
          batch.v <- rep(1, ncol(data.m))
      }
    if ((ncol(data.m) != length(batch.v)) | any(is.na(batch.v))) {
          stop("The length of batch.v does not equal to the number of column in
               data or there is NA in bathc.v, which is not allowed.")
      }
    cycleGene.l <- lapply(cycleGene.l, intersect, rownames(data.m))
    allgenes.v <- Reduce(union, cycleGene.l)
    if (sum(allgenes.v %in% rownames(data.m)) < 30) {
          stop("Less than 30 cell cycle gene found. Not enough infomation to assgin the 5 stages.")
      }
    data.m <- as.matrix(data.m[allgenes.v, ])

    cc.v <- rep(NA_character_, ncol(data.m))
    for (b in unique(batch.v)) {
        idx <- which(batch.v == b)
        dat <- data.m[, idx]
        score.m <- t(scale(t(scale(do.call(cbind, lapply(seq_along(cycleGene.l), function(i) {
            gene <- cycleGene.l[[i]]
            gene.sum <- rowSums(dat[gene, ])
            if (sum(gene.sum > 0) > 3) {
                gene <- gene[which(gene.sum > 0)]
            } else {
                stop(paste0("Less than 3 ", names(cycleGene.l)[i], " stage genes expressed in bacth ", b, ". \n Consider changing batch setting."))
            }
            mean.v <- colMeans(dat[gene, ])
            cor.v <- cor(as.matrix(t(dat[gene, ])), as.vector(mean.v))
            if (sum(cor.v > corThres) < 3) {
                  warning(paste0("Batch ", b, " phase ", names(cycleGene.l)[i], " too few genes (<3)."))
              }
            message(paste0("Batch ", b, " phase ", names(cycleGene.l)[i], " gene:", sum(cor.v > corThres), "\n"))
            return(colMeans(dat[gene[cor.v > corThres], ]))
        }))))))
        cc.v[idx] <- apply(score.m, MARGIN = 1, function(s) {
            o.v <- order(s, decreasing = TRUE)
            eta <- o.v[1]
            eta_ <- o.v[2]
            if (((abs(eta - eta_) > 1) & (abs(eta - eta_) < (length(cycleGene.l) - 1))) | ((s[eta] - s[eta_]) < tolerance)) {
                  return(NA_character_)
              }
            return(names(cycleGene.l)[eta])
        })
    }
    cc.v <- factor(cc.v, levels = names(cycleGene.l))
    return(cc.v)
}


#' @importFrom AnnotationDbi mapIds
#' @importMethodsFrom AnnotationDbi colnames get ncol nrow
#' @importMethodsFrom GenomicRanges intersect union
.getSYMBOL <- function(gname, species = c("mouse", "human"), AnnotationDb = NULL) {
    species <- match.arg(species)
    AnnotationDb <- .getAnnotationDB(AnnotationDb, species)
    SYMBOL <- suppressMessages(toupper(AnnotationDbi::mapIds(AnnotationDb, keys = gname, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")))
    return(SYMBOL)
}

#' @importMethodsFrom IRanges as.matrix "colnames<-"  diff lapply rownames "rownames<-"  table tolower toupper unique which
.estimate_Schwabe_stage <- function(data.m, cycleGene.l = NULL, gname = NULL, gname.type = c("ENSEMBL", "SYMBOL"), species = c("mouse", "human"), AnnotationDb = NULL, ...) {
    species <- match.arg(species)
    gname.type <- match.arg(gname.type)
    if (is.null(gname)) {
        message("No gname input. Rownames of x will be used.")
        gname <- rownames(data.m)
    } else {
        if (nrow(data.m) != length(gname)) stop("gname does not match nrow x!")
    	  rownames(data.m) <- gname
    }
    if (is.null(cycleGene.l)) {
        if (species == "mouse") {
            if (gname.type == "SYMBOL") {
                rownames(data.m) <- toupper(rownames(data.m))
            } else {
                rownames(data.m) <- .getSYMBOL(gname, species = "mouse", AnnotationDb = AnnotationDb)
            }
        } else {
            if (gname.type == "ENSEMBL") {
                  rownames(data.m) <- .getSYMBOL(gname, species = "human", AnnotationDb = AnnotationDb)
              }
        }
        data_env <- new.env(parent = emptyenv())
        data("RevelioGeneList", package = "tricycle", envir = data_env)
        cycleGene.l <- data_env[["RevelioGeneList"]]
    }

    return(.CCStage(data.m, cycleGene.l = cycleGene.l, ...))
}


#' @export
#' @rdname estimate_Schwabe_stage
#' @importMethodsFrom SummarizedExperiment  assay colData  order sort
#' @importMethodsFrom SingleCellExperiment altExp cbind rbind reducedDim "reducedDim<-"  reducedDimNames
#' 
estimate_Schwabe_stage <- function(x, exprs_values = "logcounts", batch.v = NULL,
                                   altexp = NULL, cycleGene.l = NULL, gname = NULL,
                                   gname.type = c("ENSEMBL", "SYMBOL"),
                                   species = c("mouse", "human"), AnnotationDb = NULL,
                                   corThres = 0.2, tolerance = 0.3) {
	message("This function is a re-implementation of Schwabe et al. 2020. If you want to use tricycle method, please run estimate_cycle_position!")
  if (is(x, "SingleCellExperiment")) {
    if (!is.null(batch.v)) {
      if ((length(batch.v) == 1) & all(batch.v %in% names(colData(x)))) {
        batch.v <- colData(x)[, batch.v]
      }
    }
    if (!is.null(altexp)) {
      y <- altExp(x, altexp)
    } else {
      y <- x
    }
    x$CCStage <- .estimate_Schwabe_stage(assay(y, exprs_values), batch.v = batch.v,
                                         cycleGene.l = cycleGene.l, gname = gname,
                                         gname.type = gname.type, species = species,
                                         AnnotationDb = AnnotationDb, corThres = corThres,
                                         tolerance = tolerance)
    out <- x
  } else if (is(x, "SummarizedExperiment")) {
    if (!is.null(batch.v)) {
      if ((length(batch.v) == 1) & all(batch.v %in% names(colData(x)))) {
        batch.v <- colData(x)[, batch.v]
      }
    }
    x$CCStage <- .estimate_Schwabe_stage(assay(x, exprs_values), batch.v = batch.v,
                                         cycleGene.l = cycleGene.l, gname = gname,
                                         gname.type = gname.type, species = species,
                                         AnnotationDb = AnnotationDb, corThres = corThres,
                                         tolerance = tolerance)
    out <- x
  } else {
    out <- .estimate_Schwabe_stage(x, batch.v = batch.v,
                                   cycleGene.l = cycleGene.l, gname = gname,
                                   gname.type = gname.type, species = species,
                                   AnnotationDb = AnnotationDb, corThres = corThres,
                                   tolerance = tolerance)
  }
  return(out)
}

