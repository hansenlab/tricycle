#' Run PCA on Gene Ontology cell cycle genes
#'
#' Run PCA on Gene Ontology cell cycle genes abd get a new \linkS4class{SingleCellExperiment}. 
#' User could use this function to learn new reference projection matrix.
#'
#' @param sce.o A \linkS4class{SingleCellExperiment} contains library size normalized **log-expression** matrix.
#' @param gname Alternative rownames of \code{sce.o}. If provided, this will be used to map genes within Gene Ontology cell cycle gene list.
#' If not provided, the rownames of \code{sce.o} will be used instead. Default: NULL
#' @param exprs_values Integer scalar or string indicating which assay of \code{sce.o} contains the **log-expression** values, which will be used to run PCA. 
#' Default: 'logcounts'
#' @param gname.type The type of gene names as in \code{gname} or rownames of \code{sce.o}. It can be either 'ENSEMBL' or 'SYMBOL'. Default: 'ENSEMBL'
#' @param species The type of species in \code{sce.o}. It can be either 'mouse' or 'human'. The corresponding AnnotationDb object will be
#'  \code{\link[org.Mm.eg.db]{org.Mm.eg.db}} and \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}. Default: 'mouse'
#' @param ntop The number of genes with highest variance to use when calculating PCA, as in \code{\link[scater]{calculatePCA}}. Default: 500
#' @param ncomponents The number of component components to obtain, as in \code{\link[scater]{calculatePCA}}. Default: 20
#' @param runSeuratBy If the value is not NULL, the subsetted data (only GO cell cycle genes) will be integrated by Seurat. See \code{\link[Seurat]{IntegrateData}}.
#' This value should be set as a column name in \code{colData(sce.o)}. Such as 'sample' or 'batch'. Default: NULL
#' @param nfeatures The number of highly variable features within each sample/batch, which will only be used when integrating with Seurat.
#'  See \code{\link[Seurat]{FindVariableFeatures}}. Default: 500
#' @param anchor.features The number of genes used in anchor finding during Seurat integration. See \code{\link[Seurat]{FindIntegrationAnchors}}. Default: 2000
#' @param name String specifying the name to be used to store the result in the \code{\link[SingleCellExperiment]{reducedDims}} of the output. Default: 'PCA'
#'
#' @details
#' The function require an output of a \linkS4class{SingleCellExperiment} object which contains the library size normalized **log-expression** matrix. The full dataset will 
#' be subsetted to genes in the Gene Ontology cell cycle gene list (GO:0007049). The corresponding AnnotationDb object will be
#'  \code{\link[org.Mm.eg.db]{org.Mm.eg.db}} and \code{\link[org.Hs.eg.db]{org.Hs.eg.db}} for mouse and human respectively. If \code{runSeuratBy} is set, the data will be 
#'  integrated to remove batch effect between samples/batches by Seurat.
#'  
#' User can use this function to make new reference projection matrix by getting the 'rotation' attribute in PCA results. Such as 
#' \code{attr(reducedDim(sce.o, 'PCA'), 'rotation')[, 1:2]}. See examples for more details.
#'
#' @return
#' A subset \linkS4class{SingleCellExperiment} object with only GO cell cycle genes will be return.
#' The PCA resulting will be save in reducedDims with chosen name \code{\link[SingleCellExperiment]{reducedDims}(..., name)}. 
#' If Seurat integration is performed, another reducedDims with name 'matched.'+\code{name} will also be included in the \linkS4class{SingleCellExperiment}.
#'
#' @name pcaGoCC
#' 
#' @author Shijie C. Zheng
#'
#' @examples
#' gocc_sce.o <- pcaGoCC(example_sce)
#' new.ref <- attr(reducedDim(gocc_sce.o, 'PCA'), 'rotation')[, seq_len(2)]
#' example_sce <- inferCCTime(example_sce)  ### Use internal NeuroRef to project and infer CC time
#' 
#' ### Use new reference to project and infer CC time
#' new_sce <- inferCCTime(example_sce, ref.m  = new.ref, dimred = 'ccProjection2') 
#' plot(example_sce$CCTime, new_sce$CCTime)
NULL


#' @importFrom  org.Hs.eg.db org.Hs.eg.db
#' @importFrom  org.Mm.eg.db org.Mm.eg.db
#' @importFrom AnnotationDbi select
#' @importFrom scater runPCA calculatePCA
#' @importFrom SummarizedExperiment colData
.pcaGoCC <- function(sce.o, gname = NULL, exprs_values = "logcounts", gname.type = c("ENSEMBL", "SYMBOL"), species = c("mouse", "human"), ntop = 500, ncomponents = 20,  
    runSeuratBy = NULL, nfeatures = 500, anchor.features = 2000, name = "PCA") {
    species <- match.arg(species)
    gname.type <- match.arg(gname.type)
    if (species == "mouse") {
        AnnotationDb <- org.Mm.eg.db::org.Mm.eg.db
    } else {
        AnnotationDb <- org.Hs.eg.db::org.Hs.eg.db
    }
    if (is.null(gname)) {
        message("No gname input. Rownames of sce.o will be used.")
        gname <- rownames(sce.o)
    } else {
        if (nrow(sce.o) != length(gname)) 
            stop("gname does not match nrow sce.o!")
    }
    
    cycle.anno <- AnnotationDbi::select(AnnotationDb, keytype = "GOALL", keys = "GO:0007049", columns = gname.type)[, gname.type]
    gene.idx <- which(gname %in% cycle.anno)
    if (length(gene.idx) < 100) 
        stop("Less than 100 Gene Ontology cell cycle genes found in your data. Check you data or gname input.")
    message(paste0(length(gene.idx), " out of ", length(cycle.anno), " Gene Ontology cell cycle genes found in your data."))
    GO.o <- sce.o[gene.idx, ]
    
    GO.o <- runPCA(GO.o, exprs_values = exprs_values, name = name, ntop = ntop, ncomponents = ncomponents)
    
    ### merge samples
    if (!is.null(runSeuratBy)) {
        if (!(runSeuratBy %in% names(colData(GO.o)))) 
            stop(paste0(runSeuratBy, " does not exist in colData(sce.o)."))
        corrected.m <- .seuratIntegrate(count.m = 2^(assay(GO.o, exprs_values)) - 1, batch = colData(GO.o)[, runSeuratBy], nfeatures = nfeatures, anchor.features = anchor.features)
        reducedDim(GO.o, paste0("matched.", name)) <- calculatePCA(corrected.m, ntop = ntop, ncomponents = ncomponents)
        GO.o@metadata <- c(GO.o@metadata, list(seurat.corrected.m = corrected.m))
    }
    return(GO.o)
}

#' @importFrom Seurat CreateSeuratObject SplitObject NormalizeData FindVariableFeatures FindIntegrationAnchors IntegrateData
.seuratIntegrate <- function(count.m, batch, nfeatures, anchor.features) {
    seurat.o <- CreateSeuratObject(counts = count.m)
    seurat.o[["batch"]] <- batch
    seurat.list <- lapply(SplitObject(seurat.o, split.by = "batch"), function(x) {
        x <- FindVariableFeatures(NormalizeData(x, verbose = FALSE), nfeatures = nfeatures, verbose = FALSE)
        return(x)
    })
    seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = anchor.features)
    seurat.integrated <- IntegrateData(anchorset = seurat.anchors)
    corrected.m <- seurat.integrated@assays$integrated@data
    return(corrected.m)
}

#' @export
#' @rdname pcaGoCC
setMethod("pcaGoCC", "SingleCellExperiment", function(sce.o, exprs_values = "logcounts", 
                                                      gname.type = c("ENSEMBL", "SYMBOL"), species = c("mouse", "human"),
                                                      ntop = 500, ncomponents = 20, 
                                                      runSeuratBy = NULL, nfeatures = 500, anchor.features = 2000, name = "PCA") {
    .pcaGoCC(sce.o, exprs_values = exprs_values, gname.type = gname.type, species = species, 
             ntop = ntop, ncomponents = ncomponents, 
             runSeuratBy = runSeuratBy, nfeatures = nfeatures, anchor.features = anchor.features, name = name)
})






