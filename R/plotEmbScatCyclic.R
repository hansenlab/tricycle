#' Plot embedding with cyclic cell cycle time
#'
#' Generate scat plot of embedding with cyclic cell cycle time or other cyclic variables
#'
#' @param sce.o A \linkS4class{SingleCellExperiment} contains the embbing to be plotted against.
#' @param color_by AAA
#' @param facet_by AAA
#' @param dimred The name or index of reducedDims in  \linkS4class{SingleCellExperiment} (\code{\link[SingleCellExperiment]{reducedDims}}). Default: 1
#' @param dim DBBBBB. Default: 1:2
#' @param fig.title The title of the figure. Default: NULL
#' @param point.size  \code{\link[scattermore]{geom_scattermore}}
#' @param point.alpha  \code{\link[scattermore]{geom_scattermore}}
#' @param x_lab sss
#' @param y_lab sss
#' @param hue.colors sss
#' @param hue.n sss
#' @param plot.legend sss
#'
#' @details
#' BALABALA
#' 
#' @return
#' A ggplot object or a list of ggplot objects.
#'
#' @name plotEmbScatCyclic
#' 
#' @author Shijie C. Zheng
#'
#' @examples
#' example_sce <- inferCCTime(example_sce)
#' plotEmbScatCyclic(example_sce, point.size = 3.1, point.alpha = 0.8)
NULL


#' @importFrom scattermore geom_scattermore
#' @importFrom dplyr filter `%>%`
.plotEmbScatCyclic <- function(emb.m, color.value, color_by, facet_var = NULL, fig.title = NULL, point.size = 2.1, point.alpha = 0.6, x_lab = NULL, y_lab = NULL, hue.colors = c("#2E22EA", "#9E3DFB", 
    "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9"), hue.n = 500, plot.legend = FALSE) {
    if (is.null(fig.title)) 
        fig.title <- paste0("(n=", nrow(emb.m), ")")
    
    if (is.null(x_lab)) 
        x_lab <- colnames(emb.m)[1]
    if (is.null(y_lab)) 
        x_lab <- colnames(emb.m)[2]
    
    x_lim <- range(emb.m[, 1]) + c(diff(range(emb.m[, 1])) * c(-0.05, 0.05))
    y_lim <- range(emb.m[, 2]) + c(diff(range(emb.m[, 2])) * c(-0.05, 0.05))
    
    tmp.df <- data.frame(x = emb.m[, 1], y = emb.m[, 2], color = color.value)
    
    scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) + geom_scattermore(pointsize = point.size, alpha = point.alpha) + scale_color_gradientn(name = color_by, limits = range(0, 2 * 
        pi), breaks = seq(from = 0, to = 2 * pi, length.out = hue.n), colors = hue.colors, guide = FALSE) + labs(y = y_lab, x = x_lab, title = fig.title) + xlim(x_lim) + ylim(y_lim) + .gg_theme
    
    if (!is.null(facet_var)) {
        tmp.df$facet <- facet_var
        facet_labels <- levels(factor(tmp.df$facet))
        lp <- lapply(seq_len(nlevels(factor(tmp.df$facet))), function(idx) {
            p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) + geom_scattermore(data = tmp.df %>% dplyr::filter(facet != levels(factor(tmp.df$facet))[idx]), pointsize = point.size, alpha = 0.4, 
                color = "gray90", show.legend = FALSE) + geom_scattermore(data = tmp.df %>% dplyr::filter(facet == levels(factor(tmp.df$facet))[idx]), pointsize = point.size, alpha = point.alpha) + 
                scale_color_gradientn(name = color_by, limits = range(0, 2 * pi), breaks = seq(from = 0, to = 2 * pi, length.out = hue.n), colors = hue.colors, guide = FALSE) + labs(y = y_lab, 
                x = x_lab, title = facet_labels[idx]) + xlim(x_lim) + ylim(y_lim) + .gg_theme
            return(p)
        })
        return(c(list(scat.p), lp))
    }
    return(scat.p)
}


#' @export
#' @rdname plotEmbScatCyclic
setMethod("plotEmbScatCyclic", "SingleCellExperiment", function(sce.o, ..., color_by = "CCTime", facet_by = NULL, dimred = 1, dim = 1:2) {
    if (length(dim) != 2) 
        stop("The function can only plot 2 dims at this time. Change dim argument.")
    emb.m <- reducedDim(sce.o, dimred)[, dim]
    color.value <- colData(sce.o)[[color_by]]
    if (!is.null(facet_by)) {
        if (!(facet_by %in% names(colData(sce.o)))) 
            stop("facet_by variable does not exist in colData(sce.o).")
        facet_var <- colData(sce.o)[[facet_by]]
    } else {
        facet_var <- NULL
    }
    .plotEmbScatCyclic(emb.m = emb.m, color.value = color.value, color_by = "CCTime", facet_var = facet_var, ...)
    
})



#' 
#' @title Get the cyclic legend
#' 
#' @description BALABALA
#' 
#' @param hue.colors BAAA
#' @param hue.n BAAA
#' @param alpha AAAA
#' @param y.inner BAAA
#' @param y.outer BAAA
#' @param y.text BAAA
#' @param ymax BAAA
#' @param text.size BAAA
#' 
#' @return A ggplot object
#' 
#' @details BALA
#' 
#' @author Shijie C. Zheng
#' 
#' @name cyclic_legend
#' @rdname cyclic_legend
#' 
#' @examples
#' cyclic_legend()
NULL


#' 
#' @import ggplot2
#' @export
cyclic_legend <- function(hue.colors = c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9"), hue.n = 500, alpha = 0.6, y.inner = 1.5, y.outer = 3, y.text = 3.8, 
    ymax = 4.5, text.size = 3) {
    hues.df = data.frame(theta = seq(from = 0, to = 2 * pi, length.out = hue.n), colors = colorRampPalette(hue.colors)(hue.n))
    hue_text.df <- data.frame(theta = c(0, 0.5 * pi, pi, 1.5 * pi), label = c("0/2π", "0.5π", "π", "1.5π"), hjust = c(0.1, 0.5, 0.5, 0.5))
    legend.p <- ggplot(hues.df) + geom_rect(aes(ymin = y.inner, ymax = y.outer, xmin = theta - 0.001, xmax = theta + 0.001, color = colors, fill = colors), alpha = alpha, ) + coord_polar(theta = "x", 
        start = -pi/2, direction = -1, clip = "on") +  coord_polar(direction=-1, start=0) +
    scale_color_identity() + scale_fill_identity() + guides(fill = FALSE, color = FALSE) + theme_void() + ylim(c(0, ymax)) + geom_text(data = hue_text.df, aes(x = theta, y = y.text, label = label, 
        hjust = hjust), size = text.size) + theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))
    return(legend.p)
}


