#' Plot embedding with cyclic cell cycle position
#'
#' Generate scat plot of embedding with cyclic cell cycle position or other cyclic variables
#'
#'
#' @param sce.o A \linkS4class{SingleCellExperiment} contains the embbing to be plotted against.
#' @param color_by The name of variable in \code{colData(sce.o)} to be used to show colors. Default: "CCPosition"
#' @param facet_by he name of variable in \code{colData(sce.o)} to be used to facet scatter plots. If NULL, no faceted panles will be returned. Default: NULL
#' @param dimred The name or index of reducedDims in  \linkS4class{SingleCellExperiment} (\code{\link[SingleCellExperiment]{reducedDims}}). Default: 1
#' @param dim The indices of \code{dimred} to be plotted. At the moment, it has to be two integers.   Default: 1:2
#' @param fig.title The title of the figure. Default: NULL
#' @param point.size  The size of the point in scatter plot used by \code{\link[scattermore]{geom_scattermore}}. Default: 2.1
#' @param point.alpha  The alpha value (transparency) of the point in scatter plot used by \code{\link[scattermore]{geom_scattermore}}. Default: 0.6
#' @param x_lab Title of x-axis. If not given, the colname of \code{dimred} will be used. Default: NULL
#' @param y_lab Title of y-axis. If not given, the colname of \code{dimred} will be used. Default: NULL
#' @param hue.colors The string vector gives the cyclic colors. The first color should look very similar to the last one.
#' Default: c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9")
#' @param hue.n The number of breaks of color scheme. Default: 500
#' @param plot.legend Whether the legend should be plotted with the scatter plot. We recommend not to use this legend but use the cyclic legend produced by \code{\link{circle_scale_legend}} instead. Default: FALSE
#'
#' @details
#' This function help users plot embedding scater plot colored by cyclic variables, such as cell cycle position, which is bound between 0 - 2pi.
#' It will take a \linkS4class{SingleCellExperiment} object as input, and plot out its \code{dimred} such as PCA, UMAP, and etc with a cyclic color scheme.
#'
#' @return
#' A ggplot object or a list of ggplot objects.
#' If \code{facet_by} is not assigned, a single ggplot plot of the scatter plot will be return,
#' Otherwise, apart from the first scatter plot showing all cells together, other faceted scatter plots will also be given in the list.
#'
#' @name plot_emb_circle_scale
#'
#' @author Shijie C. Zheng
#'
#' @examples
#' neurosphere_example <- estimate_cycle_position(neurosphere_example)
#' plot_emb_circle_scale(neurosphere_example, point.size = 3.1, point.alpha = 0.8)
NULL


#' @importFrom scattermore geom_scattermore
#' @importFrom dplyr filter .data
#' @importFrom magrittr `%>%`
.plot_emb_circle_scale <- function(emb.m, color.value, color_by, facet_var = NULL, fig.title = NULL, point.size = 2.1, point.alpha = 0.6, x_lab = NULL, y_lab = NULL, hue.colors = c(
        "#2E22EA",
        "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9"
    ), hue.n = 500, plot.legend = FALSE) {
    if (is.null(fig.title)) {
          fig.title <- paste0("(n=", nrow(emb.m), ")")
      }

    if (is.null(x_lab)) {
          x_lab <- colnames(emb.m)[1]
      }
    if (is.null(y_lab)) {
          y_lab <- colnames(emb.m)[2]
      }

    x_lim <- range(emb.m[, 1]) + c(diff(range(emb.m[, 1])) * c(-0.05, 0.05))
    y_lim <- range(emb.m[, 2]) + c(diff(range(emb.m[, 2])) * c(-0.05, 0.05))

    tmp.df <- data.frame(x = emb.m[, 1], y = emb.m[, 2], color = color.value)

    scat.p <- ggplot(tmp.df, aes_string(x = "x", y = "y", color = "color")) +
        geom_scattermore(pointsize = point.size, alpha = point.alpha) +
        scale_color_gradientn(name = color_by, limits = range(0, 2 * pi), breaks = seq(from = 0, to = 2 * pi, length.out = hue.n), colors = hue.colors, guide = FALSE) +
        labs(y = y_lab, x = x_lab, title = fig.title) +
        xlim(x_lim) +
        ylim(y_lim) +
        .gg_theme

    if (!is.null(facet_var)) {
        tmp.df$facet <- facet_var
        facet_labels <- levels(factor(tmp.df$facet))
        lp <- lapply(seq_len(nlevels(factor(tmp.df$facet))), function(idx) {
            p <- ggplot(tmp.df, aes_string(x = "x", y = "y", color = "color")) +
                geom_scattermore(data = tmp.df %>% dplyr::filter(.data[["facet"]] != levels(factor(tmp.df$facet))[idx]), pointsize = point.size, alpha = 0.4, color = "gray90", show.legend = FALSE) +
                geom_scattermore(data = tmp.df %>% dplyr::filter(.data[["facet"]] == levels(factor(tmp.df$facet))[idx]), pointsize = point.size, alpha = point.alpha) +
                scale_color_gradientn(name = color_by, limits = range(0, 2 * pi), breaks = seq(from = 0, to = 2 * pi, length.out = hue.n), colors = hue.colors, guide = FALSE) +
                labs(y = y_lab, x = x_lab, title = facet_labels[idx]) +
                xlim(x_lim) +
                ylim(y_lim) +
                .gg_theme
            return(p)
        })
        return(c(list(scat.p), lp))
    }
    return(scat.p)
}


#' @export
#' @rdname plot_emb_circle_scale
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData
setMethod("plot_emb_circle_scale", "SingleCellExperiment", function(sce.o, color_by = "CCPosition", facet_by = NULL, dimred = 1, dim = seq_len(2),
    fig.title = NULL, point.size = 2.1, point.alpha = 0.6, x_lab = NULL, y_lab = NULL,
    hue.colors = c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9"),
    hue.n = 500, plot.legend = FALSE) {
    stopifnot("The function can only plot 2 dims at this time. Change dim argument." = length(dim) == 2)
  
    emb.m <- reducedDim(sce.o, dimred)[, dim]
    color.value <- colData(sce.o)[[color_by]]
    if (!is.null(facet_by)) {
        if (!(facet_by %in% names(colData(sce.o)))) {
              stop("facet_by variable does not exist in colData(sce.o).")
          }
        facet_var <- colData(sce.o)[[facet_by]]
    } else {
        facet_var <- NULL
    }
    .plot_emb_circle_scale(
        emb.m = emb.m, color.value = color.value, color_by = color_by,
        facet_var = facet_var, fig.title = fig.title, point.size = point.size,
        point.alpha = point.alpha, x_lab = x_lab, y_lab = y_lab,
        hue.colors = hue.colors,
        hue.n = hue.n, plot.legend = plot.legend
    )
})



#'
#' @title Get the cyclic legend
#'
#' @description This function is a helper function to create the cyclic ggplot color legend.
#'
#' @usage circle_scale_legend(hue.colors = c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9"),
#'  hue.n = 500, alpha = 0.6, y.inner = 1.5, y.outer = 3, y.text = 3.8, ymax = 4.5, text.size = 3)
#'
#' @param hue.colors The string vector gives the cyclic colors. The first color should look very similar to the last one.
#' Default: c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9")
#' @param hue.n The number of breaks of color scheme. Default: 500
#' @param alpha The alpha value (transparency). Default: 0.6
#' @param y.inner The radius of inner circle of the donut. Default: 1.5
#' @param y.outer The radius of outer circle of the donut. Default: 3
#' @param y.text The radius of text position. Default: 3.8
#' @param ymax The value control the border of the legend. Default: 4.5
#' @param text.size The size of the text Default: 3
#'
#' @return A ggplot object
#'
#' @details The function will make a donut shape to serve as the cyclic color legend. The arguments should match the argument used in \code{\link{plot_emb_circle_scale}}.
#'
#' @author Shijie C. Zheng
#'
#' @name circle_scale_legend
#' @aliases circle_scale_legend
#' @rdname circle_scale_legend
#'
#' @examples
#' circle_scale_legend()
NULL


#' @importFrom grDevices colorRampPalette
#' @import ggplot2
#' @export
circle_scale_legend <- function(hue.colors = c("#2E22EA", "#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9"), hue.n = 500, alpha = 0.6, y.inner = 1.5, y.outer = 3,
    y.text = 3.8, ymax = 4.5, text.size = 3) {
    hues.df <- data.frame(theta = seq(from = 0, to = 2 * pi, length.out = hue.n), colors = colorRampPalette(hue.colors)(hue.n))
    hue_text.df <- data.frame(theta = c(0, 0.5 * pi, pi, 1.5 * pi), label = c("0/2\u03C0", "0.5\u03C0", "\u03C0", "1.5\u03C0"), hjust = c(0.1, 0.5, 0.5, 0.5))
    legend.p <- ggplot(hues.df) +
        geom_rect(aes(ymin = get("y.inner"), ymax = get("y.outer"), xmin = get("theta") - 0.001, xmax = get("theta") + 0.001, color = get("colors"), fill = get("colors")), alpha = alpha) +
        coord_polar(theta = "x", start = -pi / 2, direction = -1, clip = "on") +
        scale_color_identity() +
        scale_fill_identity() +
        guides(fill = FALSE, color = FALSE) +
        theme_void() +
        ylim(c(0, ymax)) +
        geom_text(data = hue_text.df, aes_string(x = "theta", y = "y.text", label = "label", hjust = "hjust"), size = text.size) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))
    return(legend.p)
}
