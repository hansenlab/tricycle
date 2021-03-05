#' Plot cell cycle position kernel density stratified by a factor
#'
#' @description  The function will compute and plot cell cycle position kernel density.
#'
#'
#' @param theta.v The cell cycle position - a numeric vector with range between 0 to 2pi.
#' @param color_var.v A factor variable to stratify \code{theta.v}, such as samples or 'CCStage'. The length of it should equal to the length of \code{theta.v}.
#' @param color_name The name of the color_var.v to be used in the legend.
#' @param palette.v A string vector to give the color names. It should has the length of the number of levels of \code{color_var.v}. If not given, the 'Set1' palette will be used.
#'  (See \code{\link[RColorBrewer]{display.brewer.all}}) If the number of levels
#'  of \code{color_var.v} is greater than 8, only the top 8 levels of most cell
#'   will be shown. You can show them all by feeding enough colors in \code{palette.v}. Default: NULL
#' @param fig.title The title of the figure. Default: NULL
#' @param type It can be either of 'linear' or 'circular'. 'linear' corresponds to Cartesian coordinate system and 'circular' for polar coordinate system. Default: 'linear'
#' @param bw The smoothing bandwidth to be used. It is equal to the concentration parameter of the von Mises distribution. See \code{\link[circular]{density.circular}}. Default: 30
#' @param weighted Whether the density should be weighted on the percentage of each level of \code{color_var.v}. Default: FALSE
#' @param line.size The size of the line used by \code{\link[ggplot2]{geom_path}}. Default: 0.5
#' @param line.alpha The alpha value of the line used by \code{\link[ggplot2]{geom_path}}. Default: 0.5
#' @param ... Other arguments accepted by \code{\link[ggplot2]{geom_path}}.
#'
#' @return A ggplot object
#'
#' @details The function first estimates kernel density using the von Mises distribution. Then,
#' it plots out the density in the polar coordinate system or Cartesian coordinate system. Different colors
#' represents different levels of \code{color_var.v} and the dashed black line is the marginal distribution of all cells.
#'
#'
#'
#' @name plot_ccposition_den
#' @aliases plot_ccposition_den
#' @seealso
#' \code{\link{estimate_cycle_stage}}, for inferring 5 stages of cell cycle
#'
#' @author Shijie C. Zheng
#'
#' @examples
#' neurosphere_example <- estimate_cycle_position(neurosphere_example)
#' plot_ccposition_den(neurosphere_example$tricyclePosition, neurosphere_example$sample, "sample")
#'
#' neurosphere_example <- estimate_cycle_stage(neurosphere_example, gname.type = "ENSEMBL", species = "mouse")
#' plot_ccposition_den(neurosphere_example$tricyclePosition, neurosphere_example$CCStage, "CCStage")
NULL





#' @importFrom circular circular density.circular
#' @importFrom RColorBrewer brewer.pal
#' @importFrom forcats fct_relevel fct_explicit_na
#' @import ggplot2
#'
#' @export
plot_ccposition_den <- function(theta.v, color_var.v, color_name, palette.v = NULL, fig.title = NULL, type = c("linear", "circular"), bw = 30, weighted = FALSE, line.size = 0.5, line.alpha = 0.5,
    ...) {
    if ((min(theta.v) < 0) | (max(theta.v) > 2 * pi)) {
          stop("theta.v need to be between 0 - 2pi.")
      }
    if (length(theta.v) != length(color_var.v)) {
          stop("The length of theta.v and color_var.v should be the same.")
      }
    type <- match.arg(type)

    d <- density.circular(circular(theta.v), bw = bw)
    all.df <- data.frame(x = as.numeric(d$x), y = d$y)

    color_var.v <- factor(color_var.v)
    if (is.null(palette.v)) {
        if (nlevels(color_var.v) > 8) {
            warning("The number of levels of color_var.v is greater than 8.\n Only the top 8 levels of most cell will be shown.\nYou can show them all by feeding enough colors in palette.v.")
            color_var.v[!(color_var.v %in% names(sort(table(color_var.v), decreasing = TRUE))[seq_len(8)])] <- NA
        }
        if (nlevels(color_var.v) == 2) {
            palette.v <- brewer.pal(3, "Set1")[seq_len(2)]
        } else {
            palette.v <- brewer.pal(nlevels(color_var.v), "Set1")
        }
    } else {
        if (nlevels(color_var.v) > length(palette.v)) {
              warning("The number of colors in palette.v doesn't match the number of levels in color_var.v. ")
          }
    }

    if (weighted) {
        weights <- table(color_var.v) / sum(!is.na(color_var.v), na.rm = TRUE)
    } else {
        weights <- rep(1, length(color_var.v))
    }

    if (is.null(fig.title)) {
          fig.title <- paste0("(n=", length(theta.v), ")")
      }

    strati.df <- do.call(rbind, lapply(seq_len(nlevels(color_var.v)), function(idx) {
        d <- density.circular(circular(theta.v[which(color_var.v == levels(color_var.v)[idx])]), bw = bw)
        return(data.frame(x = as.numeric(d$x), y = d$y * weights[idx], color = levels(color_var.v)[idx]))
    }))
    strati.df$color <- factor(strati.df$color, levels = levels(color_var.v))

    max.v <- max(all.df$y, strati.df$y)
    if (any(is.na(strati.df$color))) {
        strati.df$color <- fct_relevel(fct_explicit_na(strati.df$color, na_level = "NA"), "NA", after = Inf)
        scale_color <- scale_color_manual(values = c(palette.v, "grey"), name = color_name, labels = paste0(levels(strati.df$color), "\nn=", table(strati.df$color)), limits = c(
            levels(strati.df$color),
            "NA"
        ))
    } else {
        scale_color <- scale_color_manual(values = palette.v, name = color_name, labels = paste0(levels(strati.df$color), "\nn=", table(strati.df$color)), limits = levels(strati.df$color))
    }

    if (type == "linear") {
        p <- ggplot(strati.df, aes_string(x = "x", y = "y")) +
            geom_path(aes_string(color = "color"), size = line.size, alpha = line.alpha, ...) +
            geom_path(data = all.df, size = line.size, alpha = line.alpha, color = "black", linetype = "dashed", ...) +
            scale_color +
            scale_x_continuous(limits = c(0, 2 * pi), breaks = c(0, pi / 2, pi, 3 * pi / 2, 2 * pi), labels = paste0(c(0, 0.5, 1, 1.5, 2), "\u03C0"), name = "\u03b8") +
            labs(title = fig.title, y = "Density") +
            .gg_theme
    } else {
        p <- ggplot(strati.df, aes(x = get("x"), y = get("y") + max.v)) +
            geom_path(aes_string(color = "color"), size = line.size, alpha = line.alpha, ...) +
            geom_path(data = all.df, size = line.size, alpha = line.alpha, color = "black", linetype = "dashed", ...) +
            scale_color +
            coord_polar(theta = "x", start = -pi / 2, direction = -1, clip = "on") +
            scale_x_continuous(limits = c(0, 2 * pi), breaks = c(0, pi / 2, pi, 3 * pi / 2), labels = paste0(c(0, 0.5, 1, 1.5), "\u03C0"), name = "") +
            scale_y_continuous(limits = c(0, max.v * 2), breaks = c(max.v, max.v * 1.5, max.v * 2), labels = c("0", format(max.v * c(0.5, 1), digits = 3)), name = "Density") +
            labs(title = fig.title) +
            .gg_theme
    }
    return(p)
}
