#' Fit periodic loess line with circular predictor
#'
#' @description  The function will fit a loess line using cell cycle position and other variables, such as expression levels of a gene
#'  or log-transformed totalUMIs numbers. The circular nature of cell cycle position is taken into account by making 3 copies inside the function.
#'  For convenience, the function will also return a scatter plot with fitted line if needed.
#'
#'
#' @param theta.v The cell cycle position - a numeric vector with range between 0 to 2pi.
#' @param y The response variable - a numeric vector with the same length as \code{theta.v}.
#' @param span The parameter \eqn{\alpha} which controls the degree of smoothing. See \code{\link[stats]{loess}}. Default: 0.3
#' @param length.out The number of data points on the fitted lines to be output in the prediction data.frame. Default: 200
#' @param plot If \code{TRUE}, a \code{ggplot} scatter plot will be included in the output list. The figure will plot y ~ theta.v 
#' with points and the fitted \code{\link[stats]{loess}} line. Default: FALSE
#' @param fig.title The title of the figure. Default: NULL
#' @param point.size  The size of the point in scatter plot used by \code{\link[scattermore]{geom_scattermore}}. Default: 2.1
#' @param point.alpha  The alpha value (transparency) of the point in scatter plot used by \code{\link[scattermore]{geom_scattermore}}. Default: 0.6
#' @param line.size  The size of the fitted line, used by \code{\link[ggplot2]{geom_path}}. Default: 0.8
#' @param line.alpha  The alpha value (transparency) of the fitted line, used by \code{\link[ggplot2]{geom_path}}. Default: 0.8
#' @param color.vars Optional. A vector of categorical variable of the same length of \code{theta.v}, and it will be used to color points in figure. Default: NULL
#' @param color.name The name of the color variables. Used as the name for legend. Default: NULL
#' @param x_lab Title of x-axis. If not given, the colname of \code{dimred} will be used. Default: "\eqn{\theta}"
#' @param y_lab Title of y-axis. If not given, the colname of \code{dimred} will be used. Default: "y"
#' @param hue.colors The string vector gives custom colors. If not given, the default \code{\link[ggplot2]{scale_color_discrete}} will be used. Default: NULL
#' @param ... Other arguments input to \code{\link[stats]{loess}}.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item fitted - The fitted vaues on the loess line. A vector of the length of y.
#'   \item residual - The residual values from the fitted loess line, i.e. y - y.fit. A vector of the length of y.
#'   \item pred.df - The prediction \code{data.frame} by uniformly sampling theta from 0 - 2*pi. Names of variables: \code{x} and \code{y}. The number of rows equals to \code{length.out}.
#'   \item loess.o - The fitted loess object.
#'   \item fig - When \code{plot} is \code{TRUE}, a \code{ggplot} scatter plot object will be returned with other items.
#' }
#'
#' @details This function fit a normal loess line, but take the circularity of cell cycle position into account by making \code{theta.v} 3 periods
#' (\code{c(theta.v - 2 * pi, theta.v, theta.v + 2 * pi)}) and repeating y 3 times. Only the fitted values corresponding to original \code{theta.v}
#' will be returned. For convenience, the function will also return a scatter plot with fitted line if needed. 
#' Or user can use \code{pred.df} to visualize the loess line themselves.
#'
#'
#' @name fit_periodic_loess
#' @aliases fit_periodic_loess
#' @seealso
#' \code{\link{estimate_cycle_position}}, for inferring cell cycle position.
#'
#' @author Shijie C. Zheng
#'
#' @examples
#' neurosphere_example <- estimate_cycle_position(neurosphere_example)
#' top2a.idx <- which(rowData(neurosphere_example)$Gene == "Top2a")
#' fit.l <- fit_periodic_loess(neurosphere_example$tricyclePosition, assay(neurosphere_example, "logcounts")[top2a.idx, ], plot = TRUE)
#' fit.l$fig
NULL

#' @importFrom stats loess predict
#' @export
fit_periodic_loess <- function(theta.v, y, span = 0.3, length.out = 200, plot = FALSE,
                          fig.title = NULL, point.size = 2.1, point.alpha = 0.6,
                          line.size = 0.8, line.alpha = 0.8,
                          color.vars = NULL, color.name = NULL, x_lab = "\u03b8",
                          y_lab = "y", hue.colors = NULL, ...) {
    stopifnot("theta.v need to be between 0 - 2pi." = (min(theta.v) >= 0) & (max(theta.v) <= 2 * pi), 
              "The length of theta.v and y should be the same." = length(theta.v) == length(y))
    x <- c(theta.v - 2 * pi, theta.v, theta.v + 2 * pi)
    y <- rep(y, 3)
    loess.o <- loess(y ~ x, span = span, ...)
    fitted.v <- loess.o$fitted[(length(theta.v) + 1):(length(theta.v) * 2)]
    residual.v <- loess.o$residuals[(length(theta.v) + 1):(length(theta.v) * 2)]
    pred.x <- seq(0, 2 * pi, length.out = length.out)
    pred.y <- predict(loess.o, newdata = data.frame(x = pred.x))
    pred.df <-  data.frame(x = pred.x, y = pred.y)
    if (plot) {
      if (is.null(fig.title)) fig.title <- paste0("(n=", length(theta.v), ")")
      if (!is.null(color.vars)) {
        color.vars <- factor(color.vars)
        stopifnot("Length of theta.v does not match color.vars" = length(color.vars) == length(theta.v))
        if (!is.null(hue.colors)) {
          stopifnot("Number of colors does not match nlevels of color.vars" = nlevels(factor(color.vars) == length(hue.colors)))
          color_scale <- scale_color_manual(values = hue.colors, name = color.name, limits = levels(color.vars))
        } else {
          color_scale <- scale_color_discrete(name = color.name, limits = levels(color.vars))
        }
        tmp.df <- data.frame(theta = theta.v, y = y, color = color.vars)
        p_aes <- aes_string(x = "theta", y = "y", color = "color")
      } else {
        tmp.df <- data.frame(theta = theta.v, y = y)
        p_aes <- aes_string(x = "theta", y = "y")
        color_scale <- NULL
      }
      fig <- ggplot(data = tmp.df) +
        geom_scattermore(mapping = p_aes, pointsize = point.size) +
        geom_path(data = pred.df, aes_string(x = "x", y = "y"), size = line.size, alpha = line.alpha) +
        color_scale +
        labs(x = x_lab, y = y_lab, title = fig.title) +
        .gg_theme
      return(list(fitted = fitted.v, residual = residual.v, pred.df = pred.df, loess.o = loess.o, fig = fig))
    }
    
    return(list(fitted = fitted.v, residual = residual.v, pred.df = pred.df, loess.o = loess.o))
}
