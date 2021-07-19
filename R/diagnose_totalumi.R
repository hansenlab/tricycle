#' Diagnostic function for UMI based datasets
#'
#' @description  The function will fit loess line for total UMIs numbers over
#' cell cycle position to diagnose non-fitting data, of which cells are not cycling.
#'
#'
#'
#' @param theta.v The cell cycle position - a numeric vector with range between 0 to 2pi.
#' @param totalumis The total UMIs number for each cell (without log2 transformation) - a numeric vector with the same length as \code{theta.v}.
#' @param span The parameter \eqn{\alpha} which controls the degree of smoothing. See \code{\link[stats]{loess}}. Default: 0.3
#' @param length.out The number of data points on the fitted lines to be output in the prediction data.frame. Default: 200
#' @param plot If \code{TRUE}, a \code{ggplot} scatter plot will be included in the output list. The figure will plot log2(totalumis) ~ theta.v 
#' with points and the fitted \code{\link[stats]{loess}} line. Default: TRUE
#' @param fig.title The title of the figure. Default: NULL
#' @param point.size  The size of the point in scatter plot used by \code{\link[scattermore]{geom_scattermore}}. Default: 2.1
#' @param point.alpha  The alpha value (transparency) of the point in scatter plot used by \code{\link[scattermore]{geom_scattermore}}. Default: 0.6
#' @param line.size  The size of the fitted line, used by \code{\link[ggplot2]{geom_path}}. Default: 0.8
#' @param line.alpha  The alpha value (transparency) of the fitted line, used by \code{\link[ggplot2]{geom_path}}. Default: 0.8
#' @param x_lab Title of x-axis. Default: "\eqn{\theta}"
#' @param y_lab Title of y-axis. Default: "log2(totalumis)"
#' @param ... Other arguments input to \code{\link[stats]{loess}}.
#'
#' @return A diagnostic message and a list with the following elements:
#' \itemize{
#'   \item fitted - The fitted vaues on the loess line. A vector of the length of y.
#'   \item residual - The residual values from the fitted loess line, i.e. y - y.fit. A vector of the length of y.
#'   \item pred.df - The prediction \code{data.frame} by uniformly sampling theta from 0 - 2*pi. Names of variables: \code{x} and \code{y}. The number of rows equals to \code{length.out}.
#'   \item loess.o - The fitted loess object.
#'   \item rsquared -  The coefficient of determination R2. Calculated as 1 - residual sum of squares / the total sum of squares.
#'   \item fig - When \code{plot} is \code{TRUE}, a \code{ggplot} scatter plot object will be returned with other items.
#' }
#'
#' @details This function fit a loess line between cell cycle position and 
#' log2 transformed total UMI number, as described in \code{\link{fit_periodic_loess}}.
#' If almost all cells are not cycling in a dataset, the estimated cell cycle positions
#' might be incorrect due to the shifted embedding center. Using the fact that the cell
#' should have highest total UMI number at the end of S phase and almost half of that 
#' highest total UMI number at M phase, we could detect those datasets which should 
#' be analysesd and intepreted carefully when using tricycle package. For such probelmatic
#' datasets, the defaul embedding center (0, 0) could lead to wrong inference. Thus, 
#' We don't rececommend using cell cycle position values if you get warnings from the 
#' \code{diagnose_totalumi} function.
#'
#'
#' @name diagnose_totalumi
#' @aliases diagnose_totalumi
#' @seealso
#' \code{\link{fit_periodic_loess}}.
#'
#' @author Shijie C. Zheng
#'
#' @examples
#' data(neurosphere_example, package = "tricycle")
#' neurosphere_example <- estimate_cycle_position(neurosphere_example)
#' diagnose.l <- diagnose_totalumi(neurosphere_example$tricyclePosition,
#'  neurosphere_example$TotalUMIs, plot = TRUE)
NULL

#' @export
#' 
diagnose_totalumi <- function(theta.v, totalumis, span = 0.3, length.out = 200, plot = TRUE,
															 fig.title = NULL, point.size = 2.1, point.alpha = 0.6,
															 line.size = 0.8, line.alpha = 0.8, x_lab = "\u03b8",
															 y_lab = "log2(totalumis)",  ...) {
	stopifnot("theta.v need to be between 0 - 2pi." = (min(theta.v) >= 0) & (max(theta.v) <= 2 * pi), 
						"The length of theta.v and totalumis should be the same." = length(theta.v) == length(totalumis))
	
	fit.l <- fit_periodic_loess(theta.v = theta.v, y = log2(totalumis + 1),
															span = span, length.out = length.out, plot = plot,
															fig.title = fig.title, point.size = point.size,
															point.alpha = point.alpha, line.size = line.size,
															line.alpha = line.alpha, color.vars = NULL,
															color.name = NULL, x_lab = x_lab, y_lab = y_lab,
															hue.colors = NULL, ... = ...)
	peak.v <- max(fit.l$pred.df$y[which((fit.l$pred.df$x > 0.75 * pi) & (fit.l$pred.df$x < 1.25 * pi))])
	valley.v <- min(fit.l$pred.df$y[which((fit.l$pred.df$x > 1.35 * pi) & (fit.l$pred.df$x < 1.85 * pi))])
	diff.v <- peak.v - valley.v
	if (diff.v < 0.4) {
		warning("The difference of predicted log2(total UMIs) between the end of S phase and M phase is less than 0.4.
						It seems that the predicted cell cycle position might not be correct due to the shift of embedding center. 
						This could be caused by that most cells in the dataset are not cycling. 
						We don't rececommend using cell cycle position values in this case.")
	} else {
		message(paste0("Message: The difference of predicted log2(total UMIs) between the end of S phase and M phase is ", sprintf("%.2f", diff.v), "(within tolerance)."))
	}
	
	if (mean((theta.v > 0.5 * pi )& (theta.v < 1.8 * pi)) > 0.5) {
		warning("The percentage of cell in positions between 0.5\u03C0 to 1.8\u03C0 is greater than 50%.
						It seems that the predicted cell cycle position might not be correct due to the shift of embedding center. 
						This could be caused by that most cells in the dataset are not cycling. 
						We don't rececommend using cell cycle position values in this case.")
	}
	
	return(fit.l)
}
