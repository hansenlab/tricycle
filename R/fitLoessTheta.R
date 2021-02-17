#' Fit loess line with cell cycle time as predictor
#' 
#' @description  The function will fit a loess line using cell cycle time and other variables, such as expression levels of a gene
#'  or log-transformed totalUMIs numbers. The circular nature of cell cycle time is taken into account by making 3 copies inside the function.
#'  
#' @param theta.v The cell cycle time - a numeric vector with range between 0 to 2pi. 
#' @param y The response variable - a numeric vector with the same length as \code{theta.v}.
#' @param span The parameter α which controls the degree of smoothing. See \code{\link[stats]{loess}}. Default: 0.3
#' @param length.out The number of data points on the fitted lines to be output in the prediction data.frame. Default: 200
#' 
#' @return A list with the following elements:
#' \itemize{
#'   \item fitted - The fitted vaues on the loess line. A vector of the length of y.
#'   \item residual - The residual values from the fitted loess line, i.e. y - y.fit. A vector of the length of y.
#'   \item pred.df - The prediction \code{data.frame} by uniformly sampling theta from 0 - 2*pi. Names of variables: \code{x} and \code{y}. The number of rows equals to \code{length.out}.
#'   \item loess.o The fitted loess object.
#' }
#'
#' @details This function fit a normal loess line, but take the circularity of cell cycle time into account by making \code{theta.v} 3 periods
#' (\code{c(theta.v - 2 * pi, theta.v, theta.v + 2 * pi)}) and repeating y 3 times. Only the fitted values corresponding to original \code{theta.v} 
#' will be returned. User can use \code{pred.df} to visualize the loess line.
#' 
#' 
#' @name fitLoessTheta
#' @seealso
#' \code{\link{inferCCTime}}, for inferring cell cycle time.
#'
#' @author Shijie C. Zheng
#'
#' @examples
#' example_sce <- inferCCTime(example_sce)
#' top2a.idx <- which(rowData(example_sce)$Gene == "Top2a")
#' fit.l <- fitLoessTheta(example_sce$CCTime, assay(example_sce, "logcounts")[top2a.idx, ])
#' plot(example_sce$CCTime, assay(example_sce, "logcounts")[top2a.idx, ])
#' lines(fit.l$pred.df)
NULL


#' @export
fitLoessTheta <- function(theta.v, y, span = 0.3, length.out = 200) {
	if ((min(theta.v) < 0) | (max(theta.v) > 2 * pi)) stop("theta.v need to be between 0 - 2pi.")
	if (length(theta.v) != length(y)) stop("The length of theta.v and y should be the same.")
	x <- c(theta.v - 2 * pi, theta.v, theta.v + 2 * pi)
	y <- rep(y, 3)
	loess.o <- loess(y ~ x, span = span)
	fitted.v <- loess.o$fitted[(length(theta.v) + 1):(length(theta.v) * 2)]
	residual.v <- loess.o$residuals[(length(theta.v) + 1):(length(theta.v) * 2)]
	pred.x <- seq(0, 2 * pi, length.out = length.out)
	pred.y <- predict(loess.o, newdata = data.frame(x = pred.x ))
	return(list(fitted = fitted.v, residual = residual.v, pred.df = data.frame(x = pred.x, y = pred.y), loess.o = loess.o))
}



