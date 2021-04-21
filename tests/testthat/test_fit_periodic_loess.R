context("fit_periodic_loess")



test_that("fit_periodic_loess is working", {
	data(neurosphere_example)
	neurosphere_example <- estimate_cycle_position(neurosphere_example)
	theta.v <- neurosphere_example$tricyclePosition
	y <- assay(neurosphere_example, "logcounts")[1, ]
	
	expect_error(fit_periodic_loess(theta.v, y[-1]))
	expect_equal(length(fit_periodic_loess(theta.v, y)$fitted), length(theta.v))
	expect_equal(length(fit_periodic_loess(theta.v, y)$residual), length(theta.v))
	expect_lte(abs(fit_periodic_loess(theta.v, y)$rsquared), 1)
	n <- sample(100:500, 1)
	expect_equal(nrow(fit_periodic_loess(theta.v, y, length.out = n)$pred.df), n)
	
})
