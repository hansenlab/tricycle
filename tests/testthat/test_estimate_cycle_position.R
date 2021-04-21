context("estimate_cycle_position")


test_that("estimate_cycle_position", {
	data(neurosphere_example)
	
	expect_s4_class(estimate_cycle_position(neurosphere_example), "SingleCellExperiment")
	expect_type(estimate_cycle_position(assay(neurosphere_example, "logcounts")), "double")
	
	n <- sample(100:400, 1)
	expect_equal(ncol(estimate_cycle_position(neurosphere_example[, seq_len(n)])), n)

	expect_error(project_cycle_space(neurosphere_example, altexp = "aa"))
	expect_error(project_cycle_space(neurosphere_example, dimred = "aa"))

	
})


