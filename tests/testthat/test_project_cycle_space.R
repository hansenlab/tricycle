context("project_cycle_space")


test_that("project_cycle_space", {
	data(neurosphere_example)
	
	expect_s4_class(project_cycle_space(neurosphere_example), "SingleCellExperiment")
	expect_type(project_cycle_space(assay(neurosphere_example, "logcounts")), "double")
	
	n <- sample(100:400, 1)
	expect_equal(ncol(project_cycle_space(neurosphere_example[, seq_len(n)])), n)
	
	expect_error(project_cycle_space(neurosphere_example, gname.type = "SYMBOL"))
	expect_message(project_cycle_space(neurosphere_example, AnnotationDb = "aa"))
	expect_error(project_cycle_space(neurosphere_example, species = "human", AnnotationDb = "aa"), "The input AnnotationDB is not an AnnotationDb object.")
	
	
})
