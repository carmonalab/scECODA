test_that("ECODA SingleCellExperiment object constructor works", {
    data(example_data)
    d <- example_data$Zhang

    se <- ecoda(data = d$cell_counts_lowresolution, metadata = d$metadata)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(ncol(assay(se, "clr")), nrow(colData(se)))
    expect_true(all(assay(se, "freq_imp") > 0))
})
