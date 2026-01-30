test_that("ECODA constructor works", {
    data(example_data)
    d <- example_data$Zhang

    obj <- ecoda(data = d$cell_counts_lowresolution, metadata = d$metadata)

    expect_s4_class(obj, "ECODA")
    expect_equal(nrow(obj@clr), nrow(d$metadata))
    expect_true(all(obj@freq_imp > 0))
})
