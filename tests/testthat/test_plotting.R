test_that("create_long_data reshapes ECODA slots correctly", {
    data(example_data)
    d <- example_data$Zhang
    obj <- ecoda(data = d$cell_counts_lowresolution, metadata = d$metadata)

    # Test with label join
    long_df <- create_long_data(obj, data_slot = "clr", label_col = "Tissue")

    expect_true("clr_abundance" %in% colnames(long_df))
    expect_true("Tissue" %in% colnames(long_df))
    expect_equal(nrow(long_df), nrow(obj@clr) * ncol(obj@clr))
})

test_that("Plotting functions return ggplot or pheatmap objects", {
    data(example_data)
    d <- example_data$Zhang
    obj <- ecoda(data = d$cell_counts_lowresolution, metadata = d$metadata)

    p_pca <- plot_pca(obj, label_col = "Tissue", anosim_parallel = 1)
    expect_s3_class(p_pca, "ggplot")

    p_box <- plot_boxplot(obj, label_col = "Tissue")
    expect_s3_class(p_box, "ggplot")
})
