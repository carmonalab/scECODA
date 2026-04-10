test_that("create_long_data reshapes ECODA slots correctly", {
    data(example_data)
    d <- example_data$Zhang
    se <- ecoda(data = d$cell_counts_lowresolution, metadata = d$metadata)

    # Test with label join
    long_df <- create_long_data(se, assay = "clr", label_col = "Tissue")

    expect_true("value" %in% colnames(long_df))
    expect_true("Tissue" %in% colnames(long_df))
    expect_equal(nrow(long_df), nrow(assay(se, "clr")) * ncol(assay(se, "clr")))
})

test_that("Plotting functions return ggplot or pheatmap objects", {
    data(example_data)
    d <- example_data$Zhang
    se <- ecoda(data = d$cell_counts_lowresolution, metadata = d$metadata)

    p_pca <- plot_pca(se, label_col = "Tissue", anosim_parallel = 1)
    expect_s3_class(p_pca, "ggplot")

    p_box <- plot_boxplot(se, label_col = "Tissue")
    expect_s3_class(p_box, "ggplot")
})
