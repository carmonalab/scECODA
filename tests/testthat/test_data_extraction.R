test_that("get_celltype_counts handles NA and dimensions correctly", {
    cell_df <- data.frame(
        Sample = c(rep("S1", 3), rep("S2", 2)),
        Type = c("B", "T", NA, "T", "T")
    )

    counts <- get_celltype_counts(cell_df, "Sample", "Type")

    # Check dimensions: 2 samples, 3 categories (B, T, and NA)
    expect_equal(nrow(counts), 2)
    expect_equal(ncol(counts), 3)
    expect_true("NA" %in% colnames(counts))
    expect_equal(counts["S1", "NA"], 1)
})

test_that("get_sample_metadata filters non-constant columns", {
    cell_df <- data.frame(
        Sample_ID = c(rep("S1", 3), rep("S2", 2)),
        Age = c(30, 30, 30, 45, 45), # Constant per sample
        Cell_ID = paste0("C", 1:5), # Varies per sample
        Condition = c("A", "A", "A", "B", "B") # Constant per sample
    )

    sample_meta <- get_sample_metadata(cell_df, "Sample_ID")

    # Should keep Age and Condition, but drop Cell_ID
    expect_true(all(c("Age", "Condition") %in% colnames(sample_meta)))
    expect_false("Cell_ID" %in% colnames(sample_meta))
    expect_equal(nrow(sample_meta), 2)
    expect_equal(rownames(sample_meta), c("S1", "S2"))
})
