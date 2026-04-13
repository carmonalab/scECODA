test_that("get_celltype_variances calculates metrics correctly", {
    # Create an ECODA object with known variance
    # Celltype 1 varies a lot, Celltypes 2, 3 and 4 vary less
    df <- data.frame(
        A = c(10, 20, 5, 400),
        B = c(11, 30, 10, 250),
        C = c(12, 40, 15, 130)
    )
    rownames(df) <- c("ct1", "ct2", "ct3", "ct4")

    se <- ecoda(data = df)

    var_df <- get_celltype_variances(se)

    expect_equal(var_df$celltype[1], "ct4")
    expect_true(var_df$Variance[1] > var_df$Variance[2])
    expect_equal(var_df$variance_exp[nrow(var_df)], 1.0) # Cumulative should reach 1
})

test_that("get_hvcs selection logic is robust", {
    df_var <- data.frame(
        celltype = paste0("C", 1:10),
        Variance = 10:1,
        variance_exp = cumsum(10:1) / sum(10:1)
    )

    # Test top_n_hvcs (Integer)
    hvcs_n <- get_hvcs(df_var, top_n_hvcs = 3)
    expect_length(hvcs_n, 3)
    expect_equal(hvcs_n, c("C1", "C2", "C3"))

    # Test variance_explained threshold
    # C1 is ~0.18, C2 is ~0.34...
    hvcs_var <- get_hvcs(df_var, variance_explained = 0.5)
    expect_true(length(hvcs_var) >= 2)

    # Test minimum return (always at least 2)
    hvcs_min <- get_hvcs(df_var[1:5, ], top_n_hvcs = 1)
    expect_length(hvcs_min, 2)
})
