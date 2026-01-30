test_that("Frequency calculation is accurate", {
    # Test that rows sum to exactly 100
    counts <- data.frame(A = c(10, 50), B = c(90, 50))
    freq <- calc_freq(counts)

    expect_equal(rowSums(freq), c(100, 100))
    expect_s3_class(freq, "data.frame")
})

test_that("Zero imputation handles counts and frequencies correctly", {
    # Test count imputation (multiplicative replacement)
    counts_df <- data.frame(A = c(10, 0), B = c(20, 10))
    imputed_counts <- impute_zeros(counts_df, is_freq = FALSE)

    expect_false(any(imputed_counts == 0))
    expect_equal(imputed_counts[2, 1], 2 / 3) # Default xmin_factor * counts_min

    # Test frequency imputation (fraction of smallest non-zero)
    freq_df <- data.frame(A = c(0.5, 0), B = c(0.5, 1.0))
    imputed_freq <- impute_zeros(freq_df, is_freq = TRUE)

    expect_false(any(imputed_freq == 0))
    expect_equal(imputed_freq[2, 1], 0.5 * (2 / 3))
})

test_that("CLR transformation is mathematically consistent", {
    # CLR should fail or warn if zeros are present (based on your main.r code)
    # clr() function in main.r assumes positive data.
    pos_data <- data.frame(A = c(10, 5), B = c(20, 25))
    clr_res <- clr(pos_data)

    expect_equal(nrow(clr_res), 2)
    expect_equal(ncol(clr_res), 2)

    # Mathematical property: row sums of CLR-transformed data should be 0
    expect_equal(rowSums(clr_res), c(0, 0), tolerance = 1e-9)
})
