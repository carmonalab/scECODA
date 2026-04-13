test_that("Frequency calculation is accurate", {
    # Test that rows sum to exactly 100
    counts <- data.frame(A = c(10, 50), B = c(90, 50))
    freq <- calc_freq(counts)

    expect_equal(colSums(freq), c(A = 100, B = 100))
    expect_s3_class(freq, "data.frame")
})

test_that("Zero imputation handles counts and frequencies correctly", {
    # Test count imputation (multiplicative replacement)
    counts <- data.frame(A = c(10, 0), B = c(20, 10))
    counts_imp <- replace_zeros(counts, rep_method = "counts")

    expect_false(any(counts_imp == 0))
    expect_equal(counts_imp[2, 1], 0.5) # Default xmin_factor * counts_min

    # Test frequency imputation (fraction of smallest non-zero)
    freq <- data.frame(A = c(0.5, 0), B = c(0.5, 1.0))
    freq_imp <- replace_zeros(freq, rep_method = "frac_min")

    expect_false(any(freq_imp == 0))
    expect_equal(freq_imp[2, 1], 0.5 * (2 / 3))
})

test_that("CLR transformation is mathematically consistent", {
    # CLR should fail or warn if zeros are present (based on your main.r code)
    # clr() function in main.r assumes positive data.
    pos_data <- data.frame(A = c(10, 5), B = c(20, 25))
    clr_res <- calc_clr(pos_data)

    expect_equal(nrow(clr_res), 2)
    expect_equal(ncol(clr_res), 2)

    # Mathematical property: row sums of CLR-transformed data should be 0
    expect_equal(colSums(clr_res), c(A = 0, B = 0), tolerance = 1e-9)
})
