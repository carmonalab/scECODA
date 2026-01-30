test_that("calculate_pseudobulk aggregates and filters correctly", {
    # 3 genes, 6 cells
    counts <- matrix(1, nrow = 3, ncol = 6)
    rownames(counts) <- c("G1", "G2", "G3")
    colnames(counts) <- paste0("Cell", 1:6)
    samples <- c("S1", "S1", "S1", "S1", "S2", "S2")

    # Basic sum
    pb <- calculate_pseudobulk(counts, samples, min_cells = 1)
    expect_equal(pb["G1", "S1"], 4)
    expect_equal(ncol(pb), 2)

    # Test min_cells filtering
    # S2 only has 2 cells, so if we set min_cells to 4, it should drop S2
    expect_message(
        pb_filtered <- calculate_pseudobulk(counts, samples, min_cells = 4),
        "Filtered out 1 sample\\(s\\) with < 4 cells"
    )
    expect_equal(ncol(pb_filtered), 1)
    expect_equal(colnames(pb_filtered), "S1")
})
