test_that("Clustering and distance metrics return valid scores", {
    # Create two distinct clusters in 2D space
    # All scores should be high and near 1
    set.seed(42)
    mat <- rbind(
        matrix(rnorm(10, mean = 0), ncol = 2),
        matrix(rnorm(10, mean = 100), ncol = 2)
    )
    labels <- c(rep("Group1", 5), rep("Group2", 5))

    # ANOSIM
    ano <- calc_anosim(mat, labels, parallel = 1)
    expect_true(ano > 0.8)

    # ARI
    ari <- calc_ari(mat, labels)
    expect_type(ari, "double")
    expect_true(ari > 0.8)

    # Modularity
    mod <- calc_modularity(mat, labels)
    expect_true(mod > 0.8)

    # Silhouette
    sil <- calc_sil(mat, labels)
    expect_true(sil > 0.8)


    # Create two mixed clusters in 2D space
    # All scores should be low and near 0
    mat <- rbind(
        matrix(rnorm(10, mean = 0), ncol = 2),
        matrix(rnorm(10, mean = 0), ncol = 2)
    )

    # ANOSIM
    ano <- calc_anosim(mat, labels, parallel = 1)
    expect_true(ano < 0.2)

    # ARI
    ari <- calc_ari(mat, labels)
    expect_type(ari, "double")
    expect_true(ari < 0.2)

    # Modularity
    mod <- calc_modularity(mat, labels)
    expect_true(mod < 0.2)

    # Silhouette
    sil <- calc_sil(mat, labels)
    expect_true(sil < 0.2)
})
