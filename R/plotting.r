# PCA ---------------------------

#' @title Plot Principal Component Analysis and Calculate Clustering Scores
#'
#' @description Performs Principal Component Analysis (PCA) on a selected data
#'   matrix from the \code{ECODA} object (default: CLR-transformed abundances,
#'   \code{clr}) and visualizes the results. It can also calculate and display
#'   several metrics to evaluate the separation of groups defined by
#'   \code{label_col}.
#'
#' @details The clustering metrics (ARI, Modularity, Silhouette, ANOSIM) assess
#'   how well the sample groupings (\code{labels}) align with the underlying
#'   data structure in the feature space defined by the PCA.
#'
#' @param ecoda_object An \link[=ECODA-class]{ECODA} object.
#' @param slot Character string (default: \code{"clr"}). The name of the data
#'   matrix slot in the \code{ECODA} object to use for PCA. Must be one of:
#'   \code{"clr"} (CLR-transformed abundances, default), \code{"clr_hvc"}
#'   (CLR-transformed abundances of only the most highly variable cell types
#'   (HVCs)), \code{"pb"} (pseudobulk gene expression), \code{"counts"} (raw
#'   counts), \code{"counts_imp"} (imputed counts), \code{"freq"} (relative
#'   frequencies), \code{"freq_imp"} (imputed frequencies), or
#'   \code{"asin_sqrt"} (arcsin-square root transformed data).
#' @param label_col Character string (optional, default: \code{NULL}). The name
#'   of a column in \code{slot(ecoda_object, "metadata")} used to color and
#'   group samples in the plot, and for calculating clustering scores.
#' @param scale. Logical (default: \code{FALSE}). A value indicating whether the
#'   variables should be scaled to have unit variance before the PCA.
#' @param knn_k Integer (optional, default: \code{NULL}). The number of nearest
#'   neighbors (\code{k}) to use for the Shared Nearest Neighbor (SNN) graph
#'   construction, required for Modularity score calculation. If \code{NULL}, it
#'   defaults to \code{max(3, round(sqrt(N)))}, where \code{N} is the number of
#'   samples.
#' @param title Character string (optional, default: \code{NULL}). The main
#'   title for the plot. If clustering scores are calculated, they are appended
#'   to this title.
#' @param title_show_n_features Logical (optional, default: \code{TRUE}) Show
#'   the number of features (cell types or genes) used.
#' @param legend_title Character string (default: \code{"Group"}). The title for
#'   the color legend in the plot when a grouping column (\code{label_col}) is
#'   provided.
#' @param show_label_samples Logical (default: \code{FALSE}). If \code{TRUE},
#'   sample names are displayed next to the points in the plot. This
#'   automatically adds \code{"text"} to the \code{geom} parameter if it is not
#'   already present.
#' @param score_digits Integer (default: \code{3}). The number of decimal places
#'   to round the clustering and ANOSIM scores appended to the plot title.
#' @param cluster_score Logical (default: \code{TRUE}). If \code{TRUE},
#'   calculates the Adjusted Rand Index (ARI) using \code{\link{calc_ari}}.
#' @param mod_score Logical (default: \code{TRUE}). If \code{TRUE}, calculates
#'   the adjusted Modularity score using \code{\link{calc_modularity}}.
#' @param sil_score Logical (default: \code{FALSE}). If \code{TRUE}, calculates
#'   the average Silhouette width using \code{\link{calc_sil}}.
#' @param anosim_score Logical (default: \code{TRUE}). If \code{TRUE},
#'   calculates the ANOSIM statistic (R) using \code{vegan::anosim}.
#' @param anosim_distance Character string (default: \code{"euclidian"}). The
#'   distance method used for the ANOSIM calculation (e.g., "euclidean",
#'   "manhattan").
#' @param anosim_permutations Integer (default: \code{99}). The number of
#'   permutations to use when calculating the ANOSIM statistic.
#' @param anosim_parallel Integer (default: \code{detectCores()}). The number of
#'   parallel processes/cores to use for the ANOSIM calculation.
#' @param ari_nclusts Integer (optional, default: \code{NULL}). The target
#'   number of clusters (\code{k}) to use for \code{hclust} and \code{pam}. If
#'   \code{NULL}, it defaults to the number of unique levels in \code{labels}.
#' @param pointsize Numeric (default: \code{3}). Size of the points in the plot.
#' @param labelsize Numeric (default: \code{4}). Size of the variable labels in
#'   the plot.
#' @param coord_equal Logical (default: \code{TRUE}). If \code{TRUE}, forces the
#'   aspect ratio of the plot to be equal.
#' @param axes Numeric vector (default: \code{c(1, 2)}). The principal
#'   components to plot (e.g., \code{c(1, 2)} for PC1 vs PC2).
#' @param invisible Character vector (default: \code{c("var", "quali")}).
#'   Elements to hide. Can include "var" (variables/cell types),
#'   "ind" (samples), or "quali" (group centroids).
#' @param geom Character string or vector (default: \code{"point"}). The
#'   geometry to be used for the plot. Allowed values are combinations of:
#'             \itemize{
#'               \item \code{"point"} to show points for individuals (samples)
#'               \item \code{"text"} to show labels for individuals (samples)
#'               \item \code{"arrow"} to show vectors for variables (features)
#'             }
#'   The default \code{"point"} plots points for both individuals and variables.
#'   Use \code{c("point", "text")} to show both points and labels for samples.
#' @param n_hv_feat_show Integer (default: \code{Inf}). Number of most highly
#'   variable features to show based on their contribution to the selected axes.
#' @param repel Logical (default: \code{TRUE}). Whether to use \code{ggrepel} to
#'   prevent label overlap for variable names.
#'
#' @return A \code{ggplot} object via \code{factoextra} visualizing the PCA
#'   results.
#'
#' @importFrom stats prcomp
#' @importFrom factoextra fviz_pca
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggtitle scale_shape_manual coord_equal
#'   scale_color_discrete
#' @importFrom parallel detectCores
#'
#' @export plot_pca
#'
#' @examples
#' data(example_data)
#' ecoda_object <- ecoda(
#'     data = example_data$Zhang$cell_counts_lowresolution,
#'     metadata = example_data$Zhang$metadata,
#' )
#'
#' plot_pca(
#'     ecoda_object,
#'     label_col = "Tissue",
#'     title = "PCA based on cell type composition",
#'     anosim_parallel = 1,
#'     n_hv_feat_show = 5 # Shows the most highly variable features (cell types)
#' )
#'
#' # Using only the most highly variable cell types
#' plot_pca(
#'     ecoda_object,
#'     slot = "clr_hvc",
#'     label_col = "Tissue",
#'     title = "PCA based on highly variable cell types",
#'     anosim_parallel = 1,
#'     n_hv_feat_show = ncol(slot(ecoda_object, "clr_hvc"))
#' )
plot_pca <- function(ecoda_object,
                     slot = c(
                         "clr", "clr_hvc", "counts", "counts_imp",
                         "freq", "freq_imp", "asin_sqrt", "pb"
                     ),
                     label_col = NULL,
                     scale. = FALSE,
                     title = NULL,
                     title_show_n_features = TRUE,
                     legend_title = "Group",
                     show_label_samples = FALSE,
                     score_digits = 3,
                     cluster_score = TRUE,
                     mod_score = TRUE,
                     sil_score = FALSE,
                     anosim_score = TRUE,
                     anosim_distance = "euclidian",
                     anosim_permutations = 99,
                     anosim_parallel = detectCores(),
                     ari_nclusts = NULL,
                     knn_k = 3,
                     pointsize = 3,
                     labelsize = 4,
                     coord_equal = TRUE,
                     axes = c(1, 2),
                     invisible = c("var", "quali"),
                     geom = "point",
                     n_hv_feat_show = Inf,
                     repel = TRUE) {
    slot <- match.arg(slot)
    feat_mat <- slot(ecoda_object, slot)

    res.pca <- prcomp(feat_mat, scale. = scale.)


    if (!is.null(label_col)) {
        labels <- slot(ecoda_object, "metadata")[[label_col]]

        if (title_show_n_features) {
            if (slot == "clr_hvc") {
                title <- paste0(
                    title,
                    "\nHVCs: ", slot(ecoda_object, "top_n_hvcs"),
                    " Variance explained: ",
                    round(
                        slot(ecoda_object, "variance_explained"),
                        score_digits
                    )
                )
            } else if (slot == "pb") {
                title <- paste0(title, "\nNumber of genes: ", ncol(feat_mat))
            } else {
                title <- paste0(
                    title, "\nNumber of cell types: ", ncol(feat_mat)
                )
            }
        }
        if (anosim_score) {
            anosim_score <- calc_anosim(
                feat_mat, labels,
                distance = anosim_distance,
                permutations = anosim_permutations,
                parallel = anosim_parallel,
                digits = score_digits
            )
            title <- paste0(title, "\nANOSIM score: ", anosim_score)
        }
        if (cluster_score) {
            cluster_score <- calc_ari(
                feat_mat, labels,
                nclusts = ari_nclusts, digits = score_digits
            )
            title <- paste0(title, "\nARI: ", cluster_score)
        }
        if (mod_score) {
            mod_score <- calc_modularity(
                feat_mat, labels, knn_k,
                digits = score_digits
            )
            title <- paste0(title, "\nModularity score: ", mod_score)
        }
        if (sil_score) {
            sil_score <- calc_sil(feat_mat, labels, digits = score_digits)
            title <- paste0(title, "\nSilhouette score: ", sil_score)
        }
    } else {
        labels <- "none"
    }

    if (!is.infinite(n_hv_feat_show) & all(invisible %in% c("var", "quali"))) {
        invisible <- "quali"
    }

    if (show_label_samples) {
        if (!"text" %in% geom) {
            geom <- c(geom, "text")
        }
    }

    p <- fviz_pca(
        res.pca,
        axes = axes,
        habillage = labels,
        pointsize = pointsize,
        labelsize = labelsize,
        invisible = invisible,
        select.var = list(contrib = n_hv_feat_show),
        repel = repel,
        geom = geom
    ) +
        ggtitle(title) +
        scale_color_discrete(name = legend_title)

    if (!is.null(label_col)) {
        p <- p + scale_shape_manual(
            values = rep(19, length(unique(labels))),
            name = legend_title
        )
    }

    if (coord_equal) {
        p <- p + coord_equal()
    }

    return(p)
}


#' @title Plot 3-dimensional interactive Principal Component Analysis plot
#'
#' @description Performs Principal Component Analysis (PCA) on a selected data
#'   matrix from the \code{ECODA} object (default: CLR-transformed abundances,
#'   \code{clr}) and visualizes the results in 3D colored by \code{label_col}.
#'
#' @param ecoda_object An \link[=ECODA-class]{ECODA} object.
#' @param slot Character string (default: \code{"clr"}). The name of the data
#'   matrix slot in the \code{ECODA} object to use for PCA. Must be one of:
#'   \code{"clr"} (CLR-transformed abundances, default), \code{"clr_hvc"}
#'   (CLR-transformed abundances of only the most highly variable cell types
#'   (HVCs)), \code{"pb"} (pseudobulk gene expression), \code{"counts"} (raw
#'   counts), \code{"counts_imp"} (imputed counts), \code{"freq"} (relative
#'   frequencies), \code{"freq_imp"} (imputed frequencies), or
#'   \code{"asin_sqrt"} (arcsin-square root transformed data).
#' @param label_col Character string (optional, default: \code{NULL}). The name
#'   of a column in \code{slot(ecoda_object, "metadata")} used to color and
#'   group samples in the plot, and for calculating clustering scores.
#' @param scale. Logical (default: \code{FALSE}). A value indicating whether the
#'   variables should be scaled to have unit variance before the PCA.
#'
#' @return An interactive 3D \code{plotly} object visualizing the PCA results.
#'
#' @importFrom stats prcomp
#' @importFrom plotly plot_ly add_markers add_paths group2NA
#'
#' @export plot_pca3d
#'
#' @examples
#' data(example_data)
#' ecoda_object <- ecoda(
#'     data = example_data$Zhang$cell_counts_lowresolution,
#'     metadata = example_data$Zhang$metadata,
#' )
#'
#' plot_pca3d(ecoda_object, label_col = "Tissue")
plot_pca3d <- function(ecoda_object,
                       slot = c(
                           "clr", "clr_hvc", "counts", "counts_imp",
                           "freq", "freq_imp", "asin_sqrt", "pb"
                       ),
                       label_col = NULL,
                       scale. = FALSE) {
    slot <- match.arg(slot)
    feat_mat <- slot(ecoda_object, slot)

    res.pca <- prcomp(feat_mat, scale. = scale.)

    if (!is.null(label_col)) {
        labels <- slot(ecoda_object, "metadata")[[label_col]]
    } else {
        labels <- "none"
    }

    df <- as.data.frame(res.pca$x)
    df$id <- seq_len(nrow(df))
    df$vs <- factor(labels)
    ms <- replicate(2, df, simplify = FALSE)
    ms[[2]]$PC3 <- min(df$PC3)
    m <- ms %>%
        bind_rows() %>%
        group2NA("id", "vs")

    p <- plot_ly(color = ~vs) %>%
        add_markers(data = df, x = ~PC1, y = ~PC2, z = ~PC3) %>%
        add_paths(data = m, x = ~PC1, y = ~PC2, z = ~PC3, opacity = 0.2)

    return(p)
}


#' Analysis of Similarities (ANOSIM) R score
#'
#' @description Calculates the ANOSIM R-statistic to test whether there is
#'   significant separation between two or more groups (defined by
#'   \code{labels}) based on the multivariate distances among samples in the
#'   feature space (\code{feat_mat}).
#'
#' @details ANOSIM compares the mean of rank dissimilarities between groups to
#'   the mean of rank dissimilarities within groups. The R-statistic ranges from
#' -1 to 1:
#' \itemize{
#'   \item An R value close to **1** indicates clear separation of groups.
#'   \item An R value close to **0** indicates that the separation is
#'         no greater than expected by chance (i.e., poor separation).
#'   \item An R value close to **-1** indicates that within-group
#'         dissimilarities are greater than between-group dissimilarities
#'         (a very rare result).
#' }
#'
#' @param feat_mat Numeric matrix or data frame. The feature matrix (e.g., CLR
#'   abundances, with samples as rows and features as columns).
#' @param labels Vector of factors or character strings. The grouping variable
#'   (the known cluster assignments) for each row in the matrix.
#' @param permutations Integer (default: \code{999}). The number of permutations
#'   to use when calculating the ANOSIM R-statistic and p-value.
#' @param parallel Integer (optional, default: \code{detectCores()}). The number
#'   of parallel processes/cores to use for the permutation testing.
#' @param digits Integer (default: \code{3}). The number of decimal places to
#'   round the R-statistic.
#'
#' @return A numeric value representing the **ANOSIM R-statistic**.
#'
#' @importFrom vegan anosim
#'
#' @export calc_anosim
#'
#' @examples
#' data(example_data)
#' ecoda_object <- ecoda(
#'     data = example_data$GongSharma_full$cell_counts_highresolution,
#'     metadata = example_data$GongSharma_full$metadata
#' )
#'
#' # Extract necessary components
#' feat_mat <- slot(ecoda_object, "clr")
#' labels <- slot(ecoda_object, "metadata")$subject.cmv
#'
#' # Run the calculation
#' calc_anosim(feat_mat, labels, parallel = 1)
calc_anosim <- function(feat_mat,
                        labels,
                        permutations = 99,
                        parallel = 1,
                        digits = 3) {
    score <- anosim(
        x = dist(feat_mat),
        grouping = labels,
        permutations = permutations,
        parallel = parallel
    )[["statistic"]]

    return(round(score, digits))
}


#' Calculate Adjusted Rand Index (ARI)
#'
#' Calculates the Adjusted Rand Index (ARI) to measure the agreement between the
#' known cluster assignments (\code{labels}) and the cluster assignments derived
#' from two unsupervised clustering methods: Hierarchical Clustering
#' (\code{hclust}) and Partitioning Around Medoids (\code{pam}).
#'
#' @param feat_mat Numeric matrix or data frame. The feature matrix (e.g., CLR
#'   abundances).
#' @param labels Vector of factors or character strings. The true or known
#'   cluster assignments for each row in the matrix.
#' @param nclusts Integer (optional, default: \code{NULL}). The target number of
#'   clusters (\code{k}) to use for \code{hclust} and \code{pam}. If
#'   \code{NULL}, it defaults to the number of unique levels in \code{labels}.
#' @param digits Integer (default: \code{3}). The number of decimal places to
#'   round the result.
#' @param return_mean Logical (default: \code{TRUE}). If \code{TRUE}, returns
#'   the mean of the ARI scores from \code{hclust} and \code{pam}. If
#'   \code{FALSE}, returns a named list with both individual scores.
#'
#' @return A numeric value (mean ARI) or a list of two ARI scores. ARI ranges
#'   from
#'         -1 (disagreement) to +1 (perfect agreement).
#'
#' @importFrom stats dist hclust cutree
#' @importFrom cluster pam
#' @importFrom mclust adjustedRandIndex
#'
#' @export calc_ari
#'
#' @examples
#' data(example_data)
#' ecoda_object <- ecoda(
#'     data = example_data$GongSharma_full$cell_counts_highresolution,
#'     metadata = example_data$GongSharma_full$metadata
#' )
#'
#' # Extract necessary components
#' feat_mat <- slot(ecoda_object, "clr")
#' labels <- slot(ecoda_object, "metadata")$subject.cmv
#'
#' # Run the calculation
#' calc_ari(feat_mat, labels)
calc_ari <- function(feat_mat,
                     labels,
                     nclusts = NULL,
                     digits = 3,
                     return_mean = TRUE) {
    results <- list()
    dist_mat <- dist(feat_mat)

    if (is.null(nclusts)) {
        nclusts <- length(unique(labels))
    }

    # Perform hierarchical clustering
    hc <- hclust(dist_mat, method = "ward.D2")
    clust_labels <- cutree(hc, k = nclusts)
    results[["hclust_accuracy"]] <- adjustedRandIndex(
        as.numeric(as.factor(labels)), clust_labels
    )

    # Perform PAM clustering
    clust_labels <- pam(dist_mat, k = nclusts)$cluster
    results[["pamclust_accuracy"]] <- adjustedRandIndex(
        as.numeric(as.factor(labels)), clust_labels
    )

    if (return_mean) {
        return(round(mean(unlist(results)), digits))
    } else {
        results[["hclust_accuracy"]] <- round(
            results[["hclust_accuracy"]], digits
        )
        results[["pamclust_accuracy"]] <- round(
            results[["pamclust_accuracy"]], digits
        )
        return(results)
    }
}


#' Calculate Adjusted Modularity Score
#'
#' Calculates the Modularity score for a given clustering (\code{labels}) based
#' on a Shared Nearest Neighbor (SNN) graph.The score is adjusted by the
#' theoretical maximum modularity for the number of groups to always be between
#' -0.5 (poor clustering) and +1 (excellent clustering)).
#'
#' @param feat_mat Numeric matrix or data frame. The feature matrix used to
#'   compute distances (e.g., CLR abundances).
#' @param labels Vector of factors or character strings. The cluster assignments
#'   for each row in \code{feat_mat}.
#' @param knn_k Integer (optional, default: \code{NULL}). The number of nearest
#'   neighbors (\code{k}) used for SNN graph construction. If \code{NULL}, it
#'   defaults to \code{max(3, round(sqrt(N)))}, where \code{N} is the number of
#'   samples.
#' @param digits Integer (default: \code{3}). The number of decimal places to
#'   round the final adjusted score.
#'
#' @return A numeric value representing the adjusted modularity score,
#'   ranging from -0.5 to 1.0.
#'
#' @references
#' Brandes, Ulrik, et al.
#' "On finding graph clusterings with maximum modularity."
#' European Symposium on Algorithms. Springer, Berlin, Heidelberg, 2007.
#'
#' @importFrom igraph modularity
#'
#' @export calc_modularity
#'
#' @examples
#' data(example_data)
#' ecoda_object <- ecoda(
#'     data = example_data$GongSharma_full$cell_counts_highresolution,
#'     metadata = example_data$GongSharma_full$metadata
#' )
#'
#' # Extract necessary components
#' feat_mat <- slot(ecoda_object, "clr")
#' labels <- slot(ecoda_object, "metadata")$subject.cmv
#'
#' # Run the calculation
#' calc_modularity(feat_mat, labels)
calc_modularity <- function(feat_mat, labels, knn_k = 3, digits = 3) {
    ngroups <- length(unique(labels))

    if (is.null(knn_k)) {
        knn_k <- max(3, round(sqrt(nrow(feat_mat))))
    }

    # Create a graph object
    dist_mat <- dist(feat_mat)
    knn <- compute_KNN_from_dist(dist_mat, knn_k)
    g <- compute_snn_graph(knn)

    # Compute modularity
    modularity_score <- modularity(g, membership = as.numeric(factor(labels)))

    # NOTE:
    # Maximum modularity depends on the number of groups:
    # max(mod) = 1 - 1 / (number of groups)
    # see Brandes, Ulrik, et al.
    # "On finding graph clusterings with maximum modularity."

    # Adjust modularity score for number of groups
    maximum_modularity_score <- 1 - (1 / ngroups)
    adjusted_modularity_score <- modularity_score / maximum_modularity_score

    return(round(adjusted_modularity_score, digits))
}


#' Compute K-Nearest Neighbors (KNN) from Distance Matrix
#'
#' Identifies the indices of the K-nearest neighbors for each sample based on a
#' pre-computed distance matrix.
#'
#' @param dist_mat A square distance matrix or an object of class \code{dist}.
#' @param knn_k Integer. The number of neighbors (\code{k}) to return for each
#'   sample.
#'
#' @return A matrix with \code{nrow(dist_mat)} rows and \code{knn_k} columns,
#'   where each row contains the indices of the nearest neighbors.
#'
#' @importFrom RANN nn2
compute_KNN_from_dist <- function(dist_mat, knn_k) {
    # dist_mat should be a square matrix or 'dist' object
    dist_mat <- as.matrix(dist_mat)
    n <- nrow(dist_mat)
    knn <- matrix(0, nrow = n, ncol = knn_k)

    for (i in 1:n) {
        # Sort distances in row i, get indices of the 2nd to (k+1)th closest
        # (index 1 is always the node itself)
        knn[i, ] <- order(dist_mat[i, ])[2:(knn_k + 1)]
    }
    return(knn)
}


#' Compute Shared Nearest Neighbor (SNN) Graph
#'
#' Constructs a graph where nodes represent samples and edges are weighted
#' by the number of shared nearest neighbors (SNN) between them.
#'
#' @param knn A matrix of nearest neighbor indices, typically generated by
#'   \code{compute_KNN_from_dist}.
#'
#' @return An \code{igraph} graph object where edge weights correspond to the
#'   count of shared neighbors between nodes.
#'
#' @importFrom igraph graph_from_adjacency_matrix
compute_snn_graph <- function(knn) {
    # Initialize adjacency matrix
    n <- nrow(knn)
    adj_matrix <- matrix(0, n, n)
    # Count shared neighbors
    for (i in seq_len(n)) {
        for (j in knn[i, ]) {
            shared_neighbors <- length(intersect(knn[i, ], knn[j, ]))
            adj_matrix[i, j] <- shared_neighbors
            adj_matrix[j, i] <- shared_neighbors # Ensure symmetry
        }
    }
    # Create graph object
    g <- igraph::graph_from_adjacency_matrix(adj_matrix,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
    )
    return(g)
}


#' Calculate Average Silhouette Width
#'
#' Calculates the average silhouette width for a given feature matrix, using
#' pre-defined cluster assignments (\code{labels}). This metric assesses the
#' quality of the clustering.
#'
#' @param feat_mat Numeric matrix or data frame. The feature matrix (e.g.,
#'   CLR-transformed abundances or PCA scores) on which the distance calculation
#'   is based.
#' @param labels Vector of factors or character strings. The cluster assignments
#'   for each row in \code{feat_mat} (i.e., the grouping variable).
#' @param digits Integer (default: \code{3}). The number of decimal places to
#'   round the final score.
#'
#' @return A numeric value representing the mean silhouette width, typically
#'   ranging from -1 (poor clustering) to +1 (excellent clustering).
#'
#' @importFrom cluster silhouette
#' @importFrom stats dist
#'
#' @export calc_sil
#'
#' @examples
#' data(example_data)
#' ecoda_object <- ecoda(
#'     data = example_data$GongSharma_full$cell_counts_highresolution,
#'     metadata = example_data$GongSharma_full$metadata
#' )
#'
#' # Extract necessary components
#' feat_mat <- slot(ecoda_object, "clr")
#' labels <- slot(ecoda_object, "metadata")$subject.cmv
#'
#' # Run the calculation
#' calc_sil(feat_mat, labels)
calc_sil <- function(feat_mat, labels, digits = 3) {
    sils <- silhouette(
        x = as.numeric(factor(labels)),
        dist = dist(feat_mat)
    ) %>%
        as.data.frame()
    score <- mean(sils[["sil_width"]])
    return(round(score, digits))
}


# Box and bar plots ---------------------------

#' Reshapes ECODA data into a long format for plotting and analysis.
#'
#' This function takes either the relative abundance (\code{freq}) or
#' CLR-transformed abundance (\code{clr}) matrix from an
#' \link[=ECODA-class]{ECODA} object, converts it from a wide (samples x cell
#' types) to a long (sample, celltype, value) format, and optionally joins it
#' with a specified column from the sample metadata.
#'
#' @param ecoda_object An initialized \link[=ECODA-class]{ECODA} object.
#' @param data_slot Character string specifying the data matrix to use. Must be
#'   either \code{"freq"} (for relative abundance) or \code{"clr"} (for
#'   CLR-transformed abundance).
#' @param label_col Character string (optional, default: \code{NULL}). The name
#'   of a column in the \code{slot(ecoda_object, "metadata")} slot to merge into
#'   the long data frame (e.g., "Disease_State" or "Batch"). If \code{NULL},
#'   only the abundance data and sample/celltype IDs are returned.
#'
#' @return A tidy, long format data frame with columns:
#'         \itemize{
#'           \item \code{sample_id}: Sample identifier.
#'           \item \code{celltype}: Cell type name.
#'           \item \code{rel_abundance} or \code{clr_abundance}:
#'                 The quantitative value depending on the chosen
#'                 \code{data_slot}.
#'           \item \code{...}: Additional column specified by
#'                 \code{label_col} (if provided).
#'         }
#'
#' @importFrom methods slot slotNames
#' @importFrom dplyr %>% left_join
#' @importFrom tidyr pivot_longer everything
#' @importFrom rlang sym
#'
#' @seealso \link[=ECODA-class]{ECODA}
create_long_data <- function(ecoda_object,
                             data_slot,
                             label_col = NULL) {
    # Ensure the data_slot is valid (freq or clr)
    if (!(data_slot %in% c("freq", "clr"))) {
        stop(
            "Invalid data_slot. Must be 'freq' (Relative Abundance) ",
            "or 'clr' (CLR Abundance)."
        )
    }

    # Determine the appropriate value column name based on the slot
    value_name <- ifelse(data_slot == "freq", "rel_abundance", "clr_abundance")

    # 1. Extract data (either @freq or @clr)
    data_matrix <- slot(ecoda_object, data_slot)
    data_df <- as.data.frame(data_matrix)
    data_df <- data.frame(
        sample_id = rownames(data_df),
        data_df
    )

    # 2. Reshape the data from wide to long format
    long_data <- data_df %>%
        pivot_longer(
            cols = -sample_id,
            names_to = "celltype",
            values_to = value_name
        )

    # 3. Prepare the metadata (only if label_col is provided)
    if (!is.null(label_col)) {
        # --- Check 1: Ensure metadata exists ---
        if (!("metadata" %in% slotNames(ecoda_object)) ||
            is.null(slot(ecoda_object, "metadata"))) {
            stop(
                "label_col was provided but ecoda_object@metadata slot ",
                "is missing or NULL."
            )
        }

        # --- Check 2: Ensure label_col exists in metadata ---
        meta_colnames <- colnames(slot(ecoda_object, "metadata"))
        if (!(label_col %in% meta_colnames)) {
            stop(
                "label_col '", label_col,
                "' not found in ecoda_object@metadata."
            )
        }

        metadata_df <- as.data.frame(slot(ecoda_object, "metadata"))
        metadata_df <- data.frame(
            sample_id = rownames(metadata_df),
            metadata_df
        )
        metadata_df <- metadata_df[, c("sample_id", label_col)]

        # 4. Join the long data with the metadata by sample_id
        plot_data <- long_data %>%
            left_join(metadata_df, by = "sample_id")
    } else {
        # If no label_col, just return the long data without joining metadata
        plot_data <- long_data
    }

    return(plot_data)
}


#' Generates a Stacked Bar Plot of Cell Type Relative Abundance.
#'
#' This function visualizes the relative abundance (composition) of cell types
#' across samples or across aggregated groups. It automatically handles data
#' preparation and ordering based on the provided parameters.
#'
#' @param ecoda_object An \link[=ECODA-class]{ECODA} object containing cell type
#'   relative frequencies in the \code{freq} slot.
#' @param label_col Character string (optional, default: \code{NULL}). The name
#'   of a column in \code{slot(ecoda_object, "metadata")} used to define
#'   grouping or groups (required if \code{plot_by = "group"}).
#' @param plot_by Character string (default: \code{"sample"}). Specifies whether
#'   to plot the relative abundance for each individual sample (\code{"sample"})
#'   or the average relative abundance aggregated by a group (\code{"group"})
#'   defined by \code{label_col}.
#' @param custom_sample_order Character vector (optional, default: \code{NULL}).
#'   A vector of sample IDs to enforce a specific order when \code{plot_by =
#'   "sample"}. If \code{NULL}, samples are ordered first by \code{label_col}
#'   (if provided) and then naturally.
#' @param title Character string (default: \code{""}). The main title for the
#'   plot.
#' @param facet_by_label_col Logical (default: \code{TRUE}). If \code{TRUE} and
#'   \code{plot_by = "sample"}, the plot will be faceted (split) horizontally by
#'   the categories in \code{label_col}.
#'
#' @return A \code{ggplot} object representing the stacked bar plot.
#'
#' @importFrom dplyr %>% group_by summarise mutate distinct arrange pull
#' @importFrom ggplot2 ggplot aes geom_col theme_minimal theme element_text labs
#'   facet_grid
#' @importFrom rlang sym
#' @importFrom gtools mixedsort
#' @importFrom stats reformulate
#'
#' @export plot_barplot
#'
#' @seealso \code{\link{create_long_data}}, \link[=ECODA-class]{ECODA}
#'
#' @examples
#' data(example_data)
#' ecoda_object <- ecoda(
#'     data = example_data$Zhang$cell_counts_lowresolution,
#'     metadata = example_data$Zhang$metadata,
#' )
#'
#' plot_barplot(ecoda_object)
#'
#' # Plotting average cell type abundance by experimental group
#' plot_barplot(
#'     ecoda_object,
#'     label_col = "Tissue",
#'     plot_by = "group",
#'     title = "Mean Relative Abundance by Condition"
#' )
#'
#' # Plotting cell type abundance for each sample separately
#' plot_barplot(
#'     ecoda_object,
#'     label_col = "Tissue",
#'     plot_by = "sample",
#'     title = "Relative Abundance for Each Sample"
#' )
plot_barplot <- function(ecoda_object,
                         label_col = NULL,
                         plot_by = c("sample", "group"),
                         custom_sample_order = NULL,
                         title = "",
                         facet_by_label_col = TRUE) {
    # Use the helper function to get the long data from @freq
    plot_data <- create_long_data(
        ecoda_object,
        data_slot = "freq",
        label_col = label_col
    )

    plot_by <- match.arg(plot_by)

    if (plot_by == "group") {
        if (is.null(label_col)) {
            stop("label_col must be provided when plot_by = 'group'")
        }

        # Plotting by group
        plot_df <- plot_data %>%
            group_by(celltype, !!sym(label_col)) %>%
            summarise(
                mean_rel_abund = mean(rel_abundance, na.rm = TRUE),
                .groups = "drop"
            ) %>%
            mutate(x_var = !!sym(label_col), y_var = mean_rel_abund)

        # Ensure group is a factor
        plot_df$x_var <- factor(plot_df$x_var)
    } else if (plot_by == "sample") {
        # Plotting by Sample
        plot_df <- plot_data %>%
            mutate(x_var = sample_id, y_var = rel_abundance)

        current_levels <- unique(plot_df$x_var)

        # --- Sample Ordering Logic ---
        if (!is.null(custom_sample_order)) {
            if (!all(current_levels %in% custom_sample_order)) {
                stop(
                    "Custom sample order is incomplete",
                    "or contains unknown sample IDs"
                )
            }
            ordered_levels <- custom_sample_order
        } else {
            # DEFAULT: Sort samples by label_col (group) if provided

            if (!is.null(label_col)) {
                # 1. Get a unique list of samples
                # and their corresponding label_col
                sample_group_map <- plot_data %>%
                    distinct(sample_id, !!sym(label_col))

                # 2. Order the samples by the group column first
                # (using mixedsort on its values) and then
                # by sample_id (using mixedsort on its values).
                ordered_levels <- sample_group_map %>%
                    mutate(
                        ordered_group = factor(
                            !!sym(label_col),
                            levels = mixedsort(unique(!!sym(label_col))),
                            ordered = TRUE
                        )
                    ) %>%
                    arrange(ordered_group, mixedsort(sample_id)) %>%
                    pull(sample_id)
            } else {
                # If no label_col, just sort samples naturally by sample_id
                ordered_levels <- mixedsort(unique(plot_data$sample_id))
            }
        }

        # Apply the ordered factor to the x_var column
        plot_df$x_var <- factor(plot_df$x_var, levels = ordered_levels)
    }

    # --- Generate the plot ---
    p <- ggplot(plot_df, aes(x = x_var, y = y_var, fill = celltype)) +
        geom_col(position = "stack") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.title = element_text(face = "bold")
        ) +
        labs(
            title = title,
            x = "",
            y = "Relative Abundance (%)",
            fill = "Cell Type"
        )

    if (!is.null(label_col) & plot_by == "sample" & facet_by_label_col) {
        p <- p +
            facet_grid(reformulate(label_col), scales = "free_x")
    }

    return(p)
}


#' Generates Boxplots for CLR-transformed Cell Type Abundances with Optional
#' Group Comparison.
#'
#' This function visualizes the distribution of CLR-transformed abundance for
#' each cell type using boxplots, optionally splitting the data by a sample
#' metadata column and performing statistical tests for comparison between
#' groups.
#'
#' \strong{Statistical Test Logic:}
#' \itemize{
#'   \item If the number of groups (\code{label_col} levels) is \strong{2}, the
#'         function uses the specified \code{stat_method}
#'         (default: "wilcox.test") and displays pairwise significance labels
#'         or stars.
#'   \item If the number of groups is \strong{3 or more}, and
#'         \code{stat_method} is "wilcox.test" or "t.test", the function
#'         automatically switches to the global non-parametric test,
#'         \strong{"kruskal.test"}, and displays the overall p-value
#'         for the comparison across all groups.
#' }
#'
#' @param ecoda_object An \link[=ECODA-class]{ECODA} object containing the
#'   CLR-transformed abundances in the \code{clr} slot.
#' @param label_col Character string (optional, default: \code{NULL}). The name
#'   of a column in \code{slot(ecoda_object, "metadata")} used to define groups
#'   for comparison. If \code{NULL}, a single boxplot is generated per cell
#'   type.
#' @param selected_celltypes Specify selected celltypes you want to plot instead
#'   of plotting boxplots for all.
#' @param title Character string (default: \code{""}). The main title for the
#'   plot.
#' @param stat_method Character string (default: \code{"wilcox.test"}). The
#'   statistical method used for comparisons between 2 groups (e.g., "t.test",
#'   "wilcox.test"). Note: This is overridden by "kruskal.test" for 3+ groups.
#' @param paired Logical (default: \code{FALSE}). If \code{TRUE}, performs a
#'   paired statistical test (e.g., paired t-test or paired Wilcoxon test). Only
#'   applicable for 2-group comparisons.
#' @param signif_label Character string (default: \code{"p.signif"}). Controls
#'   how p-values are displayed for 2-group comparisons (e.g., "p.signif" for
#'   stars, "p.format" for numeric p-value).
#'
#' @return A \code{ggplot} object (enhanced by \code{ggpubr}) representing the
#'   CLR abundance boxplot, including jittered data points and dynamic
#'   significance information (pairwise stars for 2 groups, overall p-value for
#'   3+ groups).
#'
#' @importFrom ggplot2 aes geom_jitter labs theme element_text guides
#'   position_jitterdodge
#' @importFrom ggpubr ggboxplot stat_compare_means stat_pvalue_manual
#' @importFrom stringr str_to_title
#' @importFrom rlang sym
#' @importFrom rstatix wilcox_test add_xy_position
#' @importFrom stats as.formula
#'
#' @export plot_boxplot
#'
#' @seealso \code{\link{create_long_data}}, \link[=ECODA-class]{ECODA}
#'
#' @examples
#' data(example_data)
#' ecoda_object <- ecoda(
#'     data = example_data$Zhang$cell_counts_lowresolution,
#'     metadata = example_data$Zhang$metadata,
#' )
#'
#' # 1. Boxplots for CLR abundance without grouping (no stats calculated):
#' plot_boxplot(ecoda_object)
#'
#' # 2. Boxplots grouped by 'Treatment' (2 groups) and applying Wilcoxon test:
#' plot_boxplot(
#'     ecoda_object,
#'     label_col = "Tissue",
#'     stat_method = "wilcox.test",
#'     title = "CLR Abundance by Tissue (with Wilcoxon Test)"
#' )
plot_boxplot <- function(ecoda_object,
                         label_col = NULL,
                         selected_celltypes = NULL,
                         title = "",
                         stat_method = "wilcox.test",
                         paired = FALSE,
                         signif_label = c("p.signif", "p.format")) {
    signif_label <- match.arg(signif_label)
    plot_data <- create_long_data(
        ecoda_object,
        data_slot = "clr",
        label_col = label_col
    )

    if (!is.null(selected_celltypes)) {
        if (!all(selected_celltypes %in% unique(plot_data$celltype))) {
            stop(
                "Not all selected_celltypes found in ecoda_object. ",
                "Please check colnames(ecoda_object@clr)."
            )
        }

        plot_data <- plot_data[plot_data$celltype %in% selected_celltypes, ]
    }

    # Ensure celltype is a factor for plotting
    plot_data$celltype <- factor(plot_data$celltype)

    # --- Setup for Plotting ---

    # Determine color mapping for ggboxplot
    box_color_map <- if (!is.null(label_col)) label_col else "black"

    # Generate the base boxplot
    p <- ggboxplot(plot_data,
        x = "celltype",
        y = "clr_abundance",
        xlab = "",
        ylab = "CLR Abundance",
        outlier.shape = NA,
        color = box_color_map # Color by group variable or "black"
    )

    # --- Add Jittered Points and Significance ---

    if (is.null(label_col)) {
        # SCENARIO 1: No grouping variable provided
        # (Single boxplot per cell type)

        # Add simple jitter (no dodging)
        p <- p + geom_jitter(
            width = 0.2,
            size = 1,
            alpha = 0.6
        ) +
            guides(color = "none")
    } else {
        # SCENARIO 2: Grouping variable is provided (Comparison between groups)

        # Calculate the number of boxplots for correct jitter-dodging
        nr_of_boxplots <- length(unique(plot_data[[label_col]]))
        label_col_sym <- sym(label_col)

        # Add jittered points with dodging
        p <- p + geom_jitter(
            mapping = aes(color = !!label_col_sym), # Map color to group
            position = position_jitterdodge(jitter.width = 1 / nr_of_boxplots),
            size = 1,
            alpha = 0.6
        )

        # Add significance testing
        if (stat_method == "kruskal.test") {
            # Kruskal-Wallis: Overall test (no 'group' aesthetic needed)
            p <- p + stat_compare_means(
                method = stat_method,
                label.y.npc = "top", # Place label at the top
                label.x.npc = "center",
                label = "p.format" # Show the overall p-value
            )
        } else if (nr_of_boxplots == 2) {
            # Wilcoxon or t.test: Pairwise test (requires 'group' aesthetic)
            p <- p + stat_compare_means(
                aes(group = !!label_col_sym),
                method = stat_method,
                paired = paired,
                label = signif_label, # Show significance stars
                tip.length = 0,
                hide.ns = TRUE
            )
        } else if (nr_of_boxplots > 2) {
            y_var <- colnames(plot_data)[3]
            x_group_var <- colnames(plot_data)[2]
            fill_compare_var <- colnames(plot_data)[4]

            dsub_stats <- plot_data %>%
                group_by(!!sym(x_group_var)) %>%
                wilcox_test(as.formula(paste(y_var, "~", fill_compare_var))) %>%
                add_xy_position(x = x_group_var)

            p <- p +
                # Add p-values using the generated stats table
                stat_pvalue_manual(dsub_stats,
                    label = "p.adj.signif",
                    tip.length = 0.01
                )
        }

        # Add legend title
        p <- p + labs(color = str_to_title(label_col))
    }

    # --- Final Theme and Labels ---
    p <- p +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, ),
            legend.title = element_text(face = "bold")
        ) +
        labs(title = title)

    return(p)
}


# Heatmap ---------------------------

#' Generates a Heatmap of Cell Abundance Data from an ECODA Slot.
#'
#' This function visualizes a data matrix from a specified slot (e.g.,
#' CLR-transformed, frequency, or pseudobulk data) after mean-centering. It
#' supports optional filtering to only Highly Variable Cell Types (HVCs) and
#' includes a sample annotation sidebar based on a specified metadata column.
#'
#' @param ecoda_object An \link[=ECODA-class]{ECODA} object.
#' @param slot Character string (default: \code{"clr"}). The name of the data
#'   matrix slot in the \code{ECODA} object to use for the heatmap. Must be one
#'   of: \code{"clr"} (CLR-transformed abundances, default), \code{"clr_hvc"}
#'   (CLR-transformed abundances of only the most highly variable cell types
#'   (HVCs)), \code{"pb"} (pseudobulk gene expression), \code{"counts"} (raw
#'   counts), \code{"counts_imp"} (imputed counts), \code{"freq"} (relative
#'   frequencies), \code{"freq_imp"} (imputed frequencies), or
#'   \code{"asin_sqrt"} (arcsin-square root transformed data).
#' @param label_col Character string. The name of the column in
#'   \code{slot(ecoda_object, "metadata")} to use for annotating the samples
#'   (columns) of the heatmap.
#' @param cluster_rows Logical (default: \code{TRUE}). Whether to apply
#'   hierarchical clustering to the cell types (rows).
#' @param cluster_cols Logical (default: \code{TRUE}). Whether to apply
#'   hierarchical clustering to the samples (columns).
#' @param scale Character string (default: \code{"none"}). Method for scaling
#'   the abundance values within the heatmap using \code{pheatmap}. Options
#'   include \code{"none"}, \code{"row"}, or \code{"column"}. Note: The data is
#'   internally **mean-centered** (\code{scale(center=TRUE, scale=FALSE)})
#'   across samples before being passed to \code{pheatmap}, regardless of this
#'   argument.
#' @param clustering_method Character string (default: \code{"ward.D2"}). The
#'   clustering method to use for hierarchical clustering. Options are passed
#'   directly to \code{hclust} (e.g., \code{"complete"}, \code{"average"},
#'   \code{"ward.D2"}).
#' @param angle_col Character string (default: \code{"90"}). Angle of the sample
#'   labels (columns).
#' @param ... Additional arguments passed directly to the \code{pheatmap}
#'   function.
#'
#' @return A \code{pheatmap} object, which is typically visualized automatically
#'   when called interactively.
#'
#' @importFrom dplyr %>%
#' @importFrom pheatmap pheatmap
#'
#' @export plot_heatmap
#'
#' @seealso \link[=ECODA-class]{ECODA}, \code{\link{pheatmap}}
#'
#' @examples
#' # Example for a simple dataset:
#' data(example_data)
#' ecoda_object <- ecoda(
#'     data = example_data$Zhang$cell_counts_lowresolution,
#'     metadata = example_data$Zhang$metadata,
#' )
#'
#' plot_heatmap(ecoda_object, label_col = c("Clinical.efficacy.", "Tissue"))
#'
#' plot_heatmap(
#'     ecoda_object,
#'     label_col = c("Clinical.efficacy.", "Tissue"),
#'     # Additional arguments for pheatmap:
#'     cutree_rows = 3,
#'     cutree_cols = 3
#' )
#'
#' # Example of a large cohort with 868 samples and 69 cell types
#' ecoda_object <- ecoda(
#'     data = example_data$GongSharma_full$cell_counts_highresolution,
#'     metadata = example_data$GongSharma_full$metadata
#' )
#'
#' plot_heatmap(
#'     ecoda_object,
#'     label_col = c("subject.cmv", "age_group"),
#'     cutree_rows = 3,
#'     cutree_cols = 5,
#'     show_colnames = FALSE
#' )
#'
#' # Using only the most highly variable cell types (HVCs)
#' plot_heatmap(
#'     ecoda_object,
#'     slot = "clr_hvc",
#'     label_col = c("subject.cmv", "age_group"),
#'     cutree_rows = 3,
#'     cutree_cols = 4,
#'     show_colnames = FALSE
#' )
plot_heatmap <- function(ecoda_object,
                         slot = c(
                             "clr", "clr_hvc", "counts", "counts_imp",
                             "freq", "freq_imp", "asin_sqrt", "pb"
                         ),
                         label_col,
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         scale = "none",
                         clustering_method = "ward.D2",
                         angle_col = "90",
                         ...) {
    slot <- match.arg(slot)
    df_heatmap <- slot(ecoda_object, slot)

    df_heatmap <- df_heatmap %>%
        scale(center = TRUE, scale = FALSE) %>%
        t() %>%
        as.data.frame()

    metadata <- slot(ecoda_object, "metadata")[, label_col, drop = FALSE]
    metadata[] <- lapply(metadata, as.factor)

    heatmap <- pheatmap(
        df_heatmap,
        annotation_col = metadata,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        scale = scale,
        clustering_method = clustering_method,
        angle_col = angle_col,
        ...
    )

    return(heatmap)
}


# Correlation plot ---------------------------

#' @title Plot Cell Type Correlation Matrix
#'
#' @description Calculates the pairwise Pearson correlation matrix for all cell
#'   types (columns) using the Centered Log-Ratio (CLR) transformed abundances
#'   stored in \code{slot(ecoda_object, "clr")}. It then visualizes this matrix
#'   as a heatmap using \code{corrplot::corrplot}.
#'
#' @details The function uses the CLR matrix, where high correlation between two
#'   cell types suggests they vary together across samples, indicating potential
#'   co-occurrence or co-regulation.
#'
#' @param ecoda_object An \link[=ECODA-class]{ECODA} object.
#' @param slot Character string (default: \code{"clr"}). The name of the data
#'   matrix slot in the \code{ECODA} object to use for the correlation plot Must
#'   be one of: \code{"clr"} (CLR-transformed abundances, default),
#'   \code{"clr_hvc"} (CLR-transformed abundances of only the most highly
#'   variable cell types (HVCs)), \code{"pb"} (pseudobulk gene expression),
#'   \code{"counts"} (raw counts), \code{"counts_imp"} (imputed counts),
#'   \code{"freq"} (relative frequencies), \code{"freq_imp"} (imputed
#'   frequencies), or \code{"asin_sqrt"} (arcsin-square root transformed data).
#' @param order Character string (default: \code{"hclust"}). The ordering method
#'   for the correlation matrix. Common options include:
#'              \itemize{
#'                \item \code{"original"} (no reordering)
#'                \item \code{"hclust"} (hierarchical clustering)
#'                \item \code{"FPC"} (first principal component order)
#'              }
#' @param hclust.method Character string (default: \code{"ward.D2"}). The
#'   hierarchical clustering method to use if \code{order} is set to
#'   \code{"hclust"}.
#' @param ... Additional arguments passed to \code{corrplot::corrplot} for plot
#'   customization (e.g., \code{method}, \code{type}, \code{tl.col}).
#'
#' @return A plot object generated by \code{corrplot::corrplot}, which is a base
#'   R plot or a grid object depending on the `corrplot` version and options.
#'
#' @importFrom stats cor
#' @importFrom corrplot corrplot
#'
#' @export plot_corr
#'
#' @examples
#' data(example_data)
#' # Example of a large cohort with 868 samples and 69 cell types
#' ecoda_object <- ecoda(
#'     data = example_data$GongSharma_full$cell_counts_highresolution,
#'     metadata = example_data$GongSharma_full$metadata
#' )
#' plot_corr(ecoda_object)
plot_corr <- function(ecoda_object,
                      slot = c(
                          "clr", "clr_hvc", "counts", "counts_imp",
                          "freq", "freq_imp", "asin_sqrt", "pb"
                      ),
                      order = "hclust",
                      hclust.method = "ward.D2",
                      ...) {
    slot <- match.arg(slot)
    feat_mat <- slot(ecoda_object, slot)

    cor_matrix <- cor(feat_mat)

    corrplot(cor_matrix, order = "hclust", hclust.method = "ward.D2", ...)
}
