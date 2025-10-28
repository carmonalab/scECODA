# PCA ---------------------------

#' @title Plot Principal Component Analysis and Calculate Clustering Scores
#'
#' @description
#' Performs Principal Component Analysis (PCA) on a selected data matrix from the
#' \code{ECODA} object (default: CLR-transformed abundances, \code{clr}) and
#' visualizes the results in 2D or 3D. It can also calculate and display several
#' metrics to evaluate the separation of groups defined by \code{label_col}.
#'
#' @details
#' The clustering metrics (ARI, Modularity, Silhouette, ANOSIM) assess how well
#' the sample groupings (\code{labels}) align with the underlying data structure
#' in the feature space defined by the PCA.
#'
#' @param ecoda_object An \link[=ECODA-class]{ECODA} object.
#' @param slot Character string (default: \code{"clr"}). The name of the data matrix
#'             slot in the \code{ECODA} object to use for PCA. Must be one of:
#'             \code{"clr"} (CLR-transformed abundances, default), \code{"clr_hvc"}
#'             (CLR-transformed abundances of only the most highly variable cell
#'             types (HVCs)), \code{"pb"} (pseudobulk gene expression),
#'             \code{"counts"} (raw counts), \code{"counts_imp"} (imputed counts),
#'             \code{"freq"} (relative frequencies), \code{"freq_imp"} (imputed
#'             frequencies), or \code{"asin_sqrt"} (arcsin-square root transformed data).
#' @param label_col Character string (optional, default: \code{NULL}). The name of a
#'                  column in \code{ecoda_object@metadata} used to color and group
#'                  samples in the plot, and for calculating clustering scores.
#' @param scale. Logical (default: \code{FALSE}). A value indicating whether the
#'               variables should be scaled to have unit variance before the PCA.
#' @param pca_dims Integer (optional, default: \code{NULL}). The number of principal
#'                 components (dimensions) to retain for the PCA calculation and
#'                 downstream clustering score calculations. If \code{NULL}, all
#'                 dimensions are retained.
#' @param knn_k Integer (optional, default: \code{NULL}). The number of nearest
#'              neighbors (\code{k}) to use for the Shared Nearest Neighbor (SNN)
#'              graph construction, required for Modularity score calculation. If
#'              \code{NULL}, it defaults to \code{max(3, round(sqrt(N)))}, where
#'              \code{N} is the number of samples.
#' @param title Character string (optional, default: \code{NULL}). The main title
#'              for the plot. If clustering scores are calculated, they are appended
#'              to this title.
#' @param legend_title Character string (default: \code{"Group"}). The title for the color
#'                     legend in the plot when a grouping column (\code{label_col}) is provided.
#' @param show_label_samples Logical (default: \code{FALSE}). If \code{TRUE}, sample names are
#'                           displayed next to the points in the plot. This automatically adds
#'                           \code{"text"} to the \code{geom} parameter if it is not already present.
#' @param score_digits Integer (default: \code{3}). The number of decimal places to
#'                     round the clustering and ANOSIM scores appended to the plot title.
#' @param cluster_score Logical (default: \code{TRUE}). If \code{TRUE}, calculates
#'                      the Adjusted Rand Index (ARI) using \code{\link{calc_ari}}.
#' @param mod_score Logical (default: \code{TRUE}). If \code{TRUE}, calculates the
#'                  adjusted Modularity score using \code{\link{calc_modularity}}.
#' @param sil_score Logical (default: \code{FALSE}). If \code{TRUE}, calculates the
#'                  average Silhouette width using \code{\link{calc_sil}}.
#' @param anosim_score Logical (default: \code{TRUE}). If \code{TRUE}, calculates
#'                     the ANOSIM statistic (R) using \code{vegan::anosim}.
#' @param anosim_distance Character string (default: \code{"euclidian"}). The distance
#'                        method used for the ANOSIM calculation (e.g., "euclidean",
#'                        "manhattan").
#' @param anosim_permutations Integer (default: \code{99}). The number of permutations
#'                            to use when calculating the ANOSIM statistic.
#' @param anosim_parallel Integer (default: \code{detectCores()}). The number of parallel
#'                        processes/cores to use for the ANOSIM calculation.
#' @param pointsize Numeric (default: \code{3}). Size of the points in the plot.
#' @param labelsize Numeric (default: \code{4}). Size of the variable labels in the plot.
#' @param coord_equal Logical (default: \code{TRUE}). If \code{TRUE}, forces the
#'                    aspect ratio of the plot to be equal.
#' @param axes Numeric vector (default: \code{c(1, 2)}). The principal components
#'             to plot (e.g., \code{c(1, 2)} for PC1 vs PC2).
#' @param plotly_3d Logical (default: \code{FALSE}). If \code{TRUE}, generates an
#'                  interactive 3D scatter plot using \code{plotly} (requires \code{pca_dims >= 3}).
#' @param invisible Character vector (default: \code{c("var", "quali")}). Elements to
#'                  hide in the 2D plot. Can include "var" (variables/cell types),
#'                  "ind" (samples), or "quali" (group centroids).
#' @param geom Character string or vector (default: \code{"point"}). The geometry to be used
#'             for the plot. Allowed values are combinations of:
#'             \itemize{
#'               \item \code{"point"} to show points for individuals (samples)
#'               \item \code{"text"} to show labels for individuals (samples)
#'               \item \code{"arrow"} to show vectors for variables (features)
#'             }
#'             The default \code{"point"} plots points for both individuals and variables.
#'             Use \code{c("point", "text")} to show both points and labels for samples.
#' @param n_hv_feat_show Integer (default: \code{Inf}). Number of most highly variable features
#'                  to show based on their contribution to the selected axes.
#' @param repel Logical (default: \code{TRUE}). Whether to use \code{ggrepel} to
#'              prevent label overlap for variable names.
#'
#' @return A \code{ggplot} object (2D, via \code{factoextra}) or a \code{plotly}
#'         object (3D) visualizing the PCA results.
#'
#' @importFrom stats prcomp dist hclust cutree
#' @importFrom factoextra fviz_pca
#' @importFrom plotly plot_ly add_markers add_paths group2NA
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggtitle scale_shape_manual coord_equal scale_color_discrete
#' @importFrom parallel detectCores
#'
#' @export plot_pca
plot_pca <- function(ecoda_object,
                     slot = c("clr", "clr_hvc", "counts", "counts_imp", "freq", "freq_imp", "asin_sqrt", "pb"),
                     label_col = NULL,
                     scale. = FALSE,
                     pca_dims = NULL,
                     knn_k = NULL,
                     title = NULL,
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
                     pointsize = 3,
                     labelsize = 4,
                     coord_equal = TRUE,
                     axes = c(1, 2),
                     plotly_3d = FALSE,
                     invisible = c("var", "quali"),
                     geom = "point",
                     n_hv_feat_show = Inf,
                     repel = TRUE) {
  slot <- match.arg(slot)
  feat_mat <- slot(ecoda_object, slot)

  res.pca <- prcomp(feat_mat, scale. = scale., rank. = pca_dims)


  if (!is.null(label_col)) {
    labels <- ecoda_object@metadata[[label_col]]

    if (slot == "clr_hvc") {
      title <- paste0(
        title,
        "\nHVCs: ", ecoda_object@top_n_hvcs,
        " Variance explained: ", ecoda_object@variance_explained
      )
    }
    if (anosim_score) {
      anosim_score <- calc_anosim(
        feat_mat,
        labels,
        distance = anosim_distance,
        permutations = anosim_permutations,
        parallel = anosim_parallel,
        digits = score_digits
      )
      title <- paste0(title, "\nANOSIM score: ", anosim_score)
    }
    if (cluster_score) {
      cluster_score <- calc_ari(feat_mat, labels, digits = score_digits)
      title <- paste0(title, "\nARI: ", cluster_score)
    }
    if (mod_score) {
      mod_score <- calc_modularity(feat_mat, labels, knn_k, digits = score_digits)
      title <- paste0(title, "\nModularity score: ", mod_score)
    }
    if (sil_score) {
      sil_score <- calc_sil(feat_mat, labels, digits = score_digits)
      title <- paste0(title, "\nSilhouette score: ", sil_score)
    }
  } else {
    labels <- "none"
  }

  if (plotly_3d) {
    df <- as.data.frame(res.pca$x)
    df$id <- seq_len(nrow(df))
    df$vs <- factor(labels)
    ms <- replicate(2, df, simplify = F)
    ms[[2]]$PC3 <- min(df$PC3)
    m <- ms %>%
      bind_rows() %>%
      group2NA("id", "vs")
    # Plotting with plotly
    p <- plot_ly(color = ~vs) %>%
      add_markers(data = df, x = ~PC1, y = ~PC2, z = ~PC3) %>%
      add_paths(data = m, x = ~PC1, y = ~PC2, z = ~PC3, opacity = 0.2)
  } else {
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
  }

  return(p)
}


#' Analysis of Similarities (ANOSIM) R score
#'
#' @description
#' Calculates the ANOSIM R-statistic to test whether there is significant separation
#' between two or more groups (defined by \code{labels}) based on the multivariate
#' distances among samples in the feature space (\code{feat_mat}).
#'
#' @details
#' ANOSIM compares the mean of rank dissimilarities between groups to the mean of
#' rank dissimilarities within groups. The R-statistic ranges from \$-1\$ to \$1\$:
#' \itemize{
#'   \item An R value close to **\$1\$** indicates clear separation of groups.
#'   \item An R value close to **\$0\$** indicates that the separation is no greater than
#'         expected by chance (i.e., poor separation).
#'   \item An R value close to **\$-1\$** indicates that within-group dissimilarities are
#'         greater than between-group dissimilarities (a very rare result).
#' }
#'
#' @param feat_mat Numeric matrix or data frame. The feature matrix (e.g., CLR abundances,
#'                 with samples as rows and features as columns).
#' @param labels Vector of factors or character strings. The grouping variable (the
#'               known cluster assignments) for each row in the matrix.
#' @param distance Character string (default: \code{"euclidean"}). The distance method
#'                 to use for calculating dissimilarities. Passed to the \code{distance}
#'                 argument of \code{vegan::anosim}.
#' @param permutations Integer (default: \code{999}). The number of permutations
#'                     to use when calculating the ANOSIM R-statistic and p-value.
#' @param parallel Integer (optional, default: \code{detectCores()}). The number of parallel
#'                 processes/cores to use for the permutation testing.
#' @param digits Integer (default: \code{3}). The number of decimal places to round the R-statistic.
#'
#' @return A numeric value representing the **ANOSIM R-statistic**.
#'
#' @importFrom vegan anosim
#'
#' @export calc_anosim
calc_anosim <- function(feat_mat,
                        labels,
                        distance = "euclidean",
                        permutations = 99,
                        parallel = 1,
                        digits = 3) {
  score <- anosim(
    x = feat_mat,
    grouping = labels,
    distance = distance,
    permutations = permutations,
    parallel = parallel
  )[["statistic"]]

  return(round(score, digits))
}


#' Calculate Adjusted Rand Index (ARI)
#'
#' Calculates the Adjusted Rand Index (ARI) to measure the agreement between the
#' known cluster assignments (\code{labels}) and the cluster assignments derived
#' from two unsupervised clustering methods: Hierarchical Clustering (\code{hclust})
#' and Partitioning Around Medoids (\code{pam}).
#'
#' @param feat_mat Numeric matrix or data frame. The feature matrix (e.g., CLR abundances).
#' @param labels Vector of factors or character strings. The true or known cluster
#'               assignments for each row in the matrix.
#' @param nclusts Integer (optional, default: \code{NULL}). The target number of clusters
#'                (\code{k}) to use for \code{hclust} and \code{pam}. If \code{NULL},
#'                it defaults to the number of unique levels in \code{labels}.
#' @param digits Integer (default: \code{3}). The number of decimal places to round the result.
#' @param return_mean Logical (default: \code{TRUE}). If \code{TRUE}, returns the
#'                    mean of the ARI scores from \code{hclust} and \code{pam}.
#'                    If \code{FALSE}, returns a named list with both individual scores.
#'
#' @return A numeric value (mean ARI) or a list of two ARI scores. ARI ranges from
#'         -1 (disagreement) to +1 (perfect agreement).
#'
#' @importFrom stats dist hclust cutree
#' @importFrom cluster pam
#' @importFrom mclust adjustedRandIndex
#'
#' @export calc_ari
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
  results[["hclust_accuracy"]] <- adjustedRandIndex(as.numeric(as.factor(labels)), clust_labels)

  # Perform PAM clustering
  clust_labels <- pam(feat_mat, k = nclusts)$cluster
  results[["pamclust_accuracy"]] <- adjustedRandIndex(as.numeric(as.factor(labels)), clust_labels)

  if (return_mean) {
    return(round(mean(unlist(results)), digits))
  } else {
    results[["hclust_accuracy"]] <- round(results[["hclust_accuracy"]], digits)
    results[["pamclust_accuracy"]] <- round(results[["pamclust_accuracy"]], digits)
    return(results)
  }
}


#' Calculate Adjusted Modularity Score
#'
#' Calculates the Modularity score for a given clustering (\code{labels}) based on
#' a Shared Nearest Neighbor (SNN) graph constructed from the feature matrix
#' (\code{feat_mat}). The score is adjusted by the theoretical maximum modularity
#' for the number of groups to always be between
#' -0.5 (poor clustering) and +1 (excellent clustering)).
#'
#' @param feat_mat Numeric matrix or data frame. The feature matrix used to compute
#'                 the SNN graph (e.g., CLR abundances).
#' @param labels Vector of factors or character strings. The cluster assignments for
#'               each row in \code{feat_mat}.
#' @param digits Integer (default: \code{3}). The number of decimal places to round the final adjusted score.
#' @param knn_k Integer (optional, default: \code{NULL}). The number of nearest
#'              neighbors (\code{k}) used for SNN graph construction. If \code{NULL},
#'              it defaults to \code{max(3, round(sqrt(N)))}, where \code{N} is
#'              the number of samples.
#'
#' @return A numeric value representing the adjusted modularity score.
#'
#' @importFrom igraph modularity
#'
#' @export calc_modularity
calc_modularity <- function(feat_mat,
                            labels,
                            digits = 3,
                            knn_k = NULL) {
  ngroups <- length(unique(labels))

  if (is.null(knn_k)) {
    knn_k <- max(3, round(sqrt(nrow(feat_mat))))
  }

  # Create a graph object
  g <- compute_snn_graph(feat_mat = feat_mat, knn_k = knn_k)

  # Compute modularity
  modularity_score <- modularity(g, membership = as.numeric(factor(labels)))

  # NOTE:
  # Maximum modularity depends on the number of groups: max(mod) = 1 - 1 / (number of groups)
  # see Brandes, Ulrik, et al. On finding graph clusterings with maximum modularity.

  # Adjust modularity score for number of groups
  maximum_modularity_score <- 1 - (1 / ngroups)
  adjusted_modularity_score <- modularity_score / maximum_modularity_score

  return(round(adjusted_modularity_score, digits))
}



#' Compute K-Nearest Neighbors (KNN)
#'
#' Calculates the indices of the K-nearest neighbors for each sample (row) in a
#' feature matrix using the \code{RANN} package.
#'
#' @param feat_mat Numeric matrix or data frame. The feature matrix (e.g., CLR abundances).
#' @param knn_k Integer. The number of neighbors (\code{k}) to return.
#'
#' @return A matrix where each row corresponds to a sample and contains the indices
#'         of its \code{k} nearest neighbors.
#'
#' @importFrom RANN nn2
compute_KNN <- function(feat_mat, knn_k) {
  # Compute KNN
  knn <- nn2(as.matrix(feat_mat), k = knn_k + 1)$nn.idx
  knn <- knn[, -1] # Remove self-neighbor
  return(knn)
}



#' Compute Shared Nearest Neighbor (SNN) Graph
#'
#' Constructs a graph where nodes are samples and edges are weighted by the number
#' of shared nearest neighbors (SNN) between them. This graph is used for community
#' detection or, in this context, modularity calculation.
#'
#' @param feat_mat Numeric matrix or data frame. The feature matrix (e.g., CLR abundances).
#' @param knn_k Integer. The number of nearest neighbors (\code{k}) used for SNN calculation.
#'
#' @return An \code{igraph} graph object where edge weights represent shared neighbors.
#'
#' @importFrom igraph graph_from_adjacency_matrix
compute_snn_graph <- function(feat_mat,
                              knn_k) {
  knn <- compute_KNN(feat_mat = feat_mat, knn_k = knn_k)

  # Initialize adjacency matrix
  n <- nrow(as.matrix(feat_mat))
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
  g <- graph_from_adjacency_matrix(adj_matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  return(g)
}


#' Calculate Average Silhouette Width
#'
#' Calculates the average silhouette width for a given feature matrix, using pre-defined
#' cluster assignments (\code{labels}). This metric assesses the quality of the clustering.
#'
#' @param feat_mat Numeric matrix or data frame. The feature matrix (e.g., CLR-transformed
#'                 abundances or PCA scores) on which the distance calculation is based.
#' @param labels Vector of factors or character strings. The cluster assignments for
#'               each row in \code{feat_mat} (i.e., the grouping variable).
#' @param digits Integer (default: \code{3}). The number of decimal places to round the final score.
#'
#' @return A numeric value representing the mean silhouette width, typically ranging
#'         from -1 (poor clustering) to +1 (excellent clustering).
#'
#' @importFrom cluster silhouette
#' @importFrom stats dist
#'
#' @export calc_sil
calc_sil <- function(feat_mat,
                     labels,
                     digits = 3) {
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
#' This function takes either the relative abundance (\code{freq}) or CLR-transformed
#' abundance (\code{clr}) matrix from an \link[=ECODA-class]{ECODA} object, converts it
#' from a wide (samples x cell types) to a long (sample, celltype, value) format,
#' and optionally joins it with a specified column from the sample metadata.
#'
#' @param ecoda_object An initialized \link[=ECODA-class]{ECODA} object.
#' @param data_slot Character string specifying the data matrix to use. Must be
#'                  either \code{"freq"} (for relative abundance) or \code{"clr"}
#'                  (for CLR-transformed abundance).
#' @param label_col Character string (optional, default: \code{NULL}). The name of a
#'                  column in the \code{ecoda_object@metadata} slot to merge into
#'                  the long data frame (e.g., "Disease_State" or "Batch"). If \code{NULL},
#'                  only the abundance data and sample/celltype IDs are returned.
#'
#' @return A tidy, long format data frame with columns:
#'         \itemize{
#'           \item \code{sample_id}: Sample identifier.
#'           \item \code{celltype}: Cell type name.
#'           \item \code{rel_abundance} or \code{clr_abundance}: The quantitative value
#'                 depending on the chosen \code{data_slot}.
#'           \item \code{...}: Additional column specified by \code{label_col} (if provided).
#'         }
#'
#' @importFrom methods slot slotNames
#' @importFrom dplyr %>% left_join
#' @importFrom tidyr pivot_longer everything
#' @importFrom rlang sym
#'
#' @seealso \link[=ECODA-class]{ECODA}
#'
#' @examples
#' \dontrun{
#' # Assuming 'ecoda_obj' is a created ECODA object
#'
#' # 1. Create long data with CLR abundance only
#' long_clr <- create_long_data(ecoda_obj, data_slot = "clr")
#'
#' # 2. Create long data with relative abundance and merge the 'Treatment' metadata column
#' long_freq_labeled <- create_long_data(
#'   ecoda_obj,
#'   data_slot = "freq",
#'   label_col = "Treatment"
#' )
#' }
create_long_data <- function(ecoda_object,
                             data_slot,
                             label_col = NULL) {
  # Ensure the data_slot is valid (freq or clr)
  if (!(data_slot %in% c("freq", "clr"))) {
    stop("Invalid data_slot. Must be 'freq' (Relative Abundance) or 'clr' (CLR Abundance).")
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
    if (!("metadata" %in% slotNames(ecoda_object)) || is.null(ecoda_object@metadata)) {
      stop("label_col was provided but ecoda_object@metadata slot is missing or NULL.")
    }

    # --- Check 2: Ensure label_col exists in metadata ---
    meta_colnames <- colnames(ecoda_object@metadata)
    if (!(label_col %in% meta_colnames)) {
      stop(paste0("label_col '", label_col, "' not found in ecoda_object@metadata."))
    }

    metadata_df <- as.data.frame(ecoda_object@metadata)
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
#' @param ecoda_object An \link[=ECODA-class]{ECODA} object containing cell type relative
#'                     frequencies in the \code{freq} slot.
#' @param label_col Character string (optional, default: \code{NULL}). The name of a
#'                  column in \code{ecoda_object@metadata} used to define grouping
#'                  or groups (required if \code{plot_by = "group"}).
#' @param plot_by Character string (default: \code{"sample"}). Specifies whether
#'                to plot the relative abundance for each individual sample (\code{"sample"})
#'                or the average relative abundance aggregated by a group (\code{"group"})
#'                defined by \code{label_col}.
#' @param custom_sample_order Character vector (optional, default: \code{NULL}).
#'                            A vector of sample IDs to enforce a specific order
#'                            when \code{plot_by = "sample"}. If \code{NULL}, samples
#'                            are ordered first by \code{label_col} (if provided)
#'                            and then naturally.
#' @param title Character string (default: \code{""}). The main title for the plot.
#' @param facet_by_label_col Logical (default: \code{TRUE}). If \code{TRUE} and
#'                           \code{plot_by = "sample"}, the plot will be faceted
#'                           (split) horizontally by the categories in \code{label_col}.
#'
#' @return A \code{ggplot} object representing the stacked bar plot.
#'
#' @importFrom dplyr %>% group_by summarise mutate distinct arrange pull
#' @importFrom ggplot2 ggplot aes geom_col theme_minimal theme element_text labs facet_grid
#' @importFrom rlang sym
#' @importFrom gtools mixedsort
#' @importFrom stats reformulate
#'
#' @export plot_barplot
#'
#' @seealso \code{\link{create_long_data}}, \link[=ECODA-class]{ECODA}
#'
#' @examples
#' \dontrun{
#' # Assuming 'ecoda_obj' is a created ECODA object with metadata
#'
#' # 1. Plot frequency for every sample, ordered by sample ID:
#' p1 <- plot_barplot(ecoda_obj)
#'
#' # 2. Plot frequency for every sample, faceted and ordered by 'Treatment' column:
#' p2 <- plot_barplot(ecoda_obj, label_col = "Treatment", plot_by = "sample")
#'
#' # 3. Plot average frequency aggregated by 'Treatment' group:
#' p3 <- plot_barplot(ecoda_obj, label_col = "Treatment", plot_by = "group")
#'
#' # 4. Plot with a custom order:
#' custom_order <- c("S2", "S1", "S4", "S3")
#' p4 <- plot_barplot(ecoda_obj, custom_sample_order = custom_order)
#' }
plot_barplot <- function(ecoda_object,
                         label_col = NULL,
                         plot_by = c("sample", "group"),
                         custom_sample_order = NULL,
                         title = "",
                         facet_by_label_col = TRUE) {
  # Use the helper function to get the long data from @freq
  plot_data <- create_long_data(ecoda_object, data_slot = "freq", label_col = label_col)

  plot_by <- match.arg(plot_by)

  if (plot_by == "group") {
    if (is.null(label_col)) {
      stop("label_col must be provided when plot_by = 'group'")
    }

    # Plotting by group
    plot_df <- plot_data %>%
      group_by(celltype, !!sym(label_col)) %>%
      summarise(mean_rel_abund = mean(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
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
        stop("Custom sample order is incomplete or contains unknown sample IDs.")
      }
      ordered_levels <- custom_sample_order
    } else {
      # DEFAULT: Sort samples by label_col (group) if provided, then naturally

      if (!is.null(label_col)) {
        # 1. Get a unique list of samples and their corresponding label_col values
        sample_group_map <- plot_data %>%
          distinct(sample_id, !!sym(label_col))

        # 2. Order the samples by the group column first (using mixedsort on its values)
        #    and then by sample_id (using mixedsort on its values).
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


#' Generates Boxplots for CLR-transformed Cell Type Abundances with Optional Group Comparison.
#'
#' This function visualizes the distribution of CLR-transformed abundance for each
#' cell type using boxplots, optionally splitting the data by a sample metadata
#' column and performing statistical tests for comparison between groups.
#'
#' \strong{Statistical Test Logic:}
#' \itemize{
#'   \item If the number of groups (\code{label_col} levels) is \strong{2}, the function
#'         uses the specified \code{stat_method} (default: "wilcox.test") and displays
#'         pairwise significance labels or stars.
#'   \item If the number of groups is \strong{3 or more}, and \code{stat_method} is
#'         "wilcox.test" or "t.test", the function automatically switches to the
#'         global non-parametric test, \strong{"kruskal.test"}, and displays the
#'         overall p-value for the comparison across all groups.
#' }
#'
#' @param ecoda_object An \link[=ECODA-class]{ECODA} object containing the CLR-transformed
#'                     abundances in the \code{clr} slot.
#' @param label_col Character string (optional, default: \code{NULL}). The name of a
#'                  column in \code{ecoda_object@metadata} used to define groups for
#'                  comparison. If \code{NULL}, a single boxplot is generated per cell type.
#' @param title Character string (default: \code{""}). The main title for the plot.
#' @param stat_method Character string (default: \code{"wilcox.test"}). The statistical
#'                    method used for comparisons between 2 groups (e.g., "t.test",
#'                    "wilcox.test"). Note: This is overridden by "kruskal.test" for 3+ groups.
#' @param paired Logical (default: \code{FALSE}). If \code{TRUE}, performs a paired
#'               statistical test (e.g., paired t-test or paired Wilcoxon test). Only
#'               applicable for 2-group comparisons.
#' @param signif_label Character string (default: \code{"p.signif"}). Controls
#'                     how p-values are displayed for 2-group comparisons (e.g.,
#'                     "p.signif" for stars, "p.format" for numeric p-value).
#'
#' @return A \code{ggplot} object (enhanced by \code{ggpubr}) representing the
#'         CLR abundance boxplot, including jittered data points and dynamic
#'         significance information (pairwise stars for 2 groups, overall p-value
#'         for 3+ groups).
#'
#' @importFrom ggplot2 aes geom_jitter labs theme element_text guides position_jitterdodge
#' @importFrom ggpubr ggboxplot stat_compare_means
#' @importFrom stringr str_to_title
#' @importFrom rlang sym
#'
#' @export plot_boxplot
#'
#' @seealso \code{\link{create_long_data}}, \link[=ECODA-class]{ECODA}
#'
#' @examples
#' \dontrun{
#' # Assuming 'ecoda_obj' is a created ECODA object with metadata
#'
#' # 1. Boxplots for CLR abundance without grouping (no stats calculated):
#' p1 <- plot_boxplot(ecoda_obj)
#'
#' # 2. Boxplots grouped by 'Treatment' (2 groups) and applying Wilcoxon test:
#' p2 <- plot_boxplot(
#'   ecoda_obj,
#'   label_col = "Treatment",
#'   stat_method = "wilcox.test",
#'   title = "CLR Abundance by Treatment Group"
#' )
#'
#' # 3. Boxplots grouped by 'Dose' (3+ groups) - automatically uses Kruskal-Wallis:
#' p3 <- plot_boxplot(
#'   ecoda_obj,
#'   label_col = "Dose",
#'   stat_method = "wilcox.test", # will be overridden by kruskal.test
#'   title = "CLR Abundance by Dose Group (Kruskal-Wallis)"
#' )
#' }
plot_boxplot <- function(ecoda_object,
                         label_col = NULL,
                         title = "",
                         stat_method = "wilcox.test",
                         paired = FALSE,
                         signif_label = c("p.signif", "p.format")) {
  signif_label <- match.arg(signif_label)
  plot_data <- create_long_data(ecoda_object, data_slot = "clr", label_col = label_col)

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
    # SCENARIO 1: No grouping variable provided (Single boxplot per cell type)

    # Add simple jitter (no dodging, color determined by the base 'black' mapping)
    p <- p + geom_jitter(
      width = 0.2,
      size = 1,
      alpha = 0.6
    ) +
      # SUPPRESS THE LEGEND
      guides(color = "none")
  } else {
    # SCENARIO 2: Grouping variable is provided (Comparison between groups)

    # Calculate the number of boxplots for correct jitter-dodging
    nr_of_boxplots <- length(unique(plot_data[[label_col]]))
    label_col_sym <- sym(label_col)

    # --- DYNAMIC STATISTICAL TEST SELECTION ---
    current_stat_method <- stat_method

    if (nr_of_boxplots > 2) {
      # For 3+ non-parametric groups, use Kruskal-Wallis (overall test)
      # We only override if the user used a 2-group method like wilcox.test or t.test
      if (stat_method %in% c("wilcox.test", "t.test")) {
        current_stat_method <- "kruskal.test"
        warning(paste0(
          "stat_method automatically switched from '", stat_method,
          "' to 'kruskal.test' for multi-group (N=", nr_of_boxplots, ") comparison."
        ))
      }
    }
    # ---------------------------------------------

    # Add jittered points with dodging
    p <- p + geom_jitter(
      mapping = aes(color = !!label_col_sym), # Map color to group
      position = position_jitterdodge(jitter.width = 1 / nr_of_boxplots),
      size = 1,
      alpha = 0.6
    )

    # Add significance testing: logic changes based on Kruskal-Wallis vs. Pairwise test
    if (current_stat_method == "kruskal.test") {
      # Kruskal-Wallis: Overall test (no 'group' aesthetic needed)
      p <- p + stat_compare_means(
        method = current_stat_method,
        label.y.npc = "top", # Place label at the top
        label.x.npc = "center",
        label = "p.format" # Show the overall p-value
      )
    } else {
      # Wilcoxon or t.test: Pairwise test (requires 'group' aesthetic)
      p <- p + stat_compare_means(
        aes(group = !!label_col_sym),
        method = stat_method,
        paired = paired,
        label = signif_label, # Show significance stars
        tip.length = 0,
        hide.ns = TRUE
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
#' This function visualizes a data matrix from a specified slot (e.g., CLR-transformed,
#' frequency, or pseudobulk data) after mean-centering. It supports optional
#' filtering to only Highly Variable Cell Types (HVCs) and includes a sample
#' annotation sidebar based on a specified metadata column.
#'
#' @param ecoda_object An \link[=ECODA-class]{ECODA} object.
#' @param slot Character string (default: \code{"clr"}). The name of the data matrix
#'             slot in the \code{ECODA} object to use for the heatmap. Must be one of:
#'             \code{"clr"} (CLR-transformed abundances, default), \code{"clr_hvc"}
#'             (CLR-transformed abundances of only the most highly variable cell
#'             types (HVCs)), \code{"pb"} (pseudobulk gene expression),
#'             \code{"counts"} (raw counts), \code{"counts_imp"} (imputed counts),
#'             \code{"freq"} (relative frequencies), \code{"freq_imp"} (imputed
#'             frequencies), or \code{"asin_sqrt"} (arcsin-square root transformed data).
#' @param label_col Character string. The name of the column in \code{ecoda_object@metadata}
#'                  to use for annotating the samples (columns) of the heatmap.
#' @param cluster_rows Logical (default: \code{TRUE}). Whether to apply hierarchical
#'                     clustering to the cell types (rows).
#' @param cluster_cols Logical (default: \code{TRUE}). Whether to apply hierarchical
#'                     clustering to the samples (columns).
#' @param scale Character string (default: \code{"none"}). Method for scaling the
#'              abundance values within the heatmap using \code{pheatmap}. Options
#'              include \code{"none"}, \code{"row"}, or \code{"column"}. Note: The data
#'              is internally **mean-centered** (\code{scale(center=TRUE, scale=FALSE)})
#'              across samples before being passed to \code{pheatmap}, regardless of this argument.
#' @param clustering_method Character string (default: \code{"ward.D2"}). The clustering
#'                        method to use for hierarchical clustering. Options are
#'                        passed directly to \code{hclust} (e.g., \code{"complete"},
#'                        \code{"average"}, \code{"ward.D2"}).
#' @param angle_col Character string (default: \code{"90"}). Angle of the sample
#'                  labels (columns).
#' @param ... Additional arguments passed directly to the \code{pheatmap} function.
#'
#' @return A \code{pheatmap} object, which is typically visualized automatically
#'         when called interactively.
#'
#' @importFrom dplyr %>%
#' @importFrom pheatmap pheatmap
#'
#' @export plot_heatmap
#'
#' @seealso \link[=ECODA-class]{ECODA}, \code{\link{pheatmap}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'ecoda_obj' is a created ECODA object with metadata
#'
#' # 1. Heatmap using CLR data, clustered, and annotated by 'Condition':
#' p1 <- plot_heatmap(
#'   ecoda_obj,
#'   label_col = "Condition",
#'   main = "CLR Abundance Heatmap",
#'   fontsize = 8
#' )
#'
#' # 2. Heatmap using Relative Frequency data (clr_hvc), filtered to only HVCs,
#' #    and without clustering the samples:
#' p2 <- plot_heatmap(
#'   ecoda_obj,
#'   slot = "clr_hvc",
#'   label_col = "Batch",
#'   cluster_cols = FALSE
#' )
#' }
plot_heatmap <- function(ecoda_object,
                         slot = c("clr", "clr_hvc", "counts", "counts_imp", "freq", "freq_imp", "asin_sqrt", "pb"),
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

  metadata <- ecoda_object@metadata[, label_col, drop = FALSE]
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
#' @description
#' Calculates the pairwise Pearson correlation matrix for all cell types (columns)
#' using the Centered Log-Ratio (CLR) transformed abundances stored in
#' \code{ecoda_object@clr}. It then visualizes this matrix as a heatmap using
#' \code{corrplot::corrplot}.
#'
#' @details
#' The function uses the CLR matrix, where high correlation between two cell types
#' suggests they vary together across samples, indicating potential co-occurrence
#' or co-regulation.
#'
#' @param ecoda_object An \link[=ECODA-class]{ECODA} object.
#' @param slot Character string (default: \code{"clr"}). The name of the data matrix
#'             slot in the \code{ECODA} object to use for the correlation plot Must be one of:
#'             \code{"clr"} (CLR-transformed abundances, default), \code{"clr_hvc"}
#'             (CLR-transformed abundances of only the most highly variable cell
#'             types (HVCs)), \code{"pb"} (pseudobulk gene expression),
#'             \code{"counts"} (raw counts), \code{"counts_imp"} (imputed counts),
#'             \code{"freq"} (relative frequencies), \code{"freq_imp"} (imputed
#'             frequencies), or \code{"asin_sqrt"} (arcsin-square root transformed data).
#' @param order Character string (default: \code{"hclust"}). The ordering method
#'              for the correlation matrix. Common options include:
#'              \itemize{
#'                \item \code{"original"} (no reordering)
#'                \item \code{"hclust"} (hierarchical clustering)
#'                \item \code{"FPC"} (first principal component order)
#'              }
#' @param hclust.method Character string (default: \code{"ward.D2"}). The
#'                      hierarchical clustering method to use if \code{order} is
#'                      set to \code{"hclust"}.
#' @param ... Additional arguments passed to \code{corrplot::corrplot} for plot
#'            customization (e.g., \code{method}, \code{type}, \code{tl.col}).
#'
#' @return A plot object generated by \code{corrplot::corrplot}, which is a
#'         base R plot or a grid object depending on the `corrplot` version and options.
#'
#' @importFrom stats cor
#' @importFrom corrplot corrplot
#'
#' @export plot_corr
plot_corr <- function(ecoda_object,
                      slot = c("clr", "clr_hvc", "counts", "counts_imp", "freq", "freq_imp", "asin_sqrt", "pb"),
                      order = "hclust",
                      hclust.method = "ward.D2",
                      ...) {
  slot <- match.arg(slot)
  feat_mat <- slot(ecoda_object, slot)

  cor_matrix <- cor(feat_mat)

  corrplot(cor_matrix, order = "hclust", hclust.method = "ward.D2", ...)
}
