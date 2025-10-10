library(dplyr)
library(ggplot2)
library(gtools)
library(pheatmap)
library(rlang) # for !!sym() and as_label()
library(tidyr)


#' Creates an ECODA object from various data types.
#'
#' This is a smart constructor that can accept a Seurat object,
#' a SingleCellExperiment object, a long format data frame of cell-level data,
#' or a pre-calculated cell count dataframe
#'
#' @param data The primary input, which can be a Seurat object,
#'             a SingleCellExperiment object, a long format data frame of cell-level data,
#'             or a cell count matrix/data frame.
#' @param sample_col The column name for sample IDs (if `data` is a data frame).
#' @param celltype_col The column name for cell type annotations (if `data` is a data frame).
#' @param metadata An optional data frame of pre-calculated sample-level metadata.
#' @param long_df An optional data frame of pre-calculated sample-level metadata.
#'
#' @return A new ECODA object.
#' @export create_ecoda_object
create_ecoda_object <- function(data = NULL,
                                sample_col = NULL,
                                celltype_col = NULL,
                                get_pb = FALSE,
                                variance_explained = 0.5,
                                top_n_hvcs = NULL) {
  # --- A) Handle input from a Seurat or SingleCellExperiment object ---
  if (inherits(data, "Seurat") | inherits(data, "SingleCellExperiment")) {
    if (is.null(sample_col) || is.null(celltype_col)) {
      stop("Please provide 'sample_col' and 'celltype_col'.")
    }

    # Get cell metadata from the object
    if (inherits(data, "Seurat")) {
      cell_data_df <- as.data.frame(data@meta.data)
      if (get_pb) {
        pb <- calculate_pseudobulk(
          count_matrix = data[["RNA"]]$counts,
          sample_ids = cell_data_df[[sample_col]]
        )
      }
    } else if (inherits(data, "SingleCellExperiment")) {
      cell_data_df <- as.data.frame(data@colData)
      names(cell_data_df) <- names(data@colData)
      if (get_pb) {
        pb <- calculate_pseudobulk(
          count_matrix = assay(data, "counts"),
          sample_ids = cell_data_df[[sample_col]]
        )
      }
    }

    metadata <- get_sample_metadata(cell_data_df, sample_col)

    if (get_pb) {
      pb <- deseq2_normalize(pb, metadata)
    }

    counts <- get_celltype_counts(cell_data_df, sample_col, celltype_col)
  }

  ecoda_object <- create_ecoda_object_from_counts(
    counts,
    metadata,
    variance_explained,
    top_n_hvcs
  )

  ecoda_object@pb <- pb

  return(ecoda_object)
}


create_ecoda_object_from_counts <- function(counts = NULL,
                                            metadata = NULL,
                                            variance_explained = 0.5,
                                            top_n_hvcs = NULL) {
  # Initialize the object with default values
  ecoda_object <- methods::new("ECODA")

  if (!is.null(counts)) {
    if (!is.null(metadata)) {
      # Sort by rownames
      counts <- counts[gtools::mixedsort(rownames(counts)), ]
      metadata <- metadata[gtools::mixedsort(rownames(metadata)), , drop = FALSE]

      # Ensure correct dimensions
      if (!all(rownames(counts) == rownames(metadata))) {
        stop("Rownames of cell counts do not match metadata.")
      }
      ecoda_object@metadata <- metadata
    }

    counts_imp <- counts
    if (any(counts == 0)) {
      counts_imp <- counts_imp + 1
    }
    freq <- calc_freq(counts)
    freq_imp <- calc_freq(counts_imp)
    clr_df <- clr(counts_imp)
    sdist <- clr_df %>%
      dist() %>%
      as.matrix() %>%
      as.data.frame()

    ecoda_object@counts <- counts
    ecoda_object@counts_imp <- counts_imp
    ecoda_object@freq <- freq
    ecoda_object@freq_imp <- freq_imp
    ecoda_object@clr <- clr_df
    ecoda_object@sample_distances <- sdist

    ecoda_object <- find_highly_variable_celltypes(
      ecoda_object,
      variance_explained,
      top_n_hvcs
    )
  }

  return(ecoda_object)
}


# Find highly variable cell types ####

find_highly_variable_celltypes <- function(ecoda_object,
                                           variance_explained = 0.5,
                                           top_n_hvcs = NULL) {
  df_var <- get_celltype_variances(
    ecoda_object,
    show_plot = FALSE
  )

  hvcs <- get_highly_variable_celltypes(
    df_var,
    variance_explained
  )

  variance_explained <- (sum(df_var$Variance[1:length(hvcs)]) /
    sum(df_var$Variance))

  ecoda_object@celltype_variances <- df_var
  ecoda_object@variance_explained <- variance_explained
  ecoda_object@top_n_hvcs <- length(hvcs)
  ecoda_object@highly_variable_celltypes <- hvcs


  return(ecoda_object)
}




get_celltype_variances <- function(ecoda_object,
                                   show_plot = TRUE,
                                   label_points = TRUE,
                                   plot_title = "",
                                   smooth_method = "lm",
                                   descending = TRUE) {
  df_var <- ecoda_object@clr %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "celltype",
      values_to = "values"
    ) %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarize(
      avg_clr_abundance = mean(values, na.rm = TRUE),
      Variance = var(values, na.rm = TRUE)
    )

  if (descending) {
    df_var <- df_var %>%
      dplyr::arrange(dplyr::desc(Variance))
  } else {
    df_var <- df_var %>%
      dplyr::arrange(Variance)
  }

  total_variance <- sum(df_var$Variance)

  df_var <- df_var %>%
    dplyr::mutate(
      cumulative_variance = cumsum(Variance),
      variance_exp = (cumulative_variance / total_variance)
    )

  if (show_plot) {
    p <- varmeanplot(
      data = df_var,
      plot_title = plot_title,
      smooth_method = smooth_method,
      label_points = label_points
    )
    print(p)
  }

  return(as.data.frame(df_var))
}


get_highly_variable_celltypes <- function(df_var,
                                          variance_explained = 0.5,
                                          top_n_hvcs = NULL) {
  if (!is.null(top_n_hvcs)) {
    warning("Setting top_n_hvcs overrules variance_explained parameter")
    if (top_n_hvcs < 1) {
      top_hvcs <- df_var %>%
        slice_head(prop = top_n_hvcs) %>%
        pull(celltype)
    } else {
      top_hvcs <- df_var %>%
        slice_head(n = top_n_hvcs) %>%
        pull(celltype)
    }
  } else {
    last_hvc_index <- which.max(df_var$variance_exp >= variance_explained)

    if (length(last_hvc_index) == 0 || last_hvc_index == 0) {
      warning("Variance threshold was not met for any cell type (or data is empty). Returning top 2.")
      return(df_var$celltype[1:min(2, nrow(df_var))])
    }

    top_hvcs <- df_var %>%
      dplyr::slice(1:last_hvc_index) %>%
      dplyr::pull(celltype)
  }

  # Select at least two cell types
  if (length(top_hvcs) < 2) {
    warning("The minimum number of highly variable cell types to use is 2")
    top_hvcs <- df_var %>%
      slice_head(n = 2) %>%
      pull(celltype)
  }

  return(top_hvcs)
}



plot_varmean <- function(ecoda_object,
                         plot_title = "",
                         smooth_method = "lm",
                         label_points = TRUE,
                         highlight_cts = NULL) {
  df_var <- ecoda_object@celltype_variances

  p <- ggplot(df_var, aes(x = Relative_abundance, y = Variance)) +
    geom_point() +
    geom_smooth(method = smooth_method, color = "red", fill = "#69b3a2", se = TRUE) +
    labs(title = paste(plot_title)) +
    theme_classic() +
    xlab("Mean") +
    ylab("Variance")

  if (label_points) {
    p <- p + ggrepel::geom_text_repel(data = df_var, aes(label = celltype), vjust = -0.5)
  }

  return(p)
}




# Plotting functions ####

## Box and bar plots ####

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
  data_df <- as.data.frame(data_matrix) %>%
    tibble::rownames_to_column("sample_id")

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

    metadata_df <- as.data.frame(ecoda_object@metadata) %>%
      tibble::rownames_to_column("sample_id") %>%
      dplyr::select(sample_id, !!sym(label_col))

    # 4. Join the long data with the metadata by sample_id
    plot_data <- long_data %>%
      left_join(metadata_df, by = "sample_id")
  } else {
    # If no label_col, just return the long data without joining metadata
    plot_data <- long_data
  }

  return(plot_data)
}




plot_freq_barplot <- function(ecoda_object,
                              label_col = NULL,
                              plot_by = c("sample", "condition"),
                              custom_sample_order = NULL,
                              title = "",
                              facet_by_label_col = TRUE) {
  # Use the helper function to get the long data from @freq
  plot_data <- create_long_data(ecoda_object, data_slot = "freq", label_col = label_col)

  plot_by <- match.arg(plot_by, c("sample", "condition"))

  if (plot_by == "condition") {
    if (is.null(label_col)) {
      stop("label_col must be provided when plot_by = 'condition'")
    }

    # Plotting by Condition
    plot_df <- plot_data %>%
      group_by(celltype, !!sym(label_col)) %>%
      summarise(mean_rel_abund = mean(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
      mutate(x_var = !!sym(label_col), y_var = mean_rel_abund)

    # Ensure condition is a factor
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
      # DEFAULT: Sort samples by label_col (condition) if provided, then naturally

      if (!is.null(label_col)) {
        # 1. Get a unique list of samples and their corresponding label_col values
        sample_condition_map <- plot_data %>%
          distinct(sample_id, !!sym(label_col))

        # 2. Order the samples by the condition column first (using mixedsort on its values)
        #    and then by sample_id (using mixedsort on its values).
        ordered_levels <- sample_condition_map %>%
          mutate(
            ordered_condition = factor(
              !!sym(label_col),
              levels = gtools::mixedsort(unique(!!sym(label_col))),
              ordered = TRUE
            )
          ) %>%
          arrange(ordered_condition, gtools::mixedsort(sample_id)) %>%
          pull(sample_id)
      } else {
        # If no label_col, just sort samples naturally by sample_id
        ordered_levels <- gtools::mixedsort(unique(plot_data$sample_id))
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





plot_clr_boxplot <- function(ecoda_object,
                             label_col = NULL,
                             title = "",
                             stat_method = "wilcox.test") {
  plot_data <- create_long_data(ecoda_object, data_slot = "clr", label_col = label_col)

  # Ensure celltype is a factor for plotting
  plot_data$celltype <- factor(plot_data$celltype)

  # --- Setup for Plotting ---

  # Determine color mapping for ggboxplot
  box_color_map <- if (!is.null(label_col)) label_col else "black"

  # Generate the base boxplot
  p <- ggpubr::ggboxplot(plot_data,
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
    p <- p + ggplot2::geom_jitter(
      width = 0.2,
      size = 1,
      alpha = 0.6
    ) +
      # SUPPRESS THE LEGEND
      ggplot2::guides(color = "none")
  } else {
    # SCENARIO 2: Grouping variable is provided (Comparison between groups)

    # Calculate the number of boxplots for correct jitter-dodging
    nr_of_boxplots <- length(unique(plot_data[[label_col]]))
    label_col_sym <- sym(label_col)

    # Add jittered points with dodging
    p <- p + ggplot2::geom_jitter(
      mapping = ggplot2::aes(color = !!label_col_sym), # Map color to group
      position = ggplot2::position_jitterdodge(jitter.width = 1 / nr_of_boxplots),
      size = 1,
      alpha = 0.6
    ) +
      # Add significance testing
      ggpubr::stat_compare_means(
        ggplot2::aes(group = !!label_col_sym),
        method = stat_method,
        label = signif_label, # Show significance stars
        tip.length = 0,
        hide.ns = TRUE
      )

    # Add legend title
    p <- p + ggplot2::labs(color = stringr::str_to_title(label_col))
  }

  # --- Final Theme and Labels ---
  p <- p +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, ),
      legend.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(title = title)

  return(p)
}



## Heatmap ####

plot_heatmap <- function(ecoda_object,
                         label_col,
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         scale = "none",
                         clustering_method = "ward.D2",
                         angle_col = "90",
                         ...) {
  df_heatmap <- ecoda_object@clr %>%
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








## PCA ####

# Plot 2D and 3D PCA from feature matrix and calculate silhouette and modularity score
plot_pca <- function(ecoda_object,
                     label_col = NULL,
                     scale. = FALSE,
                     pca_dims = NULL,
                     knn_k = NULL,
                     title = NULL,
                     cluster_score = TRUE,
                     mod_score = TRUE,
                     sil_score = FALSE,
                     anosim_score = TRUE,
                     pointsize = 3,
                     labelsize = 4,
                     coord_equal = TRUE,
                     axes = c(1, 2),
                     plotly_3d = FALSE,
                     invisible = c("var", "quali"),
                     n_ct_show = Inf,
                     repel = TRUE) {
  feat_mat <- ecoda_object@clr

  res.pca <- prcomp(feat_mat, scale. = scale., rank. = pca_dims)


  if (!is.null(label_col)) {
    labels <- ecoda_object@metadata[[label_col]]

    if (anosim_score) {
      anosim_score <- round(vegan::anosim(x = feat_mat, grouping = labels, distance = "euclidean")[["statistic"]], 3)
      title <- paste0(title, "\nANOSIM score: ", anosim_score)
    }
    if (cluster_score) {
      cluster_score <- calc_ari(feat_mat, labels)
      title <- paste0(title, "\nARI: ", cluster_score)
    }
    if (mod_score) {
      mod_score <- round(calc_modularity(feat_mat, labels, knn_k), 3)
      title <- paste0(title, "\nModularity score: ", mod_score)
    }
    if (sil_score) {
      sil_score <- round(calc_sil(feat_mat, labels), 3)
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
      plotly::group2NA("id", "vs")
    # Plotting with plotly
    p <- plotly::plot_ly(color = ~vs) %>%
      plotly::add_markers(data = df, x = ~PC1, y = ~PC2, z = ~PC3) %>%
      plotly::add_paths(data = m, x = ~PC1, y = ~PC2, z = ~PC3, opacity = 0.2)
  } else {
    if (!is.infinite(n_ct_show) & all(invisible %in% c("var", "quali"))) {
      invisible <- "quali"
    }

    p <- factoextra::fviz_pca(res.pca,
      axes = axes,
      habillage = labels,
      label = "var",
      pointsize = pointsize,
      labelsize = labelsize,
      invisible = invisible,
      select.var = list(contrib = n_ct_show),
      repel = repel,
      geom = "point"
    ) +
      ggtitle(title)

    if (!is.null(label_col)) {
      p <- p + scale_shape_manual(values = rep(19, length(unique(labels))))
    }

    if (coord_equal) {
      p <- p + coord_equal()
    }
  }

  return(p)
}




# Calculate average silhouette width
calc_sil <- function(feat_mat,
                     labels) {
  sils <- cluster::silhouette(
    x = as.numeric(factor(labels)),
    dist = dist(feat_mat)
  ) %>%
    as.data.frame()
  score <- mean(sils[["sil_width"]])
  return(score)
}


calc_modularity <- function(feat_mat,
                            labels,
                            knn_k = NULL) {
  ngroups <- length(unique(labels))

  if (is.null(knn_k)) {
    knn_k <- max(3, round(sqrt(nrow(feat_mat))))
  }

  # Create a graph object
  g <- compute_snn_graph(feat_mat = feat_mat, knn_k = knn_k)

  # Compute modularity
  modularity_score <- igraph::modularity(g, membership = as.numeric(factor(labels)))

  # NOTE:
  # Maximum modularity depends on the number of groups: max(mod) = 1 - 1 / (number of groups)
  # see Brandes, Ulrik, et al. On finding graph clusterings with maximum modularity.

  # Adjust modularity score for number of groups
  maximum_modularity_score <- 1 - (1 / ngroups)
  adjusted_modularity_score <- modularity_score / maximum_modularity_score

  return(adjusted_modularity_score)
}


compute_KNN <- function(feat_mat, knn_k) {
  # Compute KNN
  knn <- RANN::nn2(as.matrix(feat_mat), k = knn_k + 1)$nn.idx
  knn <- knn[, -1] # Remove self-neighbor
  return(knn)
}


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
  g <- igraph::graph_from_adjacency_matrix(adj_matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  return(g)
}




# Calculate Adjusted Rand Index
calc_ari <- function(matrix,
                     labels,
                     nclusts = NULL,
                     digits = 3,
                     return_mean = TRUE) {
  results <- list()
  dist_mat <- dist(matrix)

  if (is.null(nclusts)) {
    nclusts <- length(unique(labels))
  }

  # Perform hierarchical clustering
  hc <- stats::hclust(dist_mat, method = "ward.D2")
  clust_labels <- stats::cutree(hc, k = nclusts)
  results[["hclust_accuracy"]] <- mclust::adjustedRandIndex(as.numeric(as.factor(labels)), clust_labels)

  # Perform PAM clustering
  clust_labels <- cluster::pam(matrix, k = nclusts)$cluster
  results[["pamclust_accuracy"]] <- mclust::adjustedRandIndex(as.numeric(as.factor(labels)), clust_labels)

  if (return_mean) {
    return(round(mean(unlist(results)), digits))
  } else {
    results[["hclust_accuracy"]] <- round(results[["hclust_accuracy"]], digits)
    results[["pamclust_accuracy"]] <- round(results[["pamclust_accuracy"]], digits)
    return(results)
  }
}
