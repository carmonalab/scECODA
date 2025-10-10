# Define ECODA object to store data
methods::setClass(
  Class = "ECODA",
  slots = list(
    counts = "data.frame",
    counts_imp = "data.frame",
    freq = "data.frame",
    freq_imp = "data.frame",
    clr = "data.frame",
    celltype_variances = "data.frame",
    variance_explained = "numeric",
    top_n_hvcs = "integer",
    highly_variable_celltypes = "character",
    metadata = "data.frame",
    pb = "data.frame",
    sample_distances = "data.frame"
  ),
  # Define the default values for each slot.
  prototype = list(
    counts = NULL,
    counts_imp = NULL,
    freq = NULL,
    freq_imp = NULL,
    clr = NULL,
    celltype_variances = NULL,
    variance_explained = NULL,
    top_n_hvcs = NULL,
    highly_variable_celltypes = NULL,
    metadata = NULL,
    pb = NULL,
    sample_distances = NULL
  )
)


calc_freq <- function(df) {
  df <- t(apply(df, 1, function(row) (row / sum(row)) * 100)) %>%
    as.data.frame()
  return(df)
}



clr <- function(df) {
  geometric_mean <- apply(df, 1, function(row) exp(mean(log(row))))
  clr_df <- apply(df, 2, function(row) log(row) - log(geometric_mean)) %>%
    as.data.frame()

  return(clr_df)
}




#' Get the cell type counts from a long data frame (e.g. seurat object metadata) where each cell is a row.
#'
#' @param sample_col The column that defines the sample ID for each cell
#' @param celltype_col The column that defines the cell type annotation for each cell
#'
#' @return A data frame with samples as rows and cell types as columns,
#'         containing the count of each cell type per sample.
#' @export get_celltype_counts
get_celltype_counts <- function(cell_data_df,
                                sample_col,
                                celltype_col) {
  cellcount_df <- table(cell_data_df[[sample_col]], cell_data_df[[celltype_col]], useNA = "ifany") %>%
    as.data.frame.matrix()

  # Check if any name is NA (the special missing value)
  dim_names <- colnames(cellcount_df)
  na_index <- is.na(dim_names)

  # Replace the NA entry with the character string "NA"
  dim_names[na_index] <- "NA"

  # Assign the fixed names back to the table
  colnames(cellcount_df) <- dim_names

  return(cellcount_df)
}



#' Extracts constant metadata for each sample from a cell-level data frame.
#'
#' This function identifies and returns metadata columns that have the same value
#' for all cells within a given sample.
#'
#' @param cell_data_df A data frame containing cell-level metadata.
#' @param sample_col The column that defines the sample ID for each cell.
#' @return A new data frame containing one row per sample and only the columns
#'         that were constant for all cells within each sample.
get_sample_metadata <- function(cell_data_df,
                                sample_col) {
  # Step 1: Group the data by `sample_col` and check for uniqueness
  distinct_counts <- cell_data_df %>%
    dplyr::group_by(!!rlang::sym(sample_col)) %>%
    dplyr::summarise(across(everything(), ~ n_distinct(.x))) %>%
    dplyr::ungroup()

  # Step 2: Identify columns that are constant across all samples
  # The column is considered constant if the max distinct count within that column is 1.
  constant_cols <- distinct_counts %>%
    # Use select(where(...)) to find constant columns, excluding the grouping column (sample_col)
    dplyr::select(where(~ max(.x) == 1)) %>%
    names()

  # Ensure the sample column itself is not checked for constancy but is included in the final set
  cols_to_keep <- unique(c(sample_col, constant_cols))

  # Step 3: Select the constant columns and unique rows
  metadata <- cell_data_df %>%
    # Select only the relevant constant columns
    dplyr::select(all_of(cols_to_keep)) %>%
    dplyr::distinct()

  # Step 4: Ensure the result is a data frame using the safe subsetting operator
  # Set row names to the sample column
  # Use double-bracket subsetting for extraction, but keep metadata as a data frame.
  rownames(metadata) <- metadata[[sample_col]]

  # Remove the sample column from the final metadata table
  # Use safe subsetting to keep the data frame structure (df[, cols] returns a df)
  cols_to_return <- setdiff(colnames(metadata), sample_col)

  # Ensure the result is a data frame even if only one column is left or if no columns are left (empty df)
  metadata <- metadata[, cols_to_return, drop = FALSE]

  return(metadata)
}




#' Calculate Pseudobulk from Count Matrix
#'
#' @param count_matrix A gene x cell count matrix (can be dense matrix or sparse Matrix)
#' @param sample_ids A vector of sample identifiers, one per cell (same length as ncol(count_matrix))
#' @param min_cells Minimum number of cells required per sample (default: 1). Samples with fewer cells will be excluded.
#'
#' @return A gene x sample pseudobulk count matrix
#'
#' @examples
#' # From Seurat object
#' pb <- calculate_pseudobulk(
#'   count_matrix = seurat[["RNA"]]$counts,
#'   sample_ids = seurat$sample_id
#' )
#'
#' # From SingleCellExperiment object
#' pb <- calculate_pseudobulk(
#'   count_matrix = assay(sce, "counts"),
#'   sample_ids = colData(sce)$sample_id
#' )
#'
calculate_pseudobulk <- function(count_matrix,
                                 sample_ids,
                                 min_cells = 1) {
  # Input validation
  if (ncol(count_matrix) != length(sample_ids)) {
    stop("Length of sample_ids must equal the number of columns in count_matrix")
  }

  if (any(is.na(sample_ids))) {
    stop("sample_ids contains NA values. Please remove or impute missing values.")
  }

  # Convert sample_ids to factor for grouping
  sample_ids <- as.factor(sample_ids)

  # Count cells per sample
  cells_per_sample <- table(sample_ids)

  # Filter samples with insufficient cells
  if (min_cells > 1) {
    valid_samples <- names(cells_per_sample)[cells_per_sample >= min_cells]

    if (length(valid_samples) == 0) {
      stop(paste("No samples have >=", min_cells, "cells"))
    }

    # Subset to valid samples
    keep_cells <- sample_ids %in% valid_samples
    count_matrix <- count_matrix[, keep_cells, drop = FALSE]
    sample_ids <- droplevels(sample_ids[keep_cells])

    # Report filtered samples
    n_filtered <- length(cells_per_sample) - length(valid_samples)
    if (n_filtered > 0) {
      message(paste("Filtered out", n_filtered, "sample(s) with <", min_cells, "cells"))
    }
  }

  # Aggregate by summing across samples
  # Note: rowsum works on rows, so we transpose, aggregate, then transpose back
  pb <- rowsum(t(count_matrix), group = sample_ids)
  pb <- t(pb) # Transpose back to genes x samples format

  # Report summary
  message(paste("Pseudobulk matrix created:", nrow(pb), "genes x", ncol(pb), "samples"))

  return(pb)
}




#' DESeq2 Normalization of Pseudobulk Data
#'
#' @param pb A gene x sample pseudobulk count matrix (Genes as rows, Samples as columns).
#' @param metadata A data.frame with sample-level metadata (rownames or a column must match colnames of pb).
#' @param hvg Optional character vector of gene names to use as highly variable genes.
#'             If NULL, will select top variable genes automatically.
#' @param nvar_genes Number of top variable genes to select (default: 2000).
#'                   Only used if hvg = NULL.
#'
#' @return A normalized expression matrix (VST-transformed) with **samples as rows** and **genes as columns**.
#'
#' @examples
#' # Auto-select HVGs
#' # pb_norm <- deseq2_normalize(pb, metadata)
#'
#' # Use pre-defined HVGs
#' # my_hvgs <- c("Gene1", "Gene2", "Gene3")
#' # pb_norm <- deseq2_normalize(pb, metadata, hvg = my_hvgs)
#'
deseq2_normalize <- function(pb,
                             metadata,
                             hvg = NULL,
                             nvar_genes = 2000) {
  # Check that all samples in pb are in metadata
  if (!all(colnames(pb) %in% rownames(metadata))) {
    stop("Not all sample names in pb are found in metadata rownames")
  }

  # Reorder metadata to match pb
  metadata <- metadata[colnames(pb), , drop = FALSE]

  suppressMessages({
    suppressWarnings({
      # Create DESeq2 dataset
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = pb,
        colData = metadata,
        design = stats::formula(paste("~ 1"))
      )

      # Estimate size factors
      dds <- DESeq2::estimateSizeFactors(dds)

      # Set minimum number of counts per gene for VST
      nsub <- min(1000, sum(rowMeans(BiocGenerics::counts(dds, normalized = TRUE)) > 10))

      # Transform counts using variance stabilizing transformation (VST)
      dds <- DESeq2::vst(dds, blind = TRUE, nsub = nsub)
      pb_norm <- SummarizedExperiment::assay(dds) # Genes x Samples format

      # Select highly variable genes
      if (!is.null(hvg)) {
        # Use user-provided HVGs
        message(paste("Using", length(hvg), "user-provided highly variable genes"))

        hvg_available <- hvg[hvg %in% rownames(pb_norm)]
        hvg_missing <- setdiff(hvg, hvg_available)

        if (length(hvg_missing) > 0) {
          warning(paste(length(hvg_missing), "genes not found in the data and will be excluded"))
        }

        if (length(hvg_available) == 0) {
          stop("None of the provided HVGs are found in the data")
        }

        pb_norm <- pb_norm[hvg_available, , drop = FALSE]
      } else {
        # Auto-select top variable genes
        rv <- MatrixGenerics::rowVars(pb_norm)
        n_genes_to_select <- min(nvar_genes, length(rv))
        select <- order(rv, decreasing = TRUE)[seq_len(n_genes_to_select)]
        select <- rownames(pb_norm)[select]
        pb_norm <- pb_norm[select, , drop = FALSE]

        message(paste("Selected top", nrow(pb_norm), "highly variable genes"))
      }
    })
  })

  # --- Transpose to Samples x Genes format (REQUIRED FORMAT) ---
  # Current format: Genes x Samples (pb_norm)
  pb_norm <- t(pb_norm) # New format: Samples x Genes

  return(as.data.frame(pb_norm))
}
