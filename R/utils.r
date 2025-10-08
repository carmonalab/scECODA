# Define ECODA object to store data
methods::setClass(
  Class = "ECODA",
  slots = list(
    counts = "data.frame",
    counts_imp = "data.frame",
    freq = "data.frame",
    freq_imp = "data.frame",
    clr = "data.frame",
    ct_var = "data.frame",
    hvcs = "data.frame",
    metadata = "data.frame",
    pb = "data.frame"
  ),
  # Define the default values for each slot.
  prototype = list(
    counts = NULL,
    counts_imp = NULL,
    freq = NULL,
    freq_imp = NULL,
    clr = NULL,
    ct_var = NULL,
    hvcs = NULL,
    metadata = NULL,
    pb = NULL
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
get_celltype_counts <- function(cell_data_df, sample_col, celltype_col) {
  cellcount_df <- table(cell_data_df[[sample_col]], cell_data_df[[celltype_col]]) %>%
    as.data.frame.matrix()
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
get_sample_metadata <- function(cell_data_df, sample_col) {
  # Step 1: Group the data by `sample_col` and check for uniqueness
  distinct_counts <- cell_data_df %>%
    dplyr::group_by(!!rlang::sym(sample_col)) %>%
    dplyr::summarise(across(everything(), ~ n_distinct(.x))) %>%
    dplyr::ungroup()

  # Step 2: Identify columns that are constant across all samples
  constant_cols <- distinct_counts %>%
    dplyr::select(where(~ max(.x) == 1)) %>%
    names()

  # Step 3: Select and return the constant columns
  metadata <- cell_data_df %>%
    dplyr::select(all_of(c(sample_col, constant_cols))) %>%
    dplyr::distinct() %>%
    as.data.frame() # Convert back to data.frame

  # Set row names to the sample column for easier access
  rownames(metadata) <- metadata[[sample_col]]
  metadata[[sample_col]] <- NULL

  return(metadata)
}










get_pb_deseq2 <- function(pb, metadata, sample_col = "Sample", hvg = NULL, nvar_genes = 2000) {
  # metadata <- get_metadata(seurat, sample_col)
  # metadata[sample_col] <- gsub("-", "_", metadata[sample_col])
  pb_norm <- as.data.frame(t(DESeq2.normalize(pb, metadata, nvar_genes)))
  return(pb_norm)
}


# get_pb <- function(seurat, sample_col = "Sample", hvg = NULL) {
#   pb <- as.matrix(Seurat::AggregateExpression(seurat, group.by = sample_col, assays = "RNA")[["RNA"]])
#   colnames(pb) <- gsub("-", "_", colnames(pb))
#   if (!is.null(hvg)) {
#     pb <- pb[hvg, ]
#   }
#   return(pb)
# }
#
# get_pb_sce <- function(sce, sample_col = "Sample", hvg = NULL) {
#   # Aggregate counts across cells by sample
#   pb_sce <- aggregateAcrossCells(sce,
#     ids = colData(sce)[[sample_col]],
#     use.assay.type = "counts"
#   )
#
#   # Extract the pseudobulk count matrix
#   pb <- assay(pb_sce, "counts")
#
#   colnames(pb) <- gsub("-", "_", colnames(pb))
#   if (!is.null(hvg)) {
#     pb <- pb[hvg, ]
#   }
#   return(pb)
# }



#' Calculate Pseudobulk from Count Matrix
#'
#' @param counts_matrix A gene x cell count matrix (can be dense matrix or sparse Matrix)
#' @param sample_ids A vector of sample identifiers, one per cell (same length as ncol(counts_matrix))
#' @param min_cells Minimum number of cells required per sample (default: 1). Samples with fewer cells will be excluded.
#'
#' @return A gene x sample pseudobulk count matrix
#'
#' @examples
#' # From Seurat object
#' pb <- calculate_pseudobulk(
#'   counts_matrix = seurat[["RNA"]]$counts,
#'   sample_ids = seurat$sample_id
#' )
#'
#' # From SingleCellExperiment object
#' pb <- calculate_pseudobulk(
#'   counts_matrix = assay(sce, "counts"),
#'   sample_ids = colData(sce)$sample_id
#' )
#'
calculate_pseudobulk <- function(counts_matrix, sample_ids, min_cells = 1) {
  # Input validation
  if (ncol(counts_matrix) != length(sample_ids)) {
    stop("Length of sample_ids must equal the number of columns in counts_matrix")
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
    counts_matrix <- counts_matrix[, keep_cells, drop = FALSE]
    sample_ids <- droplevels(sample_ids[keep_cells])

    # Report filtered samples
    n_filtered <- length(cells_per_sample) - length(valid_samples)
    if (n_filtered > 0) {
      message(paste("Filtered out", n_filtered, "sample(s) with <", min_cells, "cells"))
    }
  }

  # Aggregate by summing across samples
  # Note: rowsum works on rows, so we transpose, aggregate, then transpose back
  pb <- rowsum(t(counts_matrix), group = sample_ids)
  pb <- t(pb) # Transpose back to genes x samples format

  # Report summary
  message(paste("Pseudobulk matrix created:", nrow(pb), "genes x", ncol(pb), "samples"))

  return(pb)
}






# get_metadata <- function(sc_metadata, sample_col = "Sample") {
#   metadata <-  %>%
#     dplyr::group_by(!!sym(sample_col)) %>%
#     dplyr::slice(1)
#   return(metadata)
# }


get_pb_deseq2 <- function(pb, metadata, sample_col = "Sample", hvg = NULL, nvar_genes = 2000) {
  # metadata <- get_metadata(seurat, sample_col)
  # metadata[sample_col] <- gsub("-", "_", metadata[sample_col])
  pb_norm <- as.data.frame(t(DESeq2.normalize(pb, metadata, nvar_genes)))
  return(pb_norm)
}


DESeq2.normalize <- function(pb,
                             metadata,
                             nvar_genes = 2000) {
  suppressMessages({
    suppressWarnings({
      # Normalize pseudobulk data using DESeq2
      # do formula for design with the cluster_by elements in order
      pb <- DESeq2::DESeqDataSetFromMatrix(
        countData = pb,
        colData = metadata,
        design = stats::formula(paste("~ 1"))
      )

      pb <- DESeq2::estimateSizeFactors(pb)

      # Set minimum number of counts per gene
      nsub <- min(1000, sum(rowMeans(BiocGenerics::counts(pb, normalized = TRUE)) > 10))

      # transform counts using vst
      pb <- DESeq2::vst(pb, blind = T, nsub = nsub)
      pb <- SummarizedExperiment::assay(pb)

      # get top variable genes
      rv <- MatrixGenerics::rowVars(pb)
      select <- order(rv, decreasing = TRUE)[seq_len(min(nvar_genes, length(rv)))]
      select <- row.names(pb)[select]

      pb <- pb[select[select %in% row.names(pb)], ]
    })
  })

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
