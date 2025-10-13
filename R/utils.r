# Define ECODA object to store data

#' An S4 class to represent a compositional data analysis object (ECODA).
#'
#' This class is designed to store various forms of compositional data derived
#' from cell counts (or similar proportional data), alongside associated metadata,
#' transformation results, and variability analysis outcomes.
#'
#' @slot counts Original cell count data (samples as rows, cell types as columns).
#' @slot counts_imp Imputed cell count data, typically used to handle zero counts.
#' @slot freq Relative frequency (percentage) of cell types derived from `counts`.
#' @slot freq_imp Relative frequency (percentage) derived from `counts_imp`.
#' @slot clr Centered Log-Ratio transformed data, derived from `freq_imp`.
#' @slot celltype_variances Data frame detailing the variance metrics for each cell type.
#' @slot variance_explained Numeric value indicating the total variance captured by
#'                        highly variable cell types (HVCs).
#' @slot top_n_hvcs Integer specifying the number of top highly variable cell types selected.
#' @slot highly_variable_celltypes Character vector of names of the highly variable cell types.
#' @slot metadata Data frame of sample-level metadata (samples as rows).
#' @slot pb Data frame of pseudobulk gene expression.
#' @slot sample_distances Data frame storing calculated distances between samples.
#'
#' #' @importFrom methods setClass
setClass(
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



#' Calculate relative frequencies (percentages) row-wise.
#'
#' This function takes a matrix or data frame of counts and transforms each row
#' (assumed to be a sample or observation) into relative frequencies,
#' expressing each component as a percentage of the row's total sum.
#'
#' @param df A data frame or matrix of counts (samples/observations as rows,
#'           components as columns). Must contain only numeric, non-negative values.
#' @return A data frame of the same dimensions, where each row sums to 100.
#' @export calc_freq
#' @importFrom base t apply sum as.data.frame
#' @examples
#' counts <- data.frame(A = c(10, 50), B = c(90, 50))
#' calc_freq(counts)
calc_freq <- function(df) {
  df <- t(apply(df, 1, function(row) (row / sum(row)) * 100)) %>%
    as.data.frame()
  return(df)
}



#' Perform the Centered Log-Ratio (CLR) transformation.
#'
#' The CLR transformation is a common technique used for compositional data
#' analysis, mapping the data from the Aitchison simplex to Euclidean space.
#' It is defined as the log of the ratio between each component and the geometric
#' mean of all components in that sample. **Note:** This function assumes the input
#' data is strictly positive (i.e., does not contain zeros). Use an imputation
#' or pseudocount method prior to this function if zeros are present.
#'
#' @param df A data frame or matrix of strictly positive relative frequencies
#'           or counts (samples/observations as rows, components as columns).
#' @return A data frame of the same dimensions containing the CLR-transformed values.
#' @export clr
#' @importFrom base apply exp mean log as.data.frame
#' @examples
#' freq_imp <- data.frame(A = c(10.1, 50.1), B = c(89.9, 49.9))
#' clr(freq_imp)
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
#' @importFrom base table as.data.frame.matrix
#' @examples
#' # Create example data frame
#' cell_data_df <- data.frame(
#'   Cell_ID = paste0("C", 1:10),
#'   Sample_Name = c(rep("S1", 5), rep("S2", 5)),
#'   Cluster_Annotation = factor(c(
#'     "B_Cell", "T_Cell", "B_Cell", NA, "T_Cell",
#'     "T_Cell", "Macrophage", "B_Cell", "T_Cell", "T_Cell"
#'   ))
#' )
#'
#' # Calculate cell type counts per sample
#' celltype_counts <- get_celltype_counts(
#'   cell_data_df = cell_data_df,
#'   sample_col = "Sample_Name",
#'   celltype_col = "Cluster_Annotation"
#' )
#'
#' print(celltype_counts)
#' # Note how the NA cell type count is handled and renamed to "NA".
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
#' This function identifies columns in a cell-level metadata data frame that
#' have a **constant** value for all cells belonging to the same sample.
#' It aggregates this constant information, returning a new data frame where
#' each row represents a unique sample. Columns that vary within any sample
#' are excluded from the output.
#'
#' @param cell_data_df A data frame containing cell-level metadata.
#'                     This should include the sample ID column and all
#'                     potential metadata columns.
#' @param sample_col A character string specifying the name of the column
#'                   that defines the unique sample ID for each cell (e.g., "Sample_ID").
#' @return A new data frame where:
#'         \itemize{
#'           \item Each row corresponds to a unique sample from the input data.
#'           \item The row names are set to the values of the input `sample_col`.
#'           \item Columns contain only the metadata fields that were constant
#'                 across all cells within **each** sample. The `sample_col` itself
#'                 is excluded from the final columns but used for row names.
#'         }
#' @importFrom dplyr group_by summarise across everything n_distinct ungroup select where all_of distinct
#' @importFrom rlang sym
#' @examples
#' \dontrun{
#' # Assuming you have a data frame 'cell_df'
#' cell_df <- data.frame(
#'   Cell_ID = paste0("C", 1:10),
#'   Sample_ID = c(rep("S1", 5), rep("S2", 5)),
#'   Age = c(rep(30, 5), rep(45, 5)),
#'   Gender = c(rep("M", 5), rep("F", 5)),
#'   Cell_Type = c(rep("A", 3), rep("B", 2), rep("A", 3), rep("B", 2))
#' )
#'
#' # The 'Age' and 'Gender' columns are constant within each sample (S1 and S2).
#' # The 'Cell_ID' and 'Cell_Type' columns vary within sample S1 and/or S2.
#'
#' sample_meta <- get_sample_metadata(cell_df, "Sample_ID")
#' print(sample_meta)
#' # Output will have 'Age' and 'Gender' as columns, with row names 'S1' and 'S2'.
#' }
get_sample_metadata <- function(cell_data_df,
                                sample_col) {
  # Step 1: Group the data by `sample_col` and check for uniqueness
  distinct_counts <- cell_data_df %>%
    group_by(!!sym(sample_col)) %>%
    summarise(across(everything(), ~ n_distinct(.x))) %>%
    ungroup()

  # Step 2: Identify columns that are constant across all samples
  # The column is considered constant if the max distinct count within that column is 1.
  constant_cols <- distinct_counts %>%
    select(where(~ max(.x) == 1)) %>%
    names()

  # Ensure the sample column itself is not checked for constancy but is included in the final set
  cols_to_keep <- unique(c(sample_col, constant_cols))

  # Step 3: Select the constant columns and unique rows
  metadata <- cell_data_df %>%
    # Select only the relevant constant columns
    select(all_of(cols_to_keep)) %>%
    distinct()

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
#' This function aggregates single-cell count data into a "pseudobulk" matrix
#' by summing the counts for all cells belonging to the same sample ID.
#' It is robust to both dense and sparse count matrices. It also includes
#' filtering logic to exclude samples that do not meet a minimum cell count threshold.
#'
#' @param count_matrix A gene x cell count matrix (can be a dense matrix or sparse Matrix).
#'                     Gene identifiers should be row names and cell barcodes should be column names.
#' @param sample_ids A vector of sample identifiers, one for each column (cell) in
#'                   \code{count_matrix}. The length must equal \code{ncol(count_matrix)}.
#'                   Must not contain \code{NA} values.
#' @param min_cells Minimum number of cells required per sample (default: 1).
#'                  Samples with fewer cells than this threshold will be excluded
#'                  from the final pseudobulk matrix.
#'
#' @return A gene x sample pseudobulk count matrix. The columns correspond to
#'         the unique sample IDs, and the rows correspond to the genes.
#' @export calculate_pseudobulk
#' @importFrom base ncol length any is.na as.factor table names paste droplevels message rowsum t stop
#' @examples
#' \dontrun{
#' # Assuming seurat is a loaded Seurat object
#' pb_seurat <- calculate_pseudobulk(
#'   count_matrix = seurat[["RNA"]]$counts,
#'   sample_ids = seurat$sample_id
#' )
#'
#' # Assuming sce is a loaded SingleCellExperiment object
#' pb_sce <- calculate_pseudobulk(
#'   count_matrix = assay(sce, "counts"),
#'   sample_ids = colData(sce)$sample_id,
#'   min_cells = 5 # Example of filtering samples with < 5 cells
#' )
#' }
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
#' This function normalizes a gene x sample pseudobulk count matrix using the
#' Variance Stabilizing Transformation (VST) from the \code{DESeq2} package.
#' It estimates size factors and variance for all genes, performs the VST, and then
#' subsets the results to include only highly variable genes (HVCs),
#' either specified by the user or automatically selected based on variance.
#'
#' @param pb A gene x sample pseudobulk count matrix (Genes as rows, Samples as columns).
#'           Must contain non-negative integer counts.
#' @param metadata A data.frame with sample-level metadata. Row names must exactly
#'                 match the column names of \code{pb}. The function will reorder
#'                 the metadata to match \code{pb}.
#' @param hvg Optional character vector of gene names to use as highly variable genes.
#'            If provided, only these genes will be returned after VST.
#' @param nvar_genes Number of top variable genes to select (default: 2000).
#'                   Only used if \code{hvg = NULL}; the function will select the
#'                   top \code{nvar_genes} by variance after VST.
#'
#' @return A normalized expression matrix (VST-transformed) with **samples as rows**
#'         and **genes as columns**.
#'
#' @export deseq2_normalize
#'
#' @importFrom base all colnames rownames stop paste length warning message t as.data.frame
#' @importFrom base is.null setdiff min order rowMeans sum seq_len
#' @importFrom stats formula
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors vst
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics counts
#' @importFrom MatrixGenerics rowVars
#'
#' @examples
#' \dontrun{
#' # Assuming 'pb' is the pseudobulk matrix and 'metadata' is the sample annotation
#'
#' # 1. Auto-select top 2000 most variable genes after VST
#' pb_norm_auto <- deseq2_normalize(pb, metadata)
#'
#' # 2. Use pre-defined set of highly variable genes
#' my_hvgs <- c("Gene1", "Gene2", "Gene3")
#' pb_norm_hvg <- deseq2_normalize(pb, metadata, hvg = my_hvgs)
#' }
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
      dds <- DESeqDataSetFromMatrix(
        countData = pb,
        colData = metadata,
        design = formula(paste("~ 1"))
      )

      # Estimate size factors
      dds <- estimateSizeFactors(dds)

      # Set minimum number of counts per gene for VST
      nsub <- min(1000, sum(rowMeans(counts(dds, normalized = TRUE)) > 10))

      # Transform counts using variance stabilizing transformation (VST)
      dds <- vst(dds, blind = TRUE, nsub = nsub)
      pb_norm <- assay(dds) # Genes x Samples format

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
        rv <- rowVars(pb_norm)
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
