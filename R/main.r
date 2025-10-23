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
#' @slot hvcs Character vector of names of the highly variable cell types.
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
    asin_sqrt = "data.frame",
    pb = "data.frame",
    celltype_variances = "data.frame",
    top_n_hvcs = "integer",
    hvcs = "character",
    variance_explained = "numeric",
    metadata = "data.frame",
    sample_distances = "data.frame"
  ),
  # Define the default values for each slot.
  prototype = list(
    counts = NULL,
    counts_imp = NULL,
    freq = NULL,
    freq_imp = NULL,
    clr = NULL,
    asin_sqrt = NULL,
    pb = NULL,
    celltype_variances = NULL,
    top_n_hvcs = NULL,
    hvcs = NULL,
    variance_explained = NULL,
    metadata = NULL,
    sample_distances = NULL
  )
)


#' Creates an ECODA object from various data types.
#'
#' This is a smart constructor function used to initialize an
#' \link[=ECODA-class]{ECODA} object. It handles the processing of single-cell objects
#' (\code{Seurat} or \code{SingleCellExperiment}) or raw data frames to extract
#' cell type counts, calculate sample metadata, and optionally generate
#' DESeq2-normalized pseudobulk data.
#'
#' @param data The primary input, which can be:
#'             \itemize{
#'               \item A \code{Seurat} object.
#'               \item A \code{SingleCellExperiment} object.
#'               \item A pre-calculated sample x cell type count matrix/data frame.
#'             }
#' @param sample_col The column name in the single-cell object's metadata (or cell-level
#'                   data frame) that defines the unique sample ID for each cell.
#'                   (Required if \code{data} is a single-cell object).
#' @param celltype_col The column name in the single-cell object's metadata (or cell-level
#'                     data frame) that defines the cell type annotation for each cell.
#'                     (Required if \code{data} is a single-cell object).
#' @param get_pb Logical, if \code{TRUE} (default: \code{FALSE}), the function will
#'               calculate and store DESeq2-normalized pseudobulk data in the \code{pb} slot.
#' @param variance_explained Numeric (default: 0.5). Used in subsequent steps by
#'                         \code{create_ecoda_object_helper} to determine how many
#'                         highly variable cell types (HVCs) to select.
#' @param top_n_hvcs Integer (optional). Overrides \code{variance_explained} if provided,
#'                   specifying the exact number of top HVCs to select.
#'
#' @return A new \link[=ECODA-class]{ECODA} object populated with \code{counts}, \code{metadata},
#'         and optionally \code{pb} and the initial compositional analysis results.
#'
#' @importFrom methods new
#' @importFrom utils read.csv
#'
#' @export create_ecoda_object
#'
#' @seealso \link[=ECODA-class]{ECODA}, \code{\link{calculate_pseudobulk}}, \code{\link{get_celltype_counts}}, \code{\link{get_sample_metadata}}
#'
#' @examples
#' \dontrun{
#' # Example using a Seurat object (assuming object 'pbmc' is loaded)
#' ecoda_obj <- create_ecoda_object(
#'   data = pbmc,
#'   sample_col = "orig.ident",
#'   celltype_col = "cell_type",
#'   get_pb = TRUE
#' )
#'
#' # Example using a pre-calculated count matrix
#' # counts_df <- read.csv("cell_counts.csv", row.names = 1)
#' # ecoda_obj <- create_ecoda_object(data = counts_df)
#' }
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
      if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop(
          "Package \"SummarizedExperiment\" must be installed to use this function.",
          call. = FALSE
        )
      }
      cell_data_df <- as.data.frame(SummarizedExperiment::colData(data))
      names(cell_data_df) <- names(SummarizedExperiment::colData(data))
      if (get_pb) {
        pb <- calculate_pseudobulk(
          count_matrix = SummarizedExperiment::assay(data, "counts"),
          sample_ids = cell_data_df[[sample_col]]
        )
      }
    }
  } else {
    cell_data_df <- data
  }

  counts <- get_celltype_counts(cell_data_df, sample_col, celltype_col)
  metadata <- get_sample_metadata(cell_data_df, sample_col)

  ecoda_object <- create_ecoda_object_helper(
    counts = counts,
    metadata = metadata,
    variance_explained = variance_explained,
    top_n_hvcs = top_n_hvcs
  )

  if (get_pb) {
    pb <- deseq2_normalize(pb)
    pb <- pb[mixedsort(rownames(pb)), ]
    ecoda_object@pb <- pb
  }

  return(ecoda_object)
}


#' Creates an ECODA object from pre-calculated cell type counts.
#'
#' This is the core constructor function that initializes and performs the initial
#' compositional analysis steps for an \link[=ECODA-class]{ECODA} object, assuming the
#' cell type count matrix is already available. It handles zero-imputation,
#' calculates relative frequencies, Centered Log-Ratio (CLR) transformed data,
#' sample distances, and identifies highly variable cell types (HVCs).
#'
#' @param counts A data frame or matrix of cell type counts where **rows are samples**
#'               and **columns are cell types**. Must contain non-negative integers.
#' @param freq A data frame or matrix of **relative cell type frequencies** (proportions)
#'             where **rows are samples** and **columns are cell types**. Values must be
#'             between 0 and 100 (or close to 1 across rows). If provided, this matrix
#'             is used directly for CLR transformation.
#' @param metadata An optional data frame containing sample-level metadata. Row names
#'                 must match the row names of the \code{counts} matrix.
#' @param variance_explained Numeric (default: 0.5). The proportion of total variance
#'                           that should be captured by the selected highly variable
#'                           cell types (HVCs). Used by \code{find_hvcs}.
#' @param top_n_hvcs Integer (optional). Overrides \code{variance_explained} if provided,
#'                   specifying the exact number of top HVCs to select based on variance.
#'
#' @return A fully initialized \link[=ECODA-class]{ECODA} object, populated with counts,
#'         frequency data, CLR transformed data, sample distances, and HVCs.
#'
#' @importFrom methods new
#' @importFrom gtools mixedsort
#' @importFrom dplyr %>%
#' @importFrom stats dist
#'
#' @export create_ecoda_object_helper
#'
#' @seealso \link[=ECODA-class]{ECODA}, \code{\link{calc_freq}}, \code{\link{clr}}, \code{\link{find_hvcs}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'my_counts' and 'my_metadata' are defined data frames:
#' counts_df <- data.frame(row.names = c("S1", "S2"), A = c(10, 5), B = c(0, 15))
#' meta_df <- data.frame(row.names = c("S1", "S2"), Group = c("Treated", "Control"))
#'
#' ecoda_obj <- create_ecoda_object_helper(
#'   counts = counts_df,
#'   metadata = meta_df,
#'   top_n_hvcs = 2
#' )
#' # ecoda_obj@counts will contain the original counts
#' # ecoda_obj@clr will contain the CLR transformed data
#' }
create_ecoda_object_helper <- function(counts = NULL,
                                       freq = NULL,
                                       metadata = NULL,
                                       variance_explained = 0.5,
                                       top_n_hvcs = NULL) {
  # Initialize the object with default values
  ecoda_object <- new("ECODA")

  if (!is.null(counts) && !is.null(freq)) {
    stop("Please provide only counts or freq.")
  }

  if (!is.null(counts)) {
    counts <- counts[mixedsort(rownames(counts)), ]

    counts_imp <- counts
    if (any(counts == 0)) {
      counts_imp <- counts_imp + 1
    }
    freq <- calc_freq(counts)
    freq_imp <- calc_freq(counts_imp)
    clr_df <- clr(counts_imp)

    ecoda_object@counts <- counts
    ecoda_object@counts_imp <- counts_imp
    ecoda_object@freq <- freq
    ecoda_object@freq_imp <- freq_imp
  }

  if (!is.null(freq)) {
    freq <- freq[mixedsort(rownames(freq)), ]

    row_sums <- rowSums(freq)

    # Check if all row sums are close to 100 (using a tolerance for floating point numbers)
    # If they are intended to be percentages
    if (!all(abs(row_sums - 100) < 1e-6)) {
      # Check if they are close to 1 and should be scaled to 100
      if (all(abs(row_sums - 1) < 1e-6)) {
        warning("Frequencies sum close to 1. Rescaling all rows to sum to 100.")
        # Scale by 100
        freq <- freq * 100
      } else {
        # If they sum to neither 1 nor 100, warn and scale them to 100
        warning(paste(
          "Frequencies do not sum to 100 (or 1).",
          "Each row will be scaled so the new row sum is 100."
        ))
        # Re-scale each row so it sums to 100: (freq / row_sum) * 100
        freq <- (freq / row_sums) * 100
      }
    }

    freq_imp <- freq
    if (any(freq == 0)) {
      warning("freq contains zeros. Zeros were imputed with the smallest value in freq.")
      freq_imp <- freq_imp + min(freq_imp[freq_imp > 0])
    }
    clr_df <- clr(freq_imp)

    ecoda_object@freq <- freq
    ecoda_object@freq_imp <- freq_imp
  }

  if (!is.null(metadata)) {
    # Sort by rownames
    metadata <- metadata[mixedsort(rownames(metadata)), , drop = FALSE]

    # Ensure correct rownames
    if (!all(rownames(counts) == rownames(metadata))) {
      stop("Rownames of cell counts do not match metadata.")
    }
    ecoda_object@metadata <- metadata
  }

  ecoda_object@asin_sqrt <- freq %>%
    mutate(across(everything(), ~ . / 100)) %>%
    sqrt() %>%
    asin()

  sdist <- clr_df %>%
    dist() %>%
    as.matrix() %>%
    as.data.frame()

  ecoda_object@clr <- clr_df
  ecoda_object@sample_distances <- sdist

  ecoda_object <- find_hvcs(
    ecoda_object,
    variance_explained,
    top_n_hvcs
  )

  return(ecoda_object)
}


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
#' @param cell_data_df A data frame where each row represents a single cell, and columns
#'                     contain cell-level information, including sample ID and cell type annotation.
#' @param sample_col The column that defines the sample ID for each cell
#' @param celltype_col The column that defines the cell type annotation for each cell
#'
#' @return A data frame with samples as rows and cell types as columns,
#'         containing the count of each cell type per sample.
#' @export get_celltype_counts
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
#' @export get_sample_metadata
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


# Find highly variable cell types ---------------------------

#' Identifies and stores Highly Variable Cell Types (HVCs) in an ECODA object.
#'
#' This function calculates the variance of each cell type across all samples
#' based on the Centered Log-Ratio (CLR) transformed data. It then selects the
#' subset of highly variable cell types (HVCs) that collectively explain a
#' specified proportion of the total variance, or a fixed number of top cell types.
#'
#' @param ecoda_object An initialized \link[=ECODA-class]{ECODA} object containing the
#'                     CLR-transformed data in the \code{clr} slot.
#' @param variance_explained Numeric (default: 0.5). The target cumulative
#'                           proportion of total variance to be explained by the
#'                           selected HVCs. The function stops selecting cell types
#'                           once this threshold is met.
#' @param top_n_hvcs Integer (optional). If provided, this overrides
#'                   \code{variance_explained} and selects exactly the top N cell
#'                   types with the highest variance.
#'
#' @return The updated \link[=ECODA-class]{ECODA} object with the following slots populated:
#'         \itemize{
#'           \item \code{celltype_variances}: Data frame of cell type variances.
#'           \item \code{variance_explained}: The exact variance proportion captured by the selected HVCs.
#'           \item \code{top_n_hvcs}: The number of HVCs selected.
#'           \item \code{hvcs}: Character vector of the names of the selected HVCs.
#'         }
#'
#' @export find_hvcs
#'
#' @seealso \link[=ECODA-class]{ECODA}, \code{\link{get_celltype_variances}}, \code{\link{get_hvcs}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'ecoda_obj' is a previously created ECODA object:
#'
#' # Select HVCs that explain 75% of the total variance
#' ecoda_obj <- find_hvcs(ecoda_obj, variance_explained = 0.75)
#'
#' # Select exactly the top 5 most variable cell types
#' ecoda_obj <- find_hvcs(ecoda_obj, top_n_hvcs = 5)
#' }
find_hvcs <- function(ecoda_object,
                      variance_explained = 0.5,
                      top_n_hvcs = NULL) {
  df_var <- get_celltype_variances(
    ecoda_object,
    show_plot = FALSE
  )

  hvcs <- get_hvcs(
    df_var,
    variance_explained,
    top_n_hvcs = top_n_hvcs
  )

  variance_explained <- (sum(df_var$Variance[1:length(hvcs)]) /
    sum(df_var$Variance))

  ecoda_object@celltype_variances <- df_var
  ecoda_object@variance_explained <- variance_explained
  ecoda_object@top_n_hvcs <- length(hvcs)
  ecoda_object@hvcs <- hvcs


  return(ecoda_object)
}


#' Calculates the variance of cell types across samples.
#'
#' This function takes the Centered Log-Ratio (CLR) transformed cell type
#' abundance data from an \link[=ECODA-class]{ECODA} object, calculates the mean CLR
#' abundance and variance for each cell type, and optionally generates a
#' mean-variance plot. It also calculates the cumulative variance explained
#' by the cell types when ranked by variance.
#'
#' @param ecoda_object An initialized \link[=ECODA-class]{ECODA} object containing
#'                     CLR-transformed data in the \code{clr} slot.
#' @param show_plot Logical (default: \code{TRUE}). If \code{TRUE}, a mean-variance
#'                  plot is generated using the internal \code{plot_varmean} function.
#' @param label_points Logical (default: \code{TRUE}). If \code{TRUE}, the points
#'                     on the plot will be labeled with cell type names.
#' @param plot_title Character string (default: ""). Title for the generated plot.
#' @param smooth_method Character string (default: "lm"). Smoothing method to use
#'                      for the regression line in the mean-variance plot (e.g., "lm", "loess").
#' @param descending Logical (default: \code{TRUE}). If \code{TRUE}, the returned
#'                   data frame is sorted by variance in descending order.
#'
#' @return A data frame (samples as rows, cell types as columns) containing:
#'         \itemize{
#'           \item \code{celltype}: The name of the cell type.
#'           \item \code{avg_clr_abundance}: The mean CLR value for that cell type across all samples.
#'           \item \code{Variance}: The variance of the CLR values for that cell type across all samples.
#'           \item \code{cumulative_variance}: Cumulative sum of the variance, ranked by variance.
#'           \item \code{variance_exp}: The proportion of total variance explained cumulatively.
#'         }
#'
#' @importFrom stats var
#' @importFrom dplyr %>% everything group_by summarize arrange desc mutate
#' @importFrom tidyr pivot_longer
#'
#' @export get_celltype_variances
#'
#' @seealso \link[=ECODA-class]{ECODA}, \code{\link{plot_varmean}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'ecoda_obj' is a created ECODA object
#'
#' # Calculate variances and show the plot
#' df_var <- get_celltype_variances(ecoda_obj)
#'
#' # Calculate variances without plotting, sorted by ascending variance
#' df_var_asc <- get_celltype_variances(ecoda_obj, show_plot = FALSE, descending = FALSE)
#' }
get_celltype_variances <- function(ecoda_object,
                                   show_plot = TRUE,
                                   label_points = TRUE,
                                   plot_title = "",
                                   smooth_method = "lm",
                                   descending = TRUE) {
  df_var <- ecoda_object@clr %>%
    pivot_longer(
      cols = everything(),
      names_to = "celltype",
      values_to = "values"
    ) %>%
    group_by(celltype) %>%
    summarize(
      avg_clr_abundance = mean(values, na.rm = TRUE),
      Variance = var(values, na.rm = TRUE)
    )

  if (descending) {
    df_var <- df_var %>%
      arrange(desc(Variance))
  } else {
    df_var <- df_var %>%
      arrange(Variance)
  }

  total_variance <- sum(df_var$Variance)

  df_var <- df_var %>%
    mutate(
      cumulative_variance = cumsum(Variance),
      variance_exp = (cumulative_variance / total_variance)
    )

  if (show_plot) {
    p <- plot_varmean(
      ecoda_object,
      plot_title,
      smooth_method,
      label_points
    )
    print(p)
  }

  return(as.data.frame(df_var))
}


#' Selects Highly Variable Cell Types (HVCs) based on variance or count threshold.
#'
#' This is a utility function that selects the most variable cell types from a
#' variance-ranked data frame. Selection can be controlled either by defining the
#' cumulative variance to be explained or by specifying a fixed number of cell types.
#'
#' @param df_var A data frame (typically the output of \code{\link{get_celltype_variances}})
#'               that must contain at least the columns \code{celltype} and \code{variance_exp}
#'               (cumulative variance explained) and must be sorted by variance in
#'               descending order.
#' @param variance_explained Numeric (default: 0.5). The cumulative proportion of
#'                           total variance that the selected HVCs must explain.
#'                           This parameter is ignored if \code{top_n_hvcs} is set.
#' @param top_n_hvcs Integer or Numeric (optional).
#'                   \itemize{
#'                     \item If an integer (>= 1), it specifies the exact number of
#'                           top cell types to select.
#'                     \item If a numeric value between 0 and 1, it specifies the
#'                           proportion of cell types to select (e.g., 0.1 for the top 10%).
#'                   }
#'
#' @return A character vector containing the names of the selected highly variable cell types.
#'         The function ensures that at least two cell types are always returned.
#'
#' @importFrom dplyr %>% slice pull slice_head
#' @importFrom rlang .data
#'
#' @export get_hvcs
#'
#' @seealso \code{\link{find_hvcs}}, \code{\link{get_celltype_variances}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'df_variance' is a variance data frame sorted descendingly by variance.
#' # Select cell types explaining at least 60% of variance:
#' hvcs_60 <- get_hvcs(df_variance, variance_explained = 0.6)
#'
#' # Select the top 10 cell types:
#' hvcs_top10 <- get_hvcs(df_variance, top_n_hvcs = 10)
#'
#' # Select the top 5% of cell types:
#' hvcs_prop <- get_hvcs(df_variance, top_n_hvcs = 0.05)
#' }
get_hvcs <- function(df_var,
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
      slice(1:last_hvc_index) %>%
      pull(celltype)
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



#' @title Generates a Mean-Variance Plot for CLR-transformed Cell Type Data.
#'
#' @description
#' This function visualizes the relationship between the mean abundance (CLR) and
#' the variance (CLR) for each cell type, which is typically used to identify
#' Highly Variable Cell Types (HVCs).
#'
#' @details
#' If \code{highlight_hvcs} is \code{TRUE}, cell types previously identified and
#' stored in the \code{ecoda_object@highly_variable_celltypes} slot will be
#' highlighted in red on the plot.
#'
#' @param ecoda_object An \link[=ECODA-class]{ECODA} object containing pre-calculated
#'                     cell type variances in the \code{celltype_variances} slot and
#'                     the HVC list in the \code{highly_variable_celltypes} slot.
#' @param plot_title Character string (default: ""). The title for the plot.
#' @param highlight_hvcs Logical (default: \code{TRUE}). If \code{TRUE}, the points
#'                       corresponding to the Highly Variable Cell Types (HVCs)
#'                       stored in the ECODA object are colored red.
#' @param label_points Logical (default: \code{FALSE}). If \code{TRUE}, cell type
#'                     names are added to the points using \code{ggrepel::geom_text_repel}
#'                     to minimize overlap.
#' @param plot_fit_line Logical (default: \code{FALSE}). If \code{TRUE}, a smoothing
#'                      regression line is added to the plot.
#' @param smooth_method Character string (default: "lm"). The smoothing method to
#'                      use for the regression line if \code{plot_fit_line} is \code{TRUE}
#'                      (e.g., "lm" for linear model, "loess" for local regression).
#'
#' @return A \code{ggplot} object representing the mean-variance plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_classic xlab ylab scale_color_manual
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr mutate if_else
#'
#' @export plot_varmean
#'
#' @seealso \link[=ECODA-class]{ECODA}, \code{\link{get_celltype_variances}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'ecoda_obj' is a created ECODA object with calculated HVCs:
#'
#' # 1. Generate the plot, highlighting HVCs (default)
#' plot_varmean(ecoda_obj, plot_title = "HVCs on Mean-Variance Plot")
#'
#' # 2. Generate the plot without highlighting HVCs and add a fit line
#' plot_varmean(ecoda_obj,
#'   highlight_hvcs = FALSE,
#'   plot_fit_line = TRUE,
#'   smooth_method = "loess"
#' )
#' }
plot_varmean <- function(ecoda_object,
                         plot_title = "",
                         highlight_hvcs = TRUE,
                         label_points = FALSE,
                         plot_fit_line = FALSE,
                         smooth_method = "lm") {
  df_var <- ecoda_object@celltype_variances

  # --- 1. Create a highlighting factor column ---
  if (highlight_hvcs) {
    highlight_celltypes <- ecoda_object@hvcs

    df_var <- df_var %>%
      dplyr::mutate(
        is_highlighted = dplyr::if_else(
          .data$celltype %in% highlight_celltypes,
          "HVC selected",
          "Not selected"
        )
      )
    color_map <- c("HVC selected" = "red", "Not selected" = "black")

    p <- ggplot(df_var, aes(x = avg_clr_abundance, y = Variance, color = is_highlighted)) +
      scale_color_manual(values = color_map, name = "Cell Type Group")
  } else {
    p <- ggplot(df_var, aes(x = avg_clr_abundance, y = Variance))
  }

  p <- p +
    geom_point() +
    labs(title = paste(plot_title)) +
    theme_classic() +
    xlab("Mean (CLR)") +
    ylab("Variance (CLR)")

  if (plot_fit_line) {
    p <- p +
      geom_smooth(method = smooth_method, color = "red", fill = "#69b3a2", se = TRUE)
  }

  # if (label_points) {
  #   p <- p + geom_text_repel(data = df_var, aes(label = celltype), vjust = -0.5)
  # }

  # Ensure the legend is not shown if we are not highlighting anything
  if (!highlight_hvcs) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  if (label_points) {
    # Use the color aesthetic inside ggrepel so highlighted labels match the point color
    if (highlight_hvcs) {
      p <- p + ggrepel::geom_text_repel(data = df_var, aes(label = celltype, color = is_highlighted), size = 3)
    } else {
      p <- p + ggrepel::geom_text_repel(data = df_var, aes(label = celltype), size = 3, color = "black")
    }
  }

  return(p)
}




# Get Pseudobulk and normalize ---------------------------

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
#'   count_matrix = SummarizedExperiment::assay(sce, "counts"),
#'   sample_ids = SummarizedExperiment::colData(sce)$sample_id,
#'   min_cells = 5 # Example of filtering samples with < 5 cells
#' )
#' }
calculate_pseudobulk <- function(count_matrix,
                                 sample_ids,
                                 min_cells = 10) {
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
    count_matrix <- as.matrix(count_matrix)[, keep_cells, drop = FALSE]
    sample_ids <- droplevels(sample_ids[keep_cells])

    # Report filtered samples
    n_filtered <- length(cells_per_sample) - length(valid_samples)
    if (n_filtered > 0) {
      message(paste("Filtered out", n_filtered, "sample(s) with <", min_cells, "cells"))
    }
  }

  # Aggregate by summing across samples
  pb <- t(rowsum(t(count_matrix), group = sample_ids))

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
#' **Note:** If `deseq2_design` is not specified, the VST is performed using a
#' minimal design formula (`~ 1`) for size factor and variance estimation.
#'
#' @param pb A gene x sample pseudobulk count matrix (Genes as rows, Samples as columns).
#'           Must contain non-negative integer counts.
#' @param metadata A data.frame with sample-level metadata. Row names must exactly
#'                 match the column names of \code{pb}. The function will reorder
#'                 the metadata to match \code{pb}.
#' @param deseq2_design Optional: A \code{formula} object specifying the design
#'                      used by \code{DESeq2} to estimate size factors and variance
#'                      for the VST. This formula should reference columns in the
#'                      \code{metadata} data frame. If \code{NULL} (the default),
#'                      the minimal design \code{~ 1} is used. Note: the design
#'                      will be considered for pseudobulk normalization as the
#'                      function uses vst(blind = FALSE, ...) internally.
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
#' @importFrom stats formula
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors vst
#' @importFrom BiocGenerics counts
#'
#' @examples
#' \dontrun{
#' # Assuming 'pb' is the pseudobulk matrix and 'metadata' is the sample annotation
#'
#' # 1. Auto-select top 2000 most variable genes after VST
#' pb_norm_auto <- deseq2_normalize(pb, metadata)
#'
#' # 2. Use a specific design formula to account for a 'Group' column in metadata
#' my_design <- ~Group
#' pb_norm_design <- deseq2_normalize(
#'   pb,
#'   metadata,
#'   deseq2_design = my_design
#' )
#'
#' # 3. Use pre-defined set of highly variable genes with a specific design
#' my_hvgs <- c("Gene1", "Gene2", "Gene3")
#' pb_norm_hvg <- deseq2_normalize(
#'   pb,
#'   metadata,
#'   deseq2_design = ~ Condition + Batch,
#'   hvg = my_hvgs
#' )
#' }
deseq2_normalize <- function(pb,
                             metadata = NULL,
                             deseq2_design = NULL,
                             hvg = NULL,
                             nvar_genes = 2000) {
  if (is.null(metadata)) {
    metadata <- data.frame(
      # The column name is arbitrary, but required for the design formula (~ 1)
      dummy_group = factor(rep("A", ncol(pb))),
      # The row names MUST match the sample names in the count matrix (pb)
      row.names = colnames(pb)
    )
  } else {
    # --- Existing check/reorder for when metadata IS provided ---
    if (!all(colnames(pb) %in% rownames(metadata))) {
      stop("Not all sample names in pb are found in metadata rownames")
    }
    # Reorder metadata to match pb
    metadata <- metadata[colnames(pb), , drop = FALSE]
  }

  if (is.null(deseq2_design)) {
    deseq2_design <- formula(paste("~ 1"))
  }

  suppressMessages({
    suppressWarnings({
      # Create DESeq2 dataset (metadata is guaranteed to exist here)
      dds <- DESeqDataSetFromMatrix(
        countData = pb,
        colData = metadata,
        design = deseq2_design
      )

      # Estimate size factors
      dds <- estimateSizeFactors(dds)

      # Set minimum number of counts per gene for VST
      nsub <- min(1000, sum(rowMeans(counts(dds, normalized = TRUE)) > 10))

      # Transform counts using variance stabilizing transformation (VST)
      dds <- vst(dds, blind = FALSE, nsub = nsub)
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
        rv <- apply(pb_norm, 1, var)
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


# Stops R CMD check from complaining about "no visible binding for global variable"
# when using unquoted column names inside ggplot2 or dplyr code.
utils::globalVariables(c(
  "sample_id", "celltype", "values", "Variance",
  "cumulative_variance", "rel_abundance", "mean_rel_abund",
  "ordered_group", "x_var", "y_var", "Relative_abundance"
))
