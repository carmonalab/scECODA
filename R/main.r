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





#' Creates an ECODA object from various data types.
#'
#' This is a smart constructor function used to initialize an
#' \code{\link{ECODA}} object. It handles the processing of single-cell objects
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
#'                         \code{create_ecoda_object_from_counts} to determine how many
#'                         highly variable cell types (HVCs) to select.
#' @param top_n_hvcs Integer (optional). Overrides \code{variance_explained} if provided,
#'                   specifying the exact number of top HVCs to select.
#'
#' @return A new \code{\link{ECODA}} object populated with \code{counts}, \code{metadata},
#'         and optionally \code{pb} and the initial compositional analysis results.
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom methods new
#' @importFrom utils read.csv
#'
#' @export create_ecoda_object
#'
#' @seealso \code{\link{ECODA}}, \code{\link{calculate_pseudobulk}}, \code{\link{get_celltype_counts}}, \code{\link{get_sample_metadata}}
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
      cell_data_df <- as.data.frame(data@colData)
      names(cell_data_df) <- names(data@colData)
      if (get_pb) {
        pb <- calculate_pseudobulk(
          count_matrix = assay(data, "counts"),
          sample_ids = cell_data_df[[sample_col]]
        )
      }
    }
  } else {
    cell_data_df <- data
  }

  counts <- get_celltype_counts(cell_data_df, sample_col, celltype_col)
  metadata <- get_sample_metadata(cell_data_df, sample_col)

  ecoda_object <- create_ecoda_object_from_counts(
    counts,
    metadata,
    variance_explained,
    top_n_hvcs
  )

  if (get_pb) {
    ecoda_object@pb <- deseq2_normalize(pb, metadata)
  }

  return(ecoda_object)
}





#' Creates an ECODA object from pre-calculated cell type counts.
#'
#' This is the core constructor function that initializes and performs the initial
#' compositional analysis steps for an \code{\link{ECODA}} object, assuming the
#' cell type count matrix is already available. It handles zero-imputation,
#' calculates relative frequencies, Centered Log-Ratio (CLR) transformed data,
#' sample distances, and identifies highly variable cell types (HVCs).
#'
#' @param counts A data frame or matrix of cell type counts where **rows are samples**
#'               and **columns are cell types**. Must contain non-negative integers.
#' @param metadata An optional data frame containing sample-level metadata. Row names
#'                 must match the row names of the \code{counts} matrix.
#' @param variance_explained Numeric (default: 0.5). The proportion of total variance
#'                           that should be captured by the selected highly variable
#'                           cell types (HVCs). Used by \code{find_highly_variable_celltypes}.
#' @param top_n_hvcs Integer (optional). Overrides \code{variance_explained} if provided,
#'                   specifying the exact number of top HVCs to select based on variance.
#'
#' @return A fully initialized \code{\link{ECODA}} object, populated with counts,
#'         frequency data, CLR transformed data, sample distances, and HVCs.
#'
#' @importFrom methods new
#' @importFrom gtools mixedsort
#' @importFrom dplyr %>%
#' @importFrom stats dist
#'
#' @export create_ecoda_object_from_counts
#'
#' @seealso \code{\link{ECODA}}, \code{\link{calc_freq}}, \code{\link{clr}}, \code{\link{find_highly_variable_celltypes}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'my_counts' and 'my_metadata' are defined data frames:
#' counts_df <- data.frame(row.names = c("S1", "S2"), A = c(10, 5), B = c(0, 15))
#' meta_df <- data.frame(row.names = c("S1", "S2"), Group = c("Treated", "Control"))
#'
#' ecoda_obj <- create_ecoda_object_from_counts(
#'   counts = counts_df,
#'   metadata = meta_df,
#'   top_n_hvcs = 2
#' )
#' # ecoda_obj@counts will contain the original counts
#' # ecoda_obj@clr will contain the CLR transformed data
#' }
create_ecoda_object_from_counts <- function(counts = NULL,
                                            metadata = NULL,
                                            variance_explained = 0.5,
                                            top_n_hvcs = NULL) {
  # Initialize the object with default values
  ecoda_object <- new("ECODA")

  if (!is.null(counts)) {
    if (!is.null(metadata)) {
      # Sort by rownames
      counts <- counts[mixedsort(rownames(counts)), ]
      metadata <- metadata[mixedsort(rownames(metadata)), , drop = FALSE]

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
#' @param ecoda_object An initialized \code{\link{ECODA}} object containing the
#'                     CLR-transformed data in the \code{clr} slot.
#' @param variance_explained Numeric (default: 0.5). The target cumulative
#'                           proportion of total variance to be explained by the
#'                           selected HVCs. The function stops selecting cell types
#'                           once this threshold is met.
#' @param top_n_hvcs Integer (optional). If provided, this overrides
#'                   \code{variance_explained} and selects exactly the top N cell
#'                   types with the highest variance.
#'
#' @return The updated \code{\link{ECODA}} object with the following slots populated:
#'         \itemize{
#'           \item \code{celltype_variances}: Data frame of cell type variances.
#'           \item \code{variance_explained}: The exact variance proportion captured by the selected HVCs.
#'           \item \code{top_n_hvcs}: The number of HVCs selected.
#'           \item \code{highly_variable_celltypes}: Character vector of the names of the selected HVCs.
#'         }
#'
#' @export find_highly_variable_celltypes
#'
#' @seealso \code{\link{ECODA}}, \code{\link{get_celltype_variances}}, \code{\link{get_highly_variable_celltypes}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'ecoda_obj' is a previously created ECODA object:
#'
#' # Select HVCs that explain 75% of the total variance
#' ecoda_obj <- find_highly_variable_celltypes(ecoda_obj, variance_explained = 0.75)
#'
#' # Select exactly the top 5 most variable cell types
#' ecoda_obj <- find_highly_variable_celltypes(ecoda_obj, top_n_hvcs = 5)
#' }
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



#' Calculates the variance of cell types across samples.
#'
#' This function takes the Centered Log-Ratio (CLR) transformed cell type
#' abundance data from an \code{\link{ECODA}} object, calculates the mean CLR
#' abundance and variance for each cell type, and optionally generates a
#' mean-variance plot. It also calculates the cumulative variance explained
#' by the cell types when ranked by variance.
#'
#' @param ecoda_object An initialized \code{\link{ECODA}} object containing
#'                     CLR-transformed data in the \code{clr} slot.
#' @param show_plot Logical (default: \code{TRUE}). If \code{TRUE}, a mean-variance
#'                  plot is generated using the internal \code{varmeanplot} function.
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
#' @seealso \code{\link{ECODA}}, \code{\link{varmeanplot}}
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
#' @importFrom dplyr %>% slice pull
#' @importFrom rlang .data
#'
#' @export get_highly_variable_celltypes
#'
#' @seealso \code{\link{find_highly_variable_celltypes}}, \code{\link{get_celltype_variances}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'df_variance' is a variance data frame sorted descendingly by variance.
#' # Select cell types explaining at least 60% of variance:
#' hvcs_60 <- get_highly_variable_celltypes(df_variance, variance_explained = 0.6)
#'
#' # Select the top 10 cell types:
#' hvcs_top10 <- get_highly_variable_celltypes(df_variance, top_n_hvcs = 10)
#'
#' # Select the top 5% of cell types:
#' hvcs_prop <- get_highly_variable_celltypes(df_variance, top_n_hvcs = 0.05)
#' }
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





#' Generates a Mean-Variance Plot for CLR-transformed Cell Type Data.
#'
#' This function visualizes the relationship between the mean abundance (CLR) and
#' the variance (CLR) for each cell type, which is typically used to identify
#' Highly Variable Cell Types (HVCs).
#'
#' @param ecoda_object An \code{\link{ECODA}} object containing pre-calculated
#'                     cell type variances in the \code{celltype_variances} slot.
#' @param plot_title Character string (default: ""). The title for the plot.
#' @param smooth_method Character string (default: "lm"). The smoothing method to
#'                      use for the regression line (e.g., "lm" for linear model, "loess" for local regression).
#' @param label_points Logical (default: \code{TRUE}). If \code{TRUE}, cell type
#'                     names are added to the points using \code{ggrepel::geom_text_repel}
#'                     to minimize overlap.
#' @param highlight_cts Character vector (optional, default: \code{NULL}). A vector
#'                      of cell type names to highlight on the plot (not currently used
#'                      in the function body but kept for future functionality).
#'
#' @return A \code{ggplot} object representing the mean-variance plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_classic xlab ylab
#' @importFrom ggrepel geom_text_repel
#'
#' @export plot_varmean
#'
#' @seealso \code{\link{ECODA}}, \code{\link{get_celltype_variances}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'ecoda_obj' is a created ECODA object with calculated variances:
#'
#' # Generate the standard mean-variance plot
#' plot_varmean(ecoda_obj, plot_title = "Cell Type Variance vs Mean Abundance")
#'
#' # Generate the plot without labeling points
#' plot_varmean(ecoda_obj, label_points = FALSE)
#'
#' # Use a LOESS smoothing method
#' plot_varmean(ecoda_obj, smooth_method = "loess")
#' }
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
    xlab("Mean (CLR)") +
    ylab("Variance (CLR)")

  if (label_points) {
    p <- p + geom_text_repel(data = df_var, aes(label = celltype), vjust = -0.5)
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
