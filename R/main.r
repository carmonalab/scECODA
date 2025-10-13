library(dplyr)
library(ggplot2)
library(gtools)
library(pheatmap)
library(rlang) # for !!sym() and as_label()
library(tidyr)


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
#' @importFrom base inherits is.null stop as.data.frame names
#' @importFrom Seurat Seurat
#' @importFrom SingleCellExperiment SingleCellExperiment assay
#' @importFrom SummarizedExperiment colData
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
#' @importFrom base is.null any all stop rownames
#' @importFrom methods new
#' @importFrom gtools mixedsort
#' @importFrom dplyr %>%
#' @importFrom stats dist as.matrix
#' @import calc_freq
#' @import clr
#' @import find_highly_variable_celltypes
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




# Find highly variable cell types ####


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
#' @importFrom base sum length
#' @import get_celltype_variances
#' @import get_highly_variable_celltypes
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
#' @importFrom base as.data.frame print mean var sum
#' @importFrom dplyr %>% everything group_by summarize arrange desc mutate
#' @importFrom tidyr pivot_longer
#' @import varmeanplot
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
#' @importFrom base is.null warning length which.max min
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
#' @importFrom base paste print
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




# Plotting functions ####

## Box and bar plots ####

#' Reshapes ECODA data into a long format for plotting and analysis.
#'
#' This function takes either the relative abundance (\code{freq}) or CLR-transformed
#' abundance (\code{clr}) matrix from an \code{\link{ECODA}} object, converts it
#' from a wide (samples x cell types) to a long (sample, celltype, value) format,
#' and optionally joins it with a specified column from the sample metadata.
#'
#' @param ecoda_object An initialized \code{\link{ECODA}} object.
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
#' @importFrom base as.data.frame stop ifelse slot slotNames is.null
#' @importFrom dplyr %>% left_join select
#' @importFrom tidyr pivot_longer everything
#' @importFrom tibble rownames_to_column
#' @importFrom rlang sym
#'
#' @export create_long_data
#'
#' @seealso \code{\link{ECODA}}
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
  data_df <- as.data.frame(data_matrix) %>%
    rownames_to_column("sample_id")

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
      rownames_to_column("sample_id") %>%
      select(sample_id, !!sym(label_col))

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
#' @param ecoda_object An \code{\link{ECODA}} object containing cell type relative
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
#' @importFrom base paste stop unique as.factor factor all match.arg is.null
#' @importFrom dplyr %>% group_by summarise mutate distinct arrange pull
#' @importFrom ggplot2 ggplot aes geom_col theme_minimal theme element_text labs facet_grid
#' @importFrom rlang sym
#' @importFrom gtools mixedsort
#'
#' @export plot_freq_barplot
#'
#' @seealso \code{\link{create_long_data}}, \code{\link{ECODA}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'ecoda_obj' is a created ECODA object with metadata
#'
#' # 1. Plot frequency for every sample, ordered by sample ID:
#' p1 <- plot_freq_barplot(ecoda_obj)
#'
#' # 2. Plot frequency for every sample, faceted and ordered by 'Treatment' column:
#' p2 <- plot_freq_barplot(ecoda_obj, label_col = "Treatment", plot_by = "sample")
#'
#' # 3. Plot average frequency aggregated by 'Treatment' group:
#' p3 <- plot_freq_barplot(ecoda_obj, label_col = "Treatment", plot_by = "group")
#'
#' # 4. Plot with a custom order:
#' custom_order <- c("S2", "S1", "S4", "S3")
#' p4 <- plot_freq_barplot(ecoda_obj, custom_sample_order = custom_order)
#' }
plot_freq_barplot <- function(ecoda_object,
                              label_col = NULL,
                              plot_by = c("sample", "group"),
                              custom_sample_order = NULL,
                              title = "",
                              facet_by_label_col = TRUE) {
  # Use the helper function to get the long data from @freq
  plot_data <- create_long_data(ecoda_object, data_slot = "freq", label_col = label_col)

  plot_by <- match.arg(plot_by, c("sample", "group"))

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
#' @param ecoda_object An \code{\link{ECODA}} object containing the CLR-transformed
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
#' @param signif_label Character string (default: \code{"signif_label"}). Controls
#'                     how p-values are displayed for 2-group comparisons (e.g.,
#'                     "signif_label" for stars, "p.format" for numeric p-value).
#'
#' @return A \code{ggplot} object (enhanced by \code{ggpubr}) representing the
#'         CLR abundance boxplot, including jittered data points and dynamic
#'         significance information (pairwise stars for 2 groups, overall p-value
#'         for 3+ groups).
#'
#' @importFrom base as.data.frame factor is.null unique length paste warning
#' @importFrom ggplot2 aes geom_jitter labs theme element_text guides position_jitterdodge
#' @importFrom ggpubr ggboxplot stat_compare_means
#' @importFrom stringr str_to_title
#' @importFrom rlang sym
#'
#' @export plot_clr_boxplot
#'
#' @seealso \code{\link{create_long_data}}, \code{\link{ECODA}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'ecoda_obj' is a created ECODA object with metadata
#'
#' # 1. Boxplots for CLR abundance without grouping (no stats calculated):
#' p1 <- plot_clr_boxplot(ecoda_obj)
#'
#' # 2. Boxplots grouped by 'Treatment' (2 groups) and applying Wilcoxon test:
#' p2 <- plot_clr_boxplot(
#'   ecoda_obj,
#'   label_col = "Treatment",
#'   stat_method = "wilcox.test",
#'   title = "CLR Abundance by Treatment Group"
#' )
#'
#' # 3. Boxplots grouped by 'Dose' (3+ groups) - automatically uses Kruskal-Wallis:
#' p3 <- plot_clr_boxplot(
#'   ecoda_obj,
#'   label_col = "Dose",
#'   stat_method = "wilcox.test", # will be overridden by kruskal.test
#'   title = "CLR Abundance by Dose Group (Kruskal-Wallis)"
#' )
#' }
plot_clr_boxplot <- function(ecoda_object,
                             label_col = NULL,
                             title = "",
                             stat_method = "wilcox.test",
                             paired = FALSE,
                             signif_label = "signif_label") {
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



## Heatmap ####

#' Generates a Heatmap of CLR-transformed Cell Type Abundances.
#'
#' This function visualizes the CLR-transformed abundance matrix after mean-centering,
#' allowing for clustering of both cell types and samples, and includes
#' a sample annotation sidebar based on a specified metadata column.
#' It is important to not re-scale in order to avoid amplifying tiny differences.
#'
#' @param ecoda_object An \code{\link{ECODA}} object containing the CLR-transformed
#'                     abundances in the \code{clr} slot.
#' @param label_col Character string. The name of the column in \code{ecoda_object@metadata}
#'                  to use for annotating the samples (columns) of the heatmap.
#' @param cluster_rows Logical (default: \code{TRUE}). Whether to apply hierarchical
#'                     clustering to the cell types (rows).
#' @param cluster_cols Logical (default: \code{TRUE}). Whether to apply hierarchical
#'                     clustering to the samples (columns).
#' @param scale Character string (default: \code{"none"}). Method for scaling the
#'              CLR abundance values within the heatmap. Options include "none",
#'              "row", or "column". Note: The data is mean-centered before passing
#'              to \code{pheatmap}, but further scaling is controlled here.
#' @param clustering_method Character string (default: \code{"ward.D2"}). The clustering
#'                          method to use for hierarchical clustering. Options are
#'                          passed directly to \code{hclust} (e.g., "complete", "average", "ward.D2").
#' @param angle_col Character string (default: \code{"90"}). Angle of the sample
#'                  labels (columns).
#' @param ... Additional arguments passed directly to the \code{pheatmap} function.
#'
#' @return A \code{pheatmap} object, which is typically visualized automatically
#'         when called interactively.
#'
#' @importFrom base as.data.frame scale t as.factor lapply
#' @importFrom dplyr %>%
#' @importFrom pheatmap pheatmap
#'
#' @export plot_heatmap
#'
#' @seealso \code{\link{ECODA}}, \code{\link{pheatmap}}
#'
#' @examples
#' \dontrun{
#' # Assuming 'ecoda_obj' is a created ECODA object with metadata
#'
#' # Generate a heatmap clustered by both cell types and samples,
#' # annotated by the 'Condition' column in the metadata:
#' p1 <- plot_heatmap(
#'   ecoda_obj,
#'   label_col = "Condition",
#'   main = "CLR Abundance Heatmap",
#'   fontsize = 8
#' )
#'
#' # Generate a heatmap without clustering the samples:
#' p2 <- plot_heatmap(
#'   ecoda_obj,
#'   label_col = "Batch",
#'   cluster_cols = FALSE
#' )
#' }
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

#' @title Plot Principal Component Analysis and Calculate Clustering Scores
#'
#' @description
#' Performs Principal Component Analysis (PCA) on the CLR-transformed abundance
#' matrix (\code{ecoda_object@clr}) and visualizes the results in 2D or 3D.
#' It can also calculate and display several metrics to evaluate the separation
#' of groups defined by \code{label_col}.
#'
#' @details
#' The clustering metrics (ARI, Modularity, Silhouette, ANOSIM) assess how well
#' the sample groupings (\code{labels}) align with the underlying data structure
#' in the feature space.
#'
#' @param ecoda_object An \code{\link{ECODA}} object containing the CLR-transformed
#'                     abundances in the \code{clr} slot.
#' @param label_col Character string (optional, default: \code{NULL}). The name of a
#'                  column in \code{ecoda_object@metadata} used to color and group
#'                  samples in the plot, and for calculating clustering scores.
#' @param scale. Logical (default: \code{FALSE}). A value indicating whether the
#'               variables should be scaled to have unit variance before the analysis.
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
#' @param pointsize Numeric (default: \code{3}). Size of the points in the plot.
#' @param labelsize Numeric (default: \code{4}). Size of the variable labels in the plot.
#' @param coord_equal Logical (default: \code{TRUE}). If \code{TRUE}, forces the
#'                    aspect ratio of the plot to be equal.
#' @param axes Numeric vector (default: \code{c(1, 2)}). The principal components
#'             to plot (e.g., \code{c(1, 2)} for PC1 vs PC2).
#' @param plotly_3d Logical (default: \code{FALSE}). If \code{TRUE}, generates an
#'                  interactive 3D scatter plot using \code{plotly} (requires \code{pca_dims >= 3}).
#' @param invisible Character vector (default: \code{c("var", "quali")}). Elements to
#'                    hide in the 2D plot. Can include "var" (variables/cell types),
#'                    "ind" (samples), or "quali" (group centroids).
#' @param n_ct_show Integer (default: \code{Inf}). Number of cell types (variables)
#'                  to show based on their contribution to the selected axes. Set to
#'                  \code{Inf} to show all.
#' @param repel Logical (default: \code{TRUE}). Whether to use \code{ggrepel} to
#'                prevent label overlap for variable names.
#'
#' @return A \code{ggplot} object (2D, via \code{factoextra}) or a \code{plotly}
#'         object (3D) visualizing the PCA results.
#'
#' @importFrom base as.data.frame seq_len nrow round
#' @importFrom stats prcomp dist hclust cutree
#' @importFrom utils as.data.frame
#' @importFrom vegan anosim
#' @importFrom factoextra fviz_pca
#' @importFrom plotly plot_ly add_markers add_paths
#' @importFrom dplyr bind_rows
#' @importFrom gplots group2NA
#'
#' @export plot_pca
plot_pca <- function(ecoda_object,
                     label_col = NULL,
                     scale. = FALSE,
                     pca_dims = NULL,
                     knn_k = NULL,
                     title = NULL,
                     score_digits = 3,
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
      anosim_score <- round(anosim(x = feat_mat, grouping = labels, distance = "euclidean")[["statistic"]], score_digits)
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
    if (!is.infinite(n_ct_show) & all(invisible %in% c("var", "quali"))) {
      invisible <- "quali"
    }

    p <- fviz_pca(res.pca,
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
#' @importFrom base as.data.frame as.numeric factor mean round
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
#' @importFrom base unique length round max
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
#' @importFrom base as.matrix
#'
#' @export compute_KNN
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
#' @importFrom base matrix seq_len nrow intersect
#'
#' @export compute_snn_graph
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




#' Calculate Adjusted Rand Index (ARI)
#'
#' Calculates the Adjusted Rand Index (ARI) to measure the agreement between the
#' known cluster assignments (\code{labels}) and the cluster assignments derived
#' from two unsupervised clustering methods: Hierarchical Clustering (\code{hclust})
#' and Partitioning Around Medoids (\code{pam}).
#'
#' @param matrix Numeric matrix or data frame. The feature matrix (e.g., CLR abundances).
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
#' @importFrom base as.numeric as.factor mean unlist list round
#'
#' @export calc_ari
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
  hc <- hclust(dist_mat, method = "ward.D2")
  clust_labels <- cutree(hc, k = nclusts)
  results[["hclust_accuracy"]] <- adjustedRandIndex(as.numeric(as.factor(labels)), clust_labels)

  # Perform PAM clustering
  clust_labels <- pam(matrix, k = nclusts)$cluster
  results[["pamclust_accuracy"]] <- adjustedRandIndex(as.numeric(as.factor(labels)), clust_labels)

  if (return_mean) {
    return(round(mean(unlist(results)), digits))
  } else {
    results[["hclust_accuracy"]] <- round(results[["hclust_accuracy"]], digits)
    results[["pamclust_accuracy"]] <- round(results[["pamclust_accuracy"]], digits)
    return(results)
  }
}
