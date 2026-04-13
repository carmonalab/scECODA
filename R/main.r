# Create SummarizedExperiment object ---------------------------

# Define SummarizedExperiment object for ECODA

#' Create an SummarizedExperiment object from various data types
#'
#' This is a smart constructor function used to initialize an
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#' object. It handles the processing of single-cell objects (\code{Seurat} or
#' \code{SingleCellExperiment}) or raw data frames to extract cell type counts,
#' calculate sample metadata, and optionally generate DESeq2-normalized
#' pseudobulk data.
#'
#' @param data The primary input, which can be:
#'   \itemize{
#'     \item A \code{Seurat} object.
#'     \item A \code{SingleCellExperiment} object.
#'     \item A cell type x sample matrix or data frame (counts or frequencies).
#'   }
#' @param data_is_freq Logical (default: \code{FALSE}). If \code{TRUE}, the
#'   input \code{data} is treated as relative frequencies (proportions) rather
#'   than raw counts.
#' @param metadata A data frame containing sample-level metadata. Required if
#'   \code{data} is a matrix/data frame. If \code{data} is a single-cell object,
#'   metadata is extracted automatically.
#' @param sample_col The metadata column name defining unique sample IDs.
#'   Required if \code{data} is a single-cell object.
#' @param celltype_col The metadata column name defining cell type annotations.
#'   Required if \code{data} is a single-cell object.
#' @param get_pb Logical (default: \code{FALSE}). If \code{TRUE}, calculates and
#'   stores DESeq2-normalized pseudobulk data in \code{metadata(se)$pb}
#' @param variance_explained Numeric (default: 0.5). Proportion of variance used
#'   to determine the number of highly variable cell types (HVCs).
#' @param top_n_hvcs Integer (optional). If provided, specifies the exact number
#'   of top HVCs to select, overriding \code{variance_explained}.
#' @param pseudo_count Numeric (default: 0.5). Used if \code{rep_method =
#'   "counts"}.
#' @param pseudo_frac_min Numeric (default: 2/3). The multiplier for the minimum
#'   non-zero value when \code{rep_method = "frac_min"}.
#' @param add_to Character string. Should the value be added to \code{"all"}
#'   cells or just the \code{"zeros"}?
#'
#' @return A new
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   object.
#'
#' @importFrom methods new
#' @importFrom gtools mixedsort
#' @importFrom SummarizedExperiment colData colData<- assay
#' @importFrom S4Vectors metadata metadata<- DataFrame
#'
#' @export ecoda
#'
#' @examples
#' data(example_data)
#' Zhang <- example_data$Zhang
#' counts <- Zhang$cell_counts_lowresolution
#' freq <- calc_freq(counts)
#' metadata <- Zhang$metadata
#'
#' se <- ecoda(data = counts, metadata = metadata)
#'
#' se <- ecoda(data = freq, metadata = metadata)
ecoda <- function(data = NULL,
                  data_is_freq = FALSE,
                  metadata = NULL,
                  sample_col = NULL,
                  celltype_col = NULL,
                  get_pb = FALSE,
                  variance_explained = 0.5,
                  top_n_hvcs = NULL,
                  pseudo_count = 0.5,
                  pseudo_frac_min = 2 / 3,
                  add_to = NULL) {
    # --- A) Handle input from a Seurat or SingleCellExperiment object ---
    if (inherits(data, "Seurat") | inherits(data, "SingleCellExperiment")) {
        if (is.null(sample_col) || is.null(celltype_col)) {
            stop("Please provide 'sample_col' and 'celltype_col'.")
        }

        # Get cell metadata from the object
        if (inherits(data, "Seurat")) {
            cell_data_df <- as.data.frame(slot(data, "meta.data"))
            if (get_pb) {
                if (!is.null(data[["RNA"]]$counts)) {
                    pb <- calculate_pseudobulk(
                        count_matrix = data[["RNA"]]$counts,
                        sample_ids = cell_data_df[[sample_col]]
                    )
                } else {
                    warning("Could not find RNA counts in Seurat object")
                }
            }
        } else if (inherits(data, "SingleCellExperiment")) {
            if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
                stop("Package 'SummarizedExperiment' must be installed.")
            }
            cell_data_df <- as.data.frame(colData(data))
            names(cell_data_df) <- names(colData(data))
            if (get_pb) {
                if (!is.null(assay(data, "counts"))) {
                    pb <- calculate_pseudobulk(
                        count_matrix = assay(
                            data, "counts"
                        ),
                        sample_ids = cell_data_df[[sample_col]]
                    )
                } else {
                    warning(
                        "Could not find RNA counts ",
                        "in SingleCellExperiment object"
                    )
                }
            }
        }

        data <- get_celltype_counts(cell_data_df, sample_col, celltype_col)
        metadata <- get_sample_metadata(cell_data_df, sample_col)
    }

    data <- data[, mixedsort(colnames(data))]

    if (!is.null(metadata)) {
        metadata <- metadata[mixedsort(rownames(metadata)), , drop = FALSE]

        # Ensure correct and matching colnames
        if (!all(colnames(data) == rownames(metadata))) {
            stop(
                "Colnames of cell counts (or freq) ",
                "do not match rownames of metadata."
            )
        }
    }

    se <- ecoda_helper(
        data = data,
        data_is_freq = data_is_freq,
        variance_explained = variance_explained,
        top_n_hvcs = top_n_hvcs,
        pseudo_count = pseudo_count,
        pseudo_frac_min = pseudo_frac_min,
        add_to = add_to
    )
    if (!is.null(metadata)) colData(se) <- DataFrame(metadata)

    if (get_pb && exists("pb")) {
        pb <- deseq2_normalize(pb)
        pb <- pb[, mixedsort(colnames(pb))]
        metadata(se)$pb <- pb
    }

    return(se)
}


#' Core constructor for SummarizedExperiment objects from count/frequency
#' matrices
#'
#' This is the internal engine that initializes an
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#' object. It performs zero-imputation, Centered Log-Ratio (CLR) transformation,
#' calculates sample distances, and identifies highly variable cell types
#' (HVCs).
#'
#' @param data A matrix or data frame where **columns are samples** and
#'   **rows are cell types**.
#' @param data_is_freq Logical. If \code{TRUE}, \code{data} is treated as
#'   proportions (0-1 or 0-100). If \code{FALSE}, treated as raw integer counts.
#' @param variance_explained Numeric (default: 0.5). Target variance for HVCs.
#' @param top_n_hvcs Integer (optional). Number of top HVCs to select.
#' @param pseudo_count Numeric (default: 0.5). Used if \code{rep_method =
#'   "counts"}.
#' @param pseudo_frac_min Numeric (default: 2/3). The multiplier for the minimum
#'   non-zero value when \code{rep_method = "frac_min"}.
#' @param add_to Character string. Should the value be added to \code{"all"}
#'   cells or just the \code{"zeros"}?
#'
#' @return A fully initialized
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   object.
#'
#' @importFrom methods new
#' @importFrom dplyr %>% mutate across everything
#' @importFrom stats dist
#' @importFrom SummarizedExperiment SummarizedExperiment assay assay<-
#' @importFrom S4Vectors metadata metadata<-
ecoda_helper <- function(data = NULL,
                         data_is_freq,
                         variance_explained = 0.5,
                         top_n_hvcs = NULL,
                         pseudo_count = 0.5,
                         pseudo_frac_min = 2 / 3,
                         add_to = c("all", "zeros")) {
    rep_method <- ifelse(data_is_freq, "frac_min", "counts")
    add_to <- match.arg(add_to)

    data_imp <- data
    if (any(data == 0)) {
        if (data_is_freq) {
            message(
                "Frequencies contain zeros. Zeros were replaced with the (",
                pseudo_frac_min, ") * smallest value."
            )
        } else {
            message(
                "Counts contain zeros. A pseudo count of +", pseudo_count,
                " was added to ", add_to, " counts."
            )
        }

        data_imp <- replace_zeros(
            df = data,
            rep_method = rep_method,
            pseudo_count = pseudo_count,
            pseudo_frac_min = pseudo_frac_min,
            add_to = add_to
        )
    }

    if (data_is_freq) {
        freq <- data
        freq_imp <- data_imp
    } else {
        freq <- calc_freq(data)
        freq_imp <- calc_freq(data_imp)
    }

    # Check if all column sums are close to 100
    # (using a tolerance for floating point numbers)
    # If they are intended to be percentages
    colw_sums <- colSums(freq)
    if (!all(abs(colw_sums - 100) < 1e-6)) {
        # Check if they are close to 1 and should be scaled to 100
        if (all(abs(colw_sums - 1) < 1e-6)) {
            message(
                "Frequencies sum close to 1. Rescaling all rows to sum to 100."
            )
            # Scale by 100
            freq <- freq * 100
            freq_imp <- freq_imp * 100
        } else {
            # If they sum to neither 1 nor 100, warn and scale them to 100
            message(
                "Frequencies do not sum to 100 (or 1). ",
                "Each row will be scaled so the new row sum is 100."
            )
            # Re-scale each row so it sums to 100: (freq / row_sum) * 100
            freq <- (freq / colw_sums) * 100
            freq_imp <- (freq_imp / colSums(freq_imp)) * 100
        }
    }

    if (data_is_freq) {
        se <- SummarizedExperiment(
            assays = list(
                freq = freq,
                freq_imp = freq_imp
            )
        )
    } else {
        se <- SummarizedExperiment(
            assays = list(
                counts = data,
                counts_imp = data_imp,
                freq = freq,
                freq_imp = freq_imp
            )
        )
    }

    clr_df <- calc_clr(freq_imp)
    assay(se, "clr") <- clr_df
    metadata(se)$sample_distances <- clr_df %>%
        t() %>%
        dist()

    assay(se, "asin_sqrt") <- asin(sqrt(as.matrix(freq) / 100))

    se <- find_hvcs(
        se,
        variance_explained,
        top_n_hvcs
    )

    metadata(se)$clr_hvc <-
        calc_clr(freq_imp[metadata(se)$hvcs, , drop = FALSE])

    return(se)
}


#' Calculate relative frequencies (percentages) column-wise.
#'
#' This function takes a matrix or data frame of counts and transforms each
#' column (assumed to be a sample or observation) into relative frequencies,
#' expressing each component as a percentage of the column's total sum.
#'
#' @param df A data frame or matrix of counts (samples/observations as columns,
#'   components as rows). Must contain only numeric, non-negative values.
#' @return A data frame of the same dimensions, where each column sums to 100.
#' @export calc_freq
#' @examples
#' counts <- data.frame(A = c(10, 90), B = c(50, 50))
#' calc_freq(counts)
calc_freq <- function(df) {
    df <- sweep(df, 2, colSums(df), "/") * 100
    return(as.data.frame(df))
}


#' Perform the Centered Log-Ratio (CLR) transformation.
#'
#' The CLR transformation is a common technique used for compositional data
#' analysis, mapping the data from the simplex to Euclidean space. It is defined
#' as the log of the ratio between each component and the geometric mean of all
#' components in that sample. **Note:** This function assumes the input data is
#' strictly positive (i.e., does not contain zeros). Use an imputation or
#' pseudocount method prior to this function if zeros are present.
#'
#' @param df A data frame or matrix of strictly positive relative frequencies or
#'   counts (samples/observations as columns, components as rows).
#' @return A data frame of the same dimensions containing the CLR-transformed
#'   values.
#' @export calc_clr
#' @examples
#' freq_imp <- data.frame(A = c(10.1, 89.9), B = c(50.1, 49.9))
#' calc_clr(freq_imp)
calc_clr <- function(df) {
    log_df <- log(df)
    clr_df <- sweep(log_df, 2, colMeans(log_df), "-")

    return(as.data.frame(clr_df))
}


#' Replace zero values in count or frequency data
#'
#' This function replaces zero values in a data frame or matrix to facilitate
#' downstream transformations that require strictly positive values, such as the
#' Centered Log-Ratio (CLR) transformation.
#'
#' @details The function provides two methods for handling zeros:
#' \itemize{
#'   \item \strong{counts:} Adds a fixed value (\code{pseudo_count}).
#'   \item \strong{frac_min:} Replaces zeros with a fraction
#'   (\code{pseudo_frac_min}) of the smallest observed non-zero value in
#'   the dataset.
#' }
#'
#' @param df A data frame or matrix where zeros need to be replaced.
#' @param rep_method Character string specifying the replacement strategy. One
#'   of \code{"counts"}, or \code{"frac_min"}.
#' @param pseudo_count Numeric (default: 0.5). Used if \code{rep_method =
#'   "counts"}.
#' @param pseudo_frac_min Numeric (default: 2/3). The multiplier for the minimum
#'   non-zero value when \code{rep_method = "frac_min"}.
#' @param add_to Character string. Should the value be added to \code{"all"}
#'   cells or just the \code{"zeros"}?
#'
#' @return A data frame or matrix of the same dimensions as \code{df} with all
#'   zero values replaced by the calculated replacement value.
#'
#' @export replace_zeros
#'
#' @examples
#' # Replace zeros in a count matrix
#' counts_df <- data.frame(A = c(10, 0, 5), B = c(20, 10, 0))
#' replace_zeros(counts_df, rep_method = "counts", add_to = "all")
#'
#' # Replace zeros in a frequency matrix
#' freq_df <- data.frame(A = c(0.5, 0.5), B = c(0.2, 0.8), C = c(0.0, 1.0))
#' replace_zeros(freq_df, rep_method = "frac_min", add_to = "zeros")
replace_zeros <- function(df,
                          rep_method = c("counts", "frac_min"),
                          pseudo_count = 0.5,
                          pseudo_frac_min = 2 / 3,
                          add_to = NULL) {
    rep_method <- match.arg(rep_method)
    if (is.null(add_to)) {
        add_to <- if (rep_method == "counts") "all" else "zeros"
    } else {
        add_to <- match.arg(add_to, c("all", "zeros"))
    }

    imp_val <- switch(rep_method,
        "counts"   = pseudo_count,
        "frac_min" = min(df[df > 0], na.rm = TRUE) * pseudo_frac_min
    )

    if (add_to == "all") {
        df <- df + imp_val
    } else {
        df[df == 0] <- imp_val
    }

    return(df)
}



#' Get the cell type counts from a long data frame (e.g. seurat object metadata)
#' where each cell is a row.
#'
#' @param cell_data_df A data frame where each row represents a single cell, and
#'   columns contain cell-level information, including sample ID and cell type
#'   annotation.
#' @param sample_col The column that defines the sample ID for each cell
#' @param celltype_col The column that defines the cell type annotation for each
#'   cell
#'
#' @return A data frame with cell types as rows and samples as columns,
#'   containing the count of each cell type per sample.
#' @export get_celltype_counts
#' @examples
#' # Create example data frame
#' cell_data_df <- data.frame(
#'     Cell_ID = paste0("C", seq(10)),
#'     Sample_Name = c(rep("S1", 5), rep("S2", 5)),
#'     Cluster_Annotation = factor(c(
#'         "B_Cell", "T_Cell", "B_Cell", NA, "T_Cell",
#'         "T_Cell", "Macrophage", "B_Cell", "T_Cell", "T_Cell"
#'     ))
#' )
#' #
#' # Calculate cell type counts per sample
#' celltype_counts <- get_celltype_counts(
#'     cell_data_df = cell_data_df,
#'     sample_col = "Sample_Name",
#'     celltype_col = "Cluster_Annotation"
#' )
#' #
#' print(celltype_counts)
#' # Note how the NA cell type count is handled and renamed to "NA".
get_celltype_counts <- function(cell_data_df,
                                sample_col,
                                celltype_col) {
    cellcount_df <- table(
        cell_data_df[[celltype_col]],
        cell_data_df[[sample_col]],
        useNA = "ifany"
    ) %>%
        as.data.frame.matrix()

    return(cellcount_df)
}


#' Extracts constant metadata for each sample from a cell-level data frame.
#'
#' This function identifies columns in a cell-level metadata data frame that
#' have a **constant** value for all cells belonging to the same sample. It
#' aggregates this constant information, returning a new data frame where each
#' row represents a unique sample. Columns that vary within any sample are
#' excluded from the output.
#'
#' @param cell_data_df A data frame containing cell-level metadata. This should
#'   include the sample ID column and all potential metadata columns.
#' @param sample_col A character string specifying the name of the column that
#'   defines the unique sample ID for each cell (e.g., "Sample_ID").
#' @return A new data frame where:
#'         \itemize{
#'           \item Each row corresponds to a unique sample from the input data.
#'           \item The row names are set to the values of the input
#'                 `sample_col`.
#'           \item Columns contain only the metadata fields that were constant
#'                 across all cells within **each** sample. The `sample_col`
#'                 itself is excluded from the final columns but used for
#'                 row names.
#'         }
#' @importFrom dplyr group_by summarise across everything n_distinct ungroup
#'   select where all_of distinct %>%
#' @importFrom rlang sym
#' @export get_sample_metadata
#' @examples
#' # Assuming you have a data frame 'cell_df'
#' cell_df <- data.frame(
#'     Cell_ID = paste0("C", seq(10)),
#'     Sample_ID = c(rep("S1", 5), rep("S2", 5)),
#'     Age = c(rep(30, 5), rep(45, 5)),
#'     Gender = c(rep("M", 5), rep("F", 5)),
#'     Cell_Type = c(rep("A", 3), rep("B", 2), rep("A", 3), rep("B", 2))
#' )
#'
#' # The 'Age' and 'Gender' columns are constant within each sample (S1 and S2).
#' # The 'Cell_ID' and 'Cell_Type' columns vary within sample S1 and/or S2.
#'
#' sample_meta <- get_sample_metadata(cell_df, "Sample_ID")
#' print(sample_meta)
#' # Output will have 'Age' and 'Gender' as columns,
#' # with row names 'S1' and 'S2'.
get_sample_metadata <- function(cell_data_df,
                                sample_col) {
    # Step 1: Group the data by `sample_col` and check for uniqueness
    distinct_counts <- cell_data_df %>%
        group_by(!!sym(sample_col)) %>%
        summarise(across(everything(), ~ n_distinct(.x))) %>%
        ungroup()

    # Step 2: Identify columns that are constant across all samples
    constant_cols <- distinct_counts %>%
        select(where(~ max(.x) == 1)) %>%
        names()

    cols_to_keep <- unique(c(sample_col, constant_cols))

    # Step 3: Select the constant columns and unique rows
    metadata <- cell_data_df %>%
        # Select only the relevant constant columns
        select(all_of(cols_to_keep)) %>%
        distinct()

    rownames(metadata) <- metadata[[sample_col]]

    # Remove the sample column from the final metadata table
    cols_to_return <- setdiff(colnames(metadata), sample_col)

    metadata <- metadata[, cols_to_return, drop = FALSE] %>%
        as.data.frame()

    return(metadata)
}


# Find highly variable cell types ---------------------------

#' Identifies and stores Highly Variable Cell Types (HVCs) in an
#' SummarizedExperiment object.
#'
#' This function calculates the variance of each cell type across all samples
#' based on the Centered Log-Ratio (CLR) transformed data. It then selects the
#' subset of highly variable cell types (HVCs) that collectively explain a
#' specified proportion of the total variance, or a fixed number of top cell
#' types.
#'
#' @param se An initialized
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   object containing the CLR-transformed data in the \code{clr} assay.
#' @param variance_explained Numeric (default: 0.5). The target cumulative
#'   proportion of total variance to be explained by the selected HVCs. The
#'   function stops selecting cell types once this threshold is met.
#' @param top_n_hvcs Integer (optional). If provided, this overrides
#'   \code{variance_explained} and selects exactly the top N cell types with the
#'   highest variance.
#'
#' @importFrom S4Vectors metadata metadata<-
#'
#' @return The updated
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   object with the following metadata populated:
#'         \itemize{
#'           \item \code{celltype_variances}: Data frame of cell type variances.
#'           \item \code{variance_explained}: The exact variance proportion
#'                 captured by the selected HVCs.
#'           \item \code{top_n_hvcs}: The number of HVCs selected.
#'           \item \code{hvcs}: Character vector of the names
#'                 of the selected HVCs.
#'         }
#'
#' @export find_hvcs
#'
#' @seealso
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment},
#'   \code{\link{get_celltype_variances}}, \code{\link{get_hvcs}}
#'
#' @examples
#' data(example_data)
#' se <- ecoda(
#'     data = example_data$GongSharma_full$cell_counts_highresolution,
#'     metadata = example_data$GongSharma_full$metadata
#' )
#'
#' # Select HVCs that explain 75% of the total variance
#' se <- find_hvcs(se, variance_explained = 0.75)
#'
#' # Select exactly the top 5 most variable cell types
#' se <- find_hvcs(se, top_n_hvcs = 5)
find_hvcs <- function(se,
                      variance_explained = 0.5,
                      top_n_hvcs = NULL) {
    df_var <- get_celltype_variances(se)

    hvcs <- get_hvcs(
        df_var,
        variance_explained,
        top_n_hvcs = top_n_hvcs
    )

    variance_explained <- (sum(df_var$Variance[seq_along(hvcs)]) /
        sum(df_var$Variance))

    metadata(se)$celltype_variances <- df_var
    metadata(se)$variance_explained <- variance_explained
    metadata(se)$top_n_hvcs <- length(hvcs)
    metadata(se)$hvcs <- hvcs


    return(se)
}


#' Calculates the variance of cell types across samples.
#'
#' This function takes the Centered Log-Ratio (CLR) transformed cell type
#' abundance data from an
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#' object, calculates the mean CLR abundance and variance for each cell type. It
#' also calculates the cumulative variance explained by the cell types when
#' ranked by variance.
#'
#' @param se An initialized
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   object containing CLR-transformed data in the \code{clr} assay
#' @param descending Logical (default: \code{TRUE}). If \code{TRUE}, the
#'   returned data frame is sorted by variance in descending order.
#'
#' @return A data frame where each row is a cell type, containing:
#'         \itemize{
#'           \item \code{celltype}: The name of the cell type.
#'           \item \code{avg_clr_abundance}: The mean CLR value for that
#'                 cell type across all samples.
#'           \item \code{Variance}: The variance of the CLR values for that
#'                 cell type across all samples.
#'           \item \code{cumulative_variance}: Cumulative sum of the variance,
#'                 ranked by variance.
#'           \item \code{variance_exp}: The proportion of total variance
#'                 explained cumulatively.
#'         }
#'
#' @importFrom stats var
#' @importFrom dplyr %>% everything group_by summarize arrange desc mutate
#' @importFrom tidyr pivot_longer
#' @importFrom SummarizedExperiment assay
#'
#' @export get_celltype_variances
#'
#' @seealso
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment},
#'   \code{\link{plot_varmean}}
#'
#' @examples
#' data(example_data)
#' se <- ecoda(
#'     data = example_data$GongSharma_full$cell_counts_highresolution,
#'     metadata = example_data$GongSharma_full$metadata
#' )
#'
#' # Calculate variances without plotting, sorted by ascending variance
#' df_var_asc <- get_celltype_variances(
#'     se,
#'     descending = FALSE
#' )
get_celltype_variances <- function(se, descending = TRUE) {
    clr_data <- assay(se, "clr")

    df_var <- data.frame(
        celltype = rownames(clr_data),
        avg_clr_abundance = rowMeans(clr_data, na.rm = TRUE),
        Variance = apply(clr_data, 1, var, na.rm = TRUE)
    )

    if (descending) {
        df_var <- df_var[order(-df_var$Variance), , drop = FALSE]
    } else {
        df_var <- df_var[order(df_var$Variance), , drop = FALSE]
    }

    total_variance <- sum(df_var$Variance, na.rm = TRUE)
    df_var$cumulative_variance <- cumsum(df_var$Variance)
    df_var$variance_exp <- df_var$cumulative_variance / total_variance

    rownames(df_var) <- NULL

    return(df_var)
}


#' Selects Highly Variable Cell Types (HVCs) based on variance or count
#' threshold.
#'
#' This is a utility function that selects the most variable cell types from a
#' variance-ranked data frame. Selection can be controlled either by defining
#' the cumulative variance to be explained or by specifying a fixed number of
#' cell types.
#'
#' @param df_var A data frame (typically the output of
#'   \code{\link{get_celltype_variances}}) that must contain at least the
#'   columns \code{celltype} and \code{variance_exp} (cumulative variance
#'   explained) and must be sorted by variance in descending order.
#' @param variance_explained Numeric (default: 0.5). The cumulative proportion
#'   of total variance that the selected HVCs must explain. This parameter is
#'   ignored if \code{top_n_hvcs} is set.
#' @param top_n_hvcs Integer or Numeric (optional).
#'                   \itemize{
#'                     \item If an integer (>= 1), it specifies the exact
#'                           number of top cell types to select.
#'                     \item If a numeric value between 0 and 1, it specifies
#'                           the proportion of cell types to select
#'                           (e.g., 0.1 for the top 10%).
#'                   }
#'
#' @return A character vector containing the names of the selected highly
#'   variable cell types. The function ensures that at least two cell types are
#'   always returned.
#'
#' @importFrom dplyr %>% slice pull slice_head
#' @importFrom rlang .data
#'
#' @export get_hvcs
#'
#' @seealso \code{\link{find_hvcs}}, \code{\link{get_celltype_variances}}
#'
#' @examples
#' data(example_data)
#' ecoda_object <- ecoda(
#'     data = example_data$GongSharma_full$cell_counts_highresolution,
#'     metadata = example_data$GongSharma_full$metadata
#' )
#'
#' # Calculate variances
#' df_var <- get_celltype_variances(ecoda_object)
#'
#' # Select cell types explaining at least 60% of variance:
#' hvcs_60 <- get_hvcs(df_var, variance_explained = 0.6)
#'
#' # Select the top 10 cell types:
#' hvcs_top10 <- get_hvcs(df_var, top_n_hvcs = 10)
#'
#' # Select the top 5% of cell types:
#' hvcs_prop <- get_hvcs(df_var, top_n_hvcs = 0.05)
get_hvcs <- function(df_var,
                     variance_explained = 0.5,
                     top_n_hvcs = NULL) {
    if (!is.null(top_n_hvcs)) {
        message("Setting top_n_hvcs overrules variance_explained parameter")
        if (top_n_hvcs < 1) {
            top_hvcs <- df_var %>%
                slice_head(prop = top_n_hvcs) %>%
                pull(.data$celltype)
        } else {
            top_hvcs <- df_var %>%
                slice_head(n = top_n_hvcs) %>%
                pull(.data$celltype)
        }
    } else {
        last_hvc_index <- which.max(df_var$variance_exp >= variance_explained)

        if (length(last_hvc_index) == 0 || last_hvc_index == 0) {
            warning(
                "Variance threshold was not met for any cell type ",
                "(or data is empty). Returning top 2."
            )
            return(df_var$celltype[seq_len(min(2, nrow(df_var)))])
        }

        top_hvcs <- df_var %>%
            slice(seq_len(last_hvc_index)) %>%
            pull(.data$celltype)
    }

    # Select at least two cell types
    if (length(top_hvcs) < 2) {
        message("The minimum number of highly variable cell types to use is 2")
        top_hvcs <- df_var %>%
            slice_head(n = 2) %>%
            pull(.data$celltype)
    }

    return(top_hvcs)
}


#' @title Generates a Mean-Variance Plot for CLR-transformed Cell Type Data.
#'
#' @description This function visualizes the relationship between the mean
#'   abundance (CLR) and the variance (CLR) for each cell type, which is
#'   typically used to identify Highly Variable Cell Types (HVCs).
#'
#' @details If \code{highlight_hvcs} is \code{TRUE}, cell types previously
#'   identified and stored in \code{metadata(se)$hvcs} will be highlighted in
#'   red on the plot.
#'
#' @param se A
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'   object containing pre-calculated cell type variances in
#'   \code{metadata(se)$celltype_variances} and the HVC list in
#'   \code{metadata(se)$hvcs}
#' @param plot_title Character string (default: ""). The title for the plot.
#' @param highlight_hvcs Logical (default: \code{TRUE}). If \code{TRUE}, the
#'   points corresponding to the Highly Variable Cell Types (HVCs) stored in the
#'   SummarizedExperiment object are colored red.
#' @param labels Character (default: "only_hvc"). Options: "all" (label every
#'   point), "none" (no labels), or "only_hvc" (label only the highly variable
#'   cell types).
#' @param plot_fit_line Logical (default: \code{FALSE}). If \code{TRUE}, a
#'   smoothing regression line is added to the plot.
#' @param smooth_method Character string (default: "lm"). The smoothing method
#'   to use for the regression line if \code{plot_fit_line} is \code{TRUE}
#'   (e.g., "lm" for linear model, "loess" for local regression).
#'
#' @return A \code{ggplot} object representing the mean-variance plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_classic xlab
#'   ylab scale_color_manual theme
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr mutate if_else filter
#' @importFrom S4Vectors metadata
#'
#' @export plot_varmean
#'
#' @seealso
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment},
#'   \code{\link{get_celltype_variances}}
#'
#' @examples
#' data(example_data)
#' se <- ecoda(
#'     data = example_data$GongSharma_full$cell_counts_highresolution,
#'     metadata = example_data$GongSharma_full$metadata
#' )
#'
#' # 1. Generate the plot, highlighting HVCs (default)
#' plot_varmean(se, plot_title = "HVCs on Mean-Variance Plot")
#'
#' # 2. Generate the plot without highlighting HVCs and add a fit line
#' plot_varmean(se,
#'     highlight_hvcs = FALSE,
#'     plot_fit_line = TRUE,
#'     smooth_method = "loess"
#' )
plot_varmean <- function(se,
                         plot_title = "",
                         highlight_hvcs = TRUE,
                         labels = c("only_hvc", "all", "none"),
                         plot_fit_line = FALSE,
                         smooth_method = "lm") {
    labels <- match.arg(labels)
    df_var <- metadata(se)$celltype_variances
    highlight_celltypes <- metadata(se)$hvcs

    # --- 1. Create a highlighting factor column ---
    if (highlight_hvcs) {
        df_var <- df_var %>%
            mutate(
                is_highlighted = if_else(
                    .data$celltype %in% highlight_celltypes,
                    "HVC selected",
                    "Not selected"
                )
            )
        color_map <- c("HVC selected" = "red", "Not selected" = "black")

        p <- ggplot(
            df_var,
            aes(
                x = .data$avg_clr_abundance,
                y = .data$Variance,
                color = .data$is_highlighted
            )
        ) +
            scale_color_manual(values = color_map, name = "Cell Type Group")
    } else {
        p <- ggplot(
            df_var, aes(x = .data$avg_clr_abundance, y = .data$Variance)
        )
    }

    p <- p +
        geom_point() +
        labs(title = paste(plot_title)) +
        theme_classic() +
        xlab("Mean (CLR)") +
        ylab("Variance (CLR)")

    if (plot_fit_line) {
        p <- p +
            geom_smooth(
                method = smooth_method,
                color = "red",
                fill = "#69b3a2",
                se = TRUE
            )
    }

    # Ensure the legend is not shown if we are not highlighting anything
    if (!highlight_hvcs) p <- p + theme(legend.position = "none")

    if (labels != "none") {
        # Filter data for labeling based on user choice
        label_df <- df_var
        if (labels == "only_hvc") {
            label_df <- df_var %>%
                dplyr::filter(.data$celltype %in% highlight_celltypes)
        }

        # Use the color aesthetic inside ggrepel so
        # highlighted labels match the point color
        if (highlight_hvcs) {
            p <- p +
                geom_text_repel(
                    data = label_df,
                    aes(label = .data$celltype, color = .data$is_highlighted),
                    size = 3,
                    show.legend = FALSE
                )
        } else {
            p <- p +
                geom_text_repel(
                    data = label_df,
                    aes(label = .data$celltype),
                    size = 3,
                    color = "black"
                )
        }
    }

    return(p)
}


# Get Pseudobulk and normalize ---------------------------

#' Calculate Pseudobulk from Count Matrix
#'
#' This function aggregates single-cell count data into a "pseudobulk" matrix by
#' summing the counts for all cells belonging to the same sample ID. It is
#' robust to both dense and sparse count matrices. It also includes filtering
#' logic to exclude samples that do not meet a minimum cell count threshold.
#'
#' @param count_matrix A gene x cell count matrix (can be a dense matrix or
#'   sparse Matrix). Gene identifiers should be row names and cell barcodes
#'   should be column names.
#' @param sample_ids A vector of sample identifiers, one for each column (cell)
#'   in \code{count_matrix}. The length must equal \code{ncol(count_matrix)}.
#'   Must not contain \code{NA} values.
#' @param min_cells Minimum number of cells required per sample (default: 1).
#'   Samples with fewer cells than this threshold will be excluded from the
#'   final pseudobulk matrix.
#'
#' @return A gene x sample pseudobulk count matrix. The columns correspond to
#'   the unique sample IDs, and the rows correspond to the genes.
#' @export calculate_pseudobulk
#' @examples
#' nrow <- 100
#' ncol <- 100
#' vals <- nrow * ncol
#' mat <- round(matrix(exp(rpois(vals, lambda = 3)), nrow = nrow, ncol = ncol))
#' rownames(mat) <- paste0("Gene", seq_len(ncol))
#' colnames(mat) <- paste0("Cell", seq_len(nrow))
#' ids <- rep(c("Sample1", "Sample2"), each = nrow / 2)
#' pb <- calculate_pseudobulk(count_matrix = mat, sample_ids = ids)
calculate_pseudobulk <- function(count_matrix,
                                 sample_ids,
                                 min_cells = 10) {
    # Input validation
    if (ncol(count_matrix) != length(sample_ids)) {
        stop(
            "Length of sample_ids must equal ",
            "the number of columns in count_matrix"
        )
    }

    if (any(is.na(sample_ids))) {
        stop(
            "sample_ids contains NA values. ",
            "Please remove or replace missing values."
        )
    }

    # Convert sample_ids to factor for grouping
    sample_ids <- as.factor(sample_ids)

    # Count cells per sample
    cells_per_sample <- table(sample_ids)

    # Filter samples with insufficient cells
    if (min_cells > 1) {
        valid_samples <- names(cells_per_sample)[cells_per_sample >= min_cells]

        if (length(valid_samples) == 0) {
            stop("No samples have >= ", min_cells, "cells")
        }

        # Subset to valid samples
        keep_cells <- sample_ids %in% valid_samples
        count_matrix <- as.matrix(count_matrix)[, keep_cells, drop = FALSE]
        sample_ids <- droplevels(sample_ids[keep_cells])

        # Report filtered samples
        n_filtered <- length(cells_per_sample) - length(valid_samples)
        if (n_filtered > 0) {
            message(
                "Filtered out ", n_filtered,
                " sample(s) with < ", min_cells, " cells"
            )
        }
    }

    # Aggregate by summing across samples
    pb <- t(rowsum(t(count_matrix), group = sample_ids))

    return(pb)
}


#' DESeq2 Normalization of Pseudobulk Data
#'
#' This function normalizes a gene x sample pseudobulk count matrix using the
#' Variance Stabilizing Transformation (VST) from the \code{DESeq2} package. It
#' estimates size factors and variance for all genes, performs the VST, and then
#' subsets the results to include only highly variable genes (HVCs), either
#' specified by the user or automatically selected based on variance.
#'
#' **Note:** If `deseq2_design` is not specified, the VST is performed using a
#' minimal design formula (`~ 1`) for size factor and variance estimation.
#'
#' @param pb A gene x sample pseudobulk count matrix (Genes as rows, Samples as
#'   columns). Must contain non-negative integer counts.
#' @param metadata A data.frame with sample-level metadata. Row names must
#'   exactly match the column names of \code{pb}. The function will reorder the
#'   metadata to match \code{pb}.
#' @param deseq2_design Optional: A \code{formula} object specifying the design
#'   used by \code{DESeq2} to estimate size factors and variance for the VST.
#'   This formula should reference columns in the \code{metadata} data frame. If
#'   \code{NULL} (the default), the minimal design \code{~ 1} is used. Note: the
#'   design will be considered for pseudobulk normalization as the function uses
#'   vst(blind = FALSE, ...) internally.
#' @param hvg Optional character vector of gene names to use as highly variable
#'   genes. If provided, only these genes will be returned after VST.
#' @param nvar_genes Number of top variable genes to select (default: 2000).
#'   Only used if \code{hvg = NULL}; the function will select the top
#'   \code{nvar_genes} by variance after VST.
#'
#' @return A normalized expression matrix (VST-transformed) with\
#'         **samples as columns** and **genes as rows**.
#'
#' @export deseq2_normalize
#'
#' @importFrom stats formula
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors vst
#' @importFrom BiocGenerics counts
#' @importFrom SummarizedExperiment assay
#'
#' @examples
#' nrow <- 100
#' ncol <- 100
#' vals <- nrow * ncol
#' mat <- round(matrix(exp(rpois(vals, lambda = 3)), nrow = nrow, ncol = ncol))
#' rownames(mat) <- paste0("Gene", seq_len(ncol))
#' colnames(mat) <- paste0("Cell", seq_len(nrow))
#' ids <- rep(c("Sample1", "Sample2"), each = nrow / 2)
#' pb <- calculate_pseudobulk(count_matrix = mat, sample_ids = ids)
#' pb_norm <- deseq2_normalize(pb)
deseq2_normalize <- function(pb,
                             metadata = NULL,
                             deseq2_design = NULL,
                             hvg = NULL,
                             nvar_genes = 2000) {
    if (is.null(metadata)) {
        metadata <- data.frame(
            # Add dummy column required for the design formula (~ 1)
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

    if (is.null(deseq2_design)) deseq2_design <- formula(paste("~ 1"))

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
    pb_norm <- assay(dds) # Genes x Samples format

    # Select highly variable genes
    if (!is.null(hvg)) {
        # Use user-provided HVGs
        message("Using", length(hvg), "user-provided highly variable genes")

        hvg_available <- hvg[hvg %in% rownames(pb_norm)]
        hvg_missing <- setdiff(hvg, hvg_available)

        if (length(hvg_missing) > 0) {
            warning(
                length(hvg_missing),
                " genes not found in the data and will be excluded"
            )
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

        message("Selected top ", nrow(pb_norm), " highly variable genes")
    }

    return(as.data.frame(pb_norm))
}
