# Define ECODA object to store data
methods::setClass(
  Class = "ECODA",
  slots = list(
    data = "list",
    metadata = "data.frame|NULL"
  ),
  # Define the default values for each slot.
  prototype = list(
    data = list(
      cell_counts = data.frame(),
      frequency = data.frame(),
      frequency_imputed = data.frame(),
      clr = data.frame()
    ),
    metadata = NULL
  )
)




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
    dplyr::summarise(across(everything(), ~n_distinct(.x))) %>%
    dplyr::ungroup()

  # Step 2: Identify columns that are constant across all samples
  constant_cols <- distinct_counts %>%
    dplyr::select(where(~max(.x) == 1)) %>%
    names()

  # Step 3: Select and return the constant columns
  metadata <- cell_data_df %>%
    dplyr::select(all_of(c(sample_col, constant_cols))) %>%
    dplyr::distinct() %>%
    as.data.frame() # Convert back to data.frame

  # Set row names to the sample column for easier access
  rownames(metadata) <- metadata[[sample_col]]

  return(metadata)
}

