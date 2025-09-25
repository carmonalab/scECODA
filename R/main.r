library(dplyr)
library(tidyr)
library(rlang) # for !!sym() and as_label()


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
create_ecoda_object <- function(
    data,
    sample_col = NULL,
    celltype_col = NULL,
    metadata = NULL,
    long_df = FALSE
) {
  # Initialize the object with default values
  ecoda_object <- methods::new("ECODA")

  # --- A) Handle input from a Seurat or SingleCellExperiment object ---
  if (inherits(data, "Seurat") | inherits(data, "SingleCellExperiment")) {
    if (is.null(sample_col) || is.null(celltype_col)) {
      stop("For Seurat/SingleCellExperiment objects, please provide 'sample_col' and 'celltype_col'.")
    }

    # Get metadata from the object and convert to a data frame
    if (inherits(data, "Seurat")) {
      cell_data_df <- as.data.frame(data@meta.data)
    } else if (inherits(data, "SingleCellExperiment")) {
      cell_data_df <- as.data.frame(data@colData)
    } else if (is.data.frame(data) & long_df) {
      cell_data_df <- data
    }

    # Get cell type counts
    cell_counts <- get_celltype_counts(cell_data_df, sample_col, celltype_col)

    # Extract sample metadata
    if (is.null(metadata)) {
      metadata <- get_sample_metadata(cell_data_df, sample_col)
    }

    # --- B) Handle input from a pre-calculated count matrix ---
  } else if (is.matrix(data) || is.data.frame(data)) {

    cell_counts <- as.data.frame(data)

  } else {
    stop("Invalid input type. 'data' must be a Seurat or SingleCellExperiment object, a data frame, or a count matrix/data frame.")
  }

  # Sort by rownames
  cell_counts <- cell_counts[order(rownames(cell_counts)), ]
  metadata <- metadata[order(rownames(metadata)), ]

  # Ensure correct dimensions
  if (!all(rownames(cell_counts) == rownames(metadata))) {
    warning("Rownames of cell counts do not match metadata. Please check. Alternatively, supply metadata manually.")
  }

  ecoda_object@data$cell_counts <- cell_counts
  ecoda_object@metadata <- metadata

  return(ecoda_object)
}








#' Description lorem ipsum
#'
#' @param
#'
#' @return Description lorem ipsum
#' @export my_fun_lorem_ipsum

my_fun <- function(param1) {

  return()
}





#' Description lorem ipsum
#'
#' @param
#'
#' @return Description lorem ipsum
#' @export my_fun_lorem_ipsum

my_fun <- function(param1) {

  return()
}
