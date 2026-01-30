#' Example Data for scECODA
#'
#' A collection of single-cell compositional data from multiple studies,
#' structured to demonstrate the full scECODA workflow.
#'
#' @format A named \code{list} where each element corresponds to a different
#'   study (e.g., \code{$Adams}). Each study element is a \code{list} containing
#'   the following:
#' \itemize{
#'    \item \code{cell_counts_highresolution}: Data frame of cell type counts at
#'    high resolution annotation (samples as rows, cell types as columns).
#'    \item \code{cell_counts_lowresolution}: Data frame of cell type counts at
#'    low resolution annotation (samples as rows, cell types as columns).
#'    \item \code{metadata}: Data frame of sample-level metadata (samples as
#'    rows).
#'    \item \code{main_biologicalcondition_columnname}: Character string,
#'    indicating the main biological condition column name in the metadata.
#' }
#' @usage data(example_data)
#' @source Data originally derived from multiple public single-cell datasets and
#'   re-processed as described in the accompanying scECODA manuscript (Halter et
#'   al., *in preparation*).
"example_data"
