#' scECODA: Single-Cell Exploratory Compositional Data Analysis
#'
#' @description
#' The \code{scECODA} package provides a framework for
#' population-scale exploratory compositional data analysis (ECODA) of
#' single-cell datasets. It enables users to explore sample structures and
#' biological group relationships by transforming cell type counts into Centered
#' Log-Ratios (CLR). The package is designed to handle \code{Seurat} and
#' \code{SingleCellExperiment} objects, as well as raw count or frequency
#' matrices.
#'
#' @details
#' A core philosophy of \code{scECODA} is the identification of
#' "Highly Variable Cell Types" (HVCs). Similar to highly variable genes
#' in transcriptomics, HVCs serve as a proxy for biological relevance,
#' highlighting the cell populations that drive the most variation across
#' a cohort.
#'
#' The package facilitates:
#' \itemize{
#'   \item \strong{Data Transformation}: Automated calculation of CLR and
#'         replacement of zeros.
#'   \item \strong{Dimensionality Reduction}: PCA and 3D PCA specifically
#'         tuned for compositional data.
#'   \item \strong{Quantification}: Metrics like ANOSIM, ARI,
#'         and Modularity to measure group separation.
#'   \item \strong{Visualization}: Integrated plotting for mean-variance
#'         relationships, stacked bar charts, and clustered heatmaps.
#'   \item \strong{Pseudobulk Integration}: Convenient tools to generate
#'         and normalize gene-level pseudobulk data via DESeq2.
#' }
#'
#' @section Key functions:
#' \itemize{
#'   \item \code{\link{ecoda}}: The primary constructor to create an
#'         \code{ECODA} object from single-cell objects or data frames.
#'   \item \code{\link{find_hvcs}}: Identifies cell types driving the most
#'         variance across samples.
#'   \item \code{\link{plot_pca}}: Visualizes sample clusters and calculates
#'         separation scores (ANOSIM, ARI, etc.).
#'   \item \code{\link{plot_heatmap}}: Generates clustered heatmaps with
#'         metadata annotations.
#'   \item \code{\link{plot_varmean}}: Visualizes the mean-variance
#'         distribution of cell types.
#' }
#'
#' @author
#' \strong{Maintainer}: Christian Halter \email{scecoda@gmail.com}
#'
#' @seeAlso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/carmonalab/scECODA}
#'   \item Report bugs at \url{https://github.com/carmonalab/scECODA/issues}
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
