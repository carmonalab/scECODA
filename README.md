# scECODA

<p align="center">
<img width="154" height="154" alt="image" src="https://github.com/user-attachments/assets/ecd4f6c8-de4b-433c-b0f7-75aa2e37dee7" />
</p>

This R package faciliates ***Exploratory COmpositional Data Analysis (ECODA)*** and visualization for single-cell omics sample embeddings at the population level

It enables intuitive exploration of multi-sample datasets –such as large patient cohorts– and supports the unsupervised identification of samples with similar cell type compositional profiles, e.g. patient clustering. In addition, scECODA provides metrics to quantify the degree of separation between groups of samples, e.g. biological conditions, and to pinpoint the cell types or states whose change in abundance drives these differences.

### Package Installation

To install `scECODA` directly from the GitHub repository, run the following code from within R or RStudio:

``` r
install.packages("remotes")
library(remotes)

remotes::install_github("carmonalab/scECODA")
```

### Example
The following example uses 868 scRNA-seq samples from the blood of healthy donors (data from [Gong & Sharma *et al.*](https://pubmed.ncbi.nlm.nih.gov/39314416/)) with previously annotated cell types. It illustrates how samples naturally separate in an unsupervised manner by donor age and CMV infection status, and highlights the top cell types whose changes in abundance drive inter-sample variation.

```r
ecoda_object <- ecoda(
  seurat_object, # or SingleCellExperiment object or count data
  sample_col = "sample_id", # Metadata column containing sample annotation for each cell
  celltype_col = "celltype_annotations" # Metadata column containing cell type annotations
)
plot_pca(ecoda_object)
```
See also the [Tutorial](https://github.com/carmonalab/scECODA/tree/main?tab=readme-ov-file#tutorial) below.

<img width="2700" height="2100" alt="GongSharma" src="https://github.com/user-attachments/assets/aa8b34ba-722c-495d-a9f7-3aea92842652" />

Code to reproduce this figure: [https://github.com/carmonalab/scECODA/blob/main/data-raw/Create_readme_figure.rmd](https://github.com/carmonalab/scECODA/blob/main/data-raw/Create_readme_figure.rmd)


### Tutorial

Check out our step-by-step [scECODA tutorial](https://carmonalab.github.io/scECODA_demo/Tutorial.html) ([RMD](https://github.com/carmonalab/scECODA_demo/blob/master/Tutorial.rmd))


### Case studies

[**Case Study 1: Granularity Matters**](https://carmonalab.github.io/scECODA_demo/Case_Study_1.html) -
See how **fine-grained cell type annotation** can be crucial to uncover inter-sample biological variation missed by broad, low-resolution annotation or pseudobulk gene expression in these ECODA anlayses of i) blood samples from healthy individuals and ii) lung samples from patients with different pulmonary diseases. ([RMD](https://github.com/carmonalab/scECODA_demo/blob/master/Case_Study_1.rmd))

[**Case Study 2: Cell type composition vs. Pseudo-bulk gene expression**](https://carmonalab.github.io/scECODA_demo/Case_Study_2.html) -
See how scECODA's compositional analysis compares to **pseudobulk** analysis and outperforms it when differences are driven by **low-abundance cell types** in a semi-synthetic dataset. ([RMD](https://github.com/carmonalab/scECODA_demo/blob/master/Case_Study_2.rmd))


## References

Halter C, et al. 2025
