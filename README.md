# scECODA

The scECODA R package faciliates ***Exploratory COmpositional Data Analysis (ECODA)*** and visualization, especially for single-cell omics data embeddings at the sample/population level

It enables intuitive exploration of multi-sample datasets –such as large patient cohorts– and supports the unsupervised identification of samples with similar cell type compositional profiles (e.g., patient clustering or stratification). In addition, scECODA provides metrics to quantify the degree of separation between sample groups (e.g., biological conditions) and to pinpoint the cell types or cell states whose abundance changes drive these differences.

### Example
The following example uses 868 (celltype annotated) scRNA-seq samples from the blood of healthy donors (data from [Gong & Sharma *et al.*](https://pubmed.ncbi.nlm.nih.gov/39314416/)) . It illustrates how samples naturally separate in an unsupervised manner by donor age and CMV infection status, and highlights the top cell types whose changes in abundance drive inter-sample variation.

```r
ecoda_object <- create_ecoda_object(
  seurat_object, # or SCE_object or count data (see create_ecoda_object_helper)
  sample_col = "sample_id", # sample column
  celltype_col = "seurat_annotations" # cell type annotations
)
ecoda_object <- create_ecoda_object(Seurat_Object) # or SCE_object or count data
plot_pca(ecoda_object)
```

<img width="2700" height="2100" alt="GongSharma" src="https://github.com/user-attachments/assets/aa8b34ba-722c-495d-a9f7-3aea92842652" />


### Package Installation

To install `scECODA` directly from the GitHub repository, run the following code from within R or RStudio:

``` r
install.packages("remotes")
library(remotes)

remotes::install_github("carmonalab/scECODA")
```


### Tutorial

Check out our step-by-step [scECODA tutorial](https://carmonalab.github.io/scECODA_demo/Tutorial.html) ([RMD](https://github.com/carmonalab/scECODA_demo/blob/master/Tutorial.rmd))


### Case studies

[**Case Study 1: Granularity Matters**](https://carmonalab.github.io/scECODA_demo/Case_Study_1.html) -
See how **fine-grained cell type annotation** can be crucial to uncover inter-sample biological variation missed by broad, low-resolution annotation or pseudobulk gene expression in these ECODA anlayses of i) blood samples from healthy individuals and ii) lung samples from patients with different pulmonary diseases. ([RMD](https://github.com/carmonalab/scECODA_demo/blob/master/Case_Study_1.rmd))

[**Case Study 2: Cell type composition vs. Pseudo-bulk gene expression**](https://carmonalab.github.io/scECODA_demo/Case_Study_2.html) -
See how scECODA's compositional analysis compares to **pseudobulk** analysis and outperforms it when differences are driven by **low-abundance cell types** in a semi-synthetic dataset. ([RMD](https://github.com/carmonalab/scECODA_demo/blob/master/Case_Study_2.rmd))


## References

Halter C, et al. 2025
