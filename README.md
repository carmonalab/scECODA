# scECODA

The scECODA R package provides functions for sample representation and ***s***ingle-***c***ell ***E***xploratory ***CO***mpositional ***D***ata ***A***nalysis at the population-scale.

It enables users to explore their dataset by visualizing the sample structure and how groups of samples relate to each other
by providing numerous convenient plotting functions.
Additionally, it provides metrics to quantify how strongly groups of samples separate and
which cell types are driving this separation (differential abundance analysis).

The following is an example a healthy human cohort, comprised of 868 samples, by of the [Gong & Sharma *et al.*](https://pubmed.ncbi.nlm.nih.gov/39314416/) showing a few of the most highly variable cell types driving separation of samples:

<img width="2700" height="2100" alt="GongSharma" src="https://github.com/user-attachments/assets/aa8b34ba-722c-495d-a9f7-3aea92842652" />


### Package Installation

To install `scECODA` directly from the GitHub repository, run the following code from within R or RStudio:

``` r
install.packages("remotes")
library(remotes)

remotes::install_github("carmonalab/scECODA")
```


### scECODA Tutorial

Find a step-by-step tutorial for `scECODA` at: [scECODA tutorial](https://carmonalab.github.io/scECODA_demo/Tutorial.html) ([RMD](https://github.com/carmonalab/scECODA_demo/blob/master/Tutorial.rmd))


### Why scECODA?

scECODA allows for unprecedented insights in biological sample grouping structure:

[**Case Study 1: Granularity Matters**](https://carmonalab.github.io/scECODA_demo/Case_Study_1.html) -
Learn why **fine-grained cell type annotation** is crucial to uncover biological insights missed by broad, low-resolution annotation or pseudobulk gene expression. ([RMD](https://github.com/carmonalab/scECODA_demo/blob/master/Case_Study_1.rmd))

[**Case Study 2: Cell type composition vs. Gene expression**](https://carmonalab.github.io/scECODA_demo/Case_Study_2.html) -
See how scECODA's compositional analysis compares to traditional **pseudobulk** methods, and outperforms it especially when differences are driven by **low-abundance cell types**. ([RMD](https://github.com/carmonalab/scECODA_demo/blob/master/Case_Study_2.rmd))


## References

Longitudinal Multi-omic Immune Profiling Reveals Age-Related Immune Cell Dynamics in Healthy Adults. Gong Q. & Sharma M. *et al.*, bioRxiv [Preprint]. 2024 Sep 14:2024.09.10.612119. [doi: 10.1101/2024.09.10.612119](https://pubmed.ncbi.nlm.nih.gov/39314416/).
