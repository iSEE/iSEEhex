
<!-- README.md is generated from README.Rmd. Please edit that file -->

# iSEEhex

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/iSEE/iSEEhex)](https://github.com/iSEE/iSEEhex/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/iSEE/iSEEhex)](https://github.com/iSEE/iSEEhex/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check-bioc](https://github.com/iSEE/iSEEhex/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/iSEE/iSEEhex/actions)
[![Codecov test
coverage](https://codecov.io/gh/iSEE/iSEEhex/branch/main/graph/badge.svg)](https://app.codecov.io/gh/iSEE/iSEEhex?branch=main)
<!-- badges: end -->

The goal of `iSEEhex` is to provide panels summarising data points in
hexagonal bins for
*[iSEE](https://bioconductor.org/packages/3.15/iSEE)*.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `iSEEhex` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("iSEEhex")
```

And the development version from
[GitHub](https://github.com/iSEE/iSEEhex) with:

``` r
BiocManager::install("iSEE/iSEEhex")
```

## Example

For demonstration, we prepare an example
*[SingleCellExperiment](https://bioconductor.org/packages/3.15/SingleCellExperiment)*
object.

``` r
library(scRNAseq)

# Example data ----
sce <- ReprocessedAllenData(assays="tophat_counts")
class(sce)
#> [1] "SingleCellExperiment"
#> attr(,"package")
#> [1] "SingleCellExperiment"

library(scater)
sce <- logNormCounts(sce, exprs_values="tophat_counts")

sce <- runPCA(sce, ncomponents=4)
sce <- runTSNE(sce)
rowData(sce)$ave_count <- rowMeans(assay(sce, "tophat_counts"))
rowData(sce)$n_cells <- rowSums(assay(sce, "tophat_counts") > 0)
sce
#> class: SingleCellExperiment 
#> dim: 20816 379 
#> metadata(2): SuppInfo which_qc
#> assays(2): tophat_counts logcounts
#> rownames(20816): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
#> rowData names(2): ave_count n_cells
#> colnames(379): SRR2140028 SRR2140022 ... SRR2139341 SRR2139336
#> colData names(23): NREADS NALIGNED ... passes_qc_checks_s sizeFactor
#> reducedDimNames(2): PCA TSNE
#> mainExpName: endogenous
#> altExpNames(1): ERCC
```

Then, we create an *[iSEE](https://bioconductor.org/packages/3.15/iSEE)*
app that compares the `ReducedDimensionHexPlot` panel – defined in this
package – to the standard `ReducedDimensionPlot` defined in the
*[iSEE](https://bioconductor.org/packages/3.15/iSEE)* package.

``` r
library(iSEEhex)
#> Loading required package: iSEE
initialPanels <- list(
    ReducedDimensionPlot(
        ColorBy = "Feature name", ColorByFeatureName = "Cux2", PanelWidth = 6L),
    ReducedDimensionHexPlot(
        ColorBy = "Feature name", ColorByFeatureName = "Cux2", PanelWidth = 6L,
        BinResolution = 30)
)
app <- iSEE(se = sce, initial = initialPanels)
```

## Citation

Below is the citation output from using `citation('iSEEhex')` in R.
Please run this yourself to check for any updates on how to cite
**iSEEhex**.

``` r
print(citation('iSEEhex'), bibtex = TRUE)
#> 
#> kevinrue (2022). _Demonstration of a Bioconductor Package_. doi:
#> 10.18129/B9.bioc.MyBioconductorPackage (URL:
#> https://doi.org/10.18129/B9.bioc.MyBioconductorPackage),
#> https://github.com/kevinrue/MyBioconductorPackage/MyBioconductorPackage
#> - R package version 0.99.0, <URL:
#> http://www.bioconductor.org/packages/MyBioconductorPackage>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {Demonstration of a Bioconductor Package},
#>     author = {{kevinrue}},
#>     year = {2022},
#>     url = {http://www.bioconductor.org/packages/MyBioconductorPackage},
#>     note = {https://github.com/kevinrue/MyBioconductorPackage/MyBioconductorPackage - R package version 0.99.0},
#>     doi = {10.18129/B9.bioc.MyBioconductorPackage},
#>   }
#> 
#> kevinrue (2022). "Demonstration of a Bioconductor Package." _bioRxiv_.
#> doi: 10.1101/TODO (URL: https://doi.org/10.1101/TODO), <URL:
#> https://www.biorxiv.org/content/10.1101/TODO>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Demonstration of a Bioconductor Package},
#>     author = {{kevinrue}},
#>     year = {2022},
#>     journal = {bioRxiv},
#>     doi = {10.1101/TODO},
#>     url = {https://www.biorxiv.org/content/10.1101/TODO},
#>   }
```

Please note that the `iSEEhex` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `iSEEhex` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

-   Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*, and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.15/BiocCheck)*.
-   Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
-   The [documentation website](http://iSEE.github.io/iSEEhex) is
    automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
-   The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
-   The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.15/biocthis)*.

## Code of Conduct

Please note that the iSEEhex project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.
