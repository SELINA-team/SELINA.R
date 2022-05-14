
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SELINA

<!-- badges: start -->

<!-- badges: end -->

SELINA is a deep learning-based framework for single cell assignment
with multiple references. The algorithm consists of three main steps:
cell type balancing, pre-training and fine-tuning. The rare cell types
in reference data are first oversampled using SMOTE(Synthetic Minority
Oversampling Technique), and then the reference data is trained with a
supervised deep learning framework using MADA(Multi-Adversarial Domain
Adaptation). An autoencoder is subsquently used to fine-tune the
parameters of the pre-trained model. Finally, the labels from reference
data are transferred to the query data based on the fully-trained model.
Along with the annotation algorithm, we also collect 136 datasets which
were uniformly processed and curated to provide users with comprehensive
pre-trained models.

## Installation

You can install the released version of SELINA from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("SELINA")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("SELINA-team/SELINA.R")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SELINA)
## basic example code
```
