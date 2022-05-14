
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

We recommend you install SELINA.R with `devtools::install_github()` from
R:

``` r
# install.packages("devtools")
devtools::install_github("SELINA-team/SELINA.R")
```

## Usage

### Preprocess of query data

You could preprocess query data with steps in
[SELINA.py](https://github.com/SELINA-team/SELINA.py#preprocess-of-query-data)

### Pre-training of the reference data

This is a basic example which shows you how to train a MADA model with
your own data. Files used in here are included in folder `demo`. You can
check parameter details with command `?train_model`.

``` r
library(SELINA)
train_model(path_in = "demo/reference_data",
            path_out = "train_output",
            prefix = 'pre-trained')
```

In this step, two output files will be generated in the train\_output
folder. 1. `pre-trained_params.pt` : a file containing all parameters of
the trained model  
2\. `pre-trained_meta.pkl` : a file containing the cell types and genes
of the reference data

### Predict

This is a basic example which shows you how to train a MADA model with
your own data. Files used in here are included in folder `demo`. You can
check parameter details with command `?train_model`.

``` r
library(SELINA)
query_predict(query_expr = "demo/query_data/query.txt",
              model = "train_output/pre-trained_params.pt"
              path_out = 'predict_output',
              outprefix = 'demo', 
              disease = FALSE, 
              mode = 'single',
              seurat = 'query_res.rds',
              cell_cutoff = 5,
              prob_cutoff = 0.9)
```

This step will output eight files in the predict\_output folder. Note
that if the input is a cluster level matrix, then the last four files
will not be generated.  
1\. demo\_predictions.txt : predicted cell type for each cell in the
query data. 2. demo\_probability.txt : probability of cells predicted as
each of the reference cell types. 3. demo\_pred.png : umap plot with
cell type annotations image. 4. demo\_DiffGenes.tsv : matrix
representing the differentially expressed genes for each cell type, this
file can be used to validate the annotation results. 5. demo\_prob.txt :
prediction probability of each cell, each row represents one cell.
Clusters with lower prediction probability may contain novel cell types
that are not included in the reference. 6. demo\_cluster\_prob.png : box
plot indicating the prediction probability distribution of cells in each
cluster image. 7. demo\_unknown\_percent.txt : the percentage of cells
that are assigned unknown in each cluster, each row represents one
cluster. Clusters with a higher unknown percentage may contain novel
cell types that are not included in the reference. 8.
demo\_unknown\_percent.png : bar plot indicating the percentage of cells
that are assigned unknown in each cluster image.
