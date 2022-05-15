
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SELINA
<!-- badges: start -->

<!-- badges: end -->

SELINA is a deep learning-based framework for single cell assignment
with multiple references. The algorithm consists of three main steps:
cell type balancing, pre-training and fine-tuning. The rare cell types
in reference image.pngdata are first oversampled using SMOTE(Synthetic Minority
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
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("SELINA-team/SELINA.R")
```

## Usage

### Preprocess of query data

You could preprocess query data with steps in
[SELINA.py](https://github.com/SELINA-team/SELINA.py#preprocess-of-query-data). You will get rds object for each traning dataset.

### Pre-training of the reference data

Train model with `train_model`. You will get a list, which includes a training model and it's meta information.

Files used in here are included in folder `demo`. You can check
parameter details with command `?train_model`.

``` r
library(SELINA)

## Read in train data (seurat object) with meta information (`Celltype` and `Platform`). 
train_rds1 <- readRDS('demo/reference_data/train1_res.rds')
train_rds2 <- readRDS('demo/reference_data/train2_res.rds')
seuratlist = list(train_rds1, train_rds2)
model = train_model(seuratlist)

## You can save the model you trained.
save_model(model, path_out, prefix)
```

In this step, two output files will be generated in the path_out
folder.  
1\. `pre-trained_params.pt` : a file containing all parameters of the
trained model  
2\. `pre-trained_meta.rds` : a file containing the cell types and genes
of the reference data

### Predict

This is a basic example which shows you how to do annotation.  
Files used in here are included in folder `demo`. You can check
parameter details with command `?query_predict`.

``` r
library(SELINA)

## If you predict directly after training, then can pass the next load model step.
model <- read_model(path_model)

## Predict query cell with SELINA
queryObj = readRDS(path_query)
query_predict(queryObj,
              model,
              path_out,
              outprefix = 'demo', 
              disease = FALSE, 
              mode = 'single',
              cell_cutoff = 5,
              prob_cutoff = 0.9)
```

This step will output eight files in the predict\_output folder. Note
that if the input is a cluster level matrix, then the last four files
will not be generated.  
1\. `demo_predictions.txt` : predicted cell type for each cell in the
query data.  
2\. `demo_probability.txt` : probability of cells predicted as each of
the reference cell types.  
3\. `demo_pred.png` : umap plot with cell type annotations image.  
4\. `demo_DiffGenes.tsv` : matrix representing the differentially
expressed genes for each cell type, this file can be used to validate
the annotation results.  
5\. `demo_prob.txt` : prediction probability of each cell, each row
represents one cell. Clusters with lower prediction probability may
contain novel cell types that are not included in the reference.  
6\. `demo_cluster_prob.png` : box plot indicating the prediction
probability distribution of cells in each cluster image.  
7\. `demo_unknown_percent.txt` : the percentage of cells that are
assigned unknown in each cluster, each row represents one cluster.
Clusters with a higher unknown percentage may contain novel cell types
that are not included in the reference.  
8\. `demo_unknown_percent.png` : bar plot indicating the percentage of
cells that are assigned unknown in each cluster image.
