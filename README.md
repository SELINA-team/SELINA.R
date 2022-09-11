
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

- SELINA is available for macOS, Linux and Windows and has been tested on linux system with several dependency packages listed in DESCRIPTION file, all of which will be installed simultaneously using the following commands except for some python packages need to be prepared before. 
    - numpy=1.22.3
    - pandas=1.3.4
    - datatable=0.11.1
    - imbalanced-learn=0.8.1 

If you have gpu on your device and want to use it, you should install [cudatoolkit](https://developer.nvidia.com/cuda-downloads) and [cudnn](https://developer.nvidia.com/rdp/cudnn-archive) based on your system version before SELINA installation.

We recommend you install SELINA with `devtools::install_github()` from
R:

``` r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("SELINA-team/SELINA.R")
```

## Usage

### Preprocess of query data

You could preprocess query data with steps in
[SELINA.py](https://github.com/SELINA-team/SELINA.py#preprocess-of-query-data). 

### Pre-training of the reference data

Train model with `train_model`. You will get a list, which includes a training model and it's meta information. Files used in here are included in folder `demos`. You can check
parameter details with command `?train_model`.

*_NOTE:_* Please put expression and meta files in one folder,  meta file should include `Celltype` and `Platform` columns.
``` r
library(SELINA)

model <- train_model(path_in="demos/disease_data/reference_data",
                     disease=TRUE)
```
Save the model with `save_model`.
```r
save_model(model, path_out, prefix)
```

In this step, two output files will be generated in the path_out
folder.  
1\. `pre-trained_params.pt` : a file containing all parameters of the
trained model.  
2\. `pre-trained_meta.rds` : a file containing the cell types and genes
of the reference data.

### Prediction

Annotate query data with `query_predict`. You will get a list, which includes prediction results and corresponding probability for query data.

Files used in here are included in folder `demos`. You can check
parameter details with command `?query_predict`.

SELINA has trained models for 35 kinds of normal tissues and 3 kinds of disease tissues, you can load them with command `?load_selina_model`. All the tissue names are showed in the toggle list below.
<details>
  <summary>Tissue models (Click Me)</summary>
  
1.Normal
* Adrenal-Gland
* Airway-Epithelium
* Artery
* Bladder
* Blood
* Bone-Marrow
* Brain
* Breast
* Choroid
* Decidua
* Esophagus
* Eye
* Fallopian-Tube
* Gall-Bladder
* Heart
* Intestine
* Kidney
* Liver
* Lung
* Muscle
* Nose
* Ovary
* Pancreas
* Peritoneum
* Placenta
* Pleura
* Prostate
* Skin
* Spleen
* Stomach
* Testis
* Thyroid
* Ureter
* Uterus
* Visceral-Adipose

2.Disease
* AD (type II diabetes)
* T2D (non-small-cell lung carcinoma)
* NSCLC (Alzheimer’s disease)
</details>

Load the MADA model.
``` r
library(SELINA)

## If you predict directly after training, then can skip the next load model step.
# If you want to use models trained by yourself:
model <- read_model(path_model)

# If you want to load model SELINA prepared (Please make sure the input tissue name is included in our documentation, eg: Pancreas):
model <- load_selina_model(tissue)
```
Predict with SELINA.
```r
queryObj <- readRDS(path_query)
query_result <- query_predict(query_expr = queryObj,
                              model = model,
                              path_out = path_out,
                              outprefix = 'query', 
                              disease = TRUE, 
                              cell_cutoff = 5,
                              prob_cutoff = 0.9)
```
This step will output eight files in the `path_out` folder. 

- 1\. `query_predictions.txt`: predicted cell type for each cell in the
query data.  
- 2\. `query_probability.txt`: probability of cells predicted as each of
the reference cell types.  
- 3\. `query_pred.png`: umap plot with cell type annotations image.  
- 4\. `query_DiffGenes.tsv`: matrix representing the differentially
expressed genes for each cell type, this file can be used to validate
the annotation results.  

For disease mode, an extra file will be generated:  
- 5\. `query_cellsources.txt`: predicted cell source for each cell in the query data