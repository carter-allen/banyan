# banyan

The `banyan` package allows for network-based Bayesian analysis of spatially resolved single cell data. The core model is an extension of the stochastic block model (SBM) to integration of spatial and gene expression cell-cell similarity network. The package allows for robust detection of cell sub-populations and interaction among them using a Bayesian probabilistic framework. 

## Installation

The `banyan` package can be installed directly from this repository using `devtools`.

```
devtools::install_github("carter-allen/banyan")
```

## Requirements

The `banyan` package has multiple dependencies, one of which being `mlsbm`. The `mlsbm` package contains the core implementation of the Bayesian SBM written with `Rcpp`. You can install the development version of `mlsbm` from GitHub, or the stable version from CRAN using:

```
# development version
devtools::install_github("carter-allen/mlsbm")

# stable release
install.packages("mlsbm")
```

The `igraph` package is used to handle standard network routines. We also allow for integration with standard `Seurat` workflows by accommodating Seurat objects as inputs.

## Usage

Below are `R` commands for analyzing a publicly available anterior mouse brain data set published by [10X Genomics](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Mouse_Brain_Sagittal_Anterior).

```
# packages
library(SeuratData)
library(Seurat)
library(banyan)
library(patchwork)

# load data
InstallData("stxBrain")
brain <- LoadData("stxBrain","anterior1")

# run Seurat processing steps
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)

# fit banyan model
brain_fit <- fit_banyan(seurat_obj = brain, K = 6)
save(brain_fit,file = paste0(data_dir,"brain_fit_banyan.RData"))
load(paste0(data_dir,"brain_fit_banyan.RData"))

# plot labels
plot_labels(brain_fit)

# calculate scores
brain_fit <- get_scores(brain_fit)

# plot uncertainty
plot_uncertainty(brain_fit)

# plot propensity
plot_propensity(brain_fit,k = 1) + plot_propensity(brain_fit,k = 2)

# plot connectivity matrix
plot_connectivity_matrix(brain_fit)

# plot connectivity interval
plot_connectivity_intervals(brain_fit) 
```