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

