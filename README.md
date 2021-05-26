# GRANet

Prediction of Gene Regulatory Networks, and quantification of transcription factor activity, using scRNA-seq and scATAC-seq data.

GRANet combines chromatin accessibility and gene expression data to infer cell state specific regulons (cssRegulons) which are made of one transcription factor and its target genes.

## How to install GRANet

### Install Python dependencies  

GRANet uses **GRNBoost2** to compute co-expression modules from the scRNA-seq data, in order to make GRNBoost2 available for R, we suggest to install it in a conda environment.

Additionally part of the pre-processing steps on the scATAC-seq data uses the package **ArchR** which needs **macs2**, which can be install in the same environment as well.

``` r
# Create conda environment 
reticulate::conda_create(
  envname = "granet",
  packages = c("numpy", "macs2"),
  python_version = "3.7"
)

# indicate that we want to use a specific condaenv
reticulate::use_condaenv("granet")

# Install additional packages
reticulate::conda_install(
  envname = "granet",
  packages = "cytoolz",
  channel = "anaconda",
)

# Install pyscenic  with pip
reticulate::conda_install(
  envname = "granet",
  packages = "pyscenic",
  pip = TRUE
)


```

### Install GRANet

GRANet can be installed by cloning a local copy of this repository and intalling the package form source: 

``` r
install.packages("path/to/GRANet", repos = NULL, type="source")
library(GRANet)
```



