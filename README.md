# Spaniel Spatial Transcriptomics Analysis

Spaniel includes a series of tools to aid the quality control and analysis of Spatial Transcriptomics data. The package contains functions to create a either a Seurat object or SingleCellExperiment from a count matrix and spatial barcode file and provides a method of loading a histologial image into R. The spanielPlot function allows visualisation of metrics contained within the S4 object overlaid onto the image of the tissue.

# Installation

Spaniel can be installed from Bioconductor. Bioconductor version: Development (3.10)

## Install dependencies:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("Spaniel")

# View Vignette Online:

https://bioconductor.org/packages/devel/bioc/vignettes/Spaniel/inst/doc/spaniel-vignette.html




### Spaniel - with 10X import option

The development version of Spaniel (version 1.2) designed to import data from a 10X Genomics Visium experiment.  
This version will be tested and pushed to Bioconductor. In the meantime, if you would like to test the features described in this vignette you can install a development Spaniel (version 1.2) using the following command:


```{r install_dev, eval = FALSE}

devtools::install_github("RachelQueen1/Spaniel", ref = "Development" )


```

You can read more about the features in the development version here:

https://github.com/RachelQueen1/Spaniel/tree/Development


