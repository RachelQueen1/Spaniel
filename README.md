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
