# Spaniel Spatial Transcriptomics Analysis

Spaniel is an R package for analysing Spatial Transcriptomics data.

# Install dependencies:


## Mac Only

xcode-select --install


## All operating systems

### DevTools
install.packages('devtools')

### Seurat
Spaniel requires Seurat v3.0. This can be installed: <br/><br/>

```{r}
install.packages('Seurat')
```

### SingleCellExperiment

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment", version = "3.8")
```

### Scater

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scater", version = "3.8")
```








# Install Spaniel:

```{r}
devtools::install_github(repo = "RachelQueen1/Spaniel")
```

# View Vignette:

https://rachelqueen1.github.io/Spaniel/
