## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(warning=FALSE, message=FALSE, include = TRUE,
                        fig.height = 8, fig.width = 8, fig.align = "center",
                        echo=TRUE
                        )


## ----install, eval = FALSE-----------------------------------------------
#  BiocManager::install('Spaniel')

## ---- load libraries-----------------------------------------------------
library(Spaniel)
library(Seurat)
library(SingleCellExperiment)

## ----counts--------------------------------------------------------------
### read in test data
counts <- readRDS(file.path(system.file(package = "Spaniel"), 
                            "extdata/counts.rds"))


## ----colnames------------------------------------------------------------
colnames(counts)[1:10]

## ----rownames------------------------------------------------------------
rownames(counts)[1:10]

## ----barcodes------------------------------------------------------------

barcodesFile <- file.path(system.file(package = "Spaniel"), 
                            "1000L2_barcodes.txt")

## ----barcodesTop---------------------------------------------------------
barcodes <- read.csv(barcodesFile, sep = "\t", header = FALSE)
head(barcodes)

## ----createSeurat--------------------------------------------------------

seuratObj <- createSeurat(counts,
                        barcodesFile, 
                        projectName = "TestProj", 
                        sectionNumber = 1
)


## ----createSeurat_meta---------------------------------------------------
head(seuratObj[[]])

## ----createSeurat_counts-------------------------------------------------
GetAssayData(seuratObj, "counts")[1:10, 1:5]

## ----createSeurat_project------------------------------------------------
Project(seuratObj)

## ----readSCE-------------------------------------------------------------
sce <- createSCE(counts = counts,
                barcodeFile = barcodesFile, 
                projectName = "TestProj", 
                sectionNumber = 1)

## ----readSCE_ColData-----------------------------------------------------
head(colData(sce)[1:5,1:5])

## ----readSCE_counts------------------------------------------------------
counts(sce)[1:10, 1:5]

## ----readSCE_Project-----------------------------------------------------
colData(sce)$project[1]

## ---- fig.show='hold'----------------------------------------------------
### Load histological image into R

imgFile <- file.path(system.file(package = "Spaniel"), 
                        "HE_Rep1_resized.jpg")

image <- parseImage(imgFile)



## ---- qcplotting,  results = "hide"--------------------------------------

minGenes <- 280
minUMI <- 67500

filter <- seuratObj$nCount_RNA > minUMI &
            seuratObj$nFeature_RNA > minGenes 

spanielPlot(object = seuratObj, 
        grob = image, 
        plotType = "NoGenes", 
        showFilter = filter)


## ---- filter_seurat------------------------------------------------------
seuratFiltered <- subset(x = seuratObj, subset = nFeature_RNA > minGenes &
                                nCount_RNA > minUMI)


spanielPlot(object = seuratFiltered, 
        grob = image, 
        plotType = "NoGenes") 



## ---- select_spots, eval = FALSE-----------------------------------------
#  selectSpots(seuratFiltered, image)
#  

## ---- remove_spots-------------------------------------------------------

spotsToRemove <- file.path(system.file(package = "Spaniel"), 
                            "points_to_remove.txt")

seuratFiltered <- removeSpots(seuratFiltered, 
                                pointsToRemove = spotsToRemove)


spanielPlot(object = seuratFiltered, 
        grob = image, 
        plotType = "NoGenes") 




## ---- find_clusters, message=FALSE, warning=FALSE, echo=TRUE, results = "hide"----



seuratFiltered <- NormalizeData(object = seuratFiltered,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000)

seuratFiltered <- FindVariableFeatures(object = seuratFiltered,
                                        selection.method = "vst",
                                        nfeatures = 2000)

all.genes <- rownames(x = seuratFiltered)


seuratFiltered <- ScaleData(object = seuratFiltered, features = all.genes)
seuratFiltered <- RunPCA(object = seuratFiltered,
                            features = VariableFeatures(object = seuratFiltered)
                        )


seuratFiltered <- FindNeighbors(object = seuratFiltered, dims = 1:10)
seuratFiltered <- FindClusters(object = seuratFiltered,
                                resolution = c(0.4, 0.5, 0.6, 0.8))



## ---- genePlot, results = "hide"-----------------------------------------


gene = "Nrgn"
spanielPlot(object = seuratFiltered, grob = image, 
        plotType = "Gene", 
        gene = gene)



## ---- clusterPlot, warning= FALSE, message = FALSE, results = "hide"-----

spanielPlot(object = seuratFiltered, 
        grob = image, 
        plotType = "Cluster", 
        clusterRes = "RNA_snn_res.0.8"
)


## ---- markClusters-------------------------------------------------------
seuratFiltered <- markClusterCol(seuratFiltered, "res")

## ---- eval = FALSE,  echo=TRUE-------------------------------------------
#  saveRDS(seuratFiltered, "data.rds")

## ---- eval = FALSE, echo=TRUE--------------------------------------------
#  file.path(system.file(package = "Spaniel"), "extdata/SeuratData.rds" )

## ---- eval = FALSE, echo=TRUE--------------------------------------------
#  saveRDS(image, "image.rds")
#  

## ---- eval = FALSE-------------------------------------------------------
#  file.path(system.file(package = "Spaniel"), "extdata/image.rds" )

## ---- markclusterCols, eval = FALSE--------------------------------------
#  runShinySpaniel()
#  

## ---- eval = FALSE-------------------------------------------------------
#  spanielApp <- file.path(system.file(package = "Spaniel"), "ShinySpaniel" )

## ---- eval = FALSE-------------------------------------------------------
#  library(rsconnect)
#  rsconnect::deployApp(spanielApp)
#  

