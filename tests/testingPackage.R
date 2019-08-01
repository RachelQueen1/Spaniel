library(Seurat)
library(dplyr)

# rObjects = "../../Spatial_Transcriptomics_Report/Analysis_outs/SRR3382371/rObjects"
# # sObj= readRDS(file.path(rObjects, "seurObj.rds"))
# # imgFile = readRDS(file.path(rObjects, "image.rds"))
# 
# # sObjFiltered = readRDS(file.path(rObjects, "sObj_filtered.rds"))
# countsFiltered = sObjFiltered@raw.data[, colnames(sObjFiltered@scale.data)]
# saveRDS(countsFiltered, "testData/countsFiltered.rds")

counts = readRDS("testData/counts.rds")
countsFiltered = readRDS("testData/countsFiltered.rds")

imgFile = readRDS("testData/image.rds")

barcodesFile = "testData/1000L2_barcodes.txt"
#saveRDS(counts, "testData/counts.rds")

SeuratObj = readSeurat(countsFiltered,
                       barcodesFile, 
                       ProjectName = "TestProj", 
                       SectionNumber = 1
                       )

sce = readSCE(Counts = countsFiltered,
              BarcodeFile = barcodesFile, 
              ProjectName = "TestProj", 
              SectionNumber = 1)




sce = scater::calculateQCMetrics(
    sce
)

sce =  scater::normalise(sce)


saveRDS(sce, "sceData.rds")


minGenes = 2000
minUMI = 300000

# filter = SeuratObj@meta.data$nFeature_RNA > minGenes & 
#     SeuratObj@meta.data$nCount_RNA > minUMI
# 
# 
# 
SeuratFiltered <- subset(x = SeuratObj, subset = nFeature_RNA > minGenes &
                             nCount_RNA > minUMI)


SeuratFiltered <- NormalizeData(object = SeuratFiltered,
                           normalization.method = "LogNormalize",
                           scale.factor = 10000)

SeuratFiltered <- FindVariableFeatures(object = SeuratFiltered,
                                  selection.method = "vst",
                                  nfeatures = 2000)

all.genes <- rownames(x = SeuratObj)
SeuratFiltered <- ScaleData(object = SeuratFiltered, features = all.genes)
SeuratFiltered <- RunPCA(object = SeuratFiltered,
                    features = VariableFeatures(object = SeuratFiltered))


SeuratFiltered <- FindNeighbors(object = SeuratFiltered, dims = 1:10)
SeuratFiltered <- FindClusters(object = SeuratFiltered,
                          resolution = c(0.4, 0.5, 0.6, 0.8))



### Example plots


## Counts per spot 
ST_plot(Object = SeuratObj, 
        Grob = imgFile, 
        PlotType = "NoGenes", 
        ShowFilter = filter)

ST_plot(Object = sce, 
        Grob = imgFile, 
        PlotType = "NoGenes")

ST_plot(Object = sce, 
        Grob = imgFile, 
        PlotType = "CountsPerSpot")


ST_plot(Object = SeuratObj, Grob = imgFile, 
        PlotType = "CountsPerSpot", 
        ShowFilter = filter)

## Cluster plot
ST_plot(Object = SeuratFiltered, Grob = imgFile, 
        PlotType = "Cluster", 
        ClusterRes = "RNA_snn_res.0.6"
        )

ST_plot(Object = sce, Grob = imgFile, 
        PlotType = "Cluster", 
        ClusterRes = "y"
)

## Gene plot
ST_plot(Object = SeuratFiltered, Grob = imgFile, 
        PlotType = "Gene", 
        Gene = "Pax6")


ST_plot(Object = sce, Grob = imgFile, 
        PlotType = "Gene", 
        Gene = "Pax6")



sce = readSCE(Counts = counts,
              BarcodeFile = barcodesFile, 
              ProjectName = "TestProj", 
              SectionNumber = 1
)





