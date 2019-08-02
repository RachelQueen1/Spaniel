# Test Spatial plot
# ------------------------------------------------------------------------------

Object <- readRDS("../../data/SeuratData.rds")  
Grob <- readRDS("../../data/image.rds")
PlotType <- "CountsPerSpot"
Gene <- NULL 
ClusterRes <- NULL 
# CustomTitle <- NULL
# ScaleData <- T 
# ShowFilter <- NULL 
# pt.size <- 2
# pt.size.min <- 0 
# pt.size.max <- 5

test.p <- ST_plot(Object = Object, Grob, PlotType)


Object = readRDS("../../data/SeuratData.rds")  
Grob = image 
PlotType = "Cluster"
ClusterRes = "cluster_RNA_snn_res.0.6"

test.p <- ST_plot(Object = Object, Grob, PlotType, ClusterRes)
