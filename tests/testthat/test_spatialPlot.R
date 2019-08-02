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


"CountsPerSpot", 
"Cluster", 
"Gene"