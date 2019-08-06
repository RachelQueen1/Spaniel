#' @include spaniel_plot_internals.R
#' @include utilities.R

# Tests for functions used by Spatial Plot internal functions
# ------------------------------------------------------------------------------

# Test ungroupVars creates correct variables
# ------------------------------------------------------------------------------
test.vars <- ungroupVars("a", "b", "c")
test_that("ungroupVars extracts the correct variables", {
    expect_is(test.vars, "varGroup")
    expect_equal(length(test.vars), 3)
    expect_equal(test.vars[[1]], "a")
    expect_equal(test.vars[[2]], "b")
    expect_equal(test.vars[[3]], "c")
})


# Test %=% splits the list correctly
# ------------------------------------------------------------------------------
ungroupVars(plotTitle,cl,sz) %=% test.vars
test_that("%=% splits the list correctly", {
    expect_equal(plotTitle, "a")
    expect_equal(cl, "b")
    expect_equal(sz, "c")
})

# Test setVars creates the correct variables
# ------------------------------------------------------------------------------
## Test with Seurat Object
Object <- readRDS(file.path(system.file(package = "Spaniel"), 
                            "extdata/SeuratData.rds"))
pt.size <- 2
Gene <- "Mef2c"
ClusterRes <- "cluster_RNA_snn_res.0.4"

# No Genes
PlotType <- "NoGenes"
ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
    setVars(Object, PlotType, pt.size, Gene, ClusterRes)

test_that("setVars sets the correct variables", {
    expect_equal(plotTitle, "Number of Genes Per Spot")
    expect_equal(cl, "No_Of_Genes")
    expect_equal(sz, "No_Of_Genes")
    expect_equal(shp, "NULL")
    expect_equal(show_size_legend, "TRUE")
    expect_equal(colPlot, "nFeature_RNA")
})

# CountsPerSpot
PlotType <- "CountsPerSpot"
ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
    setVars(Object, PlotType, pt.size, Gene, ClusterRes)

test_that("setVars sets the correct variables: CountsPerSpot", {
    expect_equal(plotTitle, "Total Counts Per Spot")
    expect_equal(cl, "Exprs")
    expect_equal(sz, "Exprs")
    expect_equal(shp, "NULL")
    expect_equal(show_size_legend, "TRUE")
    expect_equal(colPlot, "nCount_RNA")
})

# Cluster
PlotType <- "Cluster"
ClusterRes <- "cluster_RNA_snn_res.0.6"

ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
    setVars(Object, PlotType, pt.size, Gene, ClusterRes)

test_that("setVars sets the correct variables", {
    expect_equal(plotTitle, "Spot Clusters")
    expect_equal(cl, "Cluster")
    expect_equal(sz, "2")
    expect_equal(shp, "Cluster")
    expect_equal(show_size_legend, "FALSE")
    expect_equal(colPlot, ClusterRes)
})

# Gene
PlotType <- "Gene"
ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
    setVars(Object, PlotType, pt.size, Gene, ClusterRes)

test_that("setVars sets the correct variables: Gene", {
    expect_equal(plotTitle, paste("Expression of", Gene))
    expect_equal(cl, "Exp")
    expect_equal(sz, "Exp")
    expect_equal(shp, "NULL")
    expect_equal(show_size_legend, "TRUE")
    expect_equal(colPlot, Gene)
})


## Test setVars with SCE Object
# No Genes
Object <- readRDS(file.path(system.file(package = "Spaniel"), 
                            "extdata/sceData.rds"))

PlotType <- "NoGenes"

ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
    setVars(Object, PlotType, pt.size, Gene)

test_that("setVars sets the correct variables", {
    expect_equal(colPlot, "total_features_by_counts")
})

# No Counts
PlotType <- "CountsPerSpot"

ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
    setVars(Object, PlotType, pt.size, Gene)

test_that("setVars sets the correct variables", {
    expect_equal(colPlot, "total_counts")
})


# Test makeGGDF creates the correct data frame
# ------------------------------------------------------------------------------
## Test with Seurat Object
Object <- readRDS(file.path(system.file(package = "Spaniel"), 
                            "extdata/SeuratData.rds"))


### using column from metadata
PlotType <- "CountsPerSpot"
colPlot <- colnames(Object@meta.data)[2]
cl <- "Exprs"

test.tmp <- makeGGDF(Object, PlotType, colPlot, cl)

test_that("setVars sets the correct variables", {
    expect_is(test.tmp, "data.frame")
    expect_equal(dim(test.tmp), c(248, 4))
    expect_equal(colnames(test.tmp), c(cl, "spot", "x", "y"))
})

# test that works for gene plots correctly
PlotType <- "Gene"
colPlot <- rownames(Object)[2]
cl <- "Exprs"

test.tmp <- makeGGDF(Object, PlotType, colPlot, cl)

test_that("setVars sets the correct variables", {
    expect_is(test.tmp, "data.frame")
    expect_equal(dim(test.tmp), c(248, 4))
    expect_equal(colnames(test.tmp), c(cl, "spot", "x", "y"))
})






### Test plot image function
# ------------------------------------------------------------------------------
test.md <- getMetadata(Object)
test.cood <- test.cood <- getCoordinates(test.md)
test.cood$plotCols <- seq(1, 248)
test.cood$plotSize <- seq(1, 248)

test.plot <- plotImage(Tmp = test.cood,
                       Grob = grid::roundrectGrob(),
                       Colour = "plotCols",
                       Size = "plotSize",
                       PlotTitle = "Title")
test_that("plotImages test plotting function", {
    expect_is(test.plot , c("gg","ggplot"))
    expect_equal(test.plot$data, test.cood)
    expect_equal(test.plot$labels$title, "Title")
    expect_equal(test.plot$labels$size, "plotSize")
    expect_equal(test.plot$labels$colour, "plotCols")
})







