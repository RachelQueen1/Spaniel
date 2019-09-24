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
## Test with Seurat object
object <- readRDS(file.path(system.file(package = "Spaniel"), 
                            "extdata/SeuratData.rds"))
pointSize <- 2
gene <- "Nrgn"
ClusterRes <- "cluster_RNA_snn_res.0.4"

# No genes
plotType <- "NoGenes"
ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
    setVars(object, plotType, pointSize, gene, ClusterRes)

test_that("setVars sets the correct variables", {
    expect_equal(plotTitle, "Number of Genes Per Spot")
    expect_equal(cl, "No_Of_Genes")
    expect_equal(sz, "No_Of_Genes")
    expect_equal(shp, "NULL")
    expect_equal(show_size_legend, "TRUE")
    expect_equal(colPlot, "nFeature_RNA")
})

# CountsPerSpot
plotType <- "CountsPerSpot"
ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
    setVars(object, plotType, pointSize, gene, ClusterRes)

test_that("setVars sets the correct variables: CountsPerSpot", {
    expect_equal(plotTitle, "Total Counts Per Spot")
    expect_equal(cl, "Exprs")
    expect_equal(sz, "Exprs")
    expect_equal(shp, "NULL")
    expect_equal(show_size_legend, "TRUE")
    expect_equal(colPlot, "nCount_RNA")
})

# Cluster
plotType <- "Cluster"
ClusterRes <- "cluster_RNA_snn_res.0.6"

ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
    setVars(object, plotType, pointSize, gene, ClusterRes)

test_that("setVars sets the correct variables", {
    expect_equal(plotTitle, "Spot Clusters")
    expect_equal(cl, "Cluster")
    expect_equal(sz, "2")
    expect_equal(shp, "Cluster")
    expect_equal(show_size_legend, "FALSE")
    expect_equal(colPlot, ClusterRes)
})

# gene
plotType <- "Gene"
ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
    setVars(object, plotType, pointSize, gene, ClusterRes)

test_that("setVars sets the correct variables: gene", {
    expect_equal(plotTitle, paste("Expression of", gene))
    expect_equal(cl, "Exp")
    expect_equal(sz, "Exp")
    expect_equal(shp, "NULL")
    expect_equal(show_size_legend, "TRUE")
    expect_equal(colPlot, gene)
})


## Test setVars with SCE object
# No genes
object <- readRDS(file.path(system.file(package = "Spaniel"), 
                            "extdata/sceData.rds"))

plotType <- "NoGenes"

ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
    setVars(object, plotType, pointSize, gene)

test_that("setVars sets the correct variables", {
    expect_equal(colPlot, "detected")
})

# No Counts
plotType <- "CountsPerSpot"

ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
    setVars(object, plotType, pointSize, gene)

test_that("setVars sets the correct variables", {
    expect_equal(colPlot, "sum")
})


# Test makeGGDF creates the correct data frame
# ------------------------------------------------------------------------------
## Test with Seurat object
object <- readRDS(file.path(system.file(package = "Spaniel"), 
                            "extdata/SeuratData.rds"))


### using column from metadata
plotType <- "CountsPerSpot"
colPlot <- colnames(object@meta.data)[2]
cl <- "Exprs"

testTmp <- makeGGDF(object, plotType, colPlot, cl)

test_that("setVars sets the correct variables", {
    expect_is(testTmp, "data.frame")
    expect_equal(dim(testTmp), c(248, 4))
    expect_equal(colnames(testTmp), c(cl, "spot", "x", "y"))
})

# test that works for gene plots correctly
plotType <- "Gene"
colPlot <- rownames(object)[2]
cl <- "Exprs"

testTmp <- makeGGDF(object, plotType, colPlot, cl)

test_that("setVars sets the correct variables", {
    expect_is(testTmp, "data.frame")
    expect_equal(dim(testTmp), c(248, 4))
    expect_equal(colnames(testTmp), c(cl, "spot", "x", "y"))
})






### Test plot image function
# ------------------------------------------------------------------------------
testMd <- getMetadata(object)
testCood <- getCoordinates(testMd)
testCood$plotCols <- seq(1, 248)
testCood$plotSize <- seq(1, 248)
testPlot <- plotImage(grob = grid::roundrectGrob(),
                       tmp = testCood, 
                       pointColour = "plotCols", 
                       pointSize = "plotSize",
                       plotTitle = "Title")

test_that("plotImages test plotting function", {
    expect_is(testPlot , c("gg","ggplot"))
    expect_equal(testPlot$data, testCood)
    expect_equal(testPlot$labels$title, "Title")
    expect_equal(testPlot$labels$size, "plotSize")
    expect_equal(testPlot$labels$colour, "plotCols")
})







