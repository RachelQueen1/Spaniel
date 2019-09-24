# Test Spatial plot
# ------------------------------------------------------------------------------

object = readRDS(file.path(system.file(package = "Spaniel"),
                           "extdata/SeuratData.rds"))
grob = readRDS(file.path(system.file(package = "Spaniel"),
                         "extdata/image.rds"))
plotType = "Cluster"
clusterRes = "cluster_RNA_snn_res.0.6"

test.p <- spanielPlot(object, 
                    grob, 
                    plotType, 
                    clusterRes)

test_that("spanielPlot creates a plot", {
    expect_is(test.p, "ggplot")
})
