# Test Spatial plot
# ------------------------------------------------------------------------------

Object = readRDS(file.path(system.file(package = "Spaniel"),
                           "extdata/SeuratData.rds"))
Grob = readRDS(file.path(system.file(package = "Spaniel"),
                            "extdata/image.rds"))
PlotType = "Cluster"
ClusterRes = "cluster_RNA_snn_res.0.6"

test.p <- ST_plot(Object = Object, 
                  Grob = Grob, 
                  PlotType = PlotType, 
                  ClusterRes = ClusterRes)

test_that("ST_plot creates a plot", {
    expect_is(test.p, "ggplot")
})
