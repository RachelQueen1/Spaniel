#' @include spanielPlotInternals.R
#' @include utilities.R

# Tests for that correct plots are produced
# ------------------------------------------------------------------------------


pathToTenXOuts <- file.path(system.file(package = "Spaniel"), "extdata/outs")
object <- createVisiumSCE(tenXDir=pathToTenXOuts, 
                       resolution="Low")

test_that("createVisiumSCE produces and sce", {
  expect_is(object, "SingleCellExperiment")
})



# Tests for setVars and ungroup vars
# ------------------------------------------------------------------------------

plotType = "NoGenes"
showFilter = filter
techType = "Visium"
byCoord = FALSE
ptSizeMax = 3
gene = NULL
clusterRes = NULL

plotVars <- setVars(object, 
        plotType, 
        ptSizeMax, 
        gene, 
        clusterRes)

test_that("PlotVars is the correct length", {
  expect_equal(length(plotVars), 6)
})


ungroupVars(plotTitle,cl,sz,shp, showSizeLegend, colPlot) %=% plotVars

test_that("ungrouping variables work", {
  expect_equal(plotTitle, "Number of Genes Per Spot")
  expect_equal(cl, "No_Of_Genes")
  expect_equal(sz, "No_Of_Genes")
  expect_equal(shp, "NULL")
  expect_equal(showSizeLegend, "TRUE")
  expect_equal(colPlot, "detected")
})


