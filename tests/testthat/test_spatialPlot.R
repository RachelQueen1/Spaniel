# Tests for functions used by Spatial Plot functions
# ------------------------------------------------------------------------------



# Test '%=%.varGroup' Grouping variables
# ------------------------------------------------------------------------------
x = seq(1, 3)
y <- ungroupVars(x)

y




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

## Test with Seurat Object

# Test setVars creates the correct variables
# ------------------------------------------------------------------------------
Object <- readRDS("../../data/SeuratData.rds")
pt.size <- 2
Gene <- NULL

PlotType <- "NoGenes"
ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colData) %=% 
    setVars(Object, PlotType, pt.size, Gene)

test_that("setVars sets the correct variables", {
    expect_equal(plotTitle, "Number of Genes Per Spot")
    expect_equal(cl, "No_Of_Genes")
    expect_equal(sz, "No_Of_Genes")
    expect_equal(show_size_legend, "TRUE")
})

PlotType <- "CountsPerSpot"
ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colData) %=% 
    setVars(Object, PlotType, pt.size, Gene)

test_that("setVars sets the correct variables", {
    expect_equal(plotTitle, "Total Counts Per Spot")
    expect_equal(cl, "Exprs")
    expect_equal(sz, "Exprs")
    expect_equal(shp, "NULL")
})

PlotType <- "Cluster"
ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colData) %=% 
    setVars(Object, PlotType, pt.size, Gene)

test_that("setVars sets the correct variables", {
    expect_equal(plotTitle, "Spot Clusters")
    expect_equal(cl, "Cluster")
    expect_equal(sz, "2")
    expect_equal(shp, "NULL")
})






