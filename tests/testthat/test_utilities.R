# Tests for utilities functions
# --------------------------------------------------------------------------------

context("Testing utilities")
# These tests were created to ensure that the utilities functions

seurat_obj <- readRDS("../data/SeuratData.rds")
sce_obj <- readRDS("../data/sceData.rds")

test.data1 <- getExprs(seurat_obj)
test_that("getExprs extracts a valid gene * spot matrix from a Seurat object", {
    expect_is(test.data1, "matrix")
    expect_equal(colnames(test.data1)[1], "ACAACTATGGGTTGGCGG")
    expect_equal(rownames(test.data1)[1], "Mef2c")
})


test.data1 <- getExprs(sce_obj)
test_that("getExprs extracts a valid gene * spot matrix from a sce object", {
    expect_is(test.data1, "matrix")
    expect_equal(colnames(test.data1)[1], "ACAACTATGGGTTGGCGG")
    expect_equal(rownames(test.data1)[1], "Mef2c")
})
