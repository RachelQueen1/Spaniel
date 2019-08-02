# Tests for utilities functions
# ------------------------------------------------------------------------------

context("Testing utilities")
# These tests were created to ensure that the utilities functions work correctly


seurat_obj <- readRDS("../../data/SeuratData.rds")
sce_obj <- readRDS("../../data/sceData.rds")


# Ensure that package works with current version of Seurat and Single Experiment

# Test getExprs with Seurat and sce objects
# ------------------------------------------------------------------------------
test.seurat <- getExprs(seurat_obj)
test_that("getExprs extracts a valid gene * spot matrix from a Seurat object", {
    expect_is(test.seurat, "matrix")
    expect_equal(colnames(test.seurat)[1], "ACAACTATGGGTTGGCGG")
    expect_equal(rownames(test.seurat)[1], "Mef2c")
})


test.sce <- getExprs(sce_obj)
test_that("getExprs extracts a valid gene * spot matrix from a sce object", {
    expect_is(test.sce, "matrix")
    expect_equal(colnames(test.sce)[1], "ACAACTATGGGTTGGCGG")
    expect_equal(rownames(test.sce)[1], "Mef2c")
})


# Test getMetadata with Seurat and sce objects 
# ------------------------------------------------------------------------------
test.md <- getMetadata(seurat_obj)
test_that("getMetadata extracts a valid metadata data.frame
          from a Seurat object", {
    expect_is(test.md, "data.frame")
    expect_equal(colnames(test.md)[1], "orig.ident")
    expect_equal(rownames(test.md)[1], "ACAACTATGGGTTGGCGG")
    expect_equal(dim(test.md), c(248, 10))
})


test.md <- getMetadata(sce_obj)
test_that("getExprs extracts a valid metadata data.frame", {
    expect_is(test.md, "data.frame")
    expect_equal(colnames(test.md)[1], "spot")
    expect_equal(rownames(test.md)[1], "ACAACTATGGGTTGGCGG")
    expect_equal(dim(test.md), c(248, 19))
})


# Test getCoordinates
# ------------------------------------------------------------------------------
test.cood <- getCoordinates(test.md)
test_that("getCoordinates extracts spot coordinates from metadata", {
    expect_is(test.cood, "data.frame")
    expect_equal(colnames(test.cood)[1], "spot")
    expect_equal(rownames(test.cood)[1], "ACAACTATGGGTTGGCGG")
    expect_equal(test.cood[1,2], 16)
    expect_equal(dim(test.cood), c(248, 3))
})


# Test updateMetadata
# ------------------------------------------------------------------------------
md = data.frame(test1 = seq(1,248),
                test2 = seq(248,1))
rownames(md) = colnames(seurat_obj)

test.seu <- updateMetadata(md, seurat_obj)
test_that("getCoordinates extracts spot coordinates from metadata", {
    expect_is(test.seu@meta.data, "data.frame")
    expect_equal(colnames(test.seu@meta.data)[1], "test1")
    expect_equal(rownames(test.seu@meta.data)[1], "ACAACTATGGGTTGGCGG")
    expect_equal(test.seu@meta.data[1,2], 248)
    expect_equal(dim(test.seu@meta.data), c(248, 2))
})

test.sce <- updateMetadata(md, sce_obj)
test_that("updateMetadata ", {
    expect_is(colData(test.sce) , "DataFrame")
    expect_equal(colnames(colData(test.sce))[1], "test1")
    expect_equal(rownames(colData(test.sce))[1], "ACAACTATGGGTTGGCGG")
    expect_equal(colData(test.sce)[1,2], 248)
    expect_equal(dim(colData(test.sce)), c(248, 2))
})





