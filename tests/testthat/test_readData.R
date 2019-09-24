# Tests for readData function
# ------------------------------------------------------------------------------

context("Testing readData")
# These tests were created to ensure that the readData functions works
# correctly and creates either a sce object or a Seurat object

# Test parseImage 
# ------------------------------------------------------------------------------
# Load test data and set variables for project
set.seed(1234)

counts <- readRDS(file.path(system.file(package = "Spaniel"), 
                            "extdata/counts.rds"))
barcodefile <- file.path(system.file(package = "Spaniel"), 
                         "1000L2_barcodes.txt")
projectname = "test_project"
sectionNo = "1A"

# create a test dataset
test.seurat <- createSeurat(counts,
                          barcodefile,
                          projectname,
                          sectionNo)

test_that("createSeurat creates a Seurat object containing barcodes correctly", {
    expect_is(test.seurat, "Seurat")
    expect_equal(test.seurat@meta.data[12, "x"], 7)
    expect_equal(test.seurat@meta.data[8, "y"], 15)
    expect_equal(colnames(test.seurat@meta.data)[4:6], 
                 c("spot", 
                   "x",
                   "y"))
})


