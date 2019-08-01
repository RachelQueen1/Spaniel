# Tests for addClusterCols function
# ------------------------------------------------------------------------------

context("Testing markClusterCol")
# These tests were created to ensure that the markClusterCol functions works
# correctly and marks the correct column with Cluster_ prefix

# Test markClusterCol with Seurat and sce objects 
# ------------------------------------------------------------------------------

# Create test data
set.seed(1234)
counts <- sample(seq(0, 4), 
                 625, 
                 replace = T, 
                 prob = c(0.65, 0.25, 0.04, 0.008, 0.002)) %>% 
    matrix(nrow = 25) %>% 
    data.frame()
colnames(counts) <- paste0("cell_", seq(1, 25))
rownames(counts) <- paste0("gene.", seq(1, 25))
# test metadata with 5 columns the last three columns contain the 
# clustering information
md <- counts %>% 
    t() %>% 
    data.frame() %>%
    select(c(gene.1, gene.2, gene.3, gene.4, gene.5)) %>%
    dplyr::rename(col1 = gene.1,
           col2 = gene.2, 
           res.0.6 = gene.3,
           res.0.8 = gene.4,
           res.1.0 = gene.5)


# Test with Seurat object
# ------------------------------------------------------------------------------
test.seurat <- Seurat::CreateSeuratObject(counts = counts,
                                          meta.data = md)

test.md.before <- getMetadata(test.seurat)

### prefix all columns containing patten with "cluster_"
pattern <- "res"
test.seurat <- markClusterCol(test.seurat, Pattern = pattern)
test.md.after <- getMetadata(test.seurat)

# check that no columns are marked before running markClusterCol
test.marked.before <- colnames(test.md.before) %>% 
    grepl("cluster_", .) %>% 
    sum()

#checked that columns containing cluser columns are marked
test.marked.after <- colnames(test.md.after) %>% 
    grepl("cluster_", .) %>% 
    sum()

#check that columns not containing cluster info remain unmarked
test.marked.after.other <- colnames(test.md.after)[1:5]%>% 
    grepl("cluster_", .) %>% 
    sum()

test_that("markClusterCol check that columns 
          are marked with cluster_ correctly, Seurat", {
    expect_is(test.seurat, "Seurat")
    expect_is(test.md.before, "data.frame")
    expect_is(test.md.after, "data.frame")
    expect_equal(colnames(test.md.before)[6], "res.0.6")
    expect_equal(colnames(test.md.after)[6], "cluster_res.0.6")
    expect_equal(test.marked.before, 0)
    expect_equal(test.marked.after, 3)
    expect_equal(test.marked.after.other, 0)
})


# Test with SingleCellExperiment object
# ------------------------------------------------------------------------------
test.sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), 
                                 colData = md)

test.md.before <- getMetadata(test.sce)

### prefix all columns containing patten with "cluster_"
pattern <- "res"
test.sce <- markClusterCol(test.sce, Pattern = pattern)
test.md.after <- getMetadata(test.sce)

# check that no columns are marked before running markClusterCol
test.marked.before <- colnames(test.md.before) %>% 
    grepl("cluster_", .) %>% 
    sum()

#checked that columns containing cluser columns are marked
test.marked.after <- colnames(test.md.after) %>% 
    grepl("cluster_", .) %>% 
    sum()

#check that columns not containing cluster info remain unmarked
test.marked.after.other <- colnames(test.md.after)[1:2]%>% 
    grepl("cluster_", .) %>% 
    sum()

test_that("markClusterCol check that columns 
          are marked with cluster_ correctly, SCE", {
              expect_is(test.sce, "SingleCellExperiment")
              expect_is(test.md.before, "data.frame")
              expect_is(test.md.after, "data.frame")
              expect_equal(colnames(test.md.before)[3], "res.0.6")
              expect_equal(colnames(test.md.after)[3], "cluster_res.0.6")
              expect_equal(test.marked.before, 0)
              expect_equal(test.marked.after, 3)
              expect_equal(test.marked.after.other, 0)
          })
