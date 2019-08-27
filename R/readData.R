#'@import Seurat
#'@import SingleCellExperiment
#'@import SummarizedExperiment

NULL

#' Create a Seurat Object From Spatial Transcriptomics Data
#'
#' This function converts a count matrix into a Seurat object. The barcodes for
#' each spot are added to the metadata of the Seurat object and are used in
#' plotting the data.
#'
#' @param Counts Raw count matrix or data frame where each row represents a
#'               gene and each column represents barcoded location on a spatial
#'               transcriptomics slide. The columns should be named using the
#'               spot barcode (eg "GTCCGATATGATTGCCGC")
#' @param BarcodeFile a tab seperated barcode file supplied by Spatial
#'               Trancscriptomics. The file should contains three column:
#'               The first column contains the Spatial Transcriptomics barcode,
#'               the second and third column equate to the x and y location
#' @param ProjectName The name of the project which is stored in the Seurat
#'                    Object.
#' @param SectionNumber The location of the sample on the slide
#' @return A Seurat Object
#' 
#' @usage  readSeurat(Counts, BarcodeFile, ProjectName = projectName,
#'                         SectionNumber = sectionNo)
#'  
#' 
#' @export
#' @examples
#' ## Data is taken from DOI: 10.1126/science.aaf2403
#' exampleCounts <- readRDS(file.path(system.file(package = "Spaniel"),
#'                             "extdata/counts.rds"))
#' exampleBarcodes <- file.path(system.file(package = "Spaniel"),
#'                             "1000L2_barcodes.txt")
#' SeuratObj <- readSeurat(exampleCounts,
#'                         exampleBarcodes,
#'                         ProjectName = "TestProj",
#'                         SectionNumber = 1
#'                         )


readSeurat <- function(Counts,
                        BarcodeFile,
                        ProjectName = "projectName",
                        SectionNumber = "sectionNo") {
    barcodes <- utils::read.csv(BarcodeFile, sep = "\t", header = FALSE)
    rownames(barcodes) <- barcodes$V1
    barcodes <- barcodes[colnames(Counts),]
    colnames(barcodes) <- c("spot", "x", "y")
    seuratObj <- Seurat::CreateSeuratObject(counts = Counts,
                                            project = paste0(ProjectName,
                                                                SectionNumber,
                                                                sep = "_"))
    seuratObj@meta.data[, c("spot",
                            "x", "y")] <- barcodes[, c("spot",
                                                        "x", "y")]
    return(seuratObj)
    
}


#' Create a SingleCellExperiment Object From Spatial Transcriptomics Data
#'
#' This function converts a count matrix into a SingleCellExperiment object. 
#' The barcodes for each spot are added to the coldata of the 
#' SingleCellExperiment object and are used in  plotting the data.
#'
#' @param Counts Raw count matrix or data frame where each row represents a
#'               gene and each column represents barcoded location on a spatial
#'               transcriptomics slide. The columns should be named using the
#'               spot barcode (eg "GTCCGATATGATTGCCGC")
#' @param BarcodeFile a tab seperated barcode file supplied by Spatial
#'               Trancscriptomics. The file should contains three column:
#'               The first column contains the Spatial Transcriptomics barcode,
#'               the second and third column equate to the x and y location
#' @param ProjectName The name of the project which is stored in the Seurat
#'                    Object.
#' @param SectionNumber The location of the sample on the slide
#' @usage  readSCE(Counts, BarcodeFile, ProjectName=projectName, 
#'                 SectionNumber=sectionNo)


#' @return A SingleCellExeriment Object
#' @export
#' @examples
#' ## Data is taken from DOI: 10.1126/science.aaf2403
#' exampleCounts <- readRDS(file.path(system.file(package = "Spaniel"),
#'                             "extdata/counts.rds"))
#' exampleBarcodes <- file.path(system.file(package = "Spaniel"),
#'                             "1000L2_barcodes.txt")
#' seuratOb <- readSCE(exampleCounts,
#'                         exampleBarcodes,
#'                         ProjectName = "TestProj",
#'                         SectionNumber = 1)

readSCE <- function(Counts,
                    BarcodeFile,
                    ProjectName="projectName",
                    SectionNumber="sectionNo"){
    
    barcodes <- utils::read.csv(BarcodeFile, sep="\t", header=FALSE)
    rownames(barcodes) <- barcodes$V1
    colnames(barcodes) <- c("spot", "x", "y")
    
    Counts <- Counts[, intersect(colnames(Counts), barcodes$spot)]
    
    barcodes <- barcodes[colnames(Counts), ]
    
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = as.matrix(Counts)), 
        colData = barcodes)
    SummarizedExperiment::colData(sce)$project <- paste0(ProjectName,
                                                        SectionNumber,
                                                        sep = "_")
    
    ### calculate QC metrics
    # sce <- scater::calculateQCMetrics(
    #     sce
    # )
    
    sce <- scater::addQCPerCell(
        sce
    )
    
    return(sce)
    
}
