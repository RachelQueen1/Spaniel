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
#' @export
#' @examples
#' ## Data is taken from DOI: 10.1126/science.aaf2403
#' exampleCounts <- readRDS(exampleCounts.rds)
#' exampleBarcodes <- "path to barcode"
#' seuratOb <- readSeurat(exampleCounts,
#'                        exampleBarcodes,
#'                        )


readSeurat <- function(Counts,
                       BarcodeFile,
                       ProjectName=projectName,
                       SectionNumber=sectionNo){


    barcodes <- read.csv(BarcodeFile, sep="\t", header=F)
    rownames(barcodes) <- barcodes$V1
    barcodes = barcodes[colnames(counts), ]
    colnames(barcodes) <- c("Barcode", "x", "y")
    seuratObj <- Seurat::CreateSeuratObject(raw.data=Counts,
                                           meta.data=barcodes,
                                           project=paste0(ProjectName,
                                                          SectionNumber,
                                                          sep = "_"))

    return(seuratObj)

}


### change
readSeuratV3 <- function(Counts,
                       BarcodeFile,
                       ProjectName=projectName,
                       SectionNumber=sectionNo){
    
    
    barcodes <- read.csv(BarcodeFile, sep="\t", header=F)
    rownames(barcodes) <- barcodes$V1
    barcodes = barcodes[colnames(counts), ]
    colnames(barcodes) <- c("Barcode", "x", "y")
    seuratObj <- Seurat::CreateSeuratObject(raw.data=Counts,
                                            meta.data=barcodes,
                                            project=paste0(ProjectName,
                                                           SectionNumber,
                                                           sep = "_"))
    
    return(seuratObj)
    
}


### change
readSingleCellExperiment <- function(Counts,
                         BarcodeFile,
                         ProjectName=projectName,
                         SectionNumber=sectionNo){
    
    
    barcodes <- read.csv(BarcodeFile, sep="\t", header=F)
    rownames(barcodes) <- barcodes$V1
    barcodes = barcodes[colnames(counts), ]
    colnames(barcodes) <- c("Barcode", "x", "y")
    seuratObj <- Seurat::CreateSeuratObject(raw.data=Counts,
                                            meta.data=barcodes,
                                            project=paste0(ProjectName,
                                                           SectionNumber,
                                                           sep = "_"))
    
    return(seuratObj)
    
}
