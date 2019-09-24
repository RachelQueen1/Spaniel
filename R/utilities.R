#'@importFrom magrittr %>%
#'@import methods
#'@import dplyr
#'@importFrom utils read.csv

### get expression data
getExprs <- function(object){
    if (is(object, "Seurat")){
        exprs = object@assays$RNA@scale.data}
    if (is(object, "SingleCellExperiment")){
        exprs = SingleCellExperiment::logcounts(object)}
    return(exprs)
    
}

### get ColData/Meta data
getMetadata <- function(object){
    if (is(object, "Seurat")){
        metaData <- object@meta.data
    }
    if (is(object, "SingleCellExperiment")){
        metaData <- object@colData %>% data.frame()}
    
    return(metaData)
}

### get coordinates
getCoordinates <- function(metaData){
    coordinates <- metaData[, c("spot", "x", "y")]
    return(coordinates)
}


### update metadata
updateMetadata <- function(metaData, object){
    if (is(object, "Seurat")){
        object@meta.data <- metaData}
    if (is(object, "SingleCellExperiment")){
        metaData <- metaData %>% S4Vectors::DataFrame()
        colData(object) <- metaData}
    return(object)
}


### Test that either a Seurat or SingleCellExperiment object is provided
testObject <- function(object){
    try(if(!is(object, "Seurat") & !is(object, "SingleCellExperiment")) 
        stop("object must be either a Seurat or SingleCellExperiment object")
    )
}

