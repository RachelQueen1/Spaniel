#'@importFrom magrittr %>%
#'@import methods
#'@import dplyr

### get expression data
getExprs <- function(Object){
    if (is(Object, "Seurat")){
        Exprs = Object@assays$RNA@scale.data}
    if (is(Object, "SingleCellExperiment")){
        Exprs = SingleCellExperiment::logcounts(Object)}
    return(Exprs)
    
}

### get ColData/Meta data
getMetadata <- function(Object){
    if (is(Object, "Seurat")){
        MetaData <- Object@meta.data
    }
    if (is(Object, "SingleCellExperiment")){
        MetaData <- Object@colData %>% data.frame()}
    
    return(MetaData)
}

### get coordinates
getCoordinates <- function(Metadata){
    Coordinates <- Metadata[, c("spot", "x", "y")]
    return(Coordinates)
}


### update metadata
updateMetadata <- function(MetaData, Object){
    if (is(Object, "Seurat")){
        Object@meta.data <- MetaData}
    if (is(Object, "SingleCellExperiment")){
        MetaData <- MetaData %>% DataFrame()
        colData(Object) <- MetaData}
    return(Object)
}


### Test that either a Seurat or SingleCellExperiment object is provided
testObject <- function(Object){
    try(if(!is(Object, "Seurat") & !is(Object, "SingleCellExperiment")) 
        stop("Object must be either a Seurat or SingleCellExperiment object")
    )
}

