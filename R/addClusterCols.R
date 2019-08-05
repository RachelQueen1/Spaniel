#' @include utilities.R
#'
NULL



#' markClusterCol
#' 
#' A function to mark the columns containing cluster information in the metadata
#' or colData of a Seurat or SCE object. Columns are marked with "cluster_" 
#' prefix.
#' @param Object Either a Seurat or SCE object containing clustering information
#' @param Pattern pattern indicating which columns contain cluster information
#' @examples 
#' SeuratObj <- readRDS(file.path(system.file(package = "Spaniel"), 
#'                      "extdata/SeuratData.rds"))
#' SeuratObj <- markClusterCol(SeuratObj, "res")
#' @export
markClusterCol <- function(Object, Pattern){
    MetaData <- getMetadata(Object)
    whichCols <- grep(Pattern, colnames(MetaData))
    colnames(MetaData)[whichCols] <- paste0("cluster_", 
                                            colnames(MetaData)[whichCols])
    Object<-updateMetadata(MetaData, Object)
    return(Object)
}



