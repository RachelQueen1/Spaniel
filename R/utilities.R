### get expression data
getExprs <- function(Object){
    if (class(Object) == "Seurat"){
        Exprs = Object@assays$RNA@scale.data}
    if (class(Object) == "SingleCellExperiment"){
        Exprs = SingleCellExperiment::logcounts(Object)}
    return(Exprs)
    
}

### get ColData/Meta data
getMetadata <- function(Object){
    if (class(Object) == "Seurat"){
        MetaData <- Object@meta.data
        print("got meta data")}
    if (class(Object) == "SingleCellExperiment"){
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
    if (class(Object) == "Seurat"){
        Object@meta.data <- MetaData}
    if (class(Object) == "SingleCellExperiment"){
        MetaData <- MetaData %>% DataFrame()
        colData(Object) <- MetaData}
    return(Object)
}



#### PLOT IMAGE ####
plotImage = function(Grob, Tmp, Colour, Size, ShowSizeLegend = TRUE, 
                     PlotTitle = NULL){
    p = ggplot2::ggplot(Tmp ,ggplot2::aes_string("x", "y", color = 
                                                     Colour, size = Size)) +
        ggplot2::xlim(1, 33) +
        ggplot2::ylim(1, 35) +
        ggplot2::annotation_custom(Grob, xmin = 1, xmax = 33, 
                                   ymin = 1, ymax = 35) +
        ggplot2::geom_point()  +
        ggplot2::labs(title = PlotTitle) +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank())
    NULL
    
    ### if show size false
    if (ShowSizeLegend == FALSE){
        p = p + ggplot2::guides(size=FALSE)}
    ### show plot
    p
}