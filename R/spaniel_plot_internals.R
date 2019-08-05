#'@importFrom magrittr %>%



# General Functions
# ------------------------------------------------------------------------------

# Grouping the left hand side
ungroupVars <- function(...) {
    List <- as.list(substitute(list(...)))[-1L]
    class(List) <- 'varGroup'
    return(List)
}

# Binary Operator
'%=%' <- function(l, r) {
    Envir <- as.environment(-1)
    for (i in seq_len(length(l))) {
        do.call('<-', list(l[[i]], r[[i]]), envir=Envir)
    }
}


# Set the variables for plot
# ------------------------------------------------------------------------------
setVars <- function(Object, 
                    PlotType, 
                    pt.size = NULL, 
                    Gene = NULL, 
                    ClusterRes = NULL){
    if (PlotType == "NoGenes") {
        plotTitle <- "Number of Genes Per Spot"
        cl <- "No_Of_Genes"
        sz <- "No_Of_Genes"
        shp <- "NULL"
        colPlot <- ifelse(class(Object) == "Seurat",
                              "nFeature_RNA",
                              "total_features_by_counts")
    }
    
    if (PlotType == "CountsPerSpot") {
        plotTitle <- "Total Counts Per Spot"
        cl <- "Exprs"
        sz <- "Exprs"
        shp <- "NULL"
        colPlot <- ifelse(class(Object) == "Seurat",
                              "nCount_RNA",
                              "total_counts")
    }
    
    if (PlotType == "Cluster") {
        plotTitle <- "Spot Clusters"
        cl <- "Cluster"
        sz <- pt.size
        shp <- "Cluster"
        colPlot <- ClusterRes
    }
    
    if (PlotType == "Gene") {
        plotTitle <- paste("Expression of", Gene)
        cl <- "Exp"
        sz <- "Exp"
        shp <- "NULL"
        colPlot <- Gene
    }
    
    show_size_legend <- ifelse(PlotType == "Cluster", F, T)
    
    return(c(plotTitle,
             cl,
             sz,
             shp,
             show_size_legend,
             colPlot))
}



# Create ggplot df for each plot type
# ------------------------------------------------------------------------------
### Make a generic function for all 4 plot types
makeGGDF <- function(Object, PlotType, colPlot, cl){
    
    ### Get Metadata and Coodinates
    MetaData <- getMetadata(Object)
    Coordinates <- MetaData[, c("x", "y")]
    Coordinates$spot <- rownames(MetaData)
    
    if (PlotType == "Gene"){
        try(if(!colPlot %in% rownames(Object))  
            stop(paste0(colPlot,  "not found in rownames(Object"))
        )
        # get expression data
        tmp <- getExprs(Object)
        tmp <- data.frame(toPlot = tmp[colPlot, ], 
                          spot = rownames(MetaData))
    } else {
        try(if(!colPlot %in% colnames(MetaData))  
            stop(paste0(colPlot, "not found in colames(Object"))
        )
        tmp <- data.frame(toPlot = MetaData[,colPlot],
                          spot = rownames(MetaData))
    }
    
    # join by spot
    tmp$spot <- tmp$spot %>% as.character()
    tmp <- tmp %>% 
        dplyr::inner_join(Coordinates, by = "spot") 
    # set colname for the to plot column
    colnames(tmp)[1] <- cl
    
    ## reverse the order of the y coordinates
    tmp$y <- 36 - tmp$y
    
    return(tmp)
}


# Plot image
# ------------------------------------------------------------------------------
plotImage <- function(Grob, Tmp, Colour, Size, ShowSizeLegend = TRUE, 
                     PlotTitle = NULL){
    p <- ggplot2::ggplot(Tmp ,ggplot2::aes_string("x", "y", color = 
                                                     Colour, size = Size)) +
        ggplot2::xlim(1, 33) +
        ggplot2::ylim(1, 35) +
        ggplot2::annotation_custom(Grob, xmin = 1, xmax = 33, 
                                   ymin = 1, ymax = 35) +
        ggplot2::geom_point(alpha = 0.6)  +
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
        p <- p + ggplot2::guides(size=FALSE)}
    ### show plot
    p + ggplot2::guides(color = ggplot2::guide_legend(), 
                        size = ggplot2::guide_legend())
}

