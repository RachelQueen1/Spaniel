#'@importFrom magrittr %>%



# general Functions
# ------------------------------------------------------------------------------

# Grouping the left hand side
ungroupVars <- function(...) {
    List <- as.list(substitute(list(...)))[-1L]
    class(List) <- 'varGroup'
    return(List)
}

# Binary Operator
'%=%' <- function(l, r) {
    envir <- as.environment(-1)
    for (i in seq_len(length(l))) {
        do.call('<-', list(l[[i]], r[[i]]), envir=envir)
    }
}


# Set the variables for plot
# ------------------------------------------------------------------------------
setVars <- function(object, 
                    plotType, 
                    pointpointSize = NULL, 
                    gene = NULL, 
                    clusterRes = NULL){
    if (plotType == "NoGenes") {
        plotTitle <- "Number of Genes Per Spot"
        cl <- "No_Of_Genes"
        sz <- "No_Of_Genes"
        shp <- "NULL"
        colPlot <- ifelse(is(object, "Seurat"),
                            "nFeature_RNA",
                            "detected")
    }
    
    if (plotType == "CountsPerSpot") {
        plotTitle <- "Total Counts Per Spot"
        cl <- "Exprs"
        sz <- "Exprs"
        shp <- "NULL"
        colPlot <- ifelse(is(object, "Seurat"),
                            "nCount_RNA",
                            "sum")
    }
    
    if (plotType == "Cluster") {
        plotTitle <- "Spot Clusters"
        cl <- "Cluster"
        sz <- pointpointSize
        shp <- "Cluster"
        colPlot <- clusterRes
    }
    
    if (plotType == "Gene") {
        plotTitle <- paste("Expression of", gene)
        cl <- "Exp"
        sz <- "Exp"
        shp <- "NULL"
        colPlot <- gene
    }
    
    showpointSizeLegend <- ifelse(plotType == "Cluster", FALSE, TRUE)
    
    return(c(plotTitle,
                cl,
                sz,
                shp,
                showpointSizeLegend,
                colPlot))
}



# Create ggplot df for each plot type
# ------------------------------------------------------------------------------
### Make a generic function for all 4 plot types
makeGGDF <- function(object, plotType, colPlot, cl){
    
    ### Get Metadata and Coodinates
    metaData <- getMetadata(object)
    coordinates <- metaData[, c("x", "y")]
    coordinates$spot <- rownames(metaData)
    
    if (plotType == "Gene"){
        try(if(!colPlot %in% rownames(object))  
            stop(paste0(colPlot,  "not found in rownames(object"))
        )
        # get expression data
        tmp <- getExprs(object)
        tmp <- data.frame(toPlot = tmp[colPlot, ], 
                            spot = rownames(metaData))
    } else {
        try(if(!colPlot %in% colnames(metaData))  
            stop(paste0(colPlot, "not found in colames(object"))
        )
        tmp <- data.frame(toPlot = metaData[,colPlot],
                            spot = rownames(metaData))
    }
    
    # join by spot
    tmp$spot <- tmp$spot %>% as.character()
    tmp <- tmp %>% 
        dplyr::inner_join(coordinates, by = "spot") 
    # set colname for the to plot column
    colnames(tmp)[1] <- cl
    
    ## reverse the order of the y coordinates
    tmp$y <- 36 - tmp$y
    
    return(tmp)
}


# Plot image
# ------------------------------------------------------------------------------
plotImage <- function(grob, tmp, pointColour, pointSize, plotTitle = NULL, 
                      sizeLegend = TRUE, plotType 
                        ){
    p <- ggplot2::ggplot(tmp ,ggplot2::aes_string("x", "y", 
                                                  color = pointColour, 
                                                  size = pointSize)) 
    
    
    if (plotType == "Original"){
        p <- p + ggplot2::xlim(1, 33) +
            ggplot2::ylim(1, 35) +
            ggplot2::annotation_custom(grob, xmin = 1, xmax = 33, 
                                   ymin = 1, ymax = 35)
    } else if (plotType == "Visium"){
        ### TO ADD
        p    
        
    }
    
    
    } else if (plotType == "Coordinates"){
        ### TO ADD
        p
    
    
    }
    
    ## add theme to plot     
    p<- p + ggplot2::geom_point(alpha = 0.6)  +
                ggplot2::labs(title = plotTitle) +
                ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                        axis.text.x=ggplot2::element_blank(),
                        axis.ticks.x=ggplot2::element_blank(),
                        axis.title.y=ggplot2::element_blank(),
                        axis.text.y=ggplot2::element_blank(),
                        axis.ticks.y=ggplot2::element_blank())
    NULL
    
    ### if show size false
    if (sizeLegend == FALSE){
        p <- p + ggplot2::guides(size=FALSE)}
    ### show plot
    p + ggplot2::guides(color = ggplot2::guide_legend(), 
                        size = ggplot2::guide_legend())
}

