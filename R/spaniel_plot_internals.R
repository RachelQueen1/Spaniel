# General Functions
# ------------------------------------------------------------------------------

# Grouping the left hand side
ungroupVars = function(...) {
    List = as.list(substitute(list(...)))[-1L]
    class(List) = 'varGroup'
    return(List)
}

# Binary Operator
'%=%' = function(l, r) {
    Envir = as.environment(-1)
    for (i in 1:length(l)) {
        do.call('<-', list(l[[i]], r[[i]]), envir=Envir)
    }
}

# Test that either a Seurat or SingleCellExperiment object is provided
testObject <- function(Object){
    try(if(class(Object) != "Seurat" & class(Object) != "SingleCellExperiment") 
        stop("Object must be either a Seurat or SingleCellExperiment object")
    )
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
        plotTitle = "Total Counts Per Spot"
        cl = "Exprs"
        sz = "Exprs"
        shp = "NULL"
        colPlot <- ifelse(class(Object) == "Seurat",
                              "nCount_RNA",
                              "total_counts")
    }
    
    if (PlotType == "Cluster") {
        plotTitle = "Spot Clusters"
        cl = "Cluster"
        sz = pt.size
        shp = "Cluster"
        colPlot = ClusterRes
    }
    
    if (PlotType == "Gene") {
        plotTitle = paste("Expression of", Gene)
        cl = "Exp"
        sz = "Exp"
        shp = "NULL"
        colPlot = Gene
    }
    
    show_size_legend = ifelse(PlotType == "Cluster", F, T)
    
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
make_ggdf <- function(Object, PlotType, colPlot, cl){
    
    ### Get Metadata and Coodinates
    MetaData = getMetadata(Object)
    Coordinates = MetaData[, c("x", "y")]
    Coordinates$spot = rownames(MetaData)
    
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
    return(tmp)
}

