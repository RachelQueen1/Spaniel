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
make_ggdf <- function(Object, MetaData, Coordinates, colPlot){
    tmp = data.frame(No_Of_Genes = MetaData[, colPlot],
                     spot = rownames(MetaData)) %>% 
        as.tbl %>% 
        inner_join(Coordinates)
    return(tmp)
}



#NoGenes
df_NoGenes <- function(Object, MetaData, Coordinates, colPlot){
    tmp = data.frame(No_Of_Genes = MetaData[, colPlot],
                     spot = rownames(MetaData)) %>% 
        as.tbl %>% 
        inner_join(Coordinates)
    return(tmp)
}


#CountsPerSpot
df_CountsPerSpot <- function(Object, MetaData, Coordinates, colPlot){
    tmp = data.frame(Exprs = MetaData[, colPlot],
                     spot = rownames(MetaData)) %>% 
        as.tbl %>% 
        inner_join(Coordinates)
    
    return(tmp)
}

#Cluster
df_Cluster <- function(){
    tmp = data.frame(cluster = MetaData[, colPlot], 
                     spot = rownames(MetaData)) %>%
        as.tbl %>% inner_join(Coordinates)
    return(tmp)
}




#Gene
df_Gene <- function(){
    
    return(tmp)
}



getGeneColPlot <- function(Object){
    if (class(Object) == "Seurat"){
        colPlot = "nFeature_RNA"}
    if (class(Object) == "SingleCellExperiment"){
        colPlot = "total_features_by_counts"}
    return(colPlot)}








genePlotDF <- function(Object, Coordinates, colPlot){
    
    return(tmp)
}


