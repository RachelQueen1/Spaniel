# General Functions
# ------------------------------------------------------------------------------
# Binary Operator
'%=%.varGroup' = function(l, r, ...) {
    Envir = as.environment(-1)
    for (i in 1:length(l)) {
        do.call('<-', list(l[[i]], r[[i]]), envir=Envir)
    }
}

# Grouping the left hand side
ungroupVars = function(...) {
    List = as.list(substitute(list(...)))[-1L]
    class(List) = 'varGroup'
    return(List)
}

# Gene plot Functions
# ------------------------------------------------------------------------------

setGeneVars <- function(Object){
    plotTitle = "Number of Genes Per Spot"
    cl = "No_Of_Genes"
    sz = "No_Of_Genes"
    return(c(plotTitle,cl,sz))
}

getGeneColPlot <- function(Object){
    if (class(Object) == "Seurat"){
        colPlot = "nFeature_RNA"}
    if (class(Object) == "SingleCellExperiment"){
        colPlot = "total_features_by_counts"}
    return(colPlot)}

genePlotDF <- function(Object, Coordinates, colPlot){
    tmp = data.frame(No_Of_Genes = MetaData[, colPlot],
                     spot = rownames(MetaData)) %>% 
        as.tbl %>% 
        inner_join(Coordinates)
    return(tmp)
}


