#' @include spanielPlotInternals.R
#' @include utilities.R
#' 
NULL

#' Spatial Transcriptomics Plot
#'
#' This function overlays information from a Seurat object or 
#' SingleCellExperiment object containing barcodes  onto a H & E image. 
#' There are 4 plots available showing a) the number of
#' genes detected per spot, b) the number of reads detected per spot,
#' c) clustering results, d) the gene expression of a selected gene.
#'
#' @param object Either a Seurat object (version 3) or a SingleCellExperiment 
#' object containing barcode coordinates in the metadata (Seurat) or
#' colData (SingleCellExperiment). 
#' @param grob an grob to be used as the backgound image see(parseImage)
#' @param plotType There are 5 types of plots avaiable:
#'                       1) NoGenes - This shows the number of genes per spot 
#'                       and uses information from "nFeature_RNA" column of 
#'                       Seurat object or "detected" from a 
#'                       SingleCellExperiment object.
#'                       2) CountsPerSpot - This shows the number of counts per
#'                       spot. It uses information from "nCount_RNA" column of 
#'                       Seurat object or "sum" from a 
#'                       singleCellExperiment object.
#'                       3) Cluster - This plot is designed to show clustering
#'                       results stored in the meta.data or colData of an object
#'                       4) Gene- This plot shows the expression of a single 
#'                       gene. This plot uses scaled/normalised expressin data 
#'                       from the scale.data slot of Seurat object or logcounts 
#'                       of a SingleCellExperiment object.
#'                       5) Other - A generic plot to plot any column from the
#'                       meta.data or colData of an object.
#'                                             
#'                       
#' @param gene Gene to plot
#' @param clusterRes which cluster resolution to plot 
#' @param ptSize Point size used for cluster plot default is 2
#' @param ptSizeMin Minimum point size used for QC and Gene Expression plots 
#'        default is 0
#' @param ptSizeMax Maximum point size used for QC and Gene Expression plots
#'         default is 5
#' @param customTitle Specify plot title (optional)
#' @param scaleData Show scaled data on plot (default is TRUE)
#' @param showFilter Logical filter showing pass/fail for spots
#' @return A ggplot spatial transcriptomics plot
#' @export
#' @examples
#' 
#' 
#' ## Data is taken from DOI: 10.1126/science.aaf2403
#' SeuratObj <- readRDS(file.path(system.file(package = "Spaniel"),
#'                         "extdata/SeuratData.rds"))
#' imgFile <- readRDS(file.path(system.file(package = "Spaniel"),
#'                         "extdata/image.rds"))
#' 
#' ## Counts per spot with a QC filter
#' minGenes <- 2000
#' minUMI <- 300000
#' filter <- SeuratObj$nFeature_RNA > minGenes &
#'             SeuratObj$nCount_RNA > minUMI
#' spanielPlot(object = SeuratObj, grob = imgFile,
#'         plotType = "CountsPerSpot",
#'         showFilter = filter)
#' 
#' ## Cluster plot
#' spanielPlot(object = SeuratObj, grob = imgFile,
#'         plotType = "Cluster",
#'         clusterRes = "cluster_RNA_snn_res.0.6")
#' 
#' ## Gene plot
#' spanielPlot(object = SeuratObj, grob = imgFile,
#'         plotType = "Gene",
#'         gene= "Nrgn")
#' @usage  spanielPlot(object, grob, plotType = c("NoGenes",
#'                                             "CountsPerSpot",
#'                                             "Cluster",
#'                                             "Gene"),
#'                 gene= NULL, clusterRes = NULL, customTitle = NULL,
#'                 scaleData = TRUE, showFilter = NULL, ptSize = 2,
#'                 ptSizeMin = 0, ptSizeMax = 5)




# Main Spaniel Plot Function
# ------------------------------------------------------------------------------
spanielPlot <- function (object,  
                        grob,
                        plotType = c("NoGenes", 
                                    "CountsPerSpot", 
                                    "Cluster", 
                                    "Gene"),
                        gene= NULL, 
                        clusterRes = NULL, 
                        customTitle = NULL,
                        scaleData = TRUE, 
                        showFilter = NULL, 
                        ptSize = 2,
                        ptSizeMin = 0, 
                        ptSizeMax = 5)

{
    ### Validate object is either a Seurat or SCE object
    testObject(object)
    
    ### set title, colour column, and size column according to plotType
    colPlot = NULL
    ungroupVars(plotTitle,cl,sz,shp, showSizeLegend, colPlot) %=% 
        setVars(object, plotType, ptSize, gene, clusterRes)
    
<<<<<<< Updated upstream:R/spatialPlot.R
    # convert shp NULL
    if(shp == "NULL"){shp = NULL}
    # convert size legend to logical
        showSizeLegend = ifelse(showSizeLegend == "TRUE", TRUE, FALSE)
    # convert sz to numeric
    if (!is.na(as.numeric(sz))){
        sz = as.numeric(sz)
    }
=======
    shp <- convertIfNULL(shp)
    showSizeLegend <- convertIfTRUE(showSizeLegend)
    sz <- convertSize(sz)
>>>>>>> Stashed changes:R/spanielPlot.R
    
    ### create data.frame for ggplot
    tmp <- makeGGDF(object, plotType, colPlot, cl)
    
    
    ### Update tmp optional arguments if supplied
    if (!is.null(showFilter)){
        tmp$Filter = showFilter
        cl = "Filter"
    }
    
    if (!is.null(customTitle)) {
        plotTitle = customTitle
    }
    
    ### Create plot
    p <- plotImage(grob, 
                    tmp,
                    cl, 
                    sz, 
                    plotTitle, 
                    showSizeLegend)
    
    
    ### FOR QC PLOTS
    if (!is.null(showFilter)){
        p <- p +
            ggplot2::scale_size_continuous(range=c(ptSizeMin,
                                                    ptSizeMax)) +
            ggplot2::guides(color= ggplot2::guide_legend(), 
                            size = ggplot2::guide_legend())
    }
    
    ### FOR GENE PLOTS AND NUMBER READS/COUNTS
    if (plotType != "Cluster" & is.null(showFilter)){
        p <- p +
            ggplot2::scale_colour_gradient(low="#ff3300", high="#ffff00") +
            ggplot2::scale_size_continuous(range=c(ptSizeMin,
                                                    ptSizeMax)) +
            ggplot2::guides(color= ggplot2::guide_legend(), 
                            size = ggplot2::guide_legend())
        
    }
    ### FOR CLUSTER PLOTS
    if (plotType == "Cluster"){
        p <- p + ggplot2::guides(size=FALSE) +
            ggplot2::scale_size_continuous(range=c(ptSize,
                                                    ptSize))
    }
    
    return(p)
}














