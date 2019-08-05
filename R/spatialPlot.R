#' @include spaniel_plot_internals.R
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
#' @param Object Either a Seurat object (version 3) or a SingleCellExperiment 
#' object containing barcode coordinates in the metadata (Seurat) or
#' colData (SingleCellExperiment). 
#' @param Grob an grob to be used as the backgound image see(parseImage)
#' @param PlotType There are 5 types of plots avaiable:
#'                       1) NoGenes - This shows the number of genes per spot 
#'                       and uses information from "nFeature_RNA" column of 
#'                       Seurat object or "total_features_by_counts" from a 
#'                       SingleCellExperiment object.
#'                       2) CountsPerSpot - This shows the number of counts per
#'                       spot. It uses information from "nCount_RNA" column of 
#'                       Seurat object or "total_counts" from a 
#'                       singleCellExperiment object.
#'                       3) Cluster - This plot is designed to show clustering
#'                       results stored in the meta.data or colData of an object
#'                       4) Gene - This plot shows the expression of a single 
#'                       gene. This plot uses scaled/normalised expressin data 
#'                       from the scale.data slot of Seurat object or logcounts 
#'                       of a SingleCellExperiment object.
#'                       5) Other - A generic plot to plot any column from the
#'                       meta.data or colData of an object.
#'                                             
#'                       
#' @param Gene Gene to plot
#' @param ClusterRes which cluster resolution to plot 
#' @param pt.size Point size used for cluster plot default is 2
#' @param pt.size.min Minimum point size used for QC and Gene Expression plots default is 0
#' @param pt.size.max Maximum point size used for QC and Gene Expression plots default is 5
#' @param CustomTitle Specify plot title (optional)
#' @param ScaleData Show scaled data on plot (default is TRUE)
#' @param ShowFilter Logical filter showing pass/fail for spots
#' @export
#' @examples
#' 
#' 
#' ## Data is taken from DOI: 10.1126/science.aaf2403
#' SeuratObj <- readRDS(file.path(system.file(package = "Spaniel"), 
#'                      "extdata/SeuratData.rds"))
#' imgFile <- readRDS(file.path(system.file(package = "Spaniel"), 
#'                      "extdata/image.rds"))
#' 
#' 
#' ## Counts per spot with a QC filter 
#' minGenes <- 2000
#' minUMI <- 300000
#' filter <- SeuratObj@meta.data$nFeature_RNA > minGenes & 
#'           SeuratObj@meta.data$nCount_RNA > minUMI
#' ST_plot(Object = SeuratObj, Grob = imgFile, 
#'        PlotType = "CountsPerSpot", 
#'        ShowFilter = filter)
#'
#' ## Cluster plot
#' ST_plot(Object = SeuratObj, Grob = imgFile, 
#'        PlotType = "Cluster", 
#'        ClusterRes = "cluster_RNA_snn_res.0.6")
#'
## Gene plot
#' ST_plot(Object = SeuratObj, Grob = imgFile, 
#'        PlotType = "Gene", 
#'        Gene = "Nrgn")
# Main Spaniel Plot Function
# ------------------------------------------------------------------------------
ST_plot <- function (Object,  
                     Grob,
                     PlotType = c("NoGenes", 
                                  "CountsPerSpot", 
                                  "Cluster", 
                                  "Gene"),
                     Gene = NULL, 
                     ClusterRes = NULL, 
                     CustomTitle = NULL,
                     ScaleData = TRUE, 
                     ShowFilter = NULL, 
                     pt.size = 2,
                     pt.size.min = 0, 
                     pt.size.max = 5)

{
    ### Validate Object is either a Seurat or SCE object
    testObject(Object)
    
    ### set title, colour column, and size column according to plotType
    ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
        setVars(Object, PlotType, pt.size, Gene, ClusterRes)
    
            # convert shp NULL
            if(shp == "NULL"){shp = NULL}
            # convert size legend to logical
            show_size_legend = ifelse(show_size_legend == "TRUE", T, F)
            # convert sz to numeric
            if (suppressWarnings(!is.na(as.numeric(sz)))){
                sz = as.numeric(sz)
            }
    
    ### create data.frame for ggplot
    tmp <- makeGGDF(Object, PlotType, colPlot, cl)
    
    
    ### Update tmp optional arguments if supplied
    if (!is.null(ShowFilter)){
        tmp$Filter = ShowFilter
        cl = "Filter"
    }
    
    if (!is.null(CustomTitle)) {
        plotTitle = CustomTitle
    }
    
    ### Create plot
    p <- plotImage(Grob, tmp, 
                  Colour = cl, 
                  Size = sz, 
                  PlotTitle = plotTitle, 
                  ShowSizeLegend = show_size_legend)
    
    
    ### FOR QC PLOTS
    if (!is.null(ShowFilter)){
        p <- p +
            ggplot2::scale_size_continuous(range=c(pt.size.min,
                                          pt.size.max)) +
            ggplot2::guides(color= ggplot2::guide_legend(), 
                            size = ggplot2::guide_legend())
    }
    
    ### FOR GENE PLOTS AND NUMBER READS/COUNTS
    if (PlotType != "Cluster" & is.null(ShowFilter)){
        p <- p +
            ggplot2::scale_colour_gradient(low="#ff3300", high="#ffff00") +
            ggplot2::scale_size_continuous(range=c(pt.size.min,
                                          pt.size.max)) +
            ggplot2::guides(color= ggplot2::guide_legend(), 
                            size = ggplot2::guide_legend())
        
    }
    ### FOR CLUSTER PLOTS
    if (PlotType == "Cluster"){
        p <- p + ggplot2::guides(size=FALSE) +
            ggplot2::scale_size_continuous(range=c(pt.size,
                                          pt.size))
    }
    
    return(p)
}














