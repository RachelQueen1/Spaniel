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
#' @param Grob
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
#' @param Gene 
#' @param NoClusters
#' @customT
#' @export
#' @examples
#' 
#' ## Data is taken from DOI: 10.1126/science.aaf2403
#' exampleCounts <- readRDS(exampleCounts.rds)
#' exampleBarcodes <- "path to barcode"
#' ## Counts per spot 
#' ST_plot(Object = SeuratObj, Grob = imgFile, 
#'        PlotType = "CountsPerSpot", 
#'        ShowFilter = filter)
#'
#' ## Cluster plot
#' ST_plot(Object = SeuratFiltered, Grob = imgFile, 
#'        PlotType = "Cluster", 
#'        ClusterRes = "RNA_snn_res.0.6")
#'
## Gene plot
#' ST_plot(Object = SeuratFiltered, Grob = imgFile, 
#'        PlotType = "Gene", 
#'        Gene = "Pax6")


###### FUNCTIONS ############
ST_plot <- function (Object,  Grob,
                     PlotType = c("NoGenes", 
                                  "CountsPerSpot", 
                                  "Cluster", 
                                  "Gene"),
                     Gene = NULL, 
                     ClusterRes = NULL, 
                     CustomTitle = NULL,
                     ScaleData = T, 
                     ShowFilter = NULL)
    
{
    MetaData = getMetadata(Object)
    Coordinates = MetaData[, c("x", "y")]
    Coordinates$spot = rownames(MetaData)
    show_size_legend = T
    
    if (PlotType == "NoGenes") {
        plotTitle = "Number of Genes Per Spot"
        cl = "No_Of_Genes"
        sz = "No_Of_Genes"
        
        ### get column to plot
        if (class(Object) == "Seurat"){
            colPlot = "nFeature_RNA"}
        if (class(Object) == "SingleCellExperiment"){
            colPlot = "total_features_by_counts"}
        
        tmp = data.frame(No_Of_Genes = MetaData[, colPlot],
                         spot = rownames(MetaData)) %>% 
            as.tbl %>% 
            inner_join(Coordinates)
        
        if (!is.null(ShowFilter)){
            tmp$Filter = ShowFilter
            cl = "Filter"
        }
        
    }
    if (PlotType == "CountsPerSpot") {
        cl = "Exprs"
        sz = "Exprs"
        plotTitle = "Total Counts Per Spot"
        
        ### get column to plot
        if (class(Object) == "Seurat"){
            colPlot = "nCount_RNA"}
        if (class(Object) == "SingleCellExperiment"){
            colPlot = "total_counts"}
        
        tmp = data.frame(Exprs = MetaData[, colPlot],
                         spot = rownames(MetaData)) %>% 
            as.tbl %>% 
            inner_join(Coordinates)
        
        ## add QC filter
        if (!is.null(ShowFilter)){
            tmp$Filter = ShowFilter
            cl = "Filter"
        }
        
    }
    if (PlotType == "Cluster") {
        plotTitle = "Spot Clusters"
        cl = "cluster"
        sz = 2
        show_size_legend = FALSE
        
        tmp = data.frame(cluster = MetaData[, ClusterRes], 
                         spot = rownames(MetaData)) %>%
            as.tbl %>% inner_join(Coordinates)
    }
    
    if (PlotType == "Gene") {
        cl = "Exp"
        sz = "Exp"
        plotTitle = paste("Expression of", Gene)
        Exprs = getExprs(Object)
        
        if(Gene %in% rownames(Exprs) == FALSE){
            print(paste(Gene, "not found in data. Check rownames"))
            
        }
        
        tmp = data.frame(Exp = Exprs[Gene, ])
        tmp$spot = rownames(tmp)
    }
    
    tmp = tmp %>% as.tbl %>% inner_join(Coordinates)
    tmp$y = 36 - tmp$y
    if (!is.null(CustomTitle)) {
        plotTitle = CustomTitle
    }
    print(show_size_legend)
    #tmp$y = tmp$y * -1
    p = plotImage(Grob, tmp, 
              Colour = cl, 
              Size = sz, 
              PlotTitle = plotTitle, 
              ShowSizeLegend = show_size_legend)
    
    if (PlotType != "Cluster" & is.null(ShowFilter)){
    p = p +
        scale_colour_gradient(low="#ff3300", high="#ffff00")+
        scale_size(range=c(1,5)) +
        guides(color=guide_legend(), size = guide_legend())
    }
    
    
    p
    
    
    
}



