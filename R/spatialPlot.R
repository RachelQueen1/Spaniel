# Main Spaniel Plot Function
# ------------------------------------------------------------------------------


ST_plot <- function (Object,  Grob,
                     PlotType = c("NoGenes", 
                                  "CountsPerSpot", 
                                  "Cluster", 
                                  "Gene"),
                     Gene = NULL, 
                     ClusterRes = NULL, 
                     CustomTitle = NULL,
                     ScaleData = T, 
                     ShowFilter = NULL, 
                     pt.size = 2,
                     pt.size.min = 0, 
                     pt.size.max = 5)

{
    MetaData = getMetadata(Object)
    Coordinates = MetaData[, c("x", "y")]
    Coordinates$spot = rownames(MetaData)
    show_size_legend = T
    
    if (PlotType == "NoGenes") {
        # set plot title, cluster, size and shp
        ungroupVars(plotTitle,cl,sz) %=% setGeneVars()
        shp = NULL
        ### get column to plot
        colPlot <- getGeneColPlot(Object)
        tmp <- genePlotDF(Object, Coordinates, colPlot)
        if (!is.null(ShowFilter)){
            tmp$Filter = ShowFilter
            cl = "Filter"
        }
    }
    
    
    if (PlotType == "CountsPerSpot") {
        cl = "Exprs"
        sz = "Exprs"
        shp = NULL
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
        sz = pt.size
        shp = "cluster"
        show_size_legend = FALSE
        
        tmp = data.frame(cluster = MetaData[, ClusterRes], 
                         spot = rownames(MetaData)) %>%
            as.tbl %>% inner_join(Coordinates)
    }
    
    if (PlotType == "Gene") {
        cl = "Exp"
        sz = "Exp"
        shp = NULL
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
    
    ### FOR QC PLOTS
    if (!is.null(ShowFilter)){
        p = p +
            scale_size_continuous(range=c(pt.size.min,
                                          pt.size.max)) +
            guides(color=guide_legend(), size = guide_legend())
    }
    
    ### FOR GENE PLOTS AND NUMBER READS/COUNTS
    if (PlotType != "Cluster" & is.null(ShowFilter)){
        p = p +
            scale_colour_gradient(low="#ff3300", high="#ffff00")+
            scale_size_continuous(range=c(pt.size.min,
                                          pt.size.max)) +
            guides(color=guide_legend(), size = guide_legend())
        
    }
    ### FOR CLUSTER PLOTS
    if (PlotType == "Cluster"){
        p = p + ggplot2::guides(size=FALSE) +
            scale_size_continuous(range=c(pt.size,
                                          pt.size))
    }
    
    # if (PlotType != "Cluster" & is.null(ShowFilter)){
    #     p = p +
    #         scale_colour_gradient(low="#ff3300", high="#ffff00")+
    #         #scale_size_continuous(range=c(pt.size.min,
    #         #                   pt.size.max)) +
    #         guides(color=guide_legend(), size = guide_legend())
    # }
    # if (PlotType == "Cluster"){
    #     p = p + ggplot2::guides(size=FALSE)
    # } 
    
    
    
    
    p
    
    
    
}














