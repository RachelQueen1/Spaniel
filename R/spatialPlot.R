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
                     ScaleData = T, 
                     ShowFilter = NULL, 
                     pt.size = 2,
                     pt.size.min = 0, 
                     pt.size.max = 5)

{
    ### Validate Object is either a Seurat or SCE object
    testObject(Object)
    
    ### set title, colour column, and size column according to plotType
    ungroupVars(plotTitle,cl,sz,shp, show_size_legend, colPlot) %=% 
        setVars(Object, PlotType, pt.size, Gene)
    
            # convert shp NULL
            if(shp == "NULL"){shp = NULL}
            # convert size legend to logical
            show_size_legend = ifelse(show_size_legend == "TRUE", T, F)
            # convert sz to numeric
            if (suppressWarnings(!is.na(as.numeric(sz)))){
                sz = as.numeric(sz)
            }
    
    ### create data.frame for ggplot
    tmp <- make_ggdf(Object, PlotType, colPlot, cl)
    
    
    ### Update tmp optional arguments if supplied
    if (!is.null(ShowFilter)){
        tmp$Filter = ShowFilter
        cl = "Filter"
    }
    
    if (!is.null(CustomTitle)) {
        plotTitle = CustomTitle
    }
    
    ### Create plot
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
    
    return(p)
}














