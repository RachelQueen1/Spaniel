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
    ### Validate Object
    testObject()
    
    ### Get Coodinates
    MetaData = getMetadata(Object)
    Coordinates = MetaData[, c("x", "y")]
    Coordinates$spot = rownames(MetaData)
    
    ### set title, colour column, and size column according to plotType
    ungroupVars(plotTitle,cl,sz,shp) %=% 
        setVars(PlotType, pt.size, Gene)
    
    ## convert shp to NULL
    if(shp == "NULL"){shp = NULL}
    
    ## convert size legend to logical
    show_size_legend = ifelse(show_size_legend == "TRUE", T, F)
    
    ## convert sz to numeric
    if (suppressWarnings(!is.na(as.numeric(sz)))){
        sz = as.numeric(sz)
    }
    
    
    ## create data.frame
    
    # 1 get column to plot
    colPlot <- getColPlot(Object, PlotType, ClusterRes)
    
    # 2 create tmp
    
    ## 3 add QC filter
    if (!is.null(ShowFilter)){
        tmp$Filter = ShowFilter
        cl = "Filter"
    }
    
    
    
    
    
    
    ## add QC filter
    
    
    
    
    if (PlotType == "CountsPerSpot") {
        
        
## add QC filter
if (!is.null(ShowFilter)){
    tmp$Filter = ShowFilter
    cl = "Filter"
}
    }
    
    
    if (PlotType == "Cluster") {
        
    }
    
    if (PlotType == "Gene") {
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














