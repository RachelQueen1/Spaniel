###### FUNCTIONS ############
ST_plot = function (SeuratObj = seuratObj, 
                    Coordinates = coordinates, 
                    Grob = g,
                    plot_type = c("NoGenes", "CountsPerSpot", "Cluster", "Gene"),
                    Gene = NULL, NoClusters = NULL, customTitle = NULL, ScaleData = T)
{
  show_size_legend = T
  if (plot_type == "NoGenes") {
    plotTitle = "Number of Genes Per Spot"
    cl = "No_Of_Genes"
    sz = "No_Of_Genes"
    #Counts = as.matrix(SeuratObj@assays$RNA@counts)
    #is.exp = Counts > 0
    tmp = data.frame(No_Of_Genes = SeuratObj@meta.data$nGene,
                     spot = rownames(SeuratObj@meta.data)) %>% as.tbl %>% inner_join(Coordinates)
  }
  if (plot_type == "CountsPerSpot") {
    cl = "Exprs"
    sz = "Exprs"
    if (ScaleData == F) {
      plotTitle = "Total Counts Per Spot (raw)"
      Counts = as.matrix(SeuratObj@raw.data)
    }
    if (ScaleData == T) {
      plotTitle = "Total Expression Per Spot (scaled)"
      Counts = as.matrix(SeuratObj@scale.data)
    }
    tmp = data.frame(Exprs = as.integer(colSums(Counts)),
                     spot = rownames(SeuratObj@meta.data)) %>% as.tbl %>% inner_join(Coordinates)
  }
  if (plot_type == "Cluster") {
    plotTitle = "Spot Clusters"
    cl = "cluster"
    sz = 2
    show_size_legend = FALSE
    
    tmp = data.frame(cluster = SeuratObj@meta.data[, NoClusters], spot = names(SeuratObj@ident)) %>%
      as.tbl %>% inner_join(Coordinates)
  }
  
  
  if (plot_type == "Gene") {
    cl = "Exp"
    sz = "Exp"
      plotTitle = paste("Expression of", Gene)
      tmp = data.frame(Exp = SeuratObj@scale.data[Gene, ])
      tmp$spot = rownames(tmp)
  }
  
  tmp = tmp %>% as.tbl %>% inner_join(Coordinates)
  tmp$y = 36 - tmp$y
  if (!is.null(customTitle)) {
    plotTitle = customTitle
  }
  print(show_size_legend)
  #tmp$y = tmp$y * -1
  plotImage(Grob, tmp, Colour = cl, Size = sz, PlotTitle = plotTitle, ShowSizeLegend = show_size_legend)
}




#### PLOT IMAGE ####
plotImage = function(Grob, Tmp, Colour, Size, ShowSizeLegend = TRUE, PlotTitle = NULL){
  p = ggplot2::ggplot(Tmp ,ggplot2::aes_string("x", "y", color = Colour, size = Size)) +
    ggplot2::xlim(1, 33) +
    ggplot2::ylim(1, 35) +
    ggplot2::annotation_custom(Grob, xmin = 1, xmax = 33, ymin = 1, ymax = 35) +
    ggplot2::geom_point()  +
    ggplot2::labs(title = PlotTitle) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  NULL
  
  ### if show size false
  if (ShowSizeLegend == FALSE){
    p = p + guides(size=FALSE)}
  ### show plot
  p
}