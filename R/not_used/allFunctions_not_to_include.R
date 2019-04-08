## read counts
readCounts = function(DataDir = dataDir){
  counts = read.csv(file.path(DataDir, "counts.csv"), sep = "\t", header = F)
  rownames(counts) = counts[,1]
  counts = counts[,-1]
  spotNames = read.csv(file.path(DataDir, "colnames.txt"), sep = " ",  header = F, stringsAsFactors = F ) %>%
    t() %>%
    data.frame() %>%
    .[,1] %>%
    as.character() %>%
    gsub("counts__", "", .) %>%
    gsub(".txt", "", .)
  colnames(counts) = spotNames
  counts = counts[!grepl("__", rownames(counts)), ]
  return(counts)}

# create an annotation table
annoCounts = function(Counts, Dataset = c('mmusculus_gene_ensembl', 'hsapiens_gene_ensembl') ){
  geneNames = gsub("\\..*", "", rownames(Counts))

  ##Query Biomart
  mart <- biomaRt::useMart(biomart = 'ensembl', dataset = Dataset)
  annotation = biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                              filters = "ensembl_gene_id",
                              values = geneNames,
                              mart = mart)
  rownames(annotation) = annotation$ensembl_gene_id

  #Match Biomart results with counts
  #anno$ensembl_gene_id[is.na(anno$ensembl_gene_id)]
  anno_table = data.frame(ensemblID = rownames(Counts))
  anno_table$gene_symbol = anno$external_gene_name[match( gsub("\\..*", "", anno_table$ensemblID) , anno$ensembl_gene_id)]
  anno_table$gene_symbol[is.na(anno_table$gene_symbol)] = anno_table$ensemblID[is.na(anno_table$gene_symbol)]

  Counts$gene = anno_table$gene_symbol

  tmp = Counts %>%
    group_by(gene) %>%
    summarise_all(funs(sum)) %>%
    data.frame()
  rownames(tmp) = tmp[,1]
  Counts = tmp[,-1]
  return(Counts)
}


readSeurat = function(DataDir, Genome, BarcodeFile, Project = sectionNo){
  counts = readCounts(DataDir)
  counts = annoCounts(counts, Genome)
  print(BarcodeFile)
  barcodes = read.csv(BarcodeFile, sep = "\t", header = F)
  rownames(barcodes) = barcodes$V1
  barcodes = barcodes[colnames(counts), ]
  colnames(barcodes) = c("Barcode", "x", "y")
  seuratObj = Seurat::CreateSeuratObject(raw.data = counts,
                                         meta.data = barcodes,
                                         project = Project)

  return(seuratObj)

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


ST_plot_V2 = function (SeuratObj = seuratObj, Coordinates = coordinates, Grob = g,
          plot_type = c("NoGenes", "CountsPerSpot", "Cluster", "Gene"),
          Features = NULL, customTitle = NULL, ScaleData = T)
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
    tmp = data.frame(cluster = SeuratObj@ident, spot = names(SeuratObj@ident)) %>%
      as.tbl %>% inner_join(Coordinates)
  }
  if (plot_type == "Gene") {
    cl = "Exp"
    sz = "Exp"
    if (length(Features) == 1) {
      Gene = Features
      plotTitle = paste("Expression of", Gene)
      tmp = data.frame(Exp = SeuratObj@scale.data[Gene, ])
      tmp$spot = rownames(tmp)
    }
    if (length(Features) > 1) {
      gl = gsub(", $", "", paste0(Features, collapse = "",
                                  sep = ", "))
      plotTitle = paste("Total Expression of", gl)
      Exp = colSums(data.frame(SeuratObj@assays$RNA@scale.data[Features,
                                                               ]))
      tmp = data.frame(Exp)
      tmp$spot = colnames(SeuratObj@assays$RNA@data)
    }
    tmp = tmp %>% as.tbl %>% inner_join(Coordinates)
  }
  tmp$y = 36 - tmp$y
  if (!is.null(customTitle)) {
    plotTitle = customTitle
  }
  print(show_size_legend)
  #tmp$y = tmp$y * -1
  plotImage(Grob, tmp, Colour = cl, Size = sz, PlotTitle = plotTitle, ShowSizeLegend = show_size_legend)
}
