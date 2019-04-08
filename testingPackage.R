library(Spaniel)


# rObjects = "../../Spatial_Transcriptomics_Report/Analysis_outs/SRR3382371/rObjects"
# # sObj= readRDS(file.path(rObjects, "seurObj.rds"))
# # imgFile = readRDS(file.path(rObjects, "image.rds"))
# 
# # sObjFiltered = readRDS(file.path(rObjects, "sObj_filtered.rds"))
# countsFiltered = sObjFiltered@raw.data[, colnames(sObjFiltered@scale.data)]
# saveRDS(countsFiltered, "testData/countsFiltered.rds")

counts = readRDS("testData/counts.rds")
countsFiltered = readRDS("testData/countsFiltered.rds")

imgFile = readRDS("testData/image.rds")

barcodesFile = "testData/1000L2_barcodes.txt"
#saveRDS(counts, "testData/counts.rds")

SeuratObj = readSeurat(countsFiltered,
                       barcodesFile, 
                       ProjectName = "TestProj", 
                       SectionNumber = 1
                       )

sce = readSCE(Counts = countsFiltered,
              BarcodeFile = barcodesFile, 
              ProjectName = "TestProj", 
              SectionNumber = 1)




sce = scater::calculateQCMetrics(
    sce
)

sce =  scater::normalise(sce)





minGenes = 2000
minUMI = 300000

# filter = SeuratObj@meta.data$nFeature_RNA > minGenes & 
#     SeuratObj@meta.data$nCount_RNA > minUMI
# 
# 
# 
SeuratFiltered <- subset(x = SeuratObj, subset = nFeature_RNA > minGenes &
                             nCount_RNA > minUMI)


SeuratFiltered <- NormalizeData(object = SeuratFiltered,
                           normalization.method = "LogNormalize",
                           scale.factor = 10000)

SeuratFiltered <- FindVariableFeatures(object = SeuratFiltered,
                                  selection.method = "vst",
                                  nfeatures = 2000)

all.genes <- rownames(x = SeuratObj)
SeuratFiltered <- ScaleData(object = SeuratFiltered, features = all.genes)
SeuratFiltered <- RunPCA(object = SeuratFiltered,
                    features = VariableFeatures(object = SeuratFiltered))


SeuratFiltered <- FindNeighbors(object = SeuratFiltered, dims = 1:10)
SeuratFiltered <- FindClusters(object = SeuratFiltered,
                          resolution = c(0.4, 0.5, 0.6, 0.8))


### Find markers
Markers <- FindAllMarkers(object = SeuratFiltered, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)
Markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top100 <- Markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)


###annotate top100
library(biomaRt)
annotate = function(geneNames){
    mart <- useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')

annotation = getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                   filters = "ensembl_gene_id",
                   values = geneNames,
                   mart = mart)

#set the row names to be ensembl_gene_id
#rownames(annotation) = annotation$ensembl_gene_id


#annotation = annotation[geneNames,]
return(annotation)
}


anno = annotate(substr(top100$gene, 1, 18))
mm = match((substr(top100$gene, 1, 18)), anno$ensembl_gene_id)
top100$geneShort = (substr(top100$gene, 1, 18))
top100$geneSymbol = anno$external_gene_name[mm]


##subset seurat
genesUse = rownames(SeuratFiltered)  %in%  unique(top100$gene) 
newCounts = SeuratFiltered@assays$RNA[genesUse, ] %>% data.frame()
mm = match(rownames(newCounts), top100$gene)
rownames(newCounts) = top100$geneSymbol[mm]

s = CreateSeuratObject(counts = newCounts, 
                   meta.data = SeuratFiltered@meta.data)

s@active.ident = SeuratFiltered@active.ident[genesUse]


newScaleData = SeuratFiltered@assays$RNA@scale.data[genesUse, ]
rownames(newScaleData) = rownames(s)
s@assays$RNA@scale.data = newScaleData
s = markClusterCol(s, "res")
SeuratObj = s
saveRDS(SeuratObj, "testData/SeuratData.rds")



##subset SCE
genesUse = rownames(sce)  %in%  unique(top100$gene) 
sce = sce[genesUse, ]
rownames(sce)
mm = match(rownames(sce), top100$gene)
rownames(sce) = top100$geneSymbol[mm]
saveRDS(sce, "testData/sceData.rds")



### Example plots
## Counts per spot 
ST_plot(Object = SeuratObj, 
        Grob = imgFile, 
        PlotType = "NoGenes", 
        ShowFilter = filter)

ST_plot(Object = sce, 
        Grob = imgFile, 
        PlotType = "NoGenes")

ST_plot(Object = sce, 
        Grob = imgFile, 
        PlotType = "CountsPerSpot")


ST_plot(Object = SeuratObj, Grob = imgFile, 
        PlotType = "CountsPerSpot", 
        ShowFilter = filter)

## Cluster plot
ST_plot(Object = SeuratFiltered, Grob = imgFile, 
        PlotType = "Cluster", 
        ClusterRes = "RNA_snn_res.0.6"
        )

ST_plot(Object = sce, Grob = imgFile, 
        PlotType = "Cluster", 
        ClusterRes = "y"
)

## Gene plot
gene = "ENSMUSG00000050010.8"
ST_plot(Object = SeuratFiltered, Grob = imgFile, 
        PlotType = "Gene", 
        Gene = gene)


ST_plot(Object = sce, Grob = imgFile, 
        PlotType = "Gene", 
        Gene = gene)



sce = readSCE(Counts = counts,
              BarcodeFile = barcodesFile, 
              ProjectName = "TestProj", 
              SectionNumber = 1
)





