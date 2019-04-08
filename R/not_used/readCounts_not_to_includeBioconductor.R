# create an annotation table
annoCounts <- function(Counts, Dataset = c('mmusculus_gene_ensembl',
                                           'hsapiens_gene_ensembl')){

    geneNames = rownames(Counts)

    if (Dataset == "hsapiens_gene_ensembl"){
        geneNames =  substr(geneNames, 1,15)
    }
    if (Dataset == "mmusculus_gene_ensembl"){
        geneNames =  substr(geneNames, 1,18)
    }

    ##Query Biomart
    mart <- biomaRt::useMart(biomart = 'ensembl', dataset = Dataset)
    annotation <- biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                                filters="ensembl_gene_id",
                                values=geneNames,
                                mart=mart)
    rownames(annotation) = annotation$ensembl_gene_id

    #Match Biomart results with counts
    #anno$ensembl_gene_id[is.na(anno$ensembl_gene_id)]
    anno_table = data.frame(ensemblID = geneNames)
    mm = match(anno_table$ensemblID, annotation$ensembl_gene_id)
    anno_table$gene_symbol = annotation$ensembl_gene_id[mm]


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


## Read Counts
readCounts = function(DataDir=dataDir, Genome){
    counts = read.csv(file.path(DataDir, "counts.csv"), sep = "\t", header = F)
    rownames(counts) = counts[ , 1]
    counts = counts[ , -1]
    spotNames = read.csv(file.path(DataDir, "colnames.txt"),
                         sep = " ",  header = F, stringsAsFactors = F ) %>%
        t() %>%
        data.frame() %>%
        .[ , 1] %>%
        as.character() %>%
        gsub("counts__", "", .) %>%
        gsub(".txt", "", .)

    colnames(counts) = spotNames
    counts <-counts[!grepl("__", rownames(counts)), ]
    counts <- annoCounts(counts, Genome)
    return(counts)
    }



