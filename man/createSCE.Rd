% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createObjects.R
\name{createSCE}
\alias{createSCE}
\title{Create a SingleCellExperiment Object From Spatial Transcriptomics Data}
\usage{
createSCE(counts, barcodeFile, projectName=projectName, 
                sectionNumber=sectionNo)
}
\arguments{
\item{counts}{Raw count matrix or data frame where each row represents a
gene and each column represents barcoded location on a spatial
transcriptomics slide. The columns should be named using the
spot barcode (eg "GTCCGATATGATTGCCGC")}

\item{barcodeFile}{a tab seperated barcode file supplied by Spatial
Trancscriptomics. The file should contains three column:
The first column contains the Spatial Transcriptomics barcode,
the second and third column equate to the x and y location}

\item{projectName}{The name of the project which is stored in the Seurat
Object.}

\item{sectionNumber}{The location of the sample on the slide}
}
\value{
A SingleCellExeriment Object
}
\description{
This function converts a count matrix into a SingleCellExperiment object. 
The barcodes for each spot are added to the coldata of the 
SingleCellExperiment object and are used in  plotting the data.
}
\examples{
## Data is taken from DOI: 10.1126/science.aaf2403
examplecounts <- readRDS(file.path(system.file(package = "Spaniel"),
                            "extdata/counts.rds"))
exampleBarcodes <- file.path(system.file(package = "Spaniel"),
                            "1000L2_barcodes.txt")
seuratOb <- createSCE(examplecounts,
                        exampleBarcodes,
                        projectName = "TestProj",
                        sectionNumber = 1)
}
