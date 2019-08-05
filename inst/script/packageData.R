#' #' Data is taken from DOI: 10.1126/science.aaf2403
#' count.rds 
#' To create the expression matrix, from the sequencing data, the paired FASTQ 
#' files were demultiplexed with a publically available perl script  
#' (https://github.com/tallulandrews/scRNASeqPipeline/blob/master/0_custom_undo_demultiplexing.
#'  pl) using the spatial barcodes encoded in read 1. 
#' Read 2 from successfully demultiplexed pairs were
#' trimmed for quality using Trimmomatic version 0.36 6. 
#' A reference was created using Ensembl  mouse reference genome Release M20 
#' (GRCm38.p6). 
#' The trimmed reads were aligned to this reference using STAR version 2.5.3a  
#' in single read alignment mode. The number of reads were quantified using HTSEQ 
#' version 0.6.1 8 and a count matrix was created. 
#' The matrix was then read into R:
# counts <- read.csv("counts_matrix.csv")


#' image.rds
#' The H & E images were cropped to region around the edges of the spots and 
#' resized to 1000 x 1071 pixels with a resolution of 72dpi using a photo editor.
