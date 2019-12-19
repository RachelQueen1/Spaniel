## get image dimensions
getImage <- function(spatialDir, resolution){
    if(resolution == "High"){
        img <- png::readPNG(file.path(spatialDir, "tissue_hires_image.png"))
    }else if(resolution == "Low"){
        img <- png::readPNG(file.path(spatialDir, "tissue_lowres_image.png"))
    }
    return(img)
}

## read 10X barcodes
read10XBarcodes <- function(spatialDir, imageDims, resolution){
    barcodes <- read.csv(file.path(spatialDir, "tissue_positions_list.csv"),
                         stringsAsFactors = FALSE, header = FALSE)
    colnames(barcodes) <- c("Barcode", "Section", "Spot_Y", 
                            "Spot_X", "Image_Y", "Image_X")
    #barcodes$barcode <- gsub("-.*", "", barcodes$barcode)
    scaleFactors <- jsonlite::fromJSON(txt = file.path(spatialDir, 
                                             "scalefactors_json.json"))
    
    if(resolution == "High"){
        barcodes$pixel_x <- barcodes$Image_X * scaleFactors$tissue_hires_scalef
        barcodes$pixel_y <- barcodes$Image_Y * scaleFactors$tissue_hires_scalef
    }else if(resolution == "Low"){
        barcodes$pixel_x <- barcodes$Image_X * scaleFactors$tissue_lowres_scalef
        barcodes$pixel_y <- barcodes$Image_Y * scaleFactors$tissue_lowres_scalef
    }
    ## invert Y coordinates 
    barcodes$pixel_y <- imageDims[2] - barcodes$pixel_y
    return(barcodes)
}


parseVisumImage <-  function(imgFile) {
    img <- png::readJPEG(imgFile)
    g <-grid::rasterGrob(
        img,
        interpolate = FALSE,
        width = grid::unit(1, "npc"),
        height = grid::unit(1, "npc")
    )
    return(g)
}

#' createVisiumSCE
#' 
#' A function to select to import 10X data into an SCE object
#' @param tenXDir The path to Space Ranger outs directory containing spatial 
#' directory and filtered_feature_bc_matric
#' @param resolution Resolution of the tissue image to be used for plotting. 
#' Can be either "High", or "Low". Default is "Low".
#' @examples
#' sce <- createVisiumSCE(tenXDir, resolution)
#' 
#' 
#' @export
createVisiumSCE <- function(tenXDir="../outs", 
                            resolution="Low"){

    filteredDir <- file.path(tenXDir, "filtered_feature_bc_matrix")
    
    sce <- DropletUtils::read10xCounts(filteredDir, 
                                       version = "3")
    
    
    ### parse image, create grob and get image dimensions
    spatialDir =  file.path(tenXDir, "spatial")
    img <- getImage(spatialDir, resolution)
    imgageDims <- dim(img)[1:2]
    grob <-grid::rasterGrob(
        img,
        interpolate = FALSE,
        width = grid::unit(1, "npc"),
        height = grid::unit(1, "npc")
    )
    
    ### add image dimensions to sce
    sce@metadata$ImgDims <- c(imgDims, resolution)
    
    ### add grob to sce
    sce@metadata$Grob <- grob
    
    ### read barcodes and add to sce
    barcodes <- read10XBarcodes(spatialDir, imgDims, resolution)
    
    sce@colData <- sce@colData %>% 
        data.frame() %>% 
        left_join(barcodes) %>% 
        DataFrame()
    
    ## add QC metrics to sce
    per.cell <- scater::perCellQCMetrics(sce)
    colData(sce) <- cbind(colData(sce), per.cell)
    
    
    
    return(sce)
    }








### 1 add option to read in png to parseImage function
img <- png::readPNG("tiny/spatial/tissue_hires_image.png")

g <- grid::rasterGrob(img, 
                      interpolate = FALSE, 
                      width = grid::unit(1, "npc"), 
                      height = grid::unit(1, "npc"))

toPlot <- colData(sce) %>% data.frame

ggplot(toPlot, aes(highResX, highResY))  +
    ggplot2::annotation_custom(g, xmin = 0, xmax = dim(img)[1], 
                               ymin = 0, ymax = dim(img)[2])  + xlim(0,  
                                                                     dim(img)[1]) + 
    ylim(0, dim(img)[2]) + geom_point()


scale <- "hires"
scale <- match.arg(arg = scale, choices = c('spot', 'fiducial', 'hires', 'lowres'))



dim(img)
scale.factors <- fromJSON(txt = file.path("tiny/spatial/scalefactors_json.json"))
spanielPlot(sObj, g, plotType = c("NoGenes"),
            gene= NULL, clusterRes = NULL, customTitle = NULL,
            scaleData = TRUE, showFilter = NULL, ptSize = 2,
            ptSizeMin = 0, ptSizeMax = 5)


barcodes$V5 %>% range()
