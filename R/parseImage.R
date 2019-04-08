#### PARSE IMAGE ####
#' This function parses a HE image to use as the background for plots
#' @param ImgFile Path to the image file
#' @export
#' @examples
#' imgFile = "data/HE_Rep1.jpg"
#' img = parseImage(imgFile)
parseImage =  function(ImgFile){
    img <- jpeg::readJPEG(ImgFile)
    g <- grid::rasterGrob(img, interpolate = FALSE, width=grid::unit(1,"npc"), 
                          height=grid::unit(1,"npc"))
    return(g)
}