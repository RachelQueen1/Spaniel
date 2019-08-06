#### PARSE IMAGE ####
#' This function parses a HE image to use as the background for plots
#' @param ImgFile Path to the image file
#' @return A rasterized grob
#' @export
#' @examples
#' imgFile = file.path(system.file(package = "Spaniel"),
#'                     "HE_Rep1_resized.jpg")
#' img = parseImage(imgFile)
parseImage =  function(ImgFile) {
    img <- jpeg::readJPEG(ImgFile)
    g <-grid::rasterGrob(
        img,
        interpolate = FALSE,
        width = grid::unit(1, "npc"),
        height = grid::unit(1, "npc")
    )
    return(g)
}