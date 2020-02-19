#' @include spanielPlotInternals.R
#' @include utilities.R

# Tests for that correct plots are produced
# ------------------------------------------------------------------------------


metaData <- data.frame(pixel_x = seq(0, 1000, by = 100),
                       pixel_y = seq(0, 1000, by = 100), 
                       x = seq(1, 11, by = 1), 
                       y = rep(12, 11))




# Test makeCoordinates selects correct columns (techType = "Original",
# byCoord = FALSE)
# ------------------------------------------------------------------------------
# expected Y Coodinates
z = 36 - metaData[1,"y"]

coordinates <- makeCoordinates(metaData, 
                               techType = "Original",
                               byCoord = FALSE,
                               imgDims = NULL)

test_that("makeCoordinates, orginal, false", {
  expect_equal(coordinates[1,1], 1)
  expect_equal(coordinates[1,2], z)
  expect_equal(colnames(coordinates), c("x", "y"))
})


# Test makeCoordinates selects correct columns (techType = "Original",
# byCoord = FALSE)
# ------------------------------------------------------------------------------
# expected Y Coodinates
z = 1000 - metaData[1,"pixel_y"]

coordinates <- makeCoordinates(metaData, 
                               techType = "Original",
                               byCoord = TRUE,
                               imgDims = c(1000, 1000))

test_that("makeCoordinates, orginal, true", {
  expect_equal(coordinates[1,1], 0)
  expect_equal(coordinates[1,2], z)
  expect_equal(colnames(coordinates), c("x", "y"))
})


# Test makeCoordinates selects correct columns (techType = "Visium",
# byCoord = TRUE)
# ------------------------------------------------------------------------------
# expected Y Coodinates
z = metaData[1,"pixel_y"]

coordinates <- makeCoordinates(metaData, 
                               techType = "Visium",
                               byCoord = TRUE,
                               imgDims = c(1000, 1000))

test_that("makeCoordinates, visium, true", {
  expect_equal(coordinates[1,1], 0)
  expect_equal(coordinates[1,2], z)
  expect_equal(colnames(coordinates), c("x", "y"))
})



# Test makeCoordinates selects correct columns (techType = "Visium",
# byCoord = TRUE)
# ------------------------------------------------------------------------------
# expected Y Coodinates
z = metaData[1,"pixel_y"]

coordinates <- makeCoordinates(metaData, 
                               techType = "Visium",
                               byCoord = FALSE,
                               imgDims = c(1000, 1000))

test_that("makeCoordinates, visium, false", {
  expect_equal(coordinates[1,1], 0)
  expect_equal(coordinates[1,2], z)
  expect_equal(colnames(coordinates), c("x", "y"))
})



# Test plotRange selects correct columns (techType = "Original",
# byCoord = FALSE)
# ------------------------------------------------------------------------------

imgDims <-  c(1000, 1000)

plotDims <- pointRange(techType = "Original", 
                       byCoord = FALSE, 
                       imgDims = NULL)

ungroupVars(x_min,x_max,y_min,y_max) %=% 
  plotDims


test_that("pointRange, orginal, false", {
  expect_equal(x_min, 1)
  expect_equal(x_max, 33)
  expect_equal(y_min, 1)
  expect_equal(y_max, 35)
})

# Test plotRange selects correct columns (techType = "Original",
# byCoord = TRUE)
# ------------------------------------------------------------------------------


plotDims <- pointRange(techType = "Original", 
                       byCoord = TRUE, 
                       imgDims = imgDims)

ungroupVars(x_min,x_max,y_min,y_max) %=% 
  plotDims


test_that("pointRange, orginal, true", {
  expect_equal(x_min, 0)
  expect_equal(x_max, 1000)
  expect_equal(y_min, 0)
  expect_equal(y_max, 1000)
})


# Test plotRange selects correct columns (techType = "Visium",
# byCoord = TRUE)
# ------------------------------------------------------------------------------


plotDims <- pointRange(techType = "Visium", 
                       byCoord = TRUE, 
                       imgDims = imgDims)

ungroupVars(x_min,x_max,y_min,y_max) %=% 
  plotDims


test_that("pointRange, orginal, true", {
  expect_equal(x_min, 0)
  expect_equal(x_max, 1000)
  expect_equal(y_min, 0)
  expect_equal(y_max, 1000)
})

