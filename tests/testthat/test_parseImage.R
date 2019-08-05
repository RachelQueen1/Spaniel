# Tests for parseImage function
# ------------------------------------------------------------------------------

context("Testing parseImage")
# These tests were created to ensure that the parseImage functions works
# correctly and creates a grob

# Test parseImage 
# ------------------------------------------------------------------------------

# load image
imgPath <-  file.path(system.file(package = "Spaniel"), 
                                 "HE_Rep1_resized.jpg")
test.grob <- parseImage(imgPath)

test_that("parseImage loads an image and creates a grob", {
              expect_is(test.grob, c("rastergrob","grob","gDesc"))
                        })
