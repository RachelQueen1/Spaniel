# Tests for parseImage function
# ------------------------------------------------------------------------------

context("Testing parseImage")
# These tests were created to ensure that the markClusterCol functions works
# correctly and marks the correct column with Cluster_ prefix

# Test parseImage 
# ------------------------------------------------------------------------------

# load image
imgPath <-  "../inst/HE_Rep1_resized.jpg"
test.grob <- parseImage(imgPath)

test_that("parseImage loads an image and creates a grob", {
              expect_is(test.grob, c("rastergrob","grob","gDesc"))
                        })
