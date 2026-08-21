### Unit tests for processStudy.R functions

library(RAIDS)


#############################################################################
### Tests RAIDSparam class
#############################################################################

context("RAIDSparam class results")

test_that("create a RAIDSparam class with all default parameters should return an object", {

    ## New RAIDSparam with all default values
    paramTest <- new("RAIDSparam")

    expect_true(inherits(paramTest, "RAIDSparam"))
    ## TODO ADD test for slots
    expect_null(paramTest@studyDF)
    expect_null(paramTest@studyDFSyn)
})