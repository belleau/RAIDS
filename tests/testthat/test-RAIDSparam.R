### Unit tests for processStudy.R functions

library(RAIDS)


#############################################################################
### Tests RAIDSparam class
#############################################################################

context("RAIDSparam class results")

test_that("create a RAIDSparam class with all default parameters should return an object", {

    exp_studyDF <- data.frame(study.id="NotDef",
                                study.desc="NotDef",
                                study.platform="NotDef",
                                stringsAsFactors=FALSE)
    
    exp_studyDFSyn <- data.frame(study.id="NotDef.Synthetic",
                            study.desc="NotDef synthetic data",
                            study.platform="Synthetic", stringsAsFactors=FALSE)


    exp_pedStudy <- data.frame(Name.ID=c("ProfileId"), Case.ID=c("ProfileId"),
                        Sample.Type=c("type"), Diagnosis="NotDef",
                        Source=c("NotDef"), stringsAsFactors=FALSE, 
                        row.names=c("ProfileId"))

    ## New RAIDSparam with all default values
    paramTest <- new("RAIDSparam")

    expect_true(inherits(paramTest, "RAIDSparam"))
    ## TODO ADD test for slots

    ## Test studyDF
    expect_true(is.data.frame(paramTest@studyDF))
    expect_identical(paramTest@studyDF, exp_studyDF)

    ## Test studyDFSyn
    expect_true(is.data.frame(paramTest@studyDFSyn))
    expect_identical(paramTest@studyDFSyn, exp_studyDFSyn)

    ## Test pedStudy
    expect_true(is.data.frame(paramTest@pedStudy))
    expect_identical(paramTest@pedStudy, exp_pedStudy)

    ## Test studyType
    expect_true(paramTest@studyType == "LD")

    ## Test blockTypeId
    expect_true(paramTest@blockTypeId == "GeneS.Ensembl.Hsapiens.v86")

    ## Test reference
    expect_true(paramTest@reference == "1KGc1.0")

    ## Test genome
    expect_true(paramTest@genome == "HG38")
})