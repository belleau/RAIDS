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

    exp_chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L, 181538259L,
        170805979L, 159345973L, 145138636L, 138394717L, 133797422L, 135086622L, 
        133275309L, 114364328L, 107043718L, 101991189L,  90338345L,  83257441L, 
        80373285L,  58617616L,  64444167L,  46709983L,  50818468L, 156040895L,
        57227415L, 16569L)
    names(exp_chrInfo) <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
    
    exp_paramAncestry <- list(ScanBamParam=NULL, PileupParam=NULL,
                                yieldSize=10000000)
    
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

    ## Test genoSource 
    expect_null(paramTest@genoSource)

    ## Test blockTypeId
    expect_true(paramTest@blockTypeId == "GeneS.Ensembl.Hsapiens.v86")

    ## Test reference
    expect_true(paramTest@reference == "1KGv1.0")

    ## Test genome
    expect_true(paramTest@genome == "HG38")

    ## Test chrInfo
    expect_true(all(paramTest@chrInfo == exp_chrInfo))

    ## Test paramAncestry
    expect_identical(paramTest@paramAncestry, exp_paramAncestry)

    ## Test profileFile
    expect_null(paramTest@profileFile)

        ## Validate the profileFileGeno parameter TODO
    
    ## Test pathProfileGDS 
    expect_null(paramTest@pathProfileGDS)
    
    ## Test fileReferenceGDS 
    expect_null(paramTest@fileReferenceGDS)

        ## Validate the fileReferenceAnnotGDS parameter TODO
    
    ## Test inferenceType
    expect_identical(paramTest@inferenceType, "PCAknn")

        ## Validate the sampleRef parameter TODO

    ## Test batch
    expect_true(paramTest@batch == 1L)

    ## Test prefix
    expect_equal(paramTest@prefix, "1")

    ## Test nbSim
    expect_true(paramTest@nbSim == 1L)

    ## Test offset
    expect_true(paramTest@offset == -1L)

    ## Test minCov
    expect_true(paramTest@minCov == 10L)

    ## Test minProb
    expect_true(paramTest@minProb == 0.999)

    ## Test seqError
    expect_true(paramTest@seqError == 0.001)

    ## Test seqErrorSyn
    expect_true(paramTest@seqErrorSyn == 0.001)

    ## Test pRecomb
    expect_true(paramTest@pRecomb == 0.01)

    ## Test np
    expect_true(paramTest@np == 1L)
  
    ## Test listPos
    expect_null(paramTest@listPos)

    ## Test syntheticRefDF
    expect_null(paramTest@syntheticRefDF)
    
    ## Test pruningMethod
    expect_equal(paramTest@pruningMethod, "corr")

    ## Test slideWindowMaxBP 
    expect_equal(paramTest@slideWindowMaxBP, 500000L)

    ## Test thresholdLD 
    expect_equal(paramTest@thresholdLD, sqrt(0.1))

    ## Test specificSNV
    expect_null(paramTest@specificSNV)

    ## Test genoType
    expect_equal(paramTest@genoType, "geno.ref")
  
    ## Test phaseType 
    expect_equal(paramTest@phaseType, "phase.ref")
    
    ## Test phase
    expect_false(paramTest@phase)
    
    ## Test PCAmissingRate 
    expect_equal(paramTest@PCAmissingRate, 0.025)
    
    ## Test PCAalgorithm
    expect_equal(paramTest@PCAalgorithm, "exact")

    ## Test eigenCount
    expect_equal(paramTest@eigenCount, 32L)
    
    ## Test eigenCountSyn 
    expect_equal(paramTest@eigenCountSyn, 15L)
    
    ## Test kList
    expect_equal(paramTest@kList, seq(2L, 15L, 1L))

    ## Test pcaList
    expect_equal(paramTest@pcaList, seq(2L, 15L, 1L))

    ## Test fieldPopInRef
    expect_equal(paramTest@fieldPopInRef, "superPop")

    ## Test fieldPopInfAnc
    expect_equal(paramTest@fieldPopInfAnc, "superPop")

    ## Test fieldPopInfAnc
    expect_equal(paramTest@fieldSubPop, "pop.group")

    ## Test verbose
    expect_false(paramTest@verbose)
})