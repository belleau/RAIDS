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

    ## Test fileReferenceAnnotGDS
    expect_null(paramTest@fileReferenceAnnotGDS)

    ## Test inferenceType
    expect_identical(paramTest@inferenceType, "PCAknn")

    ## Test sampleRef
    expect_null(paramTest@sampleRef)

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

test_that("create a RAIDSparam class with integer for studyDF parameter should generate an error", {

    expect_error(new("RAIDSparam", studyDF=33), 
        "got class \"numeric\", should be or extend class \"data.frame\"")
})

test_that("create a RAIDSparam class with data.frame with missing third column for studyDF parameter should generate an error", {

    message <- paste0("'studyDF' slot must be a data.frame with those 3 ", 
        "columns: \"study.id\", \"study.desc\", \"study.platform\".")
    
    expect_error(new("RAIDSparam", studyDF=data.frame("study.id"=c(1,2), 
        "study.desc"=c(1,2))), message)
})

test_that("create a RAIDSparam class with data.frame with missing first column for studyDF parameter should generate an error", {

    message <- paste0("'studyDF' slot must be a data.frame with those 3 ", 
        "columns: \"study.id\", \"study.desc\", \"study.platform\".")
    
    expect_error(new("RAIDSparam", studyDF=data.frame("study.desc"=c(1,2), 
        "study.platform"=c(1,2))), message)
})

test_that("create a RAIDSparam class with data.frame with missing second column for studyDF parameter should generate an error", {

    message <- paste0("'studyDF' slot must be a data.frame with those 3 ", 
        "columns: \"study.id\", \"study.desc\", \"study.platform\".")
    
    expect_error(new("RAIDSparam", studyDF=data.frame("study.id"=c(1,2), 
        "study.platform"=c(1,2))), message)
})

test_that("create a RAIDSparam class with character for studyDFSyn parameter should generate an error", {

    message <- paste0("got class \"character\", should be or ", 
        "extend class \"data.frame\"")
    
    expect_error(new("RAIDSparam", studyDFSyn="CANADA"), message)
})

test_that("create a RAIDSparam class with data.frame with missing first column for studyDFSyn parameter should generate an error", {

    message <- paste0("'studyDFSyn' slot must be a data.frame with those 3 ", 
        "columns: \"study.id\", \"study.desc\", \"study.platform\".")
    
    expect_error(new("RAIDSparam", studyDFSyn=data.frame("study.desc"=c(1,2), 
        "study.platform"=c(1,2))), message)
})

test_that("create a RAIDSparam class with data.frame with missing second column for studyDFSyn parameter should generate an error", {

    message <- paste0("'studyDFSyn' slot must be a data.frame with those 3 ", 
        "columns: \"study.id\", \"study.desc\", \"study.platform\".")
    
    expect_error(new("RAIDSparam", studyDFSyn=data.frame("study.id"=c(1,2), 
        "study.platform"=c(1,2))), message)
})

test_that("create a RAIDSparam class with data.frame with missing third column for studyDFSyn parameter should generate an error", {

    message <- paste0("'studyDFSyn' slot must be a data.frame with those 3 ", 
        "columns: \"study.id\", \"study.desc\", \"study.platform\".")
    
    expect_error(new("RAIDSparam", studyDFSyn=data.frame("study.id"=c(1,2), 
        "study.desc"=c(1,2))), message)
})

test_that("create a RAIDSparam class with integer for pedStudy parameter should generate an error", {

    message <- paste0("object: invalid object for slot \"pedStudy\" in class", 
        " \"RAIDSparam\": got class \"numeric\", should be or extend ", 
        "class \"data.frame\"")
    
    expect_error(new("RAIDSparam", pedStudy=44), message)
})

test_that("create a RAIDSparam class with data.frame with missing first column for pedStudy parameter should generate an error", {

    message <- paste0("'pedStudy' slot must be a data.frame with those 5 ", 
        "columns: \"Name.ID\", \"Case.ID\", \"Sample.Type\", \"Diagnosis\", ", 
        "and \"Source\".")
    
    expect_error(new("RAIDSparam", pedStudy=data.frame("Name"=c("1", "2"), 
        "Case.ID"=c(1,2), "Sample.Type"=c("cancer", "cancer"), 
        "Diagnosis"=c("cancer", "cancer"), "Source"=c("CSHL", "MIT"), 
        row.names=c("1", "2"))), message)
})

test_that("create a RAIDSparam class with data.frame with missing second column for pedStudy parameter should generate an error", {

    message <- paste0("'pedStudy' slot must be a data.frame with those 5 ", 
        "columns: \"Name.ID\", \"Case.ID\", \"Sample.Type\", \"Diagnosis\", ", 
        "and \"Source\".")
    
    expect_error(new("RAIDSparam", pedStudy=data.frame("Name.ID"=c("1", "2"), 
        "Case"=c(1,2), "Sample.Type"=c("cancer", "cancer"), 
        "Diagnosis"=c("cancer", "cancer"), "Source"=c("CSHL", "MIT"), 
        row.names=c("1", "2"))), message)
})

test_that("create a RAIDSparam class with data.frame with missing third column for pedStudy parameter should generate an error", {

    message <- paste0("'pedStudy' slot must be a data.frame with those 5 ", 
        "columns: \"Name.ID\", \"Case.ID\", \"Sample.Type\", \"Diagnosis\", ", 
        "and \"Source\".")
    
    expect_error(new("RAIDSparam", pedStudy=data.frame("Name.ID"=c("1", "2"), 
        "Case.ID"=c(1,2), "Type"=c("cancer", "cancer"), 
        "Diagnosis"=c("cancer", "cancer"), "Source"=c("CSHL", "MIT"), 
        row.names=c("1", "2"))), message)
})

test_that("create a RAIDSparam class with data.frame with missing fourth column for pedStudy parameter should generate an error", {

    message <- paste0("'pedStudy' slot must be a data.frame with those 5 ", 
        "columns: \"Name.ID\", \"Case.ID\", \"Sample.Type\", \"Diagnosis\", ", 
        "and \"Source\".")
    
    ## New RAIDSparam with wrong pedStudy parameter
    expect_error(new("RAIDSparam", pedStudy=data.frame("Name.ID"=c("1", "2"), 
        "Case.ID"=c(1,2), "Sample.Type"=c("cancer", "cancer"), 
        "Diagnos"=c("cancer", "cancer"), "Source"=c("CSHL", "MIT"), 
        row.names=c("1", "2"))), message)
})

test_that("create a RAIDSparam class with data.frame with missing fifth column for pedStudy parameter should generate an error", {

    message <- paste0("'pedStudy' slot must be a data.frame with those 5 ", 
        "columns: \"Name.ID\", \"Case.ID\", \"Sample.Type\", \"Diagnosis\", ", 
        "and \"Source\".")
    
    expect_error(new("RAIDSparam", pedStudy=data.frame("Name.ID"=c("1", "2"), 
        "Case.ID"=c(1,2), "Sample.Type"=c("cancer", "cancer"), 
        "Diagnosis"=c("cancer", "cancer"), "source"=c("CSHL", "MIT"), 
        row.names=c("1", "2"))), message)
})

test_that("create a RAIDSparam class with data.frame with wrong row names for pedStudy parameter should generate an error", {

    message <- paste0("'pedStudy' slot must be a data.frame with those 5 ", 
        "columns: \"Name.ID\", \"Case.ID\", \"Sample.Type\", \"Diagnosis\", ", 
        "and \"Source\". All row names should correspond to the ", 
        "Name.ID values.")
    
    ## New RAIDSparam with wrong pedStudy parameter
    expect_error(new("RAIDSparam", pedStudy=data.frame("Name.ID"=c("1", "2"), 
        "Case.ID"=c(1,2), "Sample.Type"=c("cancer", "cancer"), 
        "Diagnosis"=c("cancer", "cancer"), "Source"=c("CSHL", "MIT"), 
        row.names=c("11", "2"))), message)
})

test_that("create a RAIDSparam class with wrong name for studyType parameter should generate an error", {

    message <- paste0("'studyType' slot must have one character string ", 
        "within those 2 choices: \"LD\" and \"GeneAware\".")
    
    expect_error(new("RAIDSparam", studyType="CANADA"), message)
})

test_that("create a RAIDSparam class with vector of names for studyType parameter should generate an error", {

    message <- paste0("'studyType' slot must have one character string ", 
        "within those 2 choices: \"LD\" and \"GeneAware\".")
    
    ## New RAIDSparam with wrong studyType parameter
    expect_error(new("RAIDSparam", studyType=c("LD", "GeneAware")), message)
})

test_that("create a RAIDSparam class with wrong name for genoSource parameter should generate an error", {

    message <- paste0("'genoSource' slot must have one character ", 
                "string or NULL. The valid options are: \"VCF\", \"bam\", ", 
                "\"generic\", or \"snp-pileup\"")
    
    expect_error(new("RAIDSparam", genoSource="CANADA"), message)
})

test_that("create a RAIDSparam class with vector of strings for genoSource parameter should generate an error", {

    message <- paste0("'genoSource' slot must have one character ", 
                "string or NULL. The valid options are: \"VCF\", \"bam\", ", 
                "\"generic\", or \"snp-pileup\"")
    
    expect_error(new("RAIDSparam", genoSource=c("generic", "bam")), message)
})

test_that("create a RAIDSparam class with vector of strings for blockTypeId parameter should generate an error", {

    message <- paste0("'blockTypeId' slot must have one character ", 
            "string.")
    
    expect_error(new("RAIDSparam", blockTypeId=c("22", "44")), message)
})

test_that("create a RAIDSparam class with vector of strings for reference parameter should generate an error", {

    message <- paste0("'reference' slot must have one character string", 
                " within those 2 choices: \"1KGv1.0\", \"1k_hgdpV0.1\"")
    
    ## New RAIDSparam with wrong studyType parameter
    expect_error(new("RAIDSparam", reference=c("1KGv1.0", "1KGv1.0")), message)
})

test_that("create a RAIDSparam class with wrong string for reference parameter should generate an error", {

    message <- paste0("'reference' slot must have one character string", 
                " within those 2 choices: \"1KGv1.0\", \"1k_hgdpV0.1\"")
    
    expect_error(new("RAIDSparam", reference="1KGv13.0"), message)
})

test_that("create a RAIDSparam class with wrong string for genome parameter should generate an error", {

    message <- paste0("'genome' slot must be the character string \"HG38\"")
    
    expect_error(new("RAIDSparam", genome="1KG"), message)
})

test_that("create a RAIDSparam class with list with missing first entry for paramAncestry parameter should generate an error", {

    message <- paste0("'paramAncestry' slot must be a list with those ", 
                        "three entries: \"ScanBamParam\", \"PileupParam\", ", 
                        "and \"yieldSize\".")
    
    expect_error(new("RAIDSparam", paramAncestry=list(PileupParam=NULL, 
            yieldSize=10000000)), message)
})

test_that("create a RAIDSparam class with list with missing second entry for paramAncestry parameter should generate an error", {

    message <- paste0("'paramAncestry' slot must be a list with those ", 
                        "three entries: \"ScanBamParam\", \"PileupParam\", ", 
                        "and \"yieldSize\".")
    
    ## New RAIDSparam with wrong studyType parameter
    expect_error(new("RAIDSparam", paramAncestry=list(ScanBamParam=NULL, 
            yieldSize=10000000)), message)
})

test_that("create a RAIDSparam class with list with missing second entry for paramAncestry parameter should generate an error", {

    message <- paste0("'paramAncestry' slot must be a list with those ", 
                        "three entries: \"ScanBamParam\", \"PileupParam\", ", 
                        "and \"yieldSize\".")
    
    expect_error(new("RAIDSparam", paramAncestry=list(ScanBamParam=NULL, 
            yieldSize=10000000)), message)
})

test_that("create a RAIDSparam class with list with missing third entry for paramAncestry parameter should generate an error", {

    message <- paste0("'paramAncestry' slot must be a list with those ", 
                        "three entries: \"ScanBamParam\", \"PileupParam\", ", 
                        "and \"yieldSize\".")
    
    expect_error(new("RAIDSparam", paramAncestry=list(ScanBamParam=NULL, 
            PileupParam=NULL)), message)
})

test_that("create a RAIDSparam class with multiple strings for profileFile parameter should generate an error", {

    message <- paste0("'profileFile' slot must have one character string.")
    
    ## New RAIDSparam with wrong studyType parameter
    expect_error(new("RAIDSparam", profileFile=c("aa", "bb")), message)
})

test_that("create a RAIDSparam class with wrong extension for profileFile parameter when expecting bam should generate an error", {

    message <- paste0("'profileFile' slot must have one character string ", 
            "representing a file with extension '.bam' according to ",
            "'genoSource' slot.")
    
    expect_error(new("RAIDSparam", genoSource="bam", profileFile="test.csv"), message)
})

test_that("create a RAIDSparam class with wrong extension for profileFile parameter when expecting VCF should generate an error", {

    message <- paste0("'profileFile' slot must have one character string ", 
            "representing a file with extension '.vcf.gz' according to ",
            "'genoSource' slot.")
    
    expect_error(new("RAIDSparam", genoSource="VCF", profileFile="test.csv"), message)
})

test_that("create a RAIDSparam class with wrong extension for profileFile parameter when expecting generic should generate an error", {

    message <- paste0("'profileFile' slot must have one character string ", 
            "representing a file with extension '.txt.gz' according to ",
            "'genoSource' slot.")
    
    expect_error(new("RAIDSparam", genoSource="generic", profileFile="test.csv"), message)
})

test_that("create a RAIDSparam class with not existing dir for pathProfileGDS parameter should generate an error", {

    message <- paste0("'pathProfileGDS' slot must have one character string ", 
        "representing an existing directory.")
    
    expect_error(new("RAIDSparam", pathProfileGDS="./generic"), message)
})

test_that("create a RAIDSparam class with multiple strings for pathProfileGDS parameter should generate an error", {

    message <- paste0("'pathProfileGDS' slot must have one character string ", 
        "representing an existing directory.")
    
    expect_error(new("RAIDSparam", pathProfileGDS=c("./generic", "./testCanada")), message)
})

test_that("create a RAIDSparam class with multiple strings for fileReferenceGDS parameter should generate an error", {

    message <- paste0("'fileReferenceGDS' slot must have one character string", 
        " representing an existing file.")
    
    expect_error(new("RAIDSparam", fileReferenceGDS=c("generic.gds", "test.gds")), message)
})

test_that("create a RAIDSparam class with not existing file for fileReferenceGDS parameter should generate an error", {

    message <- paste0("'fileReferenceGDS' slot must have one character string", 
        " representing an existing file.")
    
    expect_error(new("RAIDSparam", fileReferenceGDS="./generic.gds"), message)
})

test_that("create a RAIDSparam class with not existing file for fileReferenceAnnotGDS parameter should generate an error", {

    message <- paste0("'fileReferenceAnnotGDS' slot must have one character ", 
        "string representing an existing file.")
    
    expect_error(new("RAIDSparam", fileReferenceAnnotGDS="./generic.gds"), message)
})

test_that("create a RAIDSparam class with multiple strings for fileReferenceAnnotGDS parameter should generate an error", {

    message <- paste0("'fileReferenceAnnotGDS' slot must have one character ", 
        "string representing an existing file.")
    
    expect_error(new("RAIDSparam", fileReferenceAnnotGDS=c("generic.gds", "test.gds")), message)
})

test_that("create a RAIDSparam class with multiple strings for inferenceType parameter should generate an error", {

    message <- paste0("'inferenceType' slot must be a single character ", 
        "string. The valid options are: 'PCAknn' and 'haploAdmixture'.")
    
    expect_error(new("RAIDSparam", inferenceType=c("PCAknn", "PCAknn")), message)
})

test_that("create a RAIDSparam class with wrong string for inferenceType parameter should generate an error", {

    message <- paste0("'inferenceType' slot must be a single character ", 
        "string. The valid options are: 'PCAknn' and 'haploAdmixture'.")
    
    expect_error(new("RAIDSparam", inferenceType="PCAknnT"), message)
})

test_that("create a RAIDSparam class with negative number for batch parameter should generate an error", {

    expect_error(new("RAIDSparam", batch=-1L), 
        "'batch' slot must have one positive integer.")
})

test_that("create a RAIDSparam class with multiple numbers for batch parameter should generate an error", {

    expect_error(new("RAIDSparam", batch=c(1L, 22L)), 
        "'batch' slot must have one positive integer.")
})

test_that("create a RAIDSparam class with multiple strings for prefix parameter should generate an error", {

    expect_error(new("RAIDSparam", prefix=c("a", "b")), 
        "'prefix' slot must have one character string.")
})

test_that("create a RAIDSparam class with multiple numbers for nbSim parameter should generate an error", {

    expect_error(new("RAIDSparam", nbSim=c(1L, 22L)), 
        "'nbSim' slot must have one positive integer.")
})

test_that("create a RAIDSparam class with multiple integers for offset parameter should generate an error", {

    expect_error(new("RAIDSparam", offset=c(0L, 1L)), 
        "'offset' slot must have one integer.")
})

test_that("create a RAIDSparam class with multiple integers for minCov parameter should generate an error", {

    expect_error(new("RAIDSparam", minCov=c(1L, 22L)), 
        "'minCov' slot must have one positive integer.")
})

test_that("create a RAIDSparam class with zero for minCov parameter should generate an error", {

    expect_error(new("RAIDSparam", minCov=0L), 
        "'minCov' slot must have one positive integer.")
})

test_that("create a RAIDSparam class with negative value for minProb parameter should generate an error", {

    expect_error(new("RAIDSparam", minProb=-0.001), 
        "'minProb' slot must have one positive numeric between 0 and 1.")
})

test_that("create a RAIDSparam class with value superior to 1 for minProb parameter should generate an error", {

    expect_error(new("RAIDSparam", minProb=1.001), 
        "'minProb' slot must have one positive numeric between 0 and 1.")
})

test_that("create a RAIDSparam class with multiple values for minProb parameter should generate an error", {

    expect_error(new("RAIDSparam", minProb=c(0.001, 0.1)), 
        "'minProb' slot must have one positive numeric between 0 and 1.")
})

test_that("create a RAIDSparam class with multiple values for seqError parameter should generate an error", {

    expect_error(new("RAIDSparam", seqError=c(0.001, 0.1)), 
        "'seqError' slot must have one positive numeric between 0 and 1.")
})

test_that("create a RAIDSparam class with value superior to 1 for seqError parameter should generate an error", {

    expect_error(new("RAIDSparam", seqError=1.001), 
        "'seqError' slot must have one positive numeric between 0 and 1.")
})

test_that("create a RAIDSparam class with negative value for seqError parameter should generate an error", {

    expect_error(new("RAIDSparam", seqError=-0.001), 
        "'seqError' slot must have one positive numeric between 0 and 1.")
})

test_that("create a RAIDSparam class with multiple values for seqErrorSyn parameter should generate an error", {

    expect_error(new("RAIDSparam", seqErrorSyn=c(0.001, 0.1)), 
        "'seqErrorSyn' slot must have one positive numeric between 0 and 1.")
})

test_that("create a RAIDSparam class with value superior to 1 for seqErrorSyn parameter should generate an error", {

    expect_error(new("RAIDSparam", seqErrorSyn=1.001), 
        "'seqErrorSyn' slot must have one positive numeric between 0 and 1.")
})

test_that("create a RAIDSparam class with negative value for seqErrorSyn parameter should generate an error", {

    expect_error(new("RAIDSparam", seqErrorSyn=-0.001), 
        "'seqErrorSyn' slot must have one positive numeric between 0 and 1.")
})

test_that("create a RAIDSparam class with negative value for pRecomb parameter should generate an error", {

    expect_error(new("RAIDSparam", pRecomb=-0.001), 
        "'pRecomb' slot must have one positive numeric between 0 and 1.")
})

test_that("create a RAIDSparam class with multiple values for pRecomb parameter should generate an error", {

    expect_error(new("RAIDSparam", pRecomb=c(0.001, 0.1)), 
        "'pRecomb' slot must have one positive numeric between 0 and 1.")
})

test_that("create a RAIDSparam class with negative value for pRecomb parameter should generate an error", {

    expect_error(new("RAIDSparam", pRecomb=-0.001), 
        "'pRecomb' slot must have one positive numeric between 0 and 1.")
})

test_that("create a RAIDSparam class with multiple values for np parameter should generate an error", {

    expect_error(new("RAIDSparam", np=c(2L, 3L)), 
        "'np' slot must have one positive integer.")
})

test_that("create a RAIDSparam class with value superior to 1 for np parameter should generate an error", {

    expect_error(new("RAIDSparam", np=-1L), 
        "'np' slot must have one positive integer.")
})

test_that("create a RAIDSparam class with integer for listPos parameter should generate an error", {

    expect_error(new("RAIDSparam", listPos=-1L), 
        "got class \"integer\", should be or extend class \"DataFrameOrNULL\"")
})

test_that("create a RAIDSparam class with snp.position column missing for listPos parameter should generate an error", {

    expect_error(new("RAIDSparam", listPos=data.frame("snp.chromosome"=c("1", "2"))), 
        paste0("'listPos' slot must be NULL or a data.frame with 2 columns", 
            " named \"snp.chromosome\" and \"snp.position\""))
})

test_that("create a RAIDSparam class with snp.chromosome column missing for listPos parameter should generate an error", {

    expect_error(new("RAIDSparam", listPos=data.frame("snp.position"=c("1", "2"))), 
        paste0("'listPos' slot must be NULL or a data.frame with 2 columns", 
            " named \"snp.chromosome\" and \"snp.position\""))
})

test_that("create a RAIDSparam class with sample.id column missing for syntheticRefDF parameter should generate an error", {

    expect_error(new("RAIDSparam", syntheticRefDF=data.frame(
        "pop.group"=c("AFR", "EUR"), "superPop"=c("YRI", "CEU"))), 
        paste0("'syntheticRefDF' slot must be NULL or a data.frame with 3 ", 
            "columns named \"sample.id\", \"pop.group\", and \"superPop\"."))
})

test_that("create a RAIDSparam class with pop.group column missing for syntheticRefDF parameter should generate an error", {

    expect_error(new("RAIDSparam", syntheticRefDF=data.frame(
        "sample.id"=c("S1", "S2"), "superPop"=c("AFR", "EUR"))), 
        paste0("'syntheticRefDF' slot must be NULL or a data.frame with 3 ", 
            "columns named \"sample.id\", \"pop.group\", and \"superPop\"."))
})

test_that("create a RAIDSparam class with superPop column missing for syntheticRefDF parameter should generate an error", {

    expect_error(new("RAIDSparam", syntheticRefDF=data.frame(
        "sample.id"=c("S1", "S2"), "pop.group"=c("YRI", "CEU"))), 
        paste0("'syntheticRefDF' slot must be NULL or a data.frame with 3 ", 
            "columns named \"sample.id\", \"pop.group\", and \"superPop\"."))
})

test_that("create a RAIDSparam class with array of strings for pruningMethod parameter should generate an error", {

    expect_error(new("RAIDSparam", pruningMethod=c("r", "corr")), 
        paste0("'pruningMethod' slot must have one character string within ", 
        "those 4 choices: \"corr\", \"r\", \"dprime\", \"composite\"."))
})

test_that("create a RAIDSparam class with wrong string for pruningMethod parameter should generate an error", {

    expect_error(new("RAIDSparam", pruningMethod="E"), 
        paste0("'pruningMethod' slot must have one character string within ", 
        "those 4 choices: \"corr\", \"r\", \"dprime\", \"composite\"."))
})

test_that("create a RAIDSparam class with negative integer for slideWindowMaxBP parameter should generate an error", {

    expect_error(new("RAIDSparam", slideWindowMaxBP=-1L), 
        "'slideWindowMaxBP' slot must have one positive integer value.")
})

test_that("create a RAIDSparam class with multiple integers for slideWindowMaxBP parameter should generate an error", {

    expect_error(new("RAIDSparam", slideWindowMaxBP=c(1L, 2L)), 
        "'slideWindowMaxBP' slot must have one positive integer value.")
})

test_that("create a RAIDSparam class with negative numeric for thresholdLD parameter should generate an error", {

    expect_error(new("RAIDSparam", thresholdLD=-1), 
        "'thresholdLD' slot must have one positive numeric value.")
})

test_that("create a RAIDSparam class with multiple numerics for thresholdLD parameter should generate an error", {

    expect_error(new("RAIDSparam", thresholdLD=c(1, 2)), 
        "'thresholdLD' slot must have one positive numeric value.")
})

test_that("create a RAIDSparam class with snp.position column missing for specificSNV parameter should generate an error", {

    expect_error(new("RAIDSparam", specificSNV=data.frame("snp.chromosome"=c(1, 2))), 
        paste0("'specificSNV' slot must be NULL or a ", 
                "data.frame with 2 columns named \"snp.chromosome\"", 
                " and \"snp.position\"."))
})

test_that("create a RAIDSparam class with snp.chromosome column missing for specificSNV parameter should generate an error", {

    expect_error(new("RAIDSparam", specificSNV=data.frame("snp.position"=c(1, 2))), 
        paste0("'specificSNV' slot must be NULL or a ", 
                "data.frame with 2 columns named \"snp.chromosome\"", 
                " and \"snp.position\"."))
})

## Validate the genoType parameter TODO


test_that("create a RAIDSparam class with multiple strings for phaseType parameter should generate an error", {

    expect_error(new("RAIDSparam", phaseType=c("phase.ref", "phase.ref")), 
        "'phaseType' slot must be one character string.")
})

test_that("create a RAIDSparam class with multiple logical for phase parameter should generate an error", {

    expect_error(new("RAIDSparam", phase=c(TRUE, FALSE)), 
        "'verbose' slot must have one logical value.")
})

test_that("create a RAIDSparam class with negative numeric for PCAmissingRate parameter should generate an error", {

    expect_error(new("RAIDSparam", PCAmissingRate=-0.011), 
        "'PCAmissingRate' slot must have one positive numeric value.")
})

test_that("create a RAIDSparam class with multiple numerics for PCAmissingRate parameter should generate an error", {

    expect_error(new("RAIDSparam", PCAmissingRate=c(0.01, 0.2)), 
        "'PCAmissingRate' slot must have one positive numeric value.")
})

test_that("create a RAIDSparam class with multiple strings for PCAalgorithm parameter should generate an error", {

    expect_error(new("RAIDSparam", PCAalgorithm=c("randomized", "exact")), 
        paste0("'PCAalgorithm' slot must have one character string. ", 
            "The valid options are: \"exact\" or \"randomized\"."))
})

test_that("create a RAIDSparam class with wrong string for PCAalgorithm parameter should generate an error", {

    expect_error(new("RAIDSparam", PCAalgorithm="test"), 
        paste0("'PCAalgorithm' slot must have one character string. ", 
            "The valid options are: \"exact\" or \"randomized\"."))
})

test_that("create a RAIDSparam class with multiple integers for eigenCount parameter should generate an error", {

    expect_error(new("RAIDSparam", eigenCount=c(1L, 2L)), 
        paste0("'eigenCount' slot must have one integer value."))
})

test_that("create a RAIDSparam class with multiple integers for eigenCountSyn parameter should generate an error", {

    expect_error(new("RAIDSparam", eigenCountSyn=c(1L, 2L)), 
        paste0("'eigenCountSyn' slot must have one integer value."))
})

test_that("create a RAIDSparam class with multiple integers for eigenCountSyn parameter should generate an error", {

    expect_error(new("RAIDSparam", eigenCountSyn=c(1L, 2L)), 
        paste0("'eigenCountSyn' slot must have one integer value."))
})

test_that("create a RAIDSparam class with negative integer for kList parameter should generate an error", {

    expect_error(new("RAIDSparam", kList=c(-1L)), 
        paste0("'kList' slot must have one or more positive integer values."))
})

test_that("create a RAIDSparam class with negative integer for pcaList parameter should generate an error", {

    expect_error(new("RAIDSparam", pcaList=c(-1L)), 
        paste0("'pcaList' slot must have one or more positive integer values."))
})

test_that("create a RAIDSparam class with multiple strings for fieldPopInRef parameter should generate an error", {

    expect_error(new("RAIDSparam", fieldPopInRef=c("a", "b")), 
        paste0("'fieldPopInRef' slot must have one character string."))
})

test_that("create a RAIDSparam class with multiple strings for fieldSubPop parameter should generate an error", {

    expect_error(new("RAIDSparam", fieldSubPop=c("a", "b")), 
        paste0("'fieldSubPop' slot must have one character string."))
})

test_that("create a RAIDSparam class with multiple logicals for verbose parameter should generate an error", {

    expect_error(new("RAIDSparam", verbose=c(FALSE, FALSE)), 
        paste0("'verbose' slot must have one logical value."))
})

test_that("create a RAIDSparam class with all studyDF setter and getter should return an object", {

    exp_studyDF <- data.frame(study.id="TEST1",
                                study.desc="Test",
                                study.platform="CSHL",
                                stringsAsFactors=FALSE)

    paramTest <- RAIDSparam()
    studyDF(paramTest) <- exp_studyDF

    expect_equal(studyDF(paramTest), exp_studyDF)
    expect_error(studyDF(paramTest) <- 33L)
    expect_equal(studyDF(paramTest), exp_studyDF)
})

test_that("create a RAIDSparam class with all studyDFSyn setter and getter should return an object", {

    exp_studyDFSyn <- data.frame(study.id=c("TEST12", "Test22"),
                                study.desc=c("Test", "Canada"),
                                study.platform=c("CSHL", "Alpha"),
                                stringsAsFactors=FALSE)

    paramTest <- RAIDSparam()
    studyDFSyn(paramTest) <- exp_studyDFSyn

    expect_equal(studyDFSyn(paramTest), exp_studyDFSyn)
    expect_error(studyDFSyn(paramTest) <- 323L)
    expect_equal(studyDFSyn(paramTest), exp_studyDFSyn)
})

test_that("create a RAIDSparam class with all pedStudy setter and getter should return an object", {

    exp_pedStudy <- data.frame(Name.ID=c("TEST12", "Test22"),
                                Case.ID=c("Case1", "Case12"),
                                Sample.Type=c("Cancer", "Cancer"),
                                Diagnosis=c("NA", "sarcoma"),
                                Source=c("NA", "CSHL"),
                                stringsAsFactors=FALSE, 
                                row.names=c("TEST12", "Test22"))

    paramTest <- RAIDSparam()
    pedStudy(paramTest) <- exp_pedStudy

    expect_equal(pedStudy(paramTest), exp_pedStudy)
    expect_error(pedStudy(paramTest) <- 323L)
    expect_equal(pedStudy(paramTest), exp_pedStudy)
})

test_that("create a RAIDSparam class with all studyType setter and getter should return an object", {

    exp_studyType <- "GeneAware"

    paramTest <- RAIDSparam()
    studyType(paramTest) <- exp_studyType

    expect_equal(studyType(paramTest), exp_studyType)
    expect_error(studyType(paramTest) <- 323L)
    expect_equal(studyType(paramTest), exp_studyType)
})

#############################################################################
### Tests RAIDSparam function
#############################################################################

context("RAIDSparam function results")

test_that("create a RAIDSparam function with all default parameters should return an object", {

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
    paramTest <- RAIDSparam()

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

    ## Test fileReferenceAnnotGDS
    expect_null(paramTest@fileReferenceAnnotGDS)

    ## Test inferenceType
    expect_identical(paramTest@inferenceType, "PCAknn")

    ## Test sampleRef
    expect_null(paramTest@sampleRef)

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



test_that("create a RAIDSparam function with all non-default parameters should return an object", {

    exp_studyDF <- data.frame(study.id="Test1",
                                study.desc="Test Desc",
                                study.platform="Test",
                                stringsAsFactors=FALSE)
    
    exp_studyDFSyn <- data.frame(study.id="Test 02",
                            study.desc="Synthetic data",
                            study.platform="Demo", stringsAsFactors=FALSE)

    exp_pedStudy <- data.frame(Name.ID=c("SampleE"), Case.ID=c("Profile11"),
                        Sample.Type=c("test"), Diagnosis="Cancer",
                        Source=c("CSHL"), stringsAsFactors=FALSE, 
                        row.names=c("SampleE"))

    exp_chrInfo <- c(2486422L, 24293529L, 198559L)
    names(exp_chrInfo) <- c(paste0("chr", 1:3))
    
    exp_listPos <- data.frame("snp.chromosome"=c("1", "1"), 
            "snp.position"=c(1,3))
    
    exp_syntheticRefDF <- data.frame(sample.id=c("Sample1", "Sample1"), 
        pop.group=c("EUR", "AFR"), superPop=c("EUR", "AFR"))
    
    exp_paramAncestry <- list(ScanBamParam=NULL, PileupParam=NULL,
                                yieldSize=10000000)
    
    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures",  "1KG_Test.gds")

    ## New RAIDSparam with all default values
    paramTest <- RAIDSparam(studyDF=data.frame(study.id="Test1",
        study.desc="Test Desc", study.platform="Test", stringsAsFactors=FALSE),
        studyDFSyn=data.frame(study.id="Test 02", study.desc="Synthetic data",
        study.platform="Demo", stringsAsFactors=FALSE), 
        pedStudy=data.frame(Name.ID=c("SampleE"), Case.ID=c("Profile11"),
        Sample.Type=c("test"), Diagnosis="Cancer", Source=c("CSHL"), 
        stringsAsFactors=FALSE, row.names=c("SampleE")), chrInfo=exp_chrInfo,
        studyType="GeneAware", genoSource="bam", blockTypeId="E", 
        reference="1k_hgdpV0.1", inferenceType="haploAdmixture", batch=12L,
        prefix="3", nbSim=11L, offset=0L, minProb=0.889, minCov=20L,
        pruningMethod="dprime", np=2L, pRecomb=0.21, seqError=0.12, 
        seqErrorSyn=0.02, pathProfileGDS=dataDir, fileReferenceGDS=fileGDS,
        fileReferenceAnnotGDS=fileGDS,
        listPos=data.frame("snp.chromosome"=c("1", "1"), "snp.position"=c(1,3)),
        syntheticRefDF=data.frame(sample.id=c("Sample1", "Sample1"), 
        pop.group=c("EUR", "AFR"), superPop=c("EUR", "AFR")), 
        slideWindowMaxBP=2000L, thresholdLD=0.2, genoType="geno.REF", 
        phaseType="phase.REF", phase=TRUE, PCAmissingRate=0.015, 
        PCAalgorithm="randomized", eigenCount=20L, eigenCountSyn=21L,
        kList=c(3L, 5L), pcaList=c(5L, 9L), fieldPopInRef="superPOP", 
        fieldPopInfAnc="supERPop", fieldSubPop="POP_GROUP", verbose=TRUE
    )

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
    expect_true(paramTest@studyType == "GeneAware")

    ## Test genoSource 
    expect_true(paramTest@genoSource == "bam")

    ## Test blockTypeId
    expect_true(paramTest@blockTypeId == "E")

    ## Test reference
    expect_true(paramTest@reference == "1k_hgdpV0.1")

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
    expect_equal(paramTest@pathProfileGDS, dataDir)
    
    ## Test fileReferenceGDS 
    expect_equal(paramTest@fileReferenceGDS, fileGDS)

    ## Test fileReferenceAnnotGDS
    expect_equal(paramTest@fileReferenceAnnotGDS, fileGDS)

    ## Test inferenceType
    expect_identical(paramTest@inferenceType, "haploAdmixture")

    ## Test sampleRef
    expect_null(paramTest@sampleRef)

    ## Test batch
    expect_true(paramTest@batch == 12L)

    ## Test prefix
    expect_equal(paramTest@prefix, "3")

    ## Test nbSim
    expect_true(paramTest@nbSim == 11L)

    ## Test offset
    expect_true(paramTest@offset == 0L)

    ## Test minCov
    expect_true(paramTest@minCov == 20L)

    ## Test minProb
    expect_true(paramTest@minProb == 0.889)

    ## Test seqError
    expect_true(paramTest@seqError == 0.12)

    ## Test seqErrorSyn
    expect_true(paramTest@seqErrorSyn == 0.02)

    ## Test pRecomb
    expect_true(paramTest@pRecomb == 0.21)

    ## Test np
    expect_true(paramTest@np == 2L)
  
    ## Test listPos
    expect_equal(paramTest@listPos, exp_listPos)

    ## Test syntheticRefDF
    expect_equal(paramTest@syntheticRefDF, exp_syntheticRefDF)
    
    ## Test pruningMethod
    expect_equal(paramTest@pruningMethod, "dprime")

    ## Test slideWindowMaxBP 
    expect_equal(paramTest@slideWindowMaxBP, 2000L)

    ## Test thresholdLD 
    expect_equal(paramTest@thresholdLD, 0.2)

    ## Test specificSNV
    expect_null(paramTest@specificSNV)

    ## Test genoType
    expect_equal(paramTest@genoType, "geno.REF")
  
    ## Test phaseType 
    expect_equal(paramTest@phaseType, "phase.REF")
    
    ## Test phase
    expect_true(paramTest@phase)
    
    ## Test PCAmissingRate 
    expect_equal(paramTest@PCAmissingRate, 0.015)
    
    ## Test PCAalgorithm
    expect_equal(paramTest@PCAalgorithm, "randomized")

    ## Test eigenCount
    expect_equal(paramTest@eigenCount, 20L)
    
    ## Test eigenCountSyn 
    expect_equal(paramTest@eigenCountSyn, 21L)
    
    ## Test kList
    expect_equal(paramTest@kList, c(3L, 5L))

    ## Test pcaList
    expect_equal(paramTest@pcaList, c(5L, 9L))

    ## Test fieldPopInRef
    expect_equal(paramTest@fieldPopInRef, "superPOP")

    ## Test fieldPopInfAnc
    expect_equal(paramTest@fieldPopInfAnc, "supERPop")

    ## Test fieldPopInfAnc
    expect_equal(paramTest@fieldSubPop, "POP_GROUP")

    ## Test verbose
    expect_true(paramTest@verbose)
})

