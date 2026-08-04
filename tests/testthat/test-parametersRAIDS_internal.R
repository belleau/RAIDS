### Unit tests for parametersRAIDS_internal.R functions

library(RAIDS)


context("validateParamRAIDS_subpart01() results")

test_that("validateParamRAIDS_subpart01() must return expected results when all inputs are valid", {

    dataDir <- test_path("fixtures")
    fileReferenceGDS <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")
    fileAnnotGDS <- file.path(dataDir,  "ex1_good_small_1KG_Annot_GDS.gds")

    ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
                Case.ID=c("TCGA-H01", "TCGA-H02"),
                Sample.Type=c("DNA", "DNA"),
                Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))
    
    studyInfo <- data.frame(study.id="Pancreatic.WES",
                    study.desc="Pancreatic study",
                    study.platform="WES", stringsAsFactors=FALSE)

    syntheticRefDF <- data.frame(sample.id=c("HG00150", "HG00138", "HG00330",
            "HG00275"), pop.group=c("GBR", "GBR","FIN", "FIN"),
            superPop=c("EUR", "EUR", "EUR", "EUR"), stringsAsFactors=FALSE)
    
    testParam <- paramRAIDS(studyDF=studyInfo, pedStudy=ped, 
        genoSource='generic', syntheticRefDF=syntheticRefDF, pruningMethod="r",
        fileReferenceGDS=fileReferenceGDS, fileReferenceAnnotGDS=fileAnnotGDS)

    expect_true(RAIDS:::validateParamRAIDS_subpart01(parameters=testParam))
})
