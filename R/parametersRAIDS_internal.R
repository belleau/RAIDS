#' @title Validate that the parametersRAIDS object is valid
#'
#' @description This function validates all the parameters present in the 
#' \code{parametersRAIDS} object.
#'
#' @param parameters a \code{parametersRAIDS} an object with all the 
#' parameters used by the RAIDS workflow.
#'
#' @return \code{TRUE} when all the parameters in the object are valid; 
#' otherwise \code{FALSE}
#'
#' @examples
#' ##TODO
#'
#'
#' @author Pascal Belleau and Astrid Deschênes
#' @encoding UTF-8
#' @keywords internal
validateParamRAIDS <- function(parameters) {
    pRAIDSNames <- c("reference", # HGDP1kg
                                "studyDF",  ## DONE AD
                                "studyDFSyn",
                                "pedStudy", ## DONE AD
                                "studyType",
                                "genoSource", ## DONE AD
                                "blockTypeId",
                                "genome",
                                "chrInfo", ## DONE AD
                                "profileFile",
                                "profileFileGeno",
                                "pathProfileGDS",
                                "fileReferenceGDS", ## DONE AD
                                "fileReferenceAnnotGDS", ## DONE AD
                                "inferenceType",
                                "sampleRef",
                                "batch",
                                "prefix",
                                "nbSim",
                                "offset",  ## DONE AD
                                "minCov",  ## DONE AD
                                "minProb",  ## DONE AD
                                "seqError",
                                "seqErrorSyn",
                                "pRecomb",
                                "np",
                                "listPos",
                                "syntheticRefDF", ##DONE AD
                                "pruningMethod",  ## DONE AD
                                "slideWindowMaxBP",
                                "thresholdLD",
                                "specificSNV", ## DONE AD
                                "genoType",
                                "phaseType",
                                "phase", ## DONE AS
                                "PCAmissingRate", ## DONE ADD
                                "PCAalgorithm", ## DONE AD
                                "eigenCount", ## DONE AD
                                "eigenCountSyn",
                                "kList", ## DONE AD
                                "pcaList", ## DONE AD
                                "fieldPopInRef", ## DONE AD
                                "fieldPopInfAnc", ## DONE AD
                                "fieldSubPop", ## DONE AD
                                "verbose") ## DONE AD
    
    ## The parameters object must be of class parametersRAIDS
    if (!inherits(parameters, "parametersRAIDS")) {
        stop("The parameters must be of class \'parametersRAIDS\'.")
    }
    
    ## All entries must be present in the parameters object
    if (!all(names(parameters) %in% pRAIDSNames)) { 
        stop("All parameters must be present: ", 
                paste0(pRAIDSNames, collapse=", "))
    }
    
    ## Validated parameters: 
    validateParamRAIDS_subpart01(parameters=parameters)
  
    validateParamRAIDS_subpart02(parameters=parameters)

    validateParamRAIDS_subpart03(parameters=parameters)
    
    validateParamRAIDS_subpart04(parameters=parameters)

    invisible(TRUE)
}

#' @title Validate the first subsection of the parametersRAIDS object
#'
#' @description This function validates a subsection of the parameters 
#' present in the \code{parametersRAIDS} object. The validated parameters are:
#' "studyDF", "pedStudy", "genoSource", "chrInfo", "syntheticRefDF", 
#' "pruningMethod", "fileReferenceGDS", and "fileReferenceAnnotGDS".
#'
#' @param parameters a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return \code{TRUE} when all the parameters tested in the object are valid; 
#' otherwise \code{FALSE}
#'
#' @examples
#' 
#' ## Path where demo files are located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#' path1KG <- file.path(dataDir, "tests")
#'
#' ## Path to existing demo Population Reference GDS file
#' fileReferenceGDS <- file.path(path1KG, "ex1_good_small_1KG.gds")
#' ## Path to existing demo Population Reference GDS annotation file
#' fileAnnotGDS <- file.path(path1KG, "ex1_good_small_1KG_Annot.gds")
#' 
#' ## PED Study
#' ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
#'             Case.ID=c("TCGA-H01", "TCGA-H02"),
#'             Sample.Type=c("DNA", "DNA"),
#'             Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))
#' 
#' ## The data frame containing the information about the study
#' ## The 3 mandatory columns: "study.id", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyInfo <- data.frame(study.id="Pancreatic.WES",
#'                 study.desc="Pancreatic study",
#'                 study.platform="WES",
#'                 stringsAsFactors=FALSE)
#' 
#' ## Profiles used for synthetic data set
#' syntheticRefDF <- data.frame(sample.id=c("HG00150", "HG00138", "HG00330",
#'         "HG00275"), pop.group=c("GBR", "GBR","FIN", "FIN"),
#'         superPop=c("EUR", "EUR", "EUR", "EUR"), stringsAsFactors=FALSE)
#' 
#' ## Create an object of class 'parametersRAIDS' with most parameters filled
#' ## default values
#' parameterAll <- paramRAIDS(studyDF=studyInfo, pedStudy=ped, 
#'     genoSource='generic', syntheticRefDF=syntheticRefDF, 
#'     fileReferenceGDS=fileReferenceGDS, fileReferenceAnnotGDS=fileAnnotGDS)
#'
#' ## Return TRUE when all tested parameters are valid
#' RAIDS:::validateParamRAIDS_subpart01(parameters=parameterAll)
#' 
#' 
#' @author Pascal Belleau and Astrid Deschênes
#' @encoding UTF-8
#' @keywords internal
validateParamRAIDS_subpart01 <- function(parameters) {

    ## The study data frame must have the mandatory columns
    validateStudyDataFrameParameter(studyDF=parameters$studyDF)
  
    ## The PED study must have the mandatory columns
    validatePEDStudyParameter(pedStudy=parameters$pedStudy)

    ## The genoSource must be a character string
    if(!(is.character(parameters$genoSource))) {
        stop("The \'genoSource\' parameter must be a character string.")
    }

    ## The chrInfo must be a vector of integer
    validatePositiveIntegerVector(parameters$chrInfo, "chrInfo")
  
    ## The syntheticRefDF must have the mandatory columns
    validateDataRefSynParameter(syntheticRefDF=parameters$syntheticRefDF)
  
    ## The fileReferenceGDS must be a character string and the file must exist
    if (!(is.character(parameters$fileReferenceGDS) && 
              (file.exists(parameters$fileReferenceGDS)))) {
        stop("The \'fileReferenceGDS\' must be a character string ",
                "representing the Reference GDS file. The file must exist.")
    }
  
    ## The fileReferenceAnnotGDS must be a character string and
    ## the file must exist
    if (!(is.character(parameters$fileReferenceAnnotGDS) &&
                (file.exists(parameters$fileReferenceAnnotGDS)))) {
        stop("The \'fileReferenceAnnotGDS\' must be a character string ",
                "representing the Reference Annotation GDS file. ",
                "The file must exist.")
    }

    if (!is.character(parameters$pruningMethod) || 
        !parameters$pruningMethod %in% c("corr", "r", "dprime", "composite")) {
            stop("The \'pruningMethod\' must be one of those values: ",
                    "\'corr\', \'r\', \'dprime\', or \'composite\'.")
    }

    invisible(TRUE)
}

#' @title Validate the second subsection of the parametersRAIDS object
#'
#' @description This function validates a subsection of the parameters 
#' present in the \code{parametersRAIDS} object. The validated parameters are:
#' "specificSNV", "phase", "PCAalgorithm", "PCAmissingRate" TODO.
#'
#' @param parameters a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return \code{TRUE} when all the parameters tested in the object are valid; 
#' otherwise \code{FALSE}
#'
#' @examples
#'
#' ## A data frame with the retained SNVs
#' specificSNV <-  data.frame("snp.chromosome"=c("chr1", "chr1"),
#'                             "snp.position"=c(12012, 14222))
#' 
#' ## Create an object of class 'parametersRAIDS' with most parameters filled
#' ## default values
#' parameterAll <- paramRAIDS(genoSource='generic', specificSNV=specificSNV,
#'         verbose=TRUE)
#'
#' ## Return TRUE when the tested parameters are valid
#' RAIDS:::validateParamRAIDS_subpart02(parameters=parameterAll)
#'
#' @author Pascal Belleau and Astrid Deschênes
#' @encoding UTF-8
#' @keywords internal
validateParamRAIDS_subpart02 <- function(parameters) {
   
    ## The specificSNV parameter should be a data frame with 2 columns minimum
    if (!(is.data.frame(parameters$specificSNV)) || 
                (!all(c("snp.chromosome", "snp.position") %in% 
                                    colnames(parameters$specificSNV)))) {
        stop("The \'specificSNV\' parameter must be a data frame ", 
            "with at least the two columns \'snp.chromosome\' and", 
            " \'snp.position\'.")
    }

    ## The phase parameter should be a logical
    validateLogical(parameters$phase, "verbose")
    
    ## The PCAmissingRate parameter should be a positive numeric or NaN
    if (!(is.numeric(parameters$PCAmissingRate) && 
                (parameters$PCAmissingRate > 0.0)) && 
                        !(is.nan(parameters$PCAmissingRate))) {
        stop("The \'PCAmissingRate\' parameter must be a positive numeric ", 
            "representing the SNPs missing rate cut-off or \'NaN\'.")
    }

    ## The PCAalgorithm parameter should be a character string
    if (!is.character(parameters$PCAalgorithm)) {
        stop("The \'PCAalgorithm\' parameter must be a character string ", 
            "representing the algorithm used to run the PCA analysis with ", 
            "the snpgdsPCA() function.")
    }

    invisible(TRUE)
}


#' @title Validate the third subsection of the parametersRAIDS object
#'
#' @description This function validates a subsection of the parameters 
#' present in the \code{parametersRAIDS} object. The validated parameters are:
#' "eigenCount", "kList", "pcaList", "fieldPopInRef", "fieldPopInfAnc", 
#' "fieldSubPop", and "verbose".
#'
#' @param parameters a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return \code{TRUE} when all the parameters tested in the object are valid; 
#' otherwise \code{FALSE}
#'
#' @examples
#' 
#' ## Create an object of class 'parametersRAIDS' with most parameters filled
#' ## default values
#' parameterAll <- paramRAIDS(genoSource='generic', verbose=TRUE)
#'
#' ## Return TRUE when the tested parameters are valid
#' RAIDS:::validateParamRAIDS_subpart03(parameters=parameterAll)
#'
#' @author Pascal Belleau and Astrid Deschênes
#' @encoding UTF-8
#' @importFrom S4Vectors isSingleNumber
#' @keywords internal
validateParamRAIDS_subpart03 <- function(parameters) {
   
    ## The eigenCount must be a single integer
    if(!(isSingleNumber(parameters$eigenCount))) {
        stop("The \'eigenCount\' parameter must be a single integer.")
    }

    ## The kList must be a vector of positive integers
    validatePositiveIntegerVector(parameters$kList, "kList")

    ## The pcaList must be a vector of positive integers
    validatePositiveIntegerVector(parameters$pcaList, "pcaList")

    ## The fieldPopInRef parameter should be a character string
    if (!is.character(parameters$fieldPopInRef)) {
        stop("The \'fieldPopInRef\' parameter must be a character string ", 
            "representing an existing column that contains the known ancestry", 
            " for the reference profiles in the Population Reference ", 
            "GDS file.")
    }

    ## The fieldPopInfAnc parameter should be a character string
    if (!is.character(parameters$fieldPopInfAnc)) {
        stop("The \'fieldPopInfAnc\' parameter must be a character string ", 
            "representing an existing column in the data frame containing the ", 
            "inferred super-population ancestry results.")
    }

    ## The fieldSubPop parameter should be a character string
    if (!is.character(parameters$fieldSubPop)) {
        stop("The \'fieldSubPop\' parameter must be a character string ", 
            "representing an existing column in the Population Reference ", 
            "GDS file.")
    }

    ## The verbose parameter should be a logical
    validateLogical(parameters$verbose, "verbose")

    invisible(TRUE)
}


#' @title Validate the fourth subsection of the parametersRAIDS object
#'
#' @description This function validates a subsection of the parameters 
#' present in the \code{parametersRAIDS} object. The validated parameters are:
#' "offset", "minCov", "minProb", ...
#'
#' @param parameters a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return \code{TRUE} when all the parameters tested in the object are valid; 
#' otherwise \code{FALSE}
#'
#' @examples
#' 
#' ## Create an object of class 'parametersRAIDS' with most parameters filled
#' ## default values
#' parameterAll <- paramRAIDS(minCov=20L, minProb=0.20, offset=-1L)
#'
#' ## Return TRUE when the tested parameters are valid
#' RAIDS:::validateParamRAIDS_subpart04(parameters=parameterAll)
#'
#' @author Pascal Belleau and Astrid Deschênes
#' @encoding UTF-8
#' @importFrom S4Vectors isSingleNumber
#' @keywords internal
validateParamRAIDS_subpart04 <- function(parameters) {

    ## Validate that offset is a single integer
    if (! isSingleNumber(parameters$offset)) {
        stop("The \'offset\' must be a single integer.")
    }

    ## The minCov parameter must be a single positive integer
    if (!(isSingleNumber(parameters$minCov) && (parameters$minCov >= 0.0))) {
        stop("The \'minCov\' must be a single numeric positive value")
    }

    ## The minProb and eProb must be single positive numeric between 0 and 1
    validateSingleRatio(parameters$minProb, "minProb")

    invisible(TRUE)
}