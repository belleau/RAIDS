#' @title Validate the parameters
#'
#' @description This function validates the parameters present in the 
#' \code{parametersRAIDS} object.
#'
#' @param parameters a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
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
                                "batch",
                                "prefix",
                                "nbSim",
                                "offset",
                                "minCov",
                                "minProb",
                                "seqError",
                                "seqErrorSyn",
                                "pRecomb",
                                "np",
                                "listPos",
                                "syntheticRefDF", ##DONE AD
                                "pruningMethod",
                                "slideWindowMaxBP",
                                "thresholdLD",
                                "specificSNV",
                                "genoType",
                                "phaseType",
                                "phase",
                                "PCAmissingRate",
                                "PCAalgorithm",
                                "eigenCount",
                                "eigenCountSyn",
                                "kList",
                                "pcaList",
                                "fieldPopInRef",
                                "fieldPopInfAnc",
                                "fieldSubPop",
                                "verbose") ## DONE AD
    
    ## All parameters must have an entry
    if (!all(names(parameters) %in% pRAIDSNames)){ 
        stop("All parameters must be present: ", 
                paste0(pRAIDSNames, collapse=", "))
    }
  
    validateParamRAIDS_subpart01(parameters=parameters)
  
    validateParamRAIDS_subpart02(parameters=parameters)
    
    invisible(TRUE)
}

#' @title Validate the first subsection fo the input parameters
#'
#' @description This function validates a subsection of the parameters 
#' present in the \code{parametersRAIDS} object.
#'
#' @param parameters a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return \code{TRUE} when all the parameters tested in the object are valid; 
#' otherwise \code{FALSE}
#'
#' @examples
#' ##TODO
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
  
    ## The fileReferenceGDS must be a character string and the file must exists
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

    invisible(TRUE)
}


#' @title Validate the second subsection fo the input parameters
#'
#' @description This function validates a subsection of the parameters 
#' present in the \code{parametersRAIDS} object.
#'
#' @param parameters a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return \code{TRUE} when all the parameters tested in the object are valid; 
#' otherwise \code{FALSE}
#'
#' @examples
#' ##TODO
#'
#'
#' @author Pascal Belleau and Astrid Deschênes
#' @encoding UTF-8
#' @keywords internal
validateParamRAIDS_subpart02 <- function(parameters) {

    validateLogical(parameters$verbose, "verbose")

    invisible(TRUE)
}

