#' @title Class Union for data.frame or NULL
#' 
#' @description A virtual class that groups data.frame and \code{NULL} 
#' together.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @name DataFrameOrNULL-class
#' @rdname DataFrameOrNULL-class
#' @exportClass DataFrameOrNULL
setClassUnion("DataFrameOrNULL", members = c("data.frame", "NULL"))


#' @title Class Union for character or NULL
#' 
#' @description A virtual class that groups character and \code{NULL} 
#' together.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @name CharacterOrNULL-class
#' @rdname CharacterOrNULL-class
#' @exportClass CharacterOrNULL
setClassUnion("CharacterOrNULL", members = c("character", "NULL"))


#' An S4 class to represent the RAIDS parameters
#'
#' @slot studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor). Default: \code{NULL}.
#' 
#' @slot studyDFSyn a \code{data.frame} containing the information about the
#' synthetic data to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor). Default: \code{NULL}.
#' 
#' @slot studyType a \code{character} string representing the type of study.
#' The possible choices are: "LD" and "GeneAware". The type of study affects 
#' the way the estimation of the allelic fraction is done. 
#' Default: \code{"LD"}.
#' 
#' @slot genoSource TODO
#' @slot blockTypeId TODO
#' @slot reference a \code{character} string with two possible values:
#' '1KGv1.0', '1k_hgdpV0.1'. It specifies the type of inference. 
#' Default: \code{"1KGc1.0"}.
#' 
#' @slot chrInfo TODO
#' @slot paramAncestry TODO
#' @slot profileFile TODO
#' @slot profileFileGeno TODO
#' @slot pathProfileGDS TODO
#' @slot fileReferenceGDS TODO
#' @slot fileReferenceAnnotGDS TODO
#' @slot inferenceType TODO
#' @slot sampleRef TODO
#' 
#' @slot batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore. 
#' Default: \code{1L}.
#' 
#' @slot prefix TODO
#' @slot nbSim TODO
#' 
#' @slot offset a single \code{integer} that is added to the SNP position to
#' switch from 0-based to 1-based coordinate when needed (or reverse).
#' Default: \code{-1L}.
#' 
#' @slot minCov TODO
#' @slot minProb TODO
#' @slot seqError TODO
#' @slot seqErrorSyn TODO
#' @slot pRecomb TODO
#' @slot np TODO
#' @slot listPos TODO
#' @slot syntheticRefDF TODO
#' 
#' @slot pruningMethod a \code{character} string representing the method that 
#' will be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The 4 possible values
#' are: "corr", "r", "dprime", and "composite". Default: \code{"corr"}.
#' 
#' @slot slideWindowMaxBP TODO
#' @slot thresholdLD TODO
#' @slot specificSNV TODO
#' @slot genoType TODO
#' @slot phaseType TODO
#' @slot phase TODO
#' @slot PCAmissingRate TODO
#' @slot PCAalgorithm TODO
#' @slot eigenCount TODO
#' @slot eigenCountSyn TODO
#' @slot kList TODO
#' @slot pcaList TODO
#' @slot fieldPopInRef TODO
#' @slot fieldPopInfAnc TODO
#' @slot fieldSubPop TODO
#' @slot verbose TODO
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @import methods
#' @name RAIDSparams-class
#' @rdname RAIDSparam-class
#' @exportClass RAIDSparam
setClass("RAIDSparam",
  slots = c(
    studyDF = "DataFrameOrNULL",
    studyDFSyn = "DataFrameOrNULL",
    studyType="CharacterOrNULL",
    genoSource="character",
    blockTypeId="character",
    reference="character",
    chrInfo="integer",
    paramAncestry="list",
    profileFile="character",
    profileFileGeno="character",
    pathProfileGDS="character",
    fileReferenceGDS="character",
    fileReferenceAnnotGDS="character",
    inferenceType="character",
    sampleRef="character",
    batch="integer",
    prefix="character",
    nbSim="integer",
    offset="integer",
    minCov="integer",
    minProb="numeric",
    seqError="numeric",
    seqErrorSyn="numeric",
    pRecomb="numeric",
    np="integer",
    listPos="integer",
    syntheticRefDF="data.frame",
    pruningMethod="character",
    slideWindowMaxBP="integer",
    thresholdLD="numeric",
    specificSNV="data.frame",
    genoType="character",
    phaseType="character",
    phase="logical",
    PCAmissingRate="numeric",
    PCAalgorithm="character",
    eigenCount="integer",
    eigenCountSyn="integer",
    kList="integer",
    pcaList="integer",
    fieldPopInRef="character",
    fieldPopInfAnc="character",
    fieldSubPop="character",
    verbose="logical"
  ),
  prototype = list(
    studyDF=NULL,
    studyDFSyn=NULL,
    studyType=NULL,
    genoSource=NULL,
    blockTypeId="GeneS.Ensembl.Hsapiens.v86",
    reference="1KGc1.0",
    chrInfo=NULL,
    paramAncestry=NULL,
    profileFile=NULL,
    profileFileGeno=NULL,
    pathProfileGDS=NULL,
    fileReferenceGDS=NULL,
    fileReferenceAnnotGDS=NULL,
    inferenceType="PCAknn",
    sampleRef=NULL,
    batch=1L,
    prefix="1",
    nbSim=1L,
    offset=-1L,
    minCov=10L,
    minProb=0.999,
    seqError=0.001,
    seqErrorSyn=0.001,
    pRecomb=0.01,
    np=1L,
    listPos=NULL,
    syntheticRefDF=NULL,
    pruningMethod="corr",
    slideWindowMaxBP=500000L,
    thresholdLD=sqrt(0.1),
    specificSNV=NULL,
    genoType="geno.ref",
    phaseType="phase.ref",
    phase=FALSE,
    PCAmissingRate=0.025,
    PCAalgorithm="exact",
    eigenCount=32L,
    eigenCountSyn=15L,
    kList=seq(2L, 15L, 1L),
    pcaList=seq(2L, 15L, 1L),
    fieldPopInRef="superPop",
    fieldPopInfAnc="superPop",
    fieldSubPop="pop.group",
    verbose=FALSE
  )
)

setValidity("RAIDSparam",
    function(object)
    {
        ## Validate the studyDF parameter
        if (!is.null(object@studyDF) &&  
            !(is.data.frame(object@studyDF) && 
                all(c("study.id", "study.desc", "study.platform") %in% 
                    colnames(object@studyDF)))) {
            return(paste0("'studyDF' slot must be NULL or a data.frame with ", 
                "those 3 columns: \"study.id\", \"study.desc\", ", 
                "\"study.platform\"."))
        }
      
        ## Validate the studyDFSyn parameter
        if (!is.null(object@studyDFSyn) && 
            !(is.data.frame(object@studyDFSyn) && 
                all(c("study.id", "study.desc", "study.platform") %in% 
                    colnames(object@studyDFSyn)))) {
            return(paste0("'studyDFSyn' slot must be NULL or a data.frame ", 
                "with those 3 columns: \"study.id\", \"study.desc\", ", 
                "\"study.platform\"."))
        }
      
        ## Validate the studyType parameter
        if (!is.null(object@studyType) && !(is.character(object@studyType) && 
                length(object@studyType) != 1 || 
                !object@studyType %in% c("LD", "GeneAware"))) {
            return(paste0("'studyType' slot must be NULL or a one character", 
                " string within those 2 choices: \"LD\" and \"GeneAware\"."))
        }
      
        ## Validate the genoSource parameter TODO
      
        ## Validate the reference parameter
        if (length(object@reference) != 1 || 
            !object@reference %in% c("1KGv1.0", "1k_hgdpV0.1")) {
            return(paste0("'reference' slot must have one character string", 
                " within those 2 choices: \"1KGv1.0\", \"1k_hgdpV0.1\"."))
        }
      
        ## Validate the batch parameter
        if (length(object@batch != 1) || object@batch < 1) {
            return("'batch' slot must have one positive integer.")
        }
      
        ## Validate the offset parameter
        if (length(object@offset != 1)) {
            return("'offset' slot must have one integer.")
        }
      
        ## Validate the pruningMethod parameter
        if (length(object@pruningMethod) != 1 || 
              !object@pruningMethod %in% c("corr", "r", "dprime",
                  "composite")) {
          return(paste0("'pruningMethod' slot must have one character string ", 
              "within those 4 choices: \"corr\", \"r\", \"dprime\", ", 
              "\"composite\"."))
        }

        TRUE
    }
)