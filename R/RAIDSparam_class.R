#' An S4 class to represent the RAIDS parameters
#'
#' @slot studyDF TODO
#' @slot studyDFSyn TODO
#' @slot studyType TODO
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
#' @slot offset TODO
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
#' @exportClass RAIDSparam
setClass("RAIDSparam",
  slots = c(
    studyDF = "character",
    studyDFSyn = "numeric",
    studyType="character",
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
    studyDF = NULL,
    studyDFSyn = NULL,
    studyType = NA_character_,
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