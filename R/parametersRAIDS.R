#' @title parameters class for RAIDS
#' function
#'
#' @description parameters class for RAIDS
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor).
#'
#' @param pedStudy a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings (no factor). The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter. Only \code{filePedRDS} or \code{pedStudy}
#' can be defined.
#'
#' @param studyType a \code{character} string representing the type of study.
#' The possible choices are: "LD" and "GeneAware". The type of study affects the
#' way the estimation of the allelic fraction is done. Default: \code{"LD"}.
#'
#' @param blockTypeID a \code{character} string corresponding to the block
#' type used to extract the block identifiers. The block type must be
#' present in the GDS Reference Annotation file.
#'
#' @param reference a \code{character} string with two possible values:
#' '1KGv1.0', '1k_hgdpV0.1'. It specifies the type of inference
#'
#' @param genome a \code{character} string with one possible value:
#' 'HG38'. It specifies the genme uses.
#'
#' @param chrInfo a \code{vector} of positive \code{integer} values
#' representing the length of the chromosomes. See 'details' section.
#'
#' @param paramAncestry a \code{list} parameters ...
#'
#' @param profileFile a \code{character} string representing the path and the
#' file name of the genotype file or the bam if genoSource is snp-pileup the
#' fine extension must be .txt.gz, if VCF the extension must be .vcf.gz
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the GDS Profile files will be created.
#' Default: \code{NULL}.
#'
#' @param fileReferenceGDS  a \code{character} string representing the file
#' name of the Population Reference GDS file. The file must exist.
#'
#' @param fileReferenceAnnotGDS a \code{character} string representing the
#' file name of the Population Reference GDS Annotation file. The file
#' must exist.
#'
#' @param inferenceType a \code{character} string with two possible values:
#' 'PCAnn', 'haploAdmixture'. It specifies the genotype ancestry inference
#'
#' @param batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore.
#'
#' @param listPos a \code{data.frame} containing 2 columns. The first column,
#' called "snp.chromosome" contains the name of the chromosome where the
#' SNV is located. The second column, called "snp.position" contains the
#' position of the SNV on the chromosome.
#'
#' @param minCov a single positive \code{integer} representing the minimum
#' required coverage. Default: \code{10L}.
#' @param minProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability that the calculated genotype call is correct.
#' Default: \code{0.999}.
#'
#' @param seqError a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error. Default: \code{0.001}.
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed. Default: \code{FALSE}.
#'
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Load the demo PCA on the synthetic profiles projected on the
#' ## demo 1KG reference PCA
#' data(demoPCASyntheticProfiles)
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo GDS file is located in this package
#' dataDir <- system.file("extdata/demoKNNSynthetic", package="RAIDS")
#' fileProfileGDS <- file.path(dataDir, "ex1.gds")
#'
#' ## Open GDS files
#' gdsProfile <- openfn.gds(fileProfileGDS)
#'
#' ## The function returns 0L when all parameters are valid
#' RAIDS:::validateComputeKNNRefSynthetic(gdsProfile=gdsProfile,
#'     listEigenvector=demoPCASyntheticProfiles,
#'     listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
#'     studyIDSyn="MyStudy", spRef=demoKnownSuperPop1KG,
#'     fieldPopInfAnc="Superpop", kList=c(10, 11, 12),
#'     pcaList=c(13, 14, 15))
#'
#' ## Close GDS file (it is important to always close the GDS files)
#' closefn.gds(gdsProfile)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
paramRAIDS <- function(studyDF=NULL,
                       pedStudy=NULL,
                       studyType=NULL,
                       genoSource=NULL,
                       blockTypeId="GeneS.Ensembl.Hsapiens.v86",
                       reference="1KGc1.0",
                       genome="HG38",
                       chrInfo=NULL,
                       paramAncestry=NULL,
                       profileFile=NULL,
                       pathProfileGDS=NULL,
                       fileReferenceGDS=NULL,
                       fileReferenceAnnotGDS=NULL,
                       inferenceType="PCAnn",
                       batch=1,
                       offset=-1,
                       minCov=10,
                       minProb=0.999,
                       seqError=0.001,
                       np=1L,
                       listPos=NULL,
                       verbose=FALSE) {

    if(is.null(studyDF)){
        studyDF <- data.frame(study.id="NotDef",
                              study.desc="NotDef",
                              study.platform="NotDef",
                              stringsAsFactors=FALSE)
    }
    if(is.null(pedStudy)){
        pedStudy <- data.frame(Name.ID=c("ProfileId"),
                               Case.ID=c("ProfileId"),
                               Sample.Type=c("type"),
                               Diagnosis="NotDef",
                               Source=c("NotDef"),
                               stringsAsFactors=FALSE)
        row.names(pedStudy) <- pedStudy$Name.ID
    }
    if(is.null(chrInfo) && genome=="HG38"){
        chrInfo <- Seqinfo::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]
    }else if(is.null(chrInfo)){
        stop("If your genome is't HG38 you need to specify chrInfo\n")
    }
    if(is.null(paramAncestry)){
        paramAncestry=list(ScanBamParam=NULL,
                           PileupParam=NULL,
                           yieldSize=10000000)
    }
    # studyType=c("LD", "GeneAware")
    # genoSource=c("snp-pileup", "generic", "VCF", "bam"),
    # studyType=c("LD", "GeneAware"), np=1L, blockTypeID=NULL,

    parameter <- list(reference="1KG",
                      studyDF=studyDF,
                      pedStudy=pedStudy,
                      studyType=studyType,
                      genoSource=genoSource,
                      blockTypeId="GeneS.Ensembl.Hsapiens.v86",
                      genome=genome,
                      chrInfo=chrInfo,
                      profileFile=profileFile,
                      pathProfileGDS=pathProfileGDS,
                      fileReferenceGDS=fileReferenceGDS,
                      fileReferenceAnnotGDS=fileReferenceAnnotGDS,
                      inferenceType=inferenceType,
                      batch=batch,
                      offset=offset,
                      minCov=minCov,
                      minProb=minProb,
                      seqError=seqError,
                      np=np,
                      listPos=listPos,
                      nReads = 22,
                      K=33,
                      verbose=verbose)
    class(parameter) <- "parametersRAIDS"
    return(parameter)
}


