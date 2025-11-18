#' @title Create the Profile GDS file(s) for one or multiple specific profiles
#' using the information from a RDS Sample description file and the 1KG
#' GDS file
#'
#' @description The function uses the information for the Reference GDS file
#' and the RDS Sample Description file to create the Profile GDS file. One
#' Profile GDS file is created per profile. One Profile GDS file will be
#' created for each entry present in the \code{listProfiles} parameter.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return The function returns a \code{parametersRAIDS} an object with all the RAIDS
#' parameters updated.
#'
#' @examples
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG.gds")
#'
#' ## The data.frame containing the information about the study
#' ## The 3 mandatory columns: "study.id", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyDF <- data.frame(study.id = "MYDATA",
#'                         study.desc = "Description",
#'                         study.platform = "PLATFORM",
#'                         stringsAsFactors = FALSE)
#'
#' ## The data.frame containing the information about the samples
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' samplePED <- data.frame(Name.ID=c("ex1"),
#'                     Case.ID=c("Patient_h11"),
#'                     Diagnosis=rep("Cancer"),
#'                     Sample.Type=c("Primary Tumor"),
#'                     Source=c("Databank B"), stringsAsFactors=FALSE,
#'                     drop=FALSE)
#' rownames(samplePED) <- samplePED$Name.ID
#'
#' pRAIDS <- RAIDS:::paramRAIDS(profileFile=file.path(dataDir, "ex1.txt.gz"),
#'     studyDF=studyDF,
#'     genoSource="snp-pileup",
#'     pedStudy=samplePED,
#'     pathProfileGDS=tempdir(),
#'     fileReferenceGDS=fileGDS
#'     )
#' ## Create the Profile GDS File for samples in 'listSamples' vector
#' ## (in this case, samples "ex1")
#' ## The Profile GDS file is created in the pathProfileGDS directory
#' pRAIDS <- RAIDS:::createProfile2(pRAIDS)
#'
#' ## The function returns OL when successful
#' result
#'
#' ## The Profile GDS file 'ex1.gds' has been created in the
#' ## specified directory
#' list.files(pRAIDS$pathProfileGDS)
#'
#' ## Remove Profile GDS file (created for demo purpose)
#' unlink(file.path(pRAIDS$pathProfileGDS, "ex1.gds"), force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @keywords internal

createProfile2 <- function(pRAIDS) {
    # profileFile, profileName,
    # filePedRDS=NULL, pedStudy=NULL, fileNameGDS,
    # batch=1, studyDF, listProfiles=NULL,
    # pathProfileGDS=NULL,
    # genoSource=c("snp-pileup", "generic", "VCF", "bam"),
    # paramProfile=list(ScanBamParam=NULL,
    #                   PileupParam=NULL,
    #                   yieldSize=10000000),
    # verbose=FALSE
    ## When filePedRDS is defined and pedStudy is null
    if (!(is.null(pRAIDS$filePedRDS)) && is.null(pRAIDS$pedStudy)) {
        ## The filePedRDS must be a character string and the file must exists
        if (!(is.character(pRAIDS$filePedRDS) && (file.exists(pRAIDS$filePedRDS)))) {
            stop("The \'filePedRDS\' must be a character string representing",
                 " the RDS Sample information file. The file must exist.")
        }
        ## Open the RDS Sample information file
        pedStudy <- readRDS(file=pRAIDS$filePedRDS)
    } else if (!(is.null(pRAIDS$filePedRDS) || is.null(pRAIDS$pedStudy))) {
        stop("Both \'filePedRDS\' and \'pedStudy\' parameters cannot be ",
             "defined at the same time.")
    } else if (is.null(pRAIDS$filePedRDS) && is.null(pRAIDS$pedStudy)) {
        stop("One of the parameter \'fineNamePED\' of \'pedStudy\' must ",
             "be defined.")
    }

    ## Validate input parameters
    validatecreateProfile( pedStudy=pRAIDS$pedStudy,
                           fileNameGDS=pRAIDS$fileReferenceGDS, batch=pRAIDS$batch, studyDF=pRAIDS$studyDF,
                           listProfiles=pRAIDS$listProfiles, pathProfileGDS=pRAIDS$pathProfileGDS,
                           genoSource=pRAIDS$genoSource, verbose=pRAIDS$verbose)

    # genoSource <- arg_match(pRAIDS$genoSource)

    ## Read the Reference GDS file
    gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)

    ## Extract the chromosome and position information for all SNPs in 1KG GDS

    listPos <- NULL
    if( is.null(pRAIDS$listPos)){
        if(pRAIDS$genoSource == "bam"){
            alDf <- read.gdsn(index.gdsn(gdsReference, "snp.allele"))
            alDf <- matrix(unlist(strsplit(alDf,"\\/")),nrow=2)
            pRAIDS$listPos <- data.frame(snp.chromosome = read.gdsn(index.gdsn(gdsReference, "snp.chromosome")),
                                  snp.position = read.gdsn(index.gdsn(gdsReference, "snp.position")),
                                  REF = alDf[1,],
                                  ALT = alDf[2,],
                                  stringsAsFactors = FALSE
            )
            # listChr <- unique(listPos$chr)
            # We can optimize
            # listPos <- lapply(listChr,
            #                 FUN=function(x, varDf){
            #                     return(varDf[which(varDf$chr == x),])
            #                 },
            #                 varDf=listPos)
            # names(listPos) <- paste0("chr", listChr)
            rm(alDf)
        } else{
            pRAIDS$listPos <- data.frame(snp.chromosome=read.gdsn(index.gdsn(node=gdsReference, "snp.chromosome")),
                                  snp.position=read.gdsn(index.gdsn(node=gdsReference, "snp.position")))
        }
        #pRAIDS$listPos <- listPos
    }
    ## Create a data.frame containing the information

    # Need to reformat for bam in varDf

    if(pRAIDS$verbose) {
        message("Start ", Sys.time())
        message("Sample info DONE ", Sys.time())
    }

    generateProfileGDS2(pRAIDS=pRAIDS)

    if(pRAIDS$verbose) {
        message("Genotype DONE ", Sys.time())
    }

    ## Close 1KG GDS file
    closefn.gds(gdsReference)

    ## Return successful code
    return(pRAIDS)
}
