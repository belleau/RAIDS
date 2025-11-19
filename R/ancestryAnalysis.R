#' @title Compute the list of pruned SNVs for a specific profile using the
#' information from the Reference GDS file and a linkage disequilibrium
#' analysis
#'
#' @description This function computes the list of pruned SNVs for a
#' specific profile. When
#' a group of SNVs are in linkage disequilibrium, only one SNV from that group
#' is retained. The linkage disequilibrium is calculated with the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The initial list of
#' SNVs that are passed to the \code{\link[SNPRelate]{snpgdsLDpruning}}()
#' function can be specified by the user.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return The function returns a \code{parametersRAIDS} an object with all the RAIDS
#' parameters updated.
#'
#' @examples
#'
#' ## Required library for GDS
#' library(gdsfmt)
#'
#' ## Path to the demo Reference GDS file is located in this package
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
#' samplePED <- data.frame(Name.ID = c("ex1", "ex2"),
#'                     Case.ID = c("Patient_h11", "Patient_h12"),
#'                     Diagnosis = rep("Cancer", 2),
#'                     Sample.Type = rep("Primary Tumor", 2),
#'                     Source = rep("Databank B", 2), stringsAsFactors = FALSE)
#' rownames(samplePED) <- samplePED$Name.ID
#'
#' ## Temporary Profile GDS file
#' pathProfileGDS <- tempdir()
#' profileFile <- file.path(pathProfileGDS, "ex1.gds")
#'
#' pRAIDS <- RAIDS:::paramRAIDS(profileFile=file.path(dataDir, "ex1.txt.gz"),
#'     studyDF=studyDF,
#'     genoSource="snp-pileup",
#'     pedStudy=samplePED,
#'     pathProfileGDS=pathProfileGDS,
#'     fileReferenceGDS=fileGDS
#'     )

#' ## Copy the Profile GDS file demo that has not been pruned yet
#' file.copy(file.path(dataDir, "ex1_demo.gds"), profileFile)
#'
#'
#' ## Compute the list of pruned SNVs for a specific profile 'ex1'
#' ## and save it in the Profile GDS file 'ex1.gds'
#' pruningProfile2(pRAIDS)
#'

#'
#' ## Check content of Profile GDS file
#' ## The 'pruned.study' entry should be present
#' content <- openfn.gds(profileFile)
#' content
#'
#' ## Close the Profile GDS file (important)
#' closefn.gds(content)
#'
#' ## Remove Profile GDS file (created for demo purpose)
#' unlink(profileFile, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @export
pruningProfile2 <- function(pRAIDS) {

    ## Validate input parameters
    # validatePruningSample(gdsReference=gdsReference, method=method,
    #                       currentProfile=currentProfile, studyID=studyID, listSNP=listSNP,
    #                       slideWindowMaxBP=slideWindowMaxBP, thresholdLD=thresholdLD, np=np,
    #                       verbose=verbose, chr=chr, superPopMinAF=superPopMinAF,
    #                       keepPrunedGDS=keepPrunedGDS, pathProfileGDS=pathProfileGDS,
    #                       keepFile=keepFile, pathPrunedGDS=pathPrunedGDS, outPrefix=outPrefix)
    if(!is.null(pRAIDS$specificSNV)){
        pRAIDS <- specificSNVKeep(pRAIDS)
    }
    ## Matches a character method against a table of candidate values
    # method <- arg_match(method)
    gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)

    ## Profile GDS file name
    fileGDSSample <-  validateProfileGDSExist(pathProfile=pRAIDS$pathProfileGDS,
                                              profile=pRAIDS$pedStudy$Name.ID[1])

    # filePruned <- file.path(pathPrunedGDS, paste0(outPrefix, ".rds"))
    # fileObj <- file.path(pathPrunedGDS, paste0(outPrefix, ".Obj.rds"))

    snp.id <- read.gdsn(node=index.gdsn(gdsReference, "snp.id"))
    sample.id <- read.gdsn(node=index.gdsn(gdsReference, "sample.id"))


    ## Open the GDS Sample file
    gdsSample <- openfn.gds(filename=fileGDSSample)

    ## Extract all study information from the GDS Sample file
    study.annot <- read.gdsn(node=index.gdsn(gdsSample, "study.annot"))

    ## Select study information associated to the current profile
    posSample <- which(study.annot$data.id == pRAIDS$pedStudy$Name.ID[1] &
                           study.annot$study.id == pRAIDS$studyDF$study.id)

    ## Check that the information is found for the specified profile and study
    if(length(posSample) != 1) {
        closefn.gds(gdsSample)
        stop("In pruningSample the profile \'", pRAIDS$pedStudy$Name.ID[1],
             "\' doesn't exists for the study \'", pRAIDS$studyDF$study.id[1], "\'\n")
    }

    ## Get the SNV genotype information for the current profile
    # "geno.ref"
    if(pRAIDS$genoType %in% ls.gdsn(gdsProfile)){
        g <- read.gdsn(index.gdsn(gdsSample, pRAIDS$genoType),
                       start=c(1, posSample), count=c(-1,1))

    }else{
        stop("The genotype type don't exists")
    }

    ## Close the Profile GDS file
    closefn.gds(gdsSample)

    listGeno <- which(g != 3)
    rm(g)

    listKeepPos <- listGeno

    if(!is.null(pRAIDS$specificSNV)){
        listKeepPos <- intersect(pRAIDS$specificSNV$snvKeep[pRAIDS$specificSNV$snvKeep > 0], listKeepPos)
    }



    if (length(listKeepPos) == 0) {
        stop("In pruningSample, the sample ", currentProfile,
             " doesn't have SNPs after filters\n")
    }
    listKeep <- snp.id[listKeepPos]

    sample.ref <- read.gdsn(index.gdsn(gdsReference, "sample.ref"))
    listSamples <- sample.id[which(sample.ref == 1)]

    ## Use a LD analysis to generate a subset of SNPs
    snpset <- runLDPruning(gds=gdsReference, method=pRAIDS$pruningMethod,
                           listSamples=listSamples, listKeep=listKeep,
                           slideWindowMaxBP=pRAIDS$slideWindowMaxBP, thresholdLD=pRAIDS$thresholdLD,
                           np=pRAIDS$np, verbose=pRAIDS$verbose)
    pruned <- unlist(snpset, use.names=FALSE)

    gdsSample <- openfn.gds(filename=fileGDSSample, readonly=FALSE)
    addGDSStudyPruning(gdsProfile=gdsSample, pruned=pruned)
    closefn.gds(gdsfile=gdsSample)
    closefn.gds(gdsReference)

    return(pRAIDS)
}
