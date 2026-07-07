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
    keepPos <- seq_len(length(snp.id))
    if("snp.KeepDefault" %in% ls.gdsn(gdsReference) ){
        keepPos <- read.gdsn(index.gdsn(gdsReference, "snp.KeepDefault"))
    }
    snp.id <- snp.id[keepPos]
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
    if(pRAIDS$genoType %in% ls.gdsn(gdsSample)){
        g <- read.gdsn(index.gdsn(gdsSample, pRAIDS$genoType),
                       start=c(1, posSample), count=c(-1,1))

    }else{
        stop("The genotype type don't exists")
    }

    ## Close the Profile GDS file
    closefn.gds(gdsSample)

    listGeno <- which(g != 3)
    rm(g)
    # This position are the position in keepPos (I suppose the position in genotype match keepPos)
    listKeepPos <- listGeno

    if(!is.null(pRAIDS$specificSNV)){
        listKeepPos <- intersect(pRAIDS$specificSNV$snvKeep[pRAIDS$specificSNV$snvKeep > 0], listKeepPos)
    }



    if (length(listKeepPos) == 0) {
        stop("In pruningSample, the sample ", pRAIDS$pedStudy$Name.ID[1],
             " doesn't have SNPs after filters\n")
    }
    # Not snp.id match keepPos
    listKeep <- snp.id[listKeepPos]

    sample.ref <- read.gdsn(index.gdsn(gdsReference, "sample.ref"))
    listSamples <- sample.id[which(sample.ref == 1)]

    ## Use a LD analysis to generate a subset of SNPs
    snpset <- runLDPruning(gds=gdsReference, method=pRAIDS$pruningMethod,
                           listSamples=listSamples, listKeep=listKeep,
                           slideWindowMaxBP=pRAIDS$slideWindowMaxBP, thresholdLD=pRAIDS$thresholdLD,
                           np=pRAIDS$np, verbose=pRAIDS$verbose)
    pruned <- unlist(snpset, use.names=FALSE)

    gdsProfile <- openfn.gds(filename=fileGDSSample, readonly=FALSE)
    addGDSStudyPruning(gdsProfile=gdsProfile, pruned=pruned)
    closefn.gds(gdsfile=gdsProfile)
    closefn.gds(gdsReference)

    return(pRAIDS)
}

#' @title Run most steps leading to the ancestry inference call on a specific
#' RNA profile
#'
#' @description This function runs most steps leading to the ancestry inference
#' call on a specific RNA profile. First, the function creates the
#' Profile GDS file for the specific profile using the information from a
#' RDS Sample description file and the Population Reference GDS file.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return a \code{list} containing 4 entries:
#' \describe{
#' \item{\code{pcaSample}}{ a \code{list} containing the information related
#' to the eigenvectors. The \code{list} contains those 3 entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{eigenvector.ref}}{ a \code{matrix} of \code{numeric} containing
#' the eigenvectors for the reference profiles.}
#' \item{\code{eigenvector}}{ a \code{matrix} of \code{numeric} containing the
#' eigenvectors for the current profile projected on the PCA from the
#' reference profiles.}
#' }
#' }
#' \item{\code{paraSample}}{ a \code{list} containing the results with
#' different \code{D} and \code{K} values that lead to optimal parameter
#' selection. The \code{list} contains those entries:
#' \describe{
#' \item{\code{dfPCA}}{ a \code{data.frame} containing statistical results
#' on all combined synthetic results done with a fixed value of \code{D} (the
#' number of dimensions). The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{median}}{ a \code{numeric} representing the median of the
#' minimum AUROC obtained (within super populations) for all combination of
#' the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{mad}}{ a \code{numeric} representing the MAD of the minimum
#' AUROC obtained (within super populations) for all combination of the fixed
#' \code{D} value and all tested \code{K} values. }
#' \item{\code{upQuartile}}{ a \code{numeric} representing the upper quartile
#' of the minimum AUROC obtained (within super populations) for all
#' combination of the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{k}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for a fixed \code{D} value. }
#' }
#' }
#' \item{\code{dfPop}}{ a \code{data.frame} containing statistical results on
#' all combined synthetic results done with different values of \code{D} (the
#' number of dimensions) and \code{K} (the number of neighbors).
#' The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{AUROC.min}}{ a \code{numeric} representing the minimum accuracy
#' obtained by grouping all the synthetic results by super-populations, for
#' the specified values of \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing the accuracy obtained
#' by grouping all the synthetic results for the specified values of \code{D}
#' and \code{K}.}
#' \item{\code{Accu.CM}}{ a \code{numeric} representing the value of accuracy
#' of the confusion matrix obtained by grouping all the synthetic results for
#' the specified values of \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{dfAUROC}}{ a \code{data.frame} the summary of the results by
#' super-population. The \code{data.frame} contains
#' those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{Call}}{ a \code{character} string representing the
#' super-population.}
#' \item{\code{L}}{ a \code{numeric} representing the lower value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing  the AUROC obtained for
#' the fixed values of super-population, \code{D} and \code{K}.}
#' \item{\code{H}}{ a \code{numeric} representing the higher value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{D}}{ a \code{numeric} representing the optimal \code{D} value
#' (the number of dimensions) for the specific profile.}
#' \item{\code{K}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for the specific profile.}
#' \item{\code{listD}}{ a \code{numeric} representing the optimal \code{D}
#' values (the number of dimensions) for the specific profile. More than one
#' \code{D} is possible.}
#' }
#' }
#' \item{\code{KNNSample}}{  a \code{data.frame} containing the inferred
#' ancestry for different values of \code{K} and \code{D}. The
#' \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry. }
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry. }
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry for the specified \code{D} and \code{K} values.}
#' }
#' }
#' \item{\code{KNNSynthetic}}{  a \code{data.frame} containing the inferred
#' ancestry for each synthetic data for different values of \code{K} and
#' \code{D}.
#' The \code{data.frame}
#' contains those columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current synthetic data.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry. }
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry. }
#' \item{\code{infer.superPop}}{ a \code{character} string representing the
#' inferred ancestry for the specified \code{D} and \code{K} values.}
#' \item{\code{ref.superPop}}{ a \code{character} string representing the known
#' ancestry from the reference}
#' }
#' }
#' \item{\code{Ancestry}}{ a \code{data.frame} containing the inferred
#' ancestry for the current profile. The \code{data.frame} contains those
#' columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry.}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry.}
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry.}
#' }
#' }
#' }
#'
#' @details
#'
#' The runExomeAncestry() function generates 3 types of files
#' in the OUTPUT directory.
#' \describe{
#' \item{Ancestry Inference}{ The ancestry inference CSV file
#' (".Ancestry.csv" file)}
#' \item{Inference Informaton}{ The inference information RDS file
#' (".infoCall.rds" file)}
#' \item{Synthetic Information}{ The parameter information RDS files
#' from the synthetic inference ("KNN.synt.*.rds" files in a sub-directory)}
#' }
#'
#' In addition, a sub-directory (named using the profile ID) is
#' also created.
#'
#' @references
#'
#' Galinsky KJ, Bhatia G, Loh PR, Georgiev S, Mukherjee S, Patterson NJ,
#' Price AL. Fast Principal-Component Analysis Reveals Convergent Evolution
#' of ADH1B in Europe and East Asia. Am J Hum Genet. 2016 Mar 3;98(3):456-72.
#' doi: 10.1016/j.ajhg.2015.12.022. Epub 2016 Feb 25.
#'
#' @examples
#'
#' ## Required library for GDS
#' library(SNPRelate)
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#'
#' #################################################################
#' ## The 1KG GDS file and the 1KG SNV Annotation GDS file
#' ## need to be located in the same directory
#' ## Note that the 1KG GDS file used for this example is a
#' ## simplified version and CANNOT be used for any real analysis
#' #################################################################
#' path1KG <- file.path(dataDir, "tests")
#'
#' fileReferenceGDS <- file.path(path1KG, "ex1_good_small_1KG.gds")
#' fileAnnotGDS <- file.path(path1KG, "ex1_good_small_1KG_Annot.gds")
#'
#' #################################################################
#' ## The Sample SNP pileup files (one per sample) need
#' ## to be located in the same directory.
#' #################################################################
#' demoProfileEx1 <- file.path(dataDir, "example", "snpPileup", "ex1.txt.gz")
#'
#' #################################################################
#' ## The path where the Profile GDS Files (one per sample)
#' ## will be created need to be specified.
#' #################################################################
#' pathProfileGDS <- file.path(tempdir(), "out.tmp")
#'
#' ####################################################################
#' ## Fix seed to ensure reproducible results
#' ####################################################################
#' set.seed(3043)
#'
#' gds1KG <- snpgdsOpen(fileReferenceGDS)
#' dataRef <- select1KGPop(gds1KG, nbProfiles=2L)
#' closefn.gds(gds1KG)
#'
#' ## Required library for this example to run correctly
#' if (requireNamespace("Seqinfo", quietly=TRUE) &&
#'      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#'
#'     ## Chromosome length information
#'     ## chr23 is chrX, chr24 is chrY and chrM is 25
#'     chrInfo <- Seqinfo::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]
#'
#'     \donttest{
#'
#'         res <- inferAncestryGeneAware(profileFile=demoProfileEx1,
#'             pathProfileGDS=pathProfileGDS,
#'             fileReferenceGDS=fileReferenceGDS,
#'             fileReferenceAnnotGDS=fileAnnotGDS,
#'             chrInfo=chrInfo,
#'             syntheticRefDF=dataRef,
#'             blockTypeID="GeneS.Ensembl.Hsapiens.v86",
#'             genoSource="snp-pileup")
#'
#'         unlink(pathProfileGDS, recursive=TRUE, force=TRUE)
#'
#'     }
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom utils write.csv
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @export
ancestryInferencePCAKNN <- function(pRAIDS=paramRAIDS()) {
    # profileFile, pathProfileGDS,
    # fileReferenceGDS, fileReferenceAnnotGDS,
    # chrInfo, syntheticRefDF=NULL,
    # genoSource=c("snp-pileup", "generic", "VCF", "bam"),
    # studyType=NULL,
    # ancestryType=NULL,
    #

    profileBaseName <- basename(pRAIDS$profileFile)
    pathGeno <- dirname(pRAIDS$profileFile)

    if(is.null(pRAIDS$studyType) || !(pRAIDS$studyType %in% c("LD", "GeneAware"))){
        stop("The study must be specify as LD or GeneAware in the studyType or paramRAIDS.\n")
    }

    # if(is.null(ancestryType)){
    #     studyType <- pRAIDS$studyType
    # }else{
    #     pRAIDS$studyType <- studyType
    # }




    #genoSource <- arg_match(genoSource)

    # if(genoSource == "bam"){
    #     stop("The bam is not release yet look to get a \'Devel\' version ",
    #             "or contact us")
    # }

    profileName <- gsub("\\.gz$", "", profileBaseName, ignore.case = TRUE)
    for(extCur in c( "\\.vcf$", "\\.txt$", "\\.bam", "\\.tsv", "\\.csv")){
        profileName <- gsub(extCur, "", profileName, ignore.case = TRUE)
    }
    #profileName <- "profile"
    pRAIDS$pedStudy$Name.ID[1] <- profileName
    pRAIDS$pedStudy$Case.ID[1] <- profileName

    if(file.exists(file.path(pRAIDS$pathProfileGDS, paste0(pRAIDS$pedStudy$Name.ID[1], ".gds")))){
        stop(paste0("The gds file for ", pRAIDS$pedStudy$Name.ID[1], " already exist."))
    }

    ## Generate the synthetic dataset if not included
    if(is.null(pRAIDS$syntheticRefDF)){
        if(pRAIDS$reference == "1KGc1.0"){
            pRAIDS$syntheticRefDF <- select1KGPopForSynthetic(fileReferenceGDS=pRAIDS$fileReferenceGDS,
                                                              nbProfiles=30L)
        }else if(pRAIDS$reference == "HGDP1kgV0.1"){
            # new synthetic profile selection function for HGDP1kgV0.1
            pRAIDS$syntheticRefDF <- selectHGDP1kgPopForSynthetic(fileReferenceGDS=pRAIDS$fileReferenceGDS,
                                                              nbProfiles=30L)
        }
    }

    # pRAIDS$profileFile <- profileFile
    # pRAIDS$pathProfileGDS <- pathProfileGDS
    # pRAIDS$fileReferenceGDS <- fileReferenceGDS
    # pRAIDS$fileReferenceAnnotGDS <- fileReferenceAnnotGDS


    ## Validate parameters
    validateRunExomeOrRNAAncestry(pedStudy=pRAIDS$pedStudy,
                                  studyDF=pRAIDS$studyDF,
                                  pathGeno=pathGeno,
                                  pathProfileGDS=pRAIDS$pathProfileGDS,
                                  pathOut="./",
                                  fileReferenceGDS=pRAIDS$fileReferenceGDS,
                                  fileReferenceAnnotGDS=pRAIDS$fileReferenceAnnotGDS,
                                  chrInfo=pRAIDS$chrInfo,
                                  syntheticRefDF=pRAIDS$syntheticRefDF,
                                  genoSource=pRAIDS$genoSource, verbose=pRAIDS$verbose)
    ##########################
    ## Create profile section
    ##########################

    ## parse the read.count if the pRAIDS$listPos is not
    ## it define it and return an update pRAIDS
    pRAIDS <- generateProfileRawGDS2(pRAIDS)

    ## Call the genotype from Ref.count and Alt.count
    generateProfileGenoCall(pRAIDS)

    ## Pruning
    pruningProfile2(pRAIDS)
    ## add genotype structure for the pruning positions
    addGenotypeProfile(pRAIDS)




    ##########################
    ## Section allelic fraction
    ##########################

    studyTypeLeg <- ifelse(pRAIDS$studyType=="LD", "DNA", "RNA")
    ## Profile GDS file name
    fileProfileGDS <-  validateProfileGDSExist(pathProfile=pRAIDS$pathProfileGDS,
                                               profile=pRAIDS$pedStudy$Name.ID[1])
    ## Open Profile GDS file
    gdsProfile <- snpgdsOpen(fileProfileGDS, readonly=FALSE)

    gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)
    gdsRefAnnot <- openfn.gds(pRAIDS$fileReferenceAnnotGDS)

    estimateAllelicFraction(gdsReference=gdsReference, gdsProfile=gdsProfile,
                            currentProfile=pRAIDS$pedStudy$Name.ID[1], studyID=pRAIDS$studyDF$study.id,
                            chrInfo=pRAIDS$chrInfo, studyType=studyTypeLeg, gdsRefAnnot=gdsRefAnnot,
                            blockID=pRAIDS$blockTypeId, verbose=pRAIDS$verbose)
    closefn.gds(gdsProfile)
    closefn.gds(gdsReference)
    closefn.gds(gdsRefAnnot)

    ##########################
    ## Section synthetic
    ##########################

    ## Add information related to the synthetic profiles in Profile GDS file
    prepSynthetic(fileProfileGDS=fileProfileGDS,
                  listSampleRef=pRAIDS$syntheticRefDF$sample.id,  profileID=pRAIDS$pedStudy$Name.ID[1],
                  studyDF=pRAIDS$studyDFSyn, prefix=pRAIDS$prefix, verbose=pRAIDS$verbose)

    if(pRAIDS$verbose){
        message("syntheticGeno start ", Sys.time())
    }

    # gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)
    # gdsRefAnnot <- openfn.gds(pRAIDS$fileReferenceAnnotGDS)
    # function(gdsReference, gdsRefAnnot, fileProfileGDS, profileID,
    #          listSampleRef, nbSim=1L, prefix="", pRecomb=0.01, minProb=0.999,
    #          seqError=0.001)
    # resG <- syntheticGeno(gdsReference=gdsReference, gdsRefAnnot=gdsRefAnnot,
    #                       fileProfileGDS=fileProfileGDS, profileID=pRAIDS$pedStudy$Name.ID[1],
    #                       listSampleRef=pRAIDS$syntheticRefDF$sample.id,
    #                       nbSim=pRAIDS$nbSim, prefix=pRAIDS$prefix,
    #                       pRecomb=pRAIDS$pRecomb, minProb=pRAIDS$minProb,
    #                       seqError=pRAIDS$seqErrorSyn)
    resG <-syntheticGeno2(pRAIDS)

    ##########################
    ## Section ancestry
    ##########################

    # get the superPop of the sample.id in the ref.
    # names(spRef) == sample.id spRef = vector of superPop
    # gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)
    # gdsRefAnnot <- openfn.gds(pRAIDS$fileReferenceAnnotGDS)

    spRef <- getRefSuperPop2(pRAIDS)
    # closefn.gds(gdsReference)
    # closefn.gds(gdsRefAnnot)
    # get a list of each pop.group in pRAIDS$syntheticRefDF
    sampleRM <- splitSelectByPop(pRAIDS$syntheticRefDF)

    # gdsProfile <- snpgdsOpen(fileProfileGDS)

    if(pRAIDS$verbose){
        message("SyntheticAncestry start ", Sys.time())
    }


    resSyn <- lapply(seq_len(nrow(sampleRM)),
                     FUN=function(x, sampleRM,
                                  pRAIDS, spRef) {
                         synthKNN <- computePoolSyntheticAncestryGr2(sampleRM=sampleRM[x,],
                                                                     spRef=spRef,
                                                                     listCatPop=unique(spRef),
                                                                     pRAIDS = pRAIDS)
                         return(synthKNN$matKNN)
                     }, sampleRM=sampleRM, pRAIDS=pRAIDS,
                     spRef=spRef)

    if(pRAIDS$verbose){
        message("SyntheticAncestry end ", Sys.time())
    }
    resSyn <- do.call(rbind, resSyn)

    gdsProfile <- snpgdsOpen(fileProfileGDS)
    gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)
    # gdsRefAnnot <- openfn.gds(pRAIDS$fileReferenceAnnotGDS)
    ## Extract the super-population information from the 1KG GDS file
    ## for profiles associated to the synthetic study
    pedSyn <- prepPedSynthetic1KG(gdsReference=gdsReference,
                        gdsSample=gdsProfile,
                        studyID=pRAIDS$studyDFSyn$study.id,
                        popName="superPop")
    ## Close Profile GDS file (important)
    closefn.gds(gdsProfile)
    closefn.gds(gdsReference)

    # bla <- list(resSyn=resSyn,
    #             pedSyn=pedSyn,
    #             spRef=spRef,
    #             pRAIDS=pRAIDS)
    # saveRDS(bla, "data/data.bck/bla.rds")
    # rm(bla)

    resCall <- computeAncestryFromSynthetic2(syntheticKNN=resSyn,
                            pedSyn=pedSyn,
                            spRef=spRef,
                            listCatPop=unique(spRef),
                            pRAIDS=pRAIDS)

    # resCall <- computeAncestryFromSynthetic(gdsReference=gdsReference,
    #                                         gdsProfile=gdsProfile,
    #                                         syntheticKNN=resSyn,
    #                                         pedSyn=pedSyn,
    #                                         currentProfile=pRAIDS$pedStudy$Name.ID[1],
    #                                         spRef=spRef,
    #                                         studyIDSyn=pRAIDS$studyDFSyn$study.id,
    #                                         np=pRAIDS$np,
    #                                         listCatPop=unique(spRef),
    #                                         fieldPopInRef=pRAIDS$fieldPopInRef,
    #                                         fieldPopInfAnc=pRAIDS$fieldPopInfAnc,
    #                                         kList=pRAIDS$kList,
    #                                         pcaList=pRAIDS$pcaList[pRAIDS$pcaList <= pRAIDS$eigenCount],
    #                                         algorithm=pRAIDS$PCAalgorithm,
    #                                         eigenCount=pRAIDS$eigenCount,
    #                                         missingRate=pRAIDS$missingRate,
    #                                         verbose=pRAIDS$verbose)
    if(pRAIDS$verbose){
        message("Ancestry end ", Sys.time())
    }


    #closefn.gds(gdsRefAnnot)

    resSyn[[paste0("ref.superPop")]] <- pedSyn[resSyn$sample.id, "superPop"]

    colnames(resSyn) <- c("sample.id", "D", "K", "infer.superPop",
                          "ref.superPop")

    res <- list(pcaSample=resCall$pcaSample, # PCA of the profile + 1KG
                paraSample=resCall$paraSample, # Result of the parameter selection
                KNNSample=resCall$KNNSample$matKNN, # KNN for the profile
                KNNSynthetic=resSyn, # KNN results for synthetic data
                Ancestry=resCall$Ancestry) # the ancestry call fo the profile
    ## Successful
    return(res)
}

#' @title Run a k-nearest neighbors analysis on a subset of the
#' synthetic dataset
#'
#' @description The function runs k-nearest neighbors analysis on a
#' subset of the synthetic data set. The function uses the 'knn' package.
#'
#' @param listEigenvector a \code{list} with 3 entries:
#' 'sample.id', 'eigenvector.ref' and 'eigenvector'. The \code{list} represents
#' the PCA done on the 1KG reference profiles and the synthetic profiles
#' projected onto it.
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations. Default:
#' \code{c("EAS", "EUR", "AFR", "AMR", "SAS")}.
#'
#' @param spRef \code{vector} of \code{character} strings representing the
#' known super population ancestry for the 1KG profiles. The 1KG profile
#' identifiers are used as names for the \code{vector}.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return a \code{list} containing 4 entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{vector} of \code{character} strings
#' representing the identifiers of the synthetic profiles analysed.}
#' \item{\code{sample1Kg}}{ a \code{vector} of \code{character} strings
#' representing the identifiers of the 1KG reference profiles used to
#' generate the synthetic profiles.}
#' \item{\code{sp}}{ a \code{vector} of \code{character} strings representing
#' the known super population ancestry of the 1KG reference profiles used
#' to generate the synthetic profiles.}
#' \item{\code{matKNN}}{ a \code{data.frame} containing the super population
#' inference for each synthetic profiles for different values of PCA
#' dimensions \code{D} and k-neighbors values \code{K}. The fourth column title
#' corresponds to the \code{fieldPopInfAnc} parameter.
#' The \code{data.frame} contains 4 columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing
#' the identifier of the synthetic profile analysed.}
#' \item{\code{D}}{ a \code{numeric} strings representing
#' the value of the PCA dimension used to infer the super population.}
#' \item{\code{K}}{ a \code{numeric} strings representing
#' the value of the k-neighbors used to infer the super population.}
#' \item{\code{fieldPopInfAnc} value}{ a \code{character} string representing
#' the inferred ancestry.}
#' }
#' }
#' }
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Load the demo PCA on the synthetic profiles projected on the
#' ## demo 1KG reference PCA
#' data(demoPCASyntheticProfiles)
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoKNNSynthetic", package="RAIDS")
#'
#' ## Open the Profile GDS file
#' gdsProfile <- snpgdsOpen(file.path(dataDir, "ex1.gds"))
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA.Synthetic"
#'
#' ## Projects synthetic profiles on 1KG PCA
#' results <- computeKNNRefSynthetic(gdsProfile=gdsProfile,
#'     listEigenvector=demoPCASyntheticProfiles,
#'     listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"), studyIDSyn=studyID,
#'     spRef=demoKnownSuperPop1KG)
#'
#' ## The inferred ancestry for the synthetic profiles for different values
#' ## of D and K
#' head(results$matKNN)
#'
#' ## Close Profile GDS file (important)
#' closefn.gds(gdsProfile)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom class knn
#' @encoding UTF-8
#' @export
computeKNNRefSyn <- function(listEigenvector,
                             spRef,
                             listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
                             pRAIDS) {

    fileGDSProfile <- file.path(pRAIDS$pathProfileGDS,
                                paste0(pRAIDS$pedStudy$Name.ID[1], ".gds"))
    gdsProfile <- openfn.gds(filename=fileGDSProfile)

    ## Validate the input parameters
    validateComputeKNNRefSynthetic(gdsProfile=gdsProfile,
                                   listEigenvector=listEigenvector,
                                   listCatPop=listCatPop, studyIDSyn=pRAIDS$studyDFSyn$study.id, spRef=spRef,
                                   fieldPopInfAnc=pRAIDS$fieldPopInfAnc, kList=pRAIDS$kList, pcaList=pRAIDS$pcaList)

    ## Get study information from the GDS Sample file
    studyAnnotAll <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))
    closefn.gds(gdsProfile)
    studyAnnot <- studyAnnotAll[which(studyAnnotAll$study.id ==
                                          pRAIDS$studyDFSyn$study.id & studyAnnotAll$data.id %in%
                                          listEigenvector$sample.id), ]
    pcaList <- pRAIDS$pcaList[pRAIDS$pcaList <= pRAIDS$eigenCountSyn]

    listMat <- lapply(pcaList,
                      FUN=function(pcaD, listEigenvector,
                                   listCatPop, pRAIDS){
                          resMat <- lapply(pRAIDS$kList,
                                           FUN=function(kV, pcaD,
                                                        listEigenvector,
                                                        listCatPop, pRAIDS){
                                               return(computeKNNProfile(listEigenvector,kV, pcaD, pRAIDS))
                                           },
                                           pcaD=pcaD,
                                           listEigenvector=listEigenvector,
                                           listCatPop=listCatPop,
                                           pRAIDS=pRAIDS)
                          resMat <- do.call(rbind, resMat)
                          return(resMat)
                      },
                      listEigenvector=listEigenvector,
                      listCatPop=listCatPop,
                      pRAIDS=pRAIDS)

    resMat <- do.call(rbind, listMat)

    listKNN <- list(sample.id=listEigenvector$sample.id,
                    sample1Kg=studyAnnot$case.id,
                    sp=spRef[studyAnnot$case.id],
                    matKNN=resMat)

    return(listKNN)
}

#' @title Run a k-nearest neighbors analysis on a subset of the
#' synthetic dataset
#'
#' @description The function runs k-nearest neighbors analysis on a
#' subset of the synthetic data set. The function uses the 'knn' package.
#'
#' @param listEigenvector a \code{list} with 3 entries:
#' 'sample.id', 'eigenvector.ref' and 'eigenvector'. The \code{list} represents
#' the PCA done on the 1KG reference profiles and the synthetic profiles
#' projected onto it.
#'
#' @param dfRef \code{data.frame} TODO of \code{character}
#' strings representing the
#' known super population ancestry for the 1KG profiles. The 1KG profile
#' identifiers are used as names for the \code{vector}.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return a \code{list} containing 4 entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{vector} of \code{character} strings
#' representing the identifiers of the synthetic profiles analysed.}
#' \item{\code{sample1Kg}}{ a \code{vector} of \code{character} strings
#' representing the identifiers of the 1KG reference profiles used to
#' generate the synthetic profiles.}
#' \item{\code{sp}}{ a \code{vector} of \code{character} strings representing
#' the known super population ancestry of the 1KG reference profiles used
#' to generate the synthetic profiles.}
#' \item{\code{matKNN}}{ a \code{data.frame} containing the super population
#' inference for each synthetic profiles for different values of PCA
#' dimensions \code{D} and k-neighbors values \code{K}. The fourth column title
#' corresponds to the \code{fieldPopInfAnc} parameter.
#' The \code{data.frame} contains 4 columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing
#' the identifier of the synthetic profile analysed.}
#' \item{\code{D}}{ a \code{numeric} strings representing
#' the value of the PCA dimension used to infer the super population.}
#' \item{\code{K}}{ a \code{numeric} strings representing
#' the value of the k-neighbors used to infer the super population.}
#' \item{\code{fieldPopInfAnc} value}{ a \code{character} string representing
#' the inferred ancestry.}
#' }
#' }
#' }
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn openfn.gds closefn.gds
#' @importFrom class knn
#' @encoding UTF-8
#' @export
computeKNNRefSynGeneric <- function(listEigenvector, dfRef, pRAIDS) {

    fileGDSProfile <- file.path(pRAIDS$pathProfileGDS,
                                    paste0(pRAIDS$pedStudy$Name.ID[1], ".gds"))
    gdsProfile <- openfn.gds(filename=fileGDSProfile)

    ## Validate the input parameters
    # validateComputeKNNRefSynthetic(gdsProfile=gdsProfile,
    # listEigenvector=listEigenvector,
    # listCatPop=listCatPop, studyIDSyn=pRAIDS$studyDFSyn$study.id, spRef=spRef,
    # fieldPopInfAnc=pRAIDS$fieldPopInfAnc, kList=pRAIDS$kList, pcaList=pRAIDS$pcaList)

    ## Get study information from the GDS Sample file
    studyAnnotAll <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))
    closefn.gds(gdsProfile)
    studyAnnot <- studyAnnotAll[which(studyAnnotAll$study.id ==
                    pRAIDS$studyDFSyn$study.id & studyAnnotAll$data.id %in%
                    listEigenvector$sample.id), ]
    pcaList <- pRAIDS$pcaList[pRAIDS$pcaList <= pRAIDS$eigenCountSyn]

    listMat <- lapply(pcaList,
                    FUN=function(pcaD, listEigenvector, dfRef, pRAIDS){
                        resMat <- lapply(pRAIDS$kList,
                            FUN=function(kV, pcaD, listEigenvector,
                                    dfRef, pRAIDS){
                                return(computeKNNProfileSubSet(listEigenvector, 
                                                kV, pcaD,dfRef, pRAIDS))
                            },
                            pcaD=pcaD,
                            listEigenvector=listEigenvector,
                            dfRef=dfRef,
                            pRAIDS=pRAIDS)
                        resMat <- do.call(rbind, resMat)
                        return(resMat)
                      },
                      listEigenvector=listEigenvector,
                      dfRef=dfRef,
                      pRAIDS=pRAIDS)

    resMat <- do.call(rbind, listMat)

    listKNN <- list(sample.id=listEigenvector$sample.id,
                    sample1Kg=studyAnnot$case.id,
                    sp=dfRef[studyAnnot$case.id,],
                    matKNN=resMat)

    return(listKNN)
}


#' @title Run a PCA analysis and a K-nearest neighbors analysis on a small set
#' of synthetic data using all 1KG profiles except the ones used to generate
#' the synthetic profiles
#'
#' @description The function runs a PCA analysis using 1 synthetic profile
#' from each sub-continental population. The reference profiles used to
#' create those synthetic profiles are first removed from the list
#' of 1KG reference profiles that generates the reference PCA. Then, the
#' retained synthetic
#' profiles are projected on the 1KG PCA space. Finally, a K-nearest neighbors
#' analysis using a range of K and D values is done.
#'
#' @param sampleRM a \code{vector} of \code{character} strings representing
#' the identifiers of the 1KG reference profiles that should not be used to
#' create the reference PCA. There should be one per sub-continental
#' population. Those profiles are
#' removed because those have been used to generate the synthetic profiles
#' that are going to be analysed here. The sub-continental
#' identifiers are used as names for the \code{vector}.
#'
#' @param spRef \code{vector} of \code{character} strings representing the
#' known super population ancestry for the 1KG profiles. The 1KG profile
#' identifiers are used as names for the \code{vector}.
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations. Default:
#' \code{("EAS", "EUR", "AFR", "AMR", "SAS")}.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object containing all RAIDS
#' parameters.
#'
#' @return a \code{list} containing the following entries:
#' \describe{
#' \item{sample.id}{ a \code{vector} of \code{character} strings representing
#' the identifiers of the synthetic profiles. }
#' \item{sample1Kg}{ a \code{vector} of \code{character} strings representing
#' the identifiers of the reference 1KG profiles used to generate the
#' synthetic profiles. }
#' \item{sp}{ a \code{vector} of \code{character} strings representing the
#' known ancestry for the reference 1KG profiles used to generate the
#' synthetic profiles. }
#' \item{matKNN}{ a \code{data.frame} containing 4 columns. The first column
#' 'sample.id' contains the name of the synthetic profile. The second column
#' 'D' represents the dimension D used to infer the ancestry. The third column
#' 'K' represents the number of neighbors K used to infer the ancestry. The
#' fourth column 'SuperPop' contains the inferred ancestry. }
#' }
#'
#' @references
#'
#' Galinsky KJ, Bhatia G, Loh PR, Georgiev S, Mukherjee S, Patterson NJ,
#' Price AL. Fast Principal-Component Analysis Reveals Convergent Evolution
#' of ADH1B in Europe and East Asia. Am J Hum Genet. 2016 Mar 3;98(3):456-72.
#' doi: 10.1016/j.ajhg.2015.12.022. Epub 2016 Feb 25.
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom rlang arg_match
#' @importFrom gdsfmt openfn.gds closefn.gds
#' @encoding UTF-8
#' @export
computePoolSyntheticAncestryGr2 <- function(sampleRM, spRef,
    listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"), pRAIDS) {

    fileGDSProfile <- file.path(pRAIDS$pathProfileGDS,
                                paste0(pRAIDS$pedStudy$Name.ID[1], ".gds"))
    gdsProfile <- openfn.gds(filename=fileGDSProfile)
    ## Validate the input parameters
    validateComputePoolSyntheticAncestryGr(gdsProfile=gdsProfile,
        sampleRM=sampleRM, spRef=spRef,
        studyIDSyn=pRAIDS$studyDFSyn$study.id,
        np=pRAIDS$np,
        listCatPop=listCatPop,
        pcaList=pRAIDS$pcaList,
        fieldPopInfAnc=pRAIDS$fieldPopInfAnc,
        kList=pRAIDS$kList,
        algorithm=pRAIDS$PCAalgorithm,
        eigenCount=pRAIDS$eigenCountSyn,
        missingRate=pRAIDS$PCAmissingRate,
        verbose=pRAIDS$verbose)
    
    if(pRAIDS$verbose){
        message("computePoolSyntheticAncestryGr2 p1 ", Sys.time())
    }
    # Get the list of syn profiles to process
    study.annot <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))
    study.annot <- study.annot[which(study.annot$study.id == pRAIDS$studyDFSyn$study.id &
                        study.annot$case.id %in% sampleRM),]
    closefn.gds(gdsProfile)
    ## Set algorithm
    #algorithm <- arg_match(algorithm)

    ## Calculate Principal Component Analysis (PCA) on SNV genotype dataset
    ## excluded the selected profiles used to generate the synthetic profiles
    if(pRAIDS$verbose){
        message("computePoolSyntheticAncestryGr2 p2 ", Sys.time())
    }
    pcaRef <- computePCARefRMMulti1(listRM=sampleRM, pRAIDS=pRAIDS)

    ## Calculate PCA on the synthetic profiles using 1KG PCA results

    # resPCA <- computePCAMultiSynthetic(gdsProfile=gdsProfile, listPCA=pca1KG,
    #                                    sampleRef=sampleRM, studyIDSyn=studyIDSyn, verbose=verbose)
    if(pRAIDS$verbose){
        message("computePoolSyntheticAncestryGr2 p3 ", Sys.time())
    }
    resPCA <- computePCAProfile(listPCA=pcaRef, profileId =study.annot$data.id , pRAIDS)
    ## Calculate the k-nearest neighbor analyses on a subset of the
    ## synthetic data set
    # synthKNN <- computeKNNRefSynthetic(gdsProfile=gdsProfile,
    #                                    listEigenvector=resPCA,
    #                                    listCatPop=listCatPop, studyIDSyn=studyIDSyn, spRef=spRef,
    #                                    fieldPopInfAnc=fieldPopInfAnc, kList=kList, pcaList=pcaList)
    if(pRAIDS$verbose){
        message("computePoolSyntheticAncestryGr2 p4 ", Sys.time())
    }
    synthKNN <- computeKNNRefSyn(listEigenvector=resPCA,
                                 spRef=spRef,
                                 listCatPop=listCatPop,
                                 pRAIDS=pRAIDS)
    if(pRAIDS$verbose){
        message("computePoolSyntheticAncestryGr2 end ", Sys.time())
    }
    return(synthKNN)
}

#' @title Select the optimal K and D parameters using the synthetic data and
#' infer the ancestry of a specific profile
#'
#' @description The function select the optimal K and D parameters for a
#' specific profile. The results on the synthetic data are used for the
#' parameter selection. Once the optimal parameters are selected, the
#' ancestry is inferred for the specific profile.
#'
#' @param syntheticKNN a \code{vector} of \code{character} strings representing
#' the name of files that contain the results of ancestry inference done on
#' the synthetic profiles for multiple values of _D_ and _K_. The files must
#' exist.
#'
#' @param pedSyn a \code{data.frame} containing the columns extracted from the
#' GDS Sample 'study.annot' node with a extra column named as the 'popName'
#' parameter that has been extracted from the 1KG GDS 'sample.annot' node.
#'
#' @param spRef a \code{vector} of \code{character} strings representing the
#' known super population ancestry for the 1KG profiles. The 1KG profile
#' identifiers are used as names for the \code{vector}.
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations. Default:
#' \code{("EAS", "EUR", "AFR", "AMR", "SAS")}.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return a \code{list} containing 4 entries:
#' \describe{
#' \item{\code{pcaSample}}{ a \code{list} containing the information related
#' to the eigenvectors. The \code{list} contains those 3 entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{eigenvector.ref}}{ a \code{matrix} of \code{numeric} containing
#' the eigenvectors for the reference profiles.}
#' \item{\code{eigenvector}}{ a \code{matrix} of \code{numeric} containing the
#' eigenvectors for the current profile projected on the PCA from the
#' reference profiles.}
#' }
#' }
#' \item{\code{paraSample}}{ a \code{list} containing the results with
#' different \code{D} and \code{K} values that lead to optimal parameter
#' selection. The \code{list} contains those entries:
#' \describe{
#' \item{\code{dfPCA}}{ a \code{data.frame} containing statistical results
#' on all combined synthetic results done with a fixed value of \code{D} (the
#' number of dimensions). The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{median}}{ a \code{numeric} representing the median of the
#' minimum AUROC obtained (within super populations) for all combination of
#' the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{mad}}{ a \code{numeric} representing the MAD of the minimum
#' AUROC obtained (within super populations) for all combination of the fixed
#' \code{D} value and all tested \code{K} values. }
#' \item{\code{upQuartile}}{ a \code{numeric} representing the upper quartile
#' of the minimum AUROC obtained (within super populations) for all
#' combination of the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{k}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for a fixed \code{D} value. }
#' }
#' }
#' \item{\code{dfPop}}{ a \code{data.frame} containing statistical results on
#' all combined synthetic results done with different values of \code{D} (the
#' number of dimensions) and \code{K} (the number of neighbors).
#' The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{AUROC.min}}{ a \code{numeric} representing the minimum accuracy
#' obtained by grouping all the synthetic results by super-populations, for
#' the specified values of \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing the accuracy obtained
#' by grouping all the synthetic results for the specified values of \code{D}
#' and \code{K}.}
#' \item{\code{Accu.CM}}{ a \code{numeric} representing the value of accuracy
#' of the confusion matrix obtained by grouping all the synthetic results for
#' the specified values of \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{dfAUROC}}{ a \code{data.frame} the summary of the results by
#' super-population. The \code{data.frame} contains
#' those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{Call}}{ a \code{character} string representing the
#' super-population.}
#' \item{\code{L}}{ a \code{numeric} representing the lower value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing  the AUROC obtained for the
#' fixed values of super-population, \code{D} and \code{K}.}
#' \item{\code{H}}{ a \code{numeric} representing the higher value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{D}}{ a \code{numeric} representing the optimal \code{D} value
#' (the number of dimensions) for the specific profile.}
#' \item{\code{K}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for the specific profile.}
#' \item{\code{listD}}{ a \code{numeric} representing the optimal \code{D}
#' values (the number of dimensions) for the specific profile. More than one
#' \code{D} is possible.}
#' }
#' }
#' \item{\code{KNNSample}}{ a \code{list} containing the inferred ancestry
#' using different \code{D} and \code{K} values. The \code{list} contains
#' those entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{matKNN}}{ a \code{data.frame} containing the inferred ancestry
#' for different values of \code{K} and \code{D}. The \code{data.frame}
#' contains those columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry. }
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry. }
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry for the specified \code{D} and \code{K} values.}
#' }
#' }
#' }
#' }
#' \item{\code{Ancestry}}{ a \code{data.frame} containing the inferred
#' ancestry for the current profile. The \code{data.frame} contains those
#' columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry.}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry.}
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry.}
#' }
#' }
#' }
#'
#'
#' @references
#'
#' Galinsky KJ, Bhatia G, Loh PR, Georgiev S, Mukherjee S, Patterson NJ,
#' Price AL. Fast Principal-Component Analysis Reveals Convergent Evolution
#' of ADH1B in Europe and East Asia. Am J Hum Genet. 2016 Mar 3;98(3):456-72.
#' doi: 10.1016/j.ajhg.2015.12.022. Epub 2016 Feb 25.
#'
#' @examples
#'
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## The Reference GDS file
#' path1KG <- system.file("extdata/tests", package="RAIDS")
#'
#' ## Open the Reference GDS file
#' gdsRef <- snpgdsOpen(file.path(path1KG, "ex1_good_small_1KG.gds"))
#'
#' ## Path to the demo synthetic results files
#' ## List of the KNN result files from PCA run on synthetic data
#' dataDirRes <- system.file("extdata/demoAncestryCall/ex1", package="RAIDS")
#' listFilesName <- dir(file.path(dataDirRes), ".rds")
#' listFiles <- file.path(file.path(dataDirRes) , listFilesName)
#' syntheticKNN <- lapply(listFiles, FUN=function(x){return(readRDS(x))})
#' syntheticKNN <- do.call(rbind, syntheticKNN)
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA.Synthetic"
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' ## Open the Profile GDS file
#' gdsProfile <- snpgdsOpen(file.path(dataDir, "ex1.gds"))
#' \dontrun{
#'     pedSyn <- RAIDS:::prepPedSynthetic1KG(gdsReference=gdsRef,
#'               gdsSample=gdsProfile, studyID=studyID, popName="superPop")
#'
#'     ## Run the ancestry inference on one profile called 'ex1'
#'     ## The values of K and D used for the inference are selected using the
#'     ## synthetic results listFiles=listFiles,
#'     resCall <- RAIDS:::computeAncestryFromSynthetic(gdsReference=gdsRef,
#'                                 gdsProfile=gdsProfile,
#'                                 syntheticKNN = syntheticKNN,
#'                                 pedSyn = pedSyn,
#'                                 currentProfile=c("ex1"),
#'                                 spRef=demoKnownSuperPop1KG,
#'                                 studyIDSyn=studyID, np=1L)
#'
#'     ## The ancestry called with the optimal D and K values
#'     resCall$Ancestry
#' }
#' ## Close the GDS files (important)
#' closefn.gds(gdsProfile)
#' closefn.gds(gdsRef)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @export
computeAncestryFromSynthetic2 <- function(syntheticKNN,
                                         pedSyn,
                                         spRef,
                                         listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
                                         pRAIDS=pRAIDS) {
    # gdsReference, gdsProfile,
    # syntheticKNN,
    # pedSyn,
    # currentProfile,
    # spRef,
    # studyIDSyn,
    # np=1L,
    # listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
    # fieldPopInRef="superPop",
    # fieldPopInfAnc="SuperPop",
    # kList=seq(2, 15, 1),
    # pcaList=seq(2, 15, 1),
    # algorithm=c("exact", "randomized"),
    # eigenCount=32L,
    # missingRate=NaN, verbose=FALSE


    ## Validate input parameters
    # validateComputeAncestryFromSynthetic(gdsReference=gdsReference,
    #                                          gdsProfile=gdsProfile, syntheticKNN=syntheticKNN,
    #                                          pedSyn=pedSyn,
    #                                          currentProfile=currentProfile, spRef=spRef, studyIDSyn=studyIDSyn,
    #                                          np=np, listCatPop=listCatPop, fieldPopInRef=fieldPopInRef,
    #                                          fieldPopInfAnc=fieldPopInfAnc, kList=kList, pcaList=pcaList,
    #                                          algorithm=algorithm, eigenCount=eigenCount, missingRate=missingRate,
    #                                          verbose=verbose)

    ## Matches a character method against a table of candidate values
    ## aalready validate
    #algorithm <- arg_match(algorithm)



    ## Compile all the inferred ancestry results for different values of
    ## D and K to select the optimal parameters
    listParaSample <- selParaPCAUpQuartile(matKNN=syntheticKNN,
                                           pedCall=pedSyn, refCall=pRAIDS$fieldPopInRef, predCall=pRAIDS$fieldPopInfAnc,
                                           listCall=listCatPop,kList=pRAIDS$kList,
                                           pcaList=pRAIDS$pcaList[pRAIDS$pcaList <= pRAIDS$eigenCountSyn])
    listPCARef <- computePCARefRMMulti1(listRM = NULL, pRAIDS = pRAIDS)
    listPCAProfile <- computePCAProfile(listPCA=listPCARef, profileId=pRAIDS$pedStudy$Name.ID[1], pRAIDS=pRAIDS)
    ## Project profile on the PCA created with the reference profiles
    # listPCAProfile <- computePCARefSample(gdsProfile=gdsProfile,
    #                                       currentProfile=currentProfile, studyIDRef="Ref.1KG", np=np,
    #                                       algorithm=algorithm, eigenCount=eigenCount, missingRate=missingRate,
    #                                       verbose=verbose)

    ## Run a k-nearest neighbors analysis on one specific profile
    listKNNSample <- computeKNNRefSample(listEigenvector=listPCAProfile,
                                         listCatPop=listCatPop, spRef=spRef,
                                         fieldPopInfAnc=pRAIDS$fieldPopInfAnc,
                                         kList=pRAIDS$kList, pcaList=pRAIDS$pcaList)

    ## The ancestry call for the current profile
    resCall <- listKNNSample$matKNN[
        which(listKNNSample$matKNN$D == listParaSample$D &
                  listKNNSample$matKNN$K == listParaSample$K ),]
    colnames(listParaSample$dfAUROC) <- c("D", "K", "Call", "L", "AUROC", "H")
    res <- list(pcaSample=listPCAProfile, # PCA of the profile + 1KG
                paraSample=listParaSample, # Result of the parameter selection
                KNNSample=listKNNSample, # KNN for the profile
                Ancestry=resCall) # the ancestry call fo the profile

    return(res)
}

#' @title Extract the sample information from the 1KG GDS file for a list
#' of profiles associated to a specific study in the Profile GDS file
#'
#' @description The function extracts the information for the profiles
#' associated to a specific study in the GDS Sample file. The information is
#' extracted from the 'study.annot' node as a 'data.frame'.
#'
#' Then, the function used the 1KG GDS file to extract specific information
#' about each sample and add it, as an extra column, to the 'data.frame'.
#'
#' As example, this function can extract the synthetic profiles
#' for a GDS Sample and the super-population of the 1KG samples used to
#' generate each synthetic profile would be added
#' as an extra column to the final 'data.frame'.
#'
#' @param pRAIDS an object of class TODO.
#'
#' @return \code{data.frame} containing the columns extracted from the
#' GDS Sample 'study.annot' node with a extra column named as the 'popName'
#' parameter that has been extracted from the 1KG GDS 'sample.annot' node.
#' Only the rows corresponding to the specified study ('studyID' parameter)
#' are returned.
#'
#'
#' @details
#'
#' As example, this function can extract the synthetic profiles
#' for a Profile GDS and the super-population of the 1KG samples used to
#' generate each synthetic profile would be added
#' as an extra column to the final 'data.frame'. In that situation, the
#' 'popName' parameter would correspond to the super-population column and the
#' 'studyID' parameter would be the name given to the synthetic dataset.
#'
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## TODO
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen
#' @encoding UTF-8
#' @keywords internal
prepPedSyntheticRef <- function(pRAIDS) {
    # gdsReference=gdsReference,
    # gdsSample=gdsProfile,
    # studyID=pRAIDS$studyDFSyn$study.id,
    # popName="superPop"
    studyID <- pRAIDS$studyDFSyn$study.id
    popName <- pRAIDS$fieldPopInRef
    fileProfileGDS <-  validateProfileGDSExist(pathProfile=pRAIDS$pathProfileGDS,
                                    profile=pRAIDS$pedStudy$Name.ID[1])
    ## Open Profile GDS file
    gdsProfile <- snpgdsOpen(fileProfileGDS, readonly=TRUE)

    gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)

    ## Extract study information from the Profile GDS file
    studyAnnot <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))

    ## Retain the information associated to the current study
    studyCur <- studyAnnot[which(studyAnnot$study.id == studyID),]
    rm(studyAnnot)
    if(studyCur$case.id[1] == studyCur$data.id[1]){
        tmp <- matrix(unlist(strsplit(x = studyCur$case.id, split = "\\.")),
                      nrow=4)
        studyCur$case.id <- tmp[3,]
    }
    ## Get the information from 1KG GDS file
    dataRef <- read.gdsn(index.gdsn(node=gdsReference, "sample.annot"))

    if(! popName %in% colnames(dataRef)) {
        stop("The population ", popName, " is not supported.")
    }

    ## Assign sample names to the information
    row.names(dataRef) <- read.gdsn(index.gdsn(node=gdsReference, "sample.id"))

    studyCur[[popName]] <- dataRef[studyCur$case.id, popName]
    rownames(studyCur) <- studyCur$data.id
    
    ## Important to close the files
    closefn.gds(gdsProfile)
    closefn.gds(gdsReference)
    
    return(studyCur)
}
