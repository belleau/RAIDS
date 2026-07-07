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
    pRAIDS <- generateProfileRawGDS2(pRAIDS=pRAIDS)
    generateProfileGenoCall(pRAIDS=pRAIDS)
    # generateProfileGDS2(pRAIDS=pRAIDS)

    if(pRAIDS$verbose) {
        message("Genotype DONE ", Sys.time())
    }

    ## Close 1KG GDS file
    closefn.gds(gdsReference)

    ## Return successful code
    return(pRAIDS)
}


#' @title Calculate Principal Component Analysis (PCA) on SNV genotype data set
#'
#' @description The functions calculates the principal component analysis (PCA)
#' for a list of pruned SNVs present in a Profile GDS file. The
#' \link[SNPRelate]{snpgdsPCA} function is used to do the calculation.
#'
#' @param gdsProfile an object of class \link[SNPRelate]{SNPGDSFileClass},
#' the opened Profile GDS file.
#'
#' @param refProfileIDs a \code{vector} of reference 1KG profile identifiers
#' that are present in the Profile GDS file.
#' Those profiles minus the one present in the \code{listRM} vector will be
#' used to run the PCA analysis.
#'
#' @param listRM a \code{vector} of \code{character} strings containing the
#' identifiers for the reference samples that need to be removed for the
#' PCA analysis.
#'
#' @param np a single positive \code{integer} representing the number of CPU
#' that will be used. Default: \code{1L}.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA. The 2 choices are "exact" (traditional exact
#' calculation) and "randomized" (fast PCA with randomized algorithm
#' introduced in Galinsky et al. 2016). Default: \code{"exact"}.
#'
#' @param eigenCount a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigenCount' <= 0, then all eigenvectors are returned.
#' Default: \code{32L}.
#'
#' @param missingRate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA} function
#' with "<= missingRate" only; if \code{NaN}, no missing threshold.
#' Default: \code{0.025}.
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed.
#'
#' @return a \code{list} containing 2 entries:
#' \describe{
#' \item{pruned}{ a \code{vector} of SNV identifiers specifying selected SNVs
#' for the PCA analysis.}
#' \item{pca.unrel}{ a \code{snpgdsPCAClass} object containing the eigenvalues
#' as generated by \link[SNPRelate]{snpgdsPCA} function.}
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
#' ## Required library
#' library(SNPRelate)
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA.Synthetic"
#'
#' ## Profiles that should be removed from the PCA analysis
#' ## Those profiles has been used to generate the synthetic data set
#' samplesRM <- c("HG00246", "HG00325", "HG00611", "HG01173", "HG02165",
#'     "HG01112", "HG01615", "HG01968", "HG02658", "HG01850", "HG02013",
#'     "HG02465", "HG02974", "HG03814", "HG03445", "HG03689", "HG03789",
#'     "NA12751", "NA19107", "NA18548", "NA19075", "NA19475", "NA19712",
#'     "NA19731", "NA20528", "NA20908")
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoKNNSynthetic", package="RAIDS")
#'
#' ## Open the Profile GDS file
#' gdsProfile <- snpgdsOpen(file.path(dataDir, "ex1.gds"))
#'
#' ## Compute PCA for the 1KG reference profiles excluding
#' ## the profiles used to generate the synthetic profiles
#' results <- RAIDS:::computePCARefRMMulti(gdsProfile=gdsProfile,
#'     refProfileIDs=names(demoKnownSuperPop1KG), listRM=samplesRM, np=1L,
#'     algorithm="exact", eigenCount=32L, missingRate=0.025, verbose=FALSE)
#'
#' ## The PCA on the pruned SNVs data set for selected profiles
#' head(results$pca.unrel$eigenvect)
#'
#' ## Close Profile GDS file (important)
#' closefn.gds(gdsProfile)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCA
#' @encoding UTF-8
#' @keywords internal
computePCARefRMMulti1 <- function(listRM, pRAIDS) {
    # gdsProfile, refProfileIDs, listRM, np=1L,
    # algorithm="exact", eigenCount=32L,
    # missingRate=0.025, verbose
    # unrelatedSamples <- refProfileIDs[which(!(refProfileIDs %in% listRM))]

    fileGDSProfile <- file.path(pRAIDS$pathProfileGDS,
                                paste0(pRAIDS$pedStudy$Name.ID[1], ".gds"))
    gdsProfile <- openfn.gds(filename=fileGDSProfile)
    gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)
    unrelatedSamples <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[
            read.gdsn(index.gdsn(gdsReference, "sample.ref")) == 1
        ]
    if(! is.null(listRM)){
        unrelatedSamples <- unrelatedSamples[which(!(unrelatedSamples %in% listRM))]
    }
    listPCA <- list()

    listPCA[["pruned"]] <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))

    ## Calculate Principal Component Analysis (PCA) on SNV genotype dataset
    listPCA[["pca.unrel"]] <- snpgdsPCA(gdsobj=gdsReference,
                                        sample.id=unrelatedSamples,
                                        snp.id=listPCA[["pruned"]],
                                        num.thread=pRAIDS$np,
                                        missing.rate=pRAIDS$PCAmissingRate,
                                        algorithm=pRAIDS$PCAalgorithm,
                                        eigen.cnt=pRAIDS$eigenCount,
                                        verbose=pRAIDS$verbose)
    listPCA[["snp.load"]] <- snpgdsPCASNPLoading(listPCA[["pca.unrel"]],
                                                 gdsobj=gdsReference,
                                                 num.thread=pRAIDS$np,
                                                 verbose=pRAIDS$verbose)
    closefn.gds(gdsReference)
    closefn.gds(gdsProfile)
    return(listPCA)
}


# computePCAMulti <- function(pRAIDS) {
#     # gdsProfile, listPCA,
#     # sampleRef, studyIDSyn, verbose=FALSE
#     ## Validate the input parameters
#     validateComputePCAMultiSynthetic(gdsProfile=gdsProfile,
#                                      listPCA=listPCA, sampleRef=sampleRef, studyIDSyn=studyIDSyn,
#                                      verbose=verbose)
#
#     ## Identify profiles from synthetic data set
#     study.annot <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))
#     study.annot <- study.annot[which(study.annot$study.id == studyIDSyn &
#                                          study.annot$case.id %in% sampleRef),]
#
#     ## SNP loading in principal component analysis
#     listPCA[["snp.load"]] <- snpgdsPCASNPLoading(listPCA[["pca.unrel"]],
#                                                  gdsobj=gdsProfile, num.thread=1,
#                                                  verbose=verbose)
#
#     ## Project synthetic profiles onto existing principal component axes
#     listPCA[["samp.load"]] <- snpgdsPCASampLoading(listPCA[["snp.load"]],
#                                                    gdsobj=gdsProfile,
#                                                    sample.id=study.annot$data.id,
#                                                    num.thread=1L, verbose=verbose)
#
#     rownames(listPCA[["pca.unrel"]]$eigenvect) <-
#         listPCA[["pca.unrel"]]$sample.id
#
#     rownames(listPCA[["samp.load"]]$eigenvect) <-
#         listPCA[["samp.load"]]$sample.id
#
#     ## Return the eigenvectors for the 1KG reference profiles
#     ## and the eigenvectors for the synthetic data set projected on the 1KG PCA
#     listRes <- list(sample.id=listPCA[["samp.load"]]$sample.id,
#                     eigenvector.ref=listPCA[["pca.unrel"]]$eigenvect,
#                     eigenvector=listPCA[["samp.load"]]$eigenvect)
#
#     return(listRes)
# }

#' @title Calculate Principal Component Analysis (PCA) on SNV genotype data set
#'
#' @description The functions calculates the principal component analysis (PCA)
#' for a list of pruned SNVs present in a Profile GDS file. The
#' \link[SNPRelate]{snpgdsPCA} function is used to do the calculation.
#'
#' @param gdsProfile an object of class \link[SNPRelate]{SNPGDSFileClass},
#' the opened Profile GDS file.
#'
#' @param refProfileIDs a \code{vector} of reference 1KG profile identifiers
#' that are present in the Profile GDS file.
#' Those profiles minus the one present in the \code{listRM} vector will be
#' used to run the PCA analysis.
#'
#' @param listRM a \code{vector} of \code{character} strings containing the
#' identifiers for the reference samples that need to be removed for the
#' PCA analysis.
#'
#' @param np a single positive \code{integer} representing the number of CPU
#' that will be used. Default: \code{1L}.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA. The 2 choices are "exact" (traditional exact
#' calculation) and "randomized" (fast PCA with randomized algorithm
#' introduced in Galinsky et al. 2016). Default: \code{"exact"}.
#'
#' @param eigenCount a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigenCount' <= 0, then all eigenvectors are returned.
#' Default: \code{32L}.
#'
#' @param missingRate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA} function
#' with "<= missingRate" only; if \code{NaN}, no missing threshold.
#' Default: \code{0.025}.
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed.
#'
#' @return a \code{list} containing 2 entries:
#' \describe{
#' \item{pruned}{ a \code{vector} of SNV identifiers specifying selected SNVs
#' for the PCA analysis.}
#' \item{pca.unrel}{ a \code{snpgdsPCAClass} object containing the eigenvalues
#' as generated by \link[SNPRelate]{snpgdsPCA} function.}
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
#' ## Required library
#' library(SNPRelate)
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA.Synthetic"
#'
#' ## Profiles that should be removed from the PCA analysis
#' ## Those profiles has been used to generate the synthetic data set
#' samplesRM <- c("HG00246", "HG00325", "HG00611", "HG01173", "HG02165",
#'     "HG01112", "HG01615", "HG01968", "HG02658", "HG01850", "HG02013",
#'     "HG02465", "HG02974", "HG03814", "HG03445", "HG03689", "HG03789",
#'     "NA12751", "NA19107", "NA18548", "NA19075", "NA19475", "NA19712",
#'     "NA19731", "NA20528", "NA20908")
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoKNNSynthetic", package="RAIDS")
#'
#' ## Open the Profile GDS file
#' gdsProfile <- snpgdsOpen(file.path(dataDir, "ex1.gds"))
#'
#' ## Compute PCA for the 1KG reference profiles excluding
#' ## the profiles used to generate the synthetic profiles
#' results <- RAIDS:::computePCARefRMMulti(gdsProfile=gdsProfile,
#'     refProfileIDs=names(demoKnownSuperPop1KG), listRM=samplesRM, np=1L,
#'     algorithm="exact", eigenCount=32L, missingRate=0.025, verbose=FALSE)
#'
#' ## The PCA on the pruned SNVs data set for selected profiles
#' head(results$pca.unrel$eigenvect)
#'
#' ## Close Profile GDS file (important)
#' closefn.gds(gdsProfile)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCA
#' @encoding UTF-8
#' @keywords internal
computePCAProfile <- function(listPCA, profileId, pRAIDS) {

    fileGDSProfile <- file.path(pRAIDS$pathProfileGDS,
                                paste0(pRAIDS$pedStudy$Name.ID[1], ".gds"))
    gdsProfile <- openfn.gds(filename=fileGDSProfile)
    listPCA[["samp.load"]] <- snpgdsPCASampLoading(listPCA[["snp.load"]],
                                                   gdsobj=gdsProfile, sample.id=profileId,
                                                   num.thread=pRAIDS$np, verbose=pRAIDS$verbose)

    rownames(listPCA[["pca.unrel"]]$eigenvect) <-
        listPCA[["pca.unrel"]]$sample.id
    rownames(listPCA[["samp.load"]]$eigenvect) <-
        listPCA[["samp.load"]]$sample.id

    listRes <- list(sample.id=listPCA[["samp.load"]]$sample.id,
                    eigenvector.ref=listPCA[["pca.unrel"]]$eigenvect,
                    eigenvector=listPCA[["samp.load"]]$eigenvect)

    closefn.gds(gdsProfile)
    return(listRes)
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
#' @param K a TODO.
#' 
#' @param D a TODO.
#'
#' @param pRAIDS a TODO
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
#' @importFrom gdsfmt read.gdsn index.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen
#' @importFrom class knn
#' @encoding UTF-8
#' @export
computeKNNProfile <- function(listEigenvector, K, D, pRAIDS) {

    # gdsProfile, listEigenvector,
    # listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
    # studyIDSyn, spRef, fieldPopInfAnc="SuperPop",
    # kList=seq(2, 15, 1), pcaList=seq(2, 15, 1)
    ## Assign default value if kList is NULL

    gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)

    spRef <- read.gdsn(index.gdsn(gdsReference, "sample.annot/superPop"))
    names(spRef) <- read.gdsn(index.gdsn(gdsReference, "sample.id"))
    listCatPop <- unique(spRef)
    closefn.gds(gdsReference)

    resMat <- lapply(seq_len(length(listEigenvector$sample.id)),
                     FUN=function(i, listEigenvector, D, K){
            eigenvect <- rbind(listEigenvector$eigenvector.ref,
                               listEigenvector$eigenvector[i,,drop=FALSE])

            pcaND <- eigenvect[ ,seq_len(D)]
            yPred <-
                knn(train=pcaND[rownames(eigenvect)[-1*nrow(eigenvect)],],
                    test=pcaND[rownames(eigenvect)[nrow(eigenvect)],,
                               drop=FALSE],
                    cl=factor(spRef[rownames(eigenvect)[-1*nrow(eigenvect)]],
                              levels=listCatPop, labels=listCatPop),
                    k=K,
                    prob=FALSE)


            df <- data.frame(sample.id = listEigenvector$sample.id[i],
                             D=D,
                             K=K,
                             stringsAsFactors = FALSE)
            df[[pRAIDS$fieldPopInfAnc]] <- listCatPop[as.integer(yPred)]
            return(df)
        },
    listEigenvector=listEigenvector,
    D=D,
    K=K
    )
    resMat <- do.call(rbind, resMat)
    return(resMat)
}


#' @title Run a k-nearest neighbors analysis on a subset of the
#' synthetic dataset
#'
#' @description The function runs k-nearest neighbors analysis on a
#' subset of the synthetic data set. The function uses the 'knn' package.
#'
#' @param gdsProfile an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, the
#' opened Profile GDS file.
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
#' @param studyIDSyn a \code{character} string corresponding to the study
#' identifier.
#' The study identifier must be present in the Profile GDS file.
#'
#' @param spRef \code{vector} of \code{character} strings representing the
#' known super population ancestry for the 1KG profiles. The 1KG profile
#' identifiers are used as names for the \code{vector}.
#'
#' @param fieldPopInfAnc a \code{character} string representing the name of
#' the column that will contain the inferred ancestry for the specified
#' data set. Default: \code{"SuperPop"}.
#'
#' @param kList  a \code{vector} of \code{integer} representing  the list of
#' values tested for the  K parameter. The K parameter represents the
#' number of neighbors used in the K-nearest neighbors analysis. If
#' \code{NULL}, the value \code{seq(2, 15, 1)} is assigned.
#' Default: \code{seq(2, 15, 1)}.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  D parameter. The D parameter represents the
#' number of dimensions used in the PCA analysis.  If \code{NULL},
#' the value \code{seq(2, 15, 1)} is assigned.
#' Default: \code{seq(2, 15, 1)}.
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
computeKNNProfileSubSet <- function(listEigenvector,K, D, dfRef, pRAIDS) {

    # gdsProfile, listEigenvector,
    # listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
    # studyIDSyn, spRef, fieldPopInfAnc="SuperPop",
    # kList=seq(2, 15, 1), pcaList=seq(2, 15, 1)
    ## Assign default value if kList is NULL

    # gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)

    # spRef <- read.gdsn(index.gdsn(gdsReference, "sample.annot/superPop"))
    # names(spRef) <- read.gdsn(index.gdsn(gdsReference, "sample.id"))
    listCatPop <- unique(dfRef$superPop)
    # closefn.gds(gdsReference)
    listEigenvector$eigenvector.ref <- listEigenvector$eigenvector.ref[row.names(listEigenvector$eigenvector.ref) %in% row.names(dfRef),]

    resMat <- lapply(seq_len(length(listEigenvector$sample.id)),
                     FUN=function(i, listEigenvector, D, K, dfRef){
                         eigenvect <- rbind(listEigenvector$eigenvector.ref,
                                            listEigenvector$eigenvector[i,,drop=FALSE])

                         pcaND <- eigenvect[ ,seq_len(D)]
                         yPred <-
                             knn(train=pcaND[rownames(eigenvect)[-1*nrow(eigenvect)],],
                                 test=pcaND[rownames(eigenvect)[nrow(eigenvect)],,
                                            drop=FALSE],
                                 cl=factor(dfRef[rownames(eigenvect)[-1*nrow(eigenvect)], "superPop"],
                                           levels=listCatPop, labels=listCatPop),
                                 k=K,
                                 prob=FALSE)


                         df <- data.frame(sample.id = listEigenvector$sample.id[i],
                                          D=D,
                                          K=K,
                                          stringsAsFactors = FALSE)
                         df[[pRAIDS$fieldPopInfAnc]] <- listCatPop[as.integer(yPred)]
                         return(df)
                     },
                     listEigenvector=listEigenvector,
                     D=D,
                     K=K,
                     dfRef=dfRef
    )
    resMat <- do.call(rbind, resMat)
    return(resMat)
}
