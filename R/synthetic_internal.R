#' @title Validate input parameters for syntheticGeno() function
#'
#' @description This function validates the input parameters for the
#' [syntheticGeno()] function.
#'
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the 1KG GDS file.
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the1 1KG SNV Annotation GDS file.
#'
#' @param fileProfileGDS a \code{character} string representing the file name of
#' the GDS Sample file containing the information about the sample.
#' The file must exist.
#'
#' @param listSampleRef a \code{vector} of \code{character} strings
#' representing the sample identifiers of the 1KG selected reference samples.
#'
#' @param profileID a \code{character} string representing the unique
#' identifier of the cancer sample.
#'
#' @param nbSim a single positive \code{integer} representing the number of
#' simulations that will be generated per sample + 1KG reference combination.
#'
#' @param prefix a \code{character} string that represent the prefix that will
#' be added to the name of the synthetic profiles generated by the function.
#'
#' @param pRecomb a single positive \code{numeric} between 0 and 1 that
#' represents the frequency of phase switching in the synthetic profiles.
#'
#' @param minProb a single positive \code{numeric} between 0 and 1 that
#' represents the probability that the genotype is correct.
#'
#' @param seqError a single positive \code{numeric} between 0 and 1
#' representing the sequencing error rate.
#'
#' @return The integer \code{0L} when the function is successful.
#'
#' @examples
#'
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gdsRef <- openfn.gds(file.path(dataDir,
#'                 "PopulationReferenceDemo.gds"), readonly=TRUE)
#'
#' ## The 1KG GDS Annotation file (opened)
#' gdsRefAnnot <- openfn.gds(file.path(dataDir,
#'     "PopulationReferenceSNVAnnotationDemo.gds"), readonly=TRUE)
#'
#' ## The GDS Sample file
#' gdsSample <- file.path(dataDir, "GDS_Sample_with_study_demo.gds")
#'
#' ## The validation should be successful
#' RAIDS:::validateSyntheticGeno(gdsReference=gdsRef, gdsRefAnnot=gdsRefAnnot,
#'      fileProfileGDS=gdsSample, profileID="A101TCGA",
#'      listSampleRef="A101TCGA", nbSim=1L, prefix="TCGA", pRecomb=0.02,
#'      minProb=0.999, seqError=0.002)
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gdsRef)
#' closefn.gds(gdsfile=gdsRefAnnot)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateSyntheticGeno <- function(gdsReference, gdsRefAnnot, fileProfileGDS,
                                        profileID, listSampleRef,
                                        nbSim, prefix, pRecomb, minProb,
                                        seqError) {

    ## The gds must be an object of class "gds.class"
    if (!inherits(gdsReference, "gds.class")) {
        stop("The \'gdsReference\' must be an object of class \'gds.class\'.")
    }

    ## The gdsRefAnnot must be an object of class "gds.class"
    if (!inherits(gdsRefAnnot, "gds.class")) {
        stop("The \'gdsRefAnnot\' must be an object of class \'gds.class\'.")
    }

    ## The fileProfileGDS must be an character string and the file must exist
    if (!(is.character(fileProfileGDS) && file.exists(fileProfileGDS))) {
        stop("The \'fileProfileGDS\' must be a character string and the file ",
                "must exist.")
    }

    ## The profileID must be a character string
    if (!(is.character(profileID) && length(profileID) == 1)) {
        stop("The \'profileID\' must be a character string.")
    }

    ## The listSampleRef must be a character string
    if(!is.character(listSampleRef)) {
        stop("The \'listSampleRef\' must be a vector of character strings.")
    }

    ## The parameter nbSim must be a single positive integer
    if(!(isSingleNumber(nbSim) && (nbSim >= 0))) {
        stop("The \'nbSim\' parameter must be a single positive ",
                "numeric value.")
    }

    ## The parameter prefix must be a single character string
    if(!(is.character(prefix) && (length(prefix) == 1))) {
        stop("The \'prefix\' parameter must be a single character ",
                "string.")
    }

    ## The parameter pRecomb must be a single positive integer
    if(!(isSingleNumber(pRecomb) && (pRecomb >= 0.0) && (pRecomb <= 1.0))) {
        stop("The \'pRecomb\' parameter must be a single positive ",
                "numeric value between 0 and 1.")
    }

    ## The parameter minProb must be a single positive integer
    if(!(isSingleNumber(minProb) && (minProb >= 0.0) && (minProb <= 1.0))) {
        stop("The \'minProb\' parameter must be a single positive ",
                "numeric value between 0 and 1.")
    }

    ## The parameter seqError must be a single positive integer
    if(!(isSingleNumber(seqError) && (seqError >= 0.0) && (seqError <= 1.0))) {
        stop("The \'seqError\' parameter must be a single positive ",
                "numeric value between 0 and 1.")
    }

    return(0L)
}


#' @title Validate input parameters for prepSynthetic() function
#'
#' @description This function validates the input parameters for the
#' [prepSynthetic()] function.
#'
#' @param fileProfileGDS a \code{character} string representing the file name
#' of the GDS Sample file containing the information about the sample
#' used to generate the synthetic profiles.
#'
#' @param listSampleRef a \code{vector} of \code{character} string
#' representing the
#' identifiers of the selected 1KG samples that will be used as reference to
#' generate the synthetic profiles.
#'
#' @param profileID a \code{character} string representing the profile
#' identifier present in the \code{fileProfileGDS} that will be used to
#' generate synthetic profiles.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 2 columns: "study.id" and "study.desc". Those 2 columns
#' must be in \code{character} strings (no factor). Other columns can be
#' present, such as "study.platform", but won't be used.
#'
#' @param nbSim a single positive \code{integer} representing the number of
#' simulations per combination of sample and 1KG reference.
#'
#' @param prefix a single \code{character} string representing the prefix that
#' is going to be added to the name of the synthetic profile. The prefix
#' enables the creation of multiple synthetic profile using the same
#' combination of sample and 1KG reference.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function.
#'
#' @return \code{0L} when successful.
#'
#' @examples
#'
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## The Profile GDS Sample
#' gdsSample <- file.path(dataDir, "GDS_Sample_with_study_demo.gds")
#'
#' ## The study data frame
#' studyDF <- data.frame(study.id="MYDATA.Synthetic",
#'     study.desc="MYDATA synthetic data", study.platform="PLATFORM",
#'     stringsAsFactors=FALSE)
#'
#' ## The validation should be successful
#' RAIDS:::validatePepSynthetic(fileProfileGDS=gdsSample,
#'      listSampleRef=c("Sample01", "Sample02"), profileID="A101TCGA",
#'      studyDF=studyDF, nbSim=1L, prefix="TCGA", verbose=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validatePepSynthetic <- function(fileProfileGDS,
                listSampleRef, profileID, studyDF, nbSim, prefix, verbose) {

    ## The fileProfileGDS must be a character string and the file must exists
    if (!(is.character(fileProfileGDS) && (file.exists(fileProfileGDS)))) {
        stop("The \'fileProfileGDS\' must be a character string representing ",
                "the GDS Sample information file. The file must exist.")
    }

    ## The listSampleRef must be character string
    if (!is.character(listSampleRef)) {
        stop("The \'listSampleRef\' must be a vector of character strings.")
    }

    ## The profileID must be a single character String
    if (!(is.character(profileID) && length(profileID) == 1)) {
        stop("The \'profileID\' must be a single character string.")
    }

    ## The study.id must have the 2 mandatory columns
    if(sum(c("study.id", "study.desc") %in% colnames(studyDF)) != 2 ) {
        stop("The \'studyDF\' data frame is incomplete. ",
                "One or more mandatory column is missing.\n")
    }

    ## The nbSim must be a single positive numeric
    if (!(isSingleNumber(nbSim) && nbSim > 0)) {
        stop("The \'nbSim\' must be a single positive integer.")
    }

    ## The prefix must be a single character String
    if (!(is.character(prefix) && length(prefix) == 1)) {
        stop("The \'prefix\' must be a single character string.")
    }

    ## The verbose must be a logical
    validateLogical(logical=verbose, "verbose")

    return(0L)
}


#' @title Validate input parameters for computeSyntheticROC() function
#'
#' @description This function validates the input parameters for the
#' [computeSyntheticROC()] function.
#'
#' @param matKNN a \code{data.frame} containing the inferred ancestry results
#' for fixed values of _D_ and _K_. On of the column names of the
#' \code{data.frame} must correspond to the \code{matKNNAncestryColumn}
#' argument.
#'
#' @param matKNNAncestryColumn  a \code{character} string
#' representing the
#' name of the column that contains the inferred ancestry for the specified
#' synthetic profiles. The column must be present in the \code{matKNN}
#' argument.
#'
#' @param pedCall a \code{data.frame} containing the information about
#' the super-population information from the 1KG GDS file
#' for profiles used to generate the synthetic profiles. The \code{data.frame}
#' must contained a column named as the \code{pedCallAncestryColumn} argument.
#'
#' @param pedCallAncestryColumn a \code{character} string representing the
#' name of the column that contains the known ancestry for the reference
#' profiles in the Reference GDS file. The column must be present in
#' the \code{pedCall} argument.
#'
#' @param listCall a \code{vector} of \code{character} strings representing
#' the list of all possible ancestry assignations.
#'
#' @return \code{0L} when successful.
#'
#' @examples
#'
#' ## Loading demo dataset containing pedigree information for synthetic
#' ## profiles and known ancestry of the profiles used to generate the
#' ## synthetic profiles
#' data(pedSynthetic)
#'
#' ## Loading demo dataset containing the inferred ancestry results
#' ## for the synthetic data
#' data(matKNNSynthetic)
#'
#' ## The inferred ancestry results for the synthetic data using
#' ## values of D=6 and K=5
#' matKNN <- matKNNSynthetic[matKNNSynthetic$K == 6 & matKNNSynthetic$D == 5, ]
#'
#' ## The validation should be successful
#' RAIDS:::validateComputeSyntheticRoc(matKNN=matKNN,
#'     matKNNAncestryColumn="SuperPop",
#'     pedCall=pedSynthetic, pedCallAncestryColumn="superPop",
#'     listCall=c("EAS", "EUR", "AFR", "AMR", "SAS"))
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateComputeSyntheticRoc <- function(matKNN, matKNNAncestryColumn, pedCall,
                                    pedCallAncestryColumn, listCall) {

    ## The matKNN must be a data.frame
    if(!is.data.frame(matKNN)) {
        stop("The \'matKNN\' must be a data frame.")
    }

    ## The matKNNAncestryColumn must be a single character String
    if (!(is.character(matKNNAncestryColumn) &&
                length(matKNNAncestryColumn) == 1)) {
        stop("The \'matKNNAncestryColumn\' must be a single character string.")
    }

    ## The matKNNAncestryColumn must be a column in the matKNN data frame
    if (!(matKNNAncestryColumn %in% colnames(matKNN))) {
        stop("The \'matKNNAncestryColumn\' must be a column in the \'matKNN\'",
                        " data frame.")
    }

    ## The pedCall must be a data.frame
    if(!is.data.frame(pedCall)) {
        stop("The \'pedCall\' must be a data frame.")
    }

    ## The pedCallAncestryColumn must be a single character String
    if (!(is.character(pedCallAncestryColumn) &&
                length(pedCallAncestryColumn) == 1)) {
        stop("The \'pedCallAncestryColumn\' must be a single character string.")
    }

    ## The pedCallAncestryColumn must be a column in the pedCall data frame
    if (!(pedCallAncestryColumn %in% colnames(pedCall))) {
        stop("The \'pedCallAncestryColumn\' must be a column in the ",
                "\'pedCall\' data frame.")
    }

    ## The listCall must be character string
    if (!is.character(listCall)) {
        stop("The \'listCall\' must be a vector of character strings.")
    }

    if(length(unique(matKNN$D)) != 1 | length(unique(matKNN$K)) != 1) {
        stop("The synthetic accuracy can only be caculated for one fixed value",
            " of D and K. The 2 data frames must be filterd to retain only",
            " one value.")
    }

    return(0L)
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
#' @param gdsReference an object of class
#' \code{\link[gdsfmt:gds.class]{gdsfmt::gds.class}}, the opened 1 KG GDS file.
#'
#' @param gdsSample an object of class
#' \code{\link[gdsfmt:gds.class]{gdsfmt::gds.class}}, the opened Profile GDS
#' file.
#'
#' @param studyID a \code{character} string representing the name of the
#' study that will be extracted from the GDS Sample 'study.annot' node.
#'
#' @param popName a \code{character} string representing the name of the
#' column from the \code{data.frame} stored in the 'sample.annot' node of the
#' 1KG GDS file. The column must be present in the \code{data.frame}.
#'
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
#' ## The open 1KG GDS file is required (this is a demo file)
#' dataDir <- system.file("extdata", package="RAIDS")
#' gds_1KG_file <- file.path(dataDir, "PopulationReferenceDemo.gds")
#' gds1KG <- openfn.gds(gds_1KG_file)
#'
#' fileSampleGDS <- file.path(dataDir, "GDS_Sample_with_study_demo.gds")
#' gdsSample <- openfn.gds(fileSampleGDS)
#'
#' ## Extract the study information for "TCGA.Synthetic" study present in the
#' ## Profile GDS file and merge column "superPop" from 1KG GDS to the
#' ## returned data.frame
#' ## This function enables to extract the super-population associated to the
#' ## 1KG samples that has been used to create the synthetic profiles
#' RAIDS:::prepPedSynthetic1KG(gdsReference=gds1KG, gdsSample=gdsSample,
#'     studyID="TCGA.Synthetic", popName="superPop")
#'
#' ## The GDS files must be closed
#' gdsfmt::closefn.gds(gds1KG)
#' gdsfmt::closefn.gds(gdsSample)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @keywords internal
prepPedSynthetic1KG <- function(gdsReference, gdsSample, studyID, popName) {

    ## Extract study information from the Profile GDS file
    studyAnnot <- read.gdsn(index.gdsn(gdsSample, "study.annot"))

    ## Retain the information associated to the current study
    studyCur <- studyAnnot[which(studyAnnot$study.id == studyID),]
    rm(studyAnnot)

    ## Get the information from 1KG GDS file
    dataRef <- read.gdsn(index.gdsn(node=gdsReference, "sample.annot"))

    if(! popName %in% colnames(dataRef)) {
        stop("The population ", popName, " is not supported.")
    }

    ## Assign sample names to the information
    row.names(dataRef) <- read.gdsn(index.gdsn(node=gdsReference, "sample.id"))

    studyCur[[popName]] <- dataRef[studyCur$case.id, popName]
    rownames(studyCur) <- studyCur$data.id

    return(studyCur)
}


#' @title Calculate the confusion matrix of the inferences for specific
#' values of D and K using the inferred ancestry results from the synthetic
#' profiles.
#'
#' @description The function calculates the confusion matrix of the inferences
#' for fixed values of _D_ and _K_ using the inferred ancestry results done
#' on the synthetic profiles.
#'
#' @param matKNN a \code{data.frame} containing the inferred ancestry results
#' for fixed values of _D_ and _K_. The \code{data.frame} must contained
#' those columns: "sample.id", "D", "K" and the fourth column name must
#' correspond to the \code{matKNNAncestryColumn} argument.
#'
#' @param matKNNAncestryColumn a \code{character} string representing the
#' name of the column that contains the inferred ancestry for the specified
#' synthetic profiles. The column must be present in the \code{matKNN}
#' argument.
#'
#' @param pedCall a \code{data.frame} containing the information about
#' the super-population information from the 1KG GDS file
#' for profiles used to generate the synthetic profiles. The \code{data.frame}
#' must contained a column named as the \code{pedCallAncestryColumn} argument.
#'
#' @param pedCallAncestryColumn a \code{character} string representing the
#' name of the column that contains the known ancestry for the reference
#' profiles in the Reference GDS file. The column must be present in
#' the \code{pedCall} argument.
#'
#' @param listCall a \code{vector} of \code{character} strings representing
#' the list of possible ancestry assignations.
#'
#' @return \code{list} containing 2 entries:
#' \describe{
#' \item{confMat}{ a \code{matrix} representing the confusion matrix }
#' \item{matAccuracy}{ a \code{data.frame} containing the statistics
#' associated to the confusion matrix}
#' }
#'
#' @examples
#'
#' ## Loading demo dataset containing pedigree information for synthetic
#' ## profiles and known ancestry of the profiles used to generate the
#' ## synthetic profiles
#' data(pedSynthetic)
#'
#' ## Loading demo dataset containing the inferred ancestry results
#' ## for the synthetic data
#' data(matKNNSynthetic)
#'
#' ## The inferred ancestry results for the synthetic data using
#' ## values of D=6 and K=5
#' matKNN <- matKNNSynthetic[matKNNSynthetic$K == 6 & matKNNSynthetic$D == 5, ]
#'
#' ## Compile the confusion matrix using the
#' ## synthetic profiles for fixed values of  D and K values
#' results <- RAIDS:::computeSyntheticConfMat(matKNN=matKNN,
#'     matKNNAncestryColumn="SuperPop",
#'     pedCall=pedSynthetic, pedCallAncestryColumn="superPop",
#'     listCall=c("EAS", "EUR", "AFR", "AMR", "SAS"))
#'
#' results$confMat
#' results$matAccuracy
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @encoding UTF-8
#' @keywords internal
computeSyntheticConfMat <- function(matKNN, matKNNAncestryColumn,
                                pedCall, pedCallAncestryColumn, listCall) {

    matAccuracy <- data.frame(pcaD=matKNN$D[1], K=matKNN$K[1],
                    Accu.CM=numeric(1), CM.CI=numeric(1), N=nrow(matKNN),
                    NBNA=length(which(is.na(matKNN[[matKNNAncestryColumn]]))))
    i <- 1
    if(length(unique(matKNN$D)) != 1 | length(unique(matKNN$K)) != 1){
        stop("Compute synthetic accuracy with different pca dimension or K\n")
    }

    listKeep <- which(!(is.na(matKNN[[matKNNAncestryColumn]])) )

    fCall <- factor(pedCall[matKNN$sample.id[listKeep], pedCallAncestryColumn],
                            levels=listCall, labels=listCall)

    fP <- factor(matKNN[[matKNNAncestryColumn]][listKeep],
                            levels = listCall, labels = listCall)

    cm <- table(fCall, fP)


    matAccuracy[i, 3] <- sum(diag(cm[rownames(cm) %in% listCall,
                                        colnames(cm) %in% listCall])) /
        nrow(pedCall[matKNN$sample.id, ][listKeep,])

    matAccuracy[i, 4] <- 1.96 * (matAccuracy[i, 3] * (1 - matAccuracy[i, 3]) /
                            nrow(pedCall[matKNN$sample.id, ][listKeep,]))^0.5

    ## Generate list that will be returned
    res <- list(confMat=cm, matAccuracy=matAccuracy)

    return(res)
}
