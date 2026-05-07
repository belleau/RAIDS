#' @title Random selection of a specific number of reference profiles in each
#' subcontinental population present in the 1KG GDS file ( same as select1KGPop
#' but the function doesn't need gds object as parameters but the file name
#' of the referenceGDS )
#'
#' @description The function randomly selects a fixed number of reference
#' for each subcontinental population present in the 1KG GDS file. When a
#' subcontinental population has less samples than the fixed number, all
#' samples from the subcontinental population are selected.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @param nbProfiles a single positive \code{integer} representing the number
#' of samples that will be selected for each subcontinental population present
#' in the 1KG GDS file. If the number of samples in a specific subcontinental
#' population is smaller than the \code{nbProfiles}, the number of samples
#' selected in this
#' subcontinental population will correspond to the size of this population.
#'
#' @return a \code{data.frame} containing those columns:
#' \describe{
#' \item{sample.id}{ a \code{character} string representing the sample
#' identifier. }
#' \item{pop.group}{ a \code{character} string representing the
#' subcontinental population assigned to the sample. }
#' \item{superPop}{ a \code{character} string representing the
#' super-population assigned to the sample. }
#' }
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## The number of samples needed by subcontinental population
#' ## The number is small for demonstration purpose
#' nbProfiles <- 5L
#' ## TODO
#' ## 1KG GDS Demo file
#' ## This file only one superpopulation (for demonstration purpose)
#' dataDir <- system.file("extdata", package="RAIDS")
#' fileGDS <- file.path(dataDir, "PopulationReferenceDemo.gds")
#'
#' ## Extract a selected number of random samples
#' ## for each subcontinental population
#' ## In the 1KG GDS Demo file, there is one subcontinental population
#' dataR <- select1KGPopForSynthetic(fileReferenceGDS=fileGDS, nbProfiles=nbProfiles)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @importFrom SNPRelate snpgdsOpen
#' @encoding UTF-8
#' @export
selectHGDP1kgForSynthetic <- function(pRAIDS, nbProfiles=30L) {

    ## The fileReferenceGDS must be a character string and the file must exists
    if (!(is.character(pRAIDS$fileReferenceGDS) && (file.exists(pRAIDS$fileReferenceGDS)))) {
        stop("The \'fileReferenceGDS\' must be a character string ",
             "representing the Reference GDS file. The file must exist.")
    }
    ## Validate that nbProfiles parameter is a single positive numeric
    if(! (isSingleNumber(nbProfiles) && nbProfiles > 0)) {
        stop("The \'nbProfiles\' parameter must be a single positive integer.")
    }

    gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)
    ## Select reference samples
    listRef <- read.gdsn(index.gdsn(gdsReference, "sample.ref"))
    listKeep <- which(listRef == 1)
    rm(listRef)

    # Extract information about the selected reference samples
    # Including all the subcontinental population classes represented
    sample.annot <- read.gdsn(index.gdsn(gdsReference,
                        "sample.annot"))[listKeep,]
    sample.id <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[listKeep]
    grAna <- list(g1 <- c("AFRag"),
                g2 <- c("EASag"),
                g3 <- c("EURag"),
                g4 <- c("SASag"),
                g5 <- c("AMRag", "MIDag", "OTHag"),
                g6 <- c("ADM"))




    ## For each subcontinental population, randomly select a fixed number of
    ## samples
    dfAll <- lapply(seq_len(length(grAna)),
                    FUN=function(i, grAna) {
        listGroup <- which(sample.annot$superPop %in% grAna[[i]])
        if(5 * nbProfiles > length(listGroup)){
            exit("nbProfiles for synthetic is too big\n")
        }
        tmp <- sample(listGroup, 5 * nbProfiles)
        return(data.frame(sample.id=sample.id[tmp],
                          popGr = rep(paste(grAna[[i]], collapse = "."), 5 * nbProfiles),
                          superPop=sample.annot$superPop[tmp],
                          stringsAsFactors=FALSE)) },
        grAna=grAna)

    df <- do.call(rbind, dfAll)


    closefn.gds(gdsReference)

    return(df)
}

#' @title Generate synthetic profiles for each cancer profile and 1KG
#' reference profile combination and add them to the Profile GDS file
#'
#' @description The functions uses one cancer profile in combination with one
#' 1KG reference profile to generate an synthetic profile that is saved in the
#' Profile GDS file.
#'
#' When more than one 1KG reference profiles are specified,
#' the function recursively generates synthetic profiles for
#' each cancer profile + 1KG reference profile combination.
#'
#' The number of
#' synthetic profiles generated by combination is specified by the number of
#' simulation requested.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return The integer \code{OL} when the function is successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#'
#' ## Profile GDS file (temporary)
#' fileNameGDS <- file.path(tempdir(), "ex1.gds")
#'
#' ## Copy the Profile GDS file demo that has been pruned and annotated
#' file.copy(file.path(dataDir, "ex1_demo_with_pruning_and_1KG_annot.gds"),
#'                  fileNameGDS)
#'
#' ## Information about the synthetic data set
#' syntheticStudyDF <- data.frame(study.id="MYDATA.Synthetic",
#'         study.desc="MYDATA synthetic data", study.platform="PLATFORM",
#'         stringsAsFactors=FALSE)
#'
#' ## Add information related to the synthetic profiles into the Profile GDS
#' prepSynthetic(fileProfileGDS=fileNameGDS,
#'         listSampleRef=c("HG00243", "HG00150"), profileID="ex1",
#'         studyDF=syntheticStudyDF, nbSim=1L, prefix="synthTest",
#'         verbose=FALSE)
#'
#' ## The 1KG files
#' gds1KG <- snpgdsOpen(file.path(dataDir,
#'                             "ex1_good_small_1KG.gds"))
#' gds1KGAnnot <- openfn.gds(file.path(dataDir,
#'                             "ex1_good_small_1KG_Annot.gds"))
#'
#' ## Generate the synthetic profiles and add them into the Profile GDS
#' syntheticGeno(gdsReference=gds1KG, gdsRefAnnot=gds1KGAnnot,
#'         fileProfileGDS=fileNameGDS, profileID="ex1",
#'         listSampleRef=c("HG00243", "HG00150"), nbSim=1,
#'         prefix="synthTest",
#'         pRecomb=0.01, minProb=0.999, seqError=0.001)
#'
#' ## Open Profile GDS file
#' profileGDS <- openfn.gds(fileNameGDS)
#'
#' tail(read.gdsn(index.gdsn(profileGDS, "sample.id")))
#'
#' ## Close GDS files (important)
#' closefn.gds(profileGDS)
#' closefn.gds(gds1KG)
#' closefn.gds(gds1KGAnnot)
#'
#' ## Remove Profile GDS file (created for demo purpose)
#' unlink(fileNameGDS, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats rmultinom
#' @encoding UTF-8
#' @export
syntheticGeno2 <- function(pRAIDS) {

    gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)
    gdsRefAnnot <- openfn.gds(pRAIDS$fileReferenceAnnotGDS)
    fileProfileGDS <-  validateProfileGDSExist(pathProfile=pRAIDS$pathProfileGDS,
                                               profile=pRAIDS$pedStudy$Name.ID[1])
    listSampleRef <- pRAIDS$syntheticRefDF$sample.id

    ## Validate the input parameters
    validateSyntheticGeno(gdsReference=gdsReference, gdsRefAnnot=gdsRefAnnot,
                          fileProfileGDS=fileProfileGDS, profileID=pRAIDS$pedStudy$Name.ID[1],
                          listSampleRef=listSampleRef, nbSim=pRAIDS$nbSim, prefix=pRAIDS$prefix,
                          pRecomb=pRAIDS$pRecomb, minProb=pRAIDS$minProb, seqError=pRAIDS$seqErrorSyn)

    ## Open the GDS Sample file
    #gdsSample <- openfn.gds(filename=fileProfileGDS, readonly=FALSE)
    gdsProfile <- snpgdsOpen(fileProfileGDS, readonly=FALSE)

    ## The name of the simulated profiles
    sampleSim <- paste(paste0(pRAIDS$prefix, ".", pRAIDS$pedStudy$Name.ID[1]),
                       paste(rep(listSampleRef,each=pRAIDS$nbSim),
                             seq_len(pRAIDS$nbSim), sep="."), sep = ".")

    sample.id <- read.gdsn(index.gdsn(gdsProfile, "sample.id"))

    if (length(which(sampleSim %in% sample.id)) > 0) {
        closefn.gds(gdsProfile)
        closefn.gds(gdsReference)
        closefn.gds(gdsRefAnnot)
        stop("Error data.id of the simulation exists change prefix\n")
    }

    ## Find information about the 1KG reference samples used to generate
    ## the simulated profiles
    sample.1kg <- read.gdsn(index.gdsn(gdsReference, "sample.id"))
    # listPosRef <- which(sample.id %in% listSampleRef)
    listPosRef.1kg <- which(sample.1kg %in% listSampleRef)

    superPop <- read.gdsn(index.gdsn(gdsReference,
                                     "sample.annot/superPop"))[listPosRef.1kg]

    # if (! all.equal(sample.id[listPosRef], sample.1kg[listPosRef.1kg])) {
    #     stop("The order between 1kg and the list of samples is not the same.\n")
    # }

    ## Get indexes of the SNV associated to the sample from the GDS Sample file
    list1KG <- read.gdsn(index.gdsn(gdsProfile, "snp.index"))

    ## Create a table with the coverage and low allelic fraction information
    infoSNV <- data.frame(count.tot=read.gdsn(index.gdsn(gdsProfile,
                                                         "Total.count"))[list1KG],
                          lap=read.gdsn(index.gdsn(gdsProfile, "lap")))

    nbSNV <- nrow(infoSNV)

    ## Define a table for each "count.tot","lap" and Freq (number of occurrence)
    ## to reduce the number of sampling calls later
    df <- as.data.frame(table(infoSNV[,c("count.tot","lap")]))
    df <- df[df$Freq > 0, ]
    df <- df[order(df$count.tot, df$lap), ]
    # order of SNV relatively to df
    listOrderSNP <- order(infoSNV$count.tot, infoSNV$lap)

    # pos in listOrderSNP of each entries of df
    # df[i, ] where i < nrow(df) link to
    # infoSNV[listOrderSNP[(posDF[i]+1):(posDF[i+1])],]
    posDF <- c(0,cumsum(df$Freq))

    block.Annot <- read.gdsn(index.gdsn(gdsRefAnnot, "block.annot"))

    posSP <- data.frame(EAS=which(block.Annot$block.id == "EAS.0.05.500k"),
                        EUR=which(block.Annot$block.id == "EUR.0.05.500k"),
                        AFR=which(block.Annot$block.id == "AFR.0.05.500k"),
                        AMR=which(block.Annot$block.id == "AMR.0.05.500k"),
                        SAS=which(block.Annot$block.id == "SAS.0.05.500k"))

    blockDF <- data.frame(EAS=read.gdsn(index.gdsn(gdsRefAnnot, "block"),
                                        start=c(1,posSP$EAS), count = c(-1,1))[list1KG],
                          EUR=read.gdsn(index.gdsn(gdsRefAnnot, "block"),
                                        start=c(1,posSP$EUR), count = c(-1,1))[list1KG],
                          AFR=read.gdsn(index.gdsn(gdsRefAnnot, "block"),
                                        start=c(1,posSP$AFR), count = c(-1,1))[list1KG],
                          AMR=read.gdsn(index.gdsn(gdsRefAnnot, "block"),
                                        start=c(1,posSP$AMR), count = c(-1,1))[list1KG],
                          SAS=read.gdsn(index.gdsn(gdsRefAnnot, "block"),
                                        start=c(1,posSP$SAS), count = c(-1,1))[list1KG])

    # For each reference simulate
    for(r in seq_len(length(listPosRef.1kg))) {

        curSynt <- listPosRef.1kg[r]
        #r.1kg <- which(sample.id[listPosRef[r]] == sample.1kg)
        # get the genotype of the sample r
        g <- read.gdsn(index.gdsn(gdsReference, "genotype"),
                       start=c(1, curSynt), count=c(-1, 1))[list1KG]

        # Order the SNV by count.tot and, lap (low allelic proportion)
        gOrder <- g[listOrderSNP]

        matSim1 <- matrix(nrow=sum(df$Freq), ncol=pRAIDS$nbSim)
        matSim2 <- matrix(nrow=sum(df$Freq), ncol=pRAIDS$nbSim)

        # Loop on the read.count and lap
        # Faster to group the read.count and lap to run rmultinom
        for(i in seq_len(nrow(df))){

            lap <- as.numeric(as.character(df$lap[i]))

            ## Number of SNV heterozygote corresponding to
            ## df$count.tot[i] and df$lap[i]
            hetero <- which(gOrder[(posDF[i]+1):(posDF[i+1])] == 1)
            nbHetero <- length(hetero)
            # Define the tree prob for the muultinomial
            p1 <- lap * (1- 3 *pRAIDS$seqErrorSyn) + (1 - lap) * pRAIDS$seqErrorSyn

            p2 <- (1 - lap) * (1- 3 *pRAIDS$seqErrorSyn) + lap * pRAIDS$seqErrorSyn

            p3 <- 2 * pRAIDS$seqErrorSyn

            tmp <- rmultinom(nbHetero * pRAIDS$nbSim,
                             as.numeric(as.character(df$count.tot[i])),  c(p1, p2, p3))
            # depth of allele 1
            matSim1[listOrderSNP[hetero + posDF[i]],] <- matrix(tmp[1,],
                                                                ncol=pRAIDS$nbSim)
            # depth of allele 2
            matSim2[listOrderSNP[hetero + posDF[i]],] <- matrix(tmp[2,],
                                                                ncol=pRAIDS$nbSim)

            # number of SNV homozygote corresponding to
            # df$count.tot[i] and df$lap[i]
            nbHomo <- df$Freq[i] - nbHetero
            homo <- which(gOrder[(posDF[i]+1):(posDF[i+1])] != 1)

            tmpHomo <- rmultinom(nbHomo * pRAIDS$nbSim,
                                 as.numeric(as.character(df$count.tot[i])),
                                 c(1- 3 * pRAIDS$seqErrorSyn, pRAIDS$seqErrorSyn, 2*pRAIDS$seqErrorSyn))

            # The order between between ref and alt is done with the phase
            # later
            matSim1[listOrderSNP[homo + posDF[i]],] <- matrix(tmpHomo[1,],
                                                              ncol=pRAIDS$nbSim)
            matSim2[listOrderSNP[homo + posDF[i]],] <- matrix(tmpHomo[2,],
                                                              ncol=pRAIDS$nbSim)
        }

        # superPop of the 1kg sample r is the same
        # for 1kg in list listPosRef.1kg and listPosRef
        curSP <- superPop[r]
        # define a negative block for SNV not in block
        # blockDF[,curSP][which(blockDF[,curSP] == 0)] <-
        #     -1*seq_len(length(which(blockDF[,curSP] == 0)))
        if(length(which(blockDF[,curSP] == 0)) > 0){
            stop("There is block set to 0\n")
        }

        listB <- unique(blockDF[,curSP])

        # block where the phase switch
        recombSwitch <- matrix(sample(x=c(0, 1), size=pRAIDS$nbSim *(length(listB)),
                                      replace=TRUE,
                                      prob=c(1-pRAIDS$pRecomb, pRAIDS$pRecomb)), ncol=pRAIDS$nbSim)

        # indice for each zone with the same phase
        blockZone <- apply(recombSwitch, 2, cumsum)

        rownames(blockZone) <- listB


        ## FOR_LOOP modification to be validated by Pascal
        ## Remove commented code and this text after validation

        # We have to manage multiple simulation which means
        # different number of zone for the different simulations
        LAPparent <- matrix(nrow = nbSNV, ncol = pRAIDS$nbSim)
        # for(i in seq_len(nbSim)){
        #     # list of zone with the same phase relatively to 1KG
        #     listZone <- unique(blockZone[,i])
        #
        #     ## matrix if the lap is the first entry in the phase or
        #     ## the second for each zone
        #     lapPos <- matrix(sample(x=c(0,1), size=1 *(length(listZone)),
        #                                 replace=TRUE), ncol=1)
        #
        #     rownames(lapPos) <- listZone
        #
        #     LAPparent[, i] <-
        #                 lapPos[as.character(blockZone[as.character(blockDF[,
        #                                                         curSP]),i]),]
        # }


        # We have to manage multiple simulations which means
        # different number of zones for the different simulations
        lapValues <- vapply(seq_len(pRAIDS$nbSim), function(i) {
            # list of zone with the same phase relatively to 1KG
            listZone <- unique(blockZone[,i])

            ## matrix if the lap is the first entry in the phase or
            ## the second for each zone
            lapPos <- matrix(sample(x=c(0,1), size=1 *(length(listZone)),
                                    replace=TRUE), ncol=1)

            rownames(lapPos) <- listZone

            return(lapPos[as.character(blockZone[as.character(blockDF[,
                                                                      curSP]),i]),])
        }, double(nbSNV))
        LAPparent[, seq_len(pRAIDS$nbSim)] <- lapValues

        phaseVal <- read.gdsn(index.gdsn(gdsRefAnnot, "phase"),
                              start=c(1,listPosRef.1kg[r]), count=c(-1,1))[list1KG]

        # mat1 is lap mat2 is 1-lap
        # LAPparent if 0 lap left and 1 lap is right

        # Ok note phaseVal must be the first allele
        tmp <- phaseVal + g * LAPparent
        refC <- matSim1 * ((tmp+1) %% 2) + matSim2 * ((tmp) %% 2)
        altC <- matSim1 * ((tmp) %% 2) + matSim2 * ((tmp+1) %% 2)
        rm(phaseVal, tmp)

        # infoSNV$count.tot
        listCount <- table(infoSNV$count.tot)
        cutOffA <- data.frame(count=unlist(vapply(as.integer(names(listCount)),
                                                  FUN=function(x, minProb, eProb){
                                                      return(max(2,qbinom(minProb, x,eProb))) },
                                                  FUN.VALUE=numeric(1), minProb=pRAIDS$minProb,
                                                  eProb= 2 * pRAIDS$seqErrorSyn )),
                              allele=unlist(vapply(as.integer(names(listCount)),
                                                   FUN=function(x, minProb, eProb){
                                                       return(max(2,qbinom(minProb, x,eProb))) },
                                                   FUN.VALUE=numeric(1), minProb=pRAIDS$minProb, eProb=pRAIDS$seqErrorSyn)))
        row.names(cutOffA) <- names(listCount)

        gSyn <- matrix(rep(-1, pRAIDS$nbSim * nrow(infoSNV)), nrow=nrow(infoSNV))

        # count total multiply by 0 if too much error
        gSyn <- gSyn +
            (infoSNV$count.tot - (refC + altC) <
                 cutOffA[as.character(infoSNV$count.tot), "count"]) *
            ((refC == 0 |  altC == 0) + # 1 if homozygot
                 (refC >= cutOffA[as.character(infoSNV$count.tot), "allele"]) *
                 (altC >= cutOffA[as.character(infoSNV$count.tot), "allele"])
             # 1 if both allele are higher than cutoff hetero
            ) * # 1 if homozygote or hetero and 0 if both > 0 both can't
            # decide if error or hetero
            (1 + (altC > 0) * (1 + (refC == 0)))
        # if altC == 0 than 1, altC > 0 and refC == 0 than 3

        ## Append the profile names of to the Profile GDS file "sample.id" node
        appendGDSSampleOnly(gds=gdsProfile,
                            listSamples=paste(paste0(pRAIDS$prefix, ".", pRAIDS$pedStudy$Name.ID[1]),
                                              paste(rep(sample.1kg[curSynt], each=pRAIDS$nbSim),
                                                    seq_len(pRAIDS$nbSim), sep="."), sep = "."))

        ## Append the genotype matrix to the GDS Sample file "genotype" node
        appendGDSgenotypeMat(gdsProfile, gSyn)
    }

    ## Close GDS Sample file
    closefn.gds(gdsProfile)
    closefn.gds(gdsReference)
    closefn.gds(gdsRefAnnot)

    return(0L)
}
