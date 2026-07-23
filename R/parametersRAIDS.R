#' @title parameters class for RAIDS
#' function
#'
#' @description parameters class for RAIDS
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor). Default: \code{NULL}.
#'
#' @param studyDFSyn a \code{data.frame} containing the information about the
#' synthetic data to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor). Default: \code{NULL}.
#'
#' @param pedStudy a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings (no factor). The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter. Only \code{filePedRDS} or \code{pedStudy}
#' can be defined. Default: \code{NULL}.
#'
#' @param studyType a \code{character} string representing the type of study.
#' The possible choices are: "LD" and "GeneAware". The type of study affects the
#' way the estimation of the allelic fraction is done. Default: \code{"LD"}.
#'
#' @param genoSource a TODO. Default: \code{NULL}.
#' 
#' @param blockTypeId a \code{character} string corresponding to the block
#' type used to extract the block identifiers. The block type must be
#' present in the GDS Reference Annotation file. 
#' Default: \code{"GeneS.Ensembl.Hsapiens.v86"}.
#'
#' @param reference a \code{character} string with two possible values:
#' '1KGv1.0', '1k_hgdpV0.1'. It specifies the type of inference. 
#' Default: \code{"1KGc1.0"}.
#'
#' @param genome a \code{character} string with one possible value:
#' 'HG38'. It specifies the genome used. Default: \code{"HG38"}.
#'
#' @param chrInfo a \code{vector} of positive \code{integer} values
#' representing the length of the chromosomes. See 'details' section. 
#' Default: \code{NULL}.
#'
#' @param paramAncestry a \code{list} parameters ...TODO. Default: \code{NULL}.
#'
#' @param profileFile a \code{character} string representing the path and the
#' file name of the genotype file or the bam if genoSource is snp-pileup the
#' fine extension must be .txt.gz, if VCF the extension must be .vcf.gz. 
#' Default: \code{NULL}.
#' 
#' @param profileFileGeno a \code{character} string TODO. 
#' Default: \code{NULL}.
#' 
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the GDS Profile files will be created.
#' Default: \code{NULL}.
#'
#' @param fileReferenceGDS  a \code{character} string representing the file
#' name of the Population Reference GDS file. The file must exist. 
#' Default: \code{NULL}.
#'
#' @param fileReferenceAnnotGDS a \code{character} string representing the
#' file name of the Population Reference GDS Annotation file. The file
#' must exist. Default: \code{NULL}.
#'
#' @param inferenceType a \code{character} string with two possible values:
#' 'PCAknn', 'haploAdmixture'. It specifies the genotype ancestry inference. 
#' Default: \code{"PCAknn"}.
#'
#' @param sampleRef a \code{array} is a subset of sample.id
#' of the gdsReference. If it NULL mean all the sample.id.
#' Only the sample in sampleRef will be used for the ancestry inference.
#' Example: remove the admix samples from the reference.
#' Default: \code{NULL}.
#' 
#' @param batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore. 
#' Default: \code{1}.
#' 
#' @param prefix a TODO. Default: \code{"1"}.
#' 
#' @param nbSim a TODO. Default: \code{1}.
#' 
#' @param offset a TODO. Default: \code{-1}.
#' 
#' @param minCov a single positive \code{integer} representing the minimum
#' required coverage. Default: \code{10L}.
#' 
#' @param minProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability that the calculated genotype call is correct.
#' Default: \code{0.999}.
#' 
#' @param seqError a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error. Default: \code{0.001}.
#'
#' @param seqErrorSyn a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error for synthetic.
#' Default: \code{0.001}.
#' 
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#' 
#' @param pRecomb a TODO. Default: \code{0.01}.
#' 
#' @param listPos a \code{data.frame} containing 2 columns. The first column,
#' called "snp.chromosome" contains the name of the chromosome where the
#' SNV is located. The second column, called "snp.position" contains the
#' position of the SNV on the chromosome. Default: \code{NULL}.
#'
#' 
#'
#' @param syntheticRefDF a \code{data.frame} containing a subset of
#' reference profiles for each sub-population present in the Reference GDS
#' file. The \code{data.frame} must have those columns:
#' \describe{
#' \item{sample.id}{ a \code{character} string representing the sample
#' identifier. }
#' \item{pop.group}{ a \code{character} string representing the
#' subcontinental population assigned to the sample. }
#' \item{superPop}{ a \code{character} string representing the
#' super-population assigned to the sample. }
#' }
#' Default: \code{NULL}.
#'
#' @param pruningMethod a \code{character} string that represents the method that will
#' be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The 4 possible values
#' are: "corr", "r", "dprime" and "composite". Default: \code{"corr"}.
#'
#' @param slideWindowMaxBP a TODO. Default: \code{500000L}.
#' 
#' @param thresholdLD a TODO. Default: \code{sqrt(0.1)}.
#' 
#' @param specificSNV a \code{data.frame} containing 2 columns. The first column,
#' called "snp.chromosome" contains the name of the chromosome where the
#' SNV is located. The second column, called "snp.position" contains the
#' position of the SNV on the chromosome. Optionel the column snvKeep add by
#' contains index in the reference or -1  if not in the reference. It is using
#' used in the pruning step. Default: \code{NULL}.
#'
#' @param genoType a TODO. Default: \code{"geno.ref"}.
#' 
#' @param phaseType a TODO. Default: \code{"phase.ref"}.
#' 
#' @param phase a \code{logical} TODO. Default: \code{FALSE}.
#' 
#' @param PCAmissingRate a TODO. Default: \code{0.025}.
#' 
#' @param PCAalgorithm a TODO. Default: \code{"exact"}.
#' 
#' @param eigenCount a TODO. Default: \code{32L}.
#' 
#' @param eigenCountSyn a TODO. Default: \code{15L}.
#' 
#' @param kList a TODO. Default: \code{seq(2,15,1)}.
#' 
#' @param pcaList a TODO. Default: \code{seq(2,15,1)}.
#' 
#' @param fieldPopInRef a TODO. Default: \code{"superPop"}.
#' 
#' @param fieldPopInfAnc a TODO. Default: \code{"superPop"}.
#' 
#' @param fieldSubPop a TODO. Default: \code{"pop.group"}.
#' 
#' @param verbose a \code{logical} indicating if message information should be
#' printed. Default: \code{FALSE}.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Load the demo PCA on the synthetic profiles projected on the
#' ## demo 1KG reference PCA
#' data(demoPCASyntheticProfiles)
#'
#'  ## TODO
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @importFrom rlang arg_match
#' @export
paramRAIDS <- function(studyDF=NULL,
                        studyDFSyn=NULL,
                        pedStudy=NULL,
                        studyType="LD",
                        genoSource=NULL,
                        blockTypeId="GeneS.Ensembl.Hsapiens.v86",
                        reference="1KGc1.0",
                        genome="HG38",
                        chrInfo=NULL,
                        paramAncestry=NULL,
                        profileFile=NULL,
                        profileFileGeno=NULL,
                        pathProfileGDS=NULL,
                        fileReferenceGDS=NULL,
                        fileReferenceAnnotGDS=NULL,
                        inferenceType="PCAknn",
                        sampleRef=NULL,
                        batch=1,
                        prefix="1",
                        nbSim=1,
                        offset=-1,
                        minCov=10,
                        minProb=0.999,
                        seqError=0.001,
                        seqErrorSyn=0.001,
                        pRecomb=0.01,
                        np=1L,
                        listPos=NULL,
                        syntheticRefDF=NULL,
                        pruningMethod=c("corr", "r", "dprime", "composite"),
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
                        kList=seq(2,15,1),
                        pcaList=seq(2,15,1),
                        fieldPopInRef="superPop",
                        fieldPopInfAnc="superPop",
                        fieldSubPop = "pop.group",
                        verbose=FALSE) {

    # listSNP=NULL, # not yet implemented
    # superPopMinAF=NULL,
    # keepPrunedGDS=TRUE,
    # pathProfileGDS=NULL,
    # keepFile=FALSE,
    # pathPrunedGDS=".", outPrefix="pruned"
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
    if(is.null(studyDFSyn)){
        studyDFSyn <- data.frame(study.id=paste0(studyDF$study.id, 
            ".Synthetic"),
            study.desc=paste0(studyDF$study.id, " synthetic data"),
            study.platform="Synthetic", stringsAsFactors=FALSE)
    }

    if(is.null(chrInfo) && genome=="HG38"){
        ## TODO PASCAL you cannot call :: inside function
        chrInfo <- Seqinfo::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[seq_len(25)]
    }else if(is.null(chrInfo)){
        stop("If your genome is't HG38 you need to specify chrInfo\n")
    }
    if(is.null(paramAncestry)){
        paramAncestry=list(ScanBamParam=NULL,
                            PileupParam=NULL,
                            yieldSize=10000000)
    }

    # studyDF,
    # currentProfile, pathProfileGDS, chrInfo, syntheticRefDF,
    # studyDFSyn, listProfileRef, studyType=c("LD", "GeneAware"),
    # np=1L, blockTypeID=NULL, verbose=FALSE
    pruningMethod <- arg_match(pruningMethod)

    parameters <- list(reference=reference, # HGDP1kg
                        studyDF=studyDF,
                        studyDFSyn=studyDFSyn,
                        pedStudy=pedStudy,
                        studyType=studyType,
                        genoSource=genoSource,
                        blockTypeId=blockTypeId,
                        genome=genome,
                        chrInfo=chrInfo,
                        profileFile=profileFile,
                        profileFileGeno=profileFileGeno,
                        pathProfileGDS=pathProfileGDS,
                        fileReferenceGDS=fileReferenceGDS,
                        fileReferenceAnnotGDS=fileReferenceAnnotGDS,
                        inferenceType=inferenceType,
                        sampleRef=sampleRef,
                        batch=batch,
                        prefix=prefix,
                        nbSim=nbSim,
                        offset=offset,
                        minCov=minCov,
                        minProb=minProb,
                        seqError=seqError,
                        seqErrorSyn=seqErrorSyn,
                        pRecomb=pRecomb,
                        np=np,
                        listPos=listPos,
                        syntheticRefDF=syntheticRefDF,
                        pruningMethod=pruningMethod,
                        slideWindowMaxBP=slideWindowMaxBP,
                        thresholdLD=thresholdLD,
                        specificSNV=specificSNV,
                        genoType=genoType,
                        phaseType=phaseType,
                        phase=phase,
                        PCAmissingRate=PCAmissingRate,
                        PCAalgorithm=PCAalgorithm,
                        eigenCount=eigenCount,
                        eigenCountSyn=eigenCountSyn,
                        kList=kList,
                        pcaList=pcaList,
                        fieldPopInRef=fieldPopInRef,
                        fieldPopInfAnc=fieldPopInfAnc,
                        fieldSubPop=fieldSubPop,
                        verbose=verbose)
    # class(parameters) <- "parametersRAIDS"
    paramters <- structure(parameters, class = "parametersRAIDS")
    return(parameters)
}

specificSNVKeep <- function(pRAIDS){
    if(!is.null(pRAIDS$specificSNV)){
        if(! "snvKeep" %in% colnames(pRAIDS$specificSNV)){
            gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)

            snp.chromosome <- read.gdsn(node=index.gdsn(gdsReference, "snp.chromosome"))
            keepPos <- seq_len(length(snp.chromosome))
            if("snp.KeepDefault" %in% ls.gdsn(gdsReference) ){
                keepPos <- read.gdsn(index.gdsn(gdsReference, "snp.KeepDefault"))
            }
            snp.chromosome <- snp.chromosome[keepPos]
            snp.position <- read.gdsn(node=index.gdsn(gdsReference, "snp.position"))[keepPos]

            z <- cbind(c(pRAIDS$specificSNV$snp.chromosome,
                         snp.chromosome,
                         pRAIDS$specificSNV$snp.chromosome),
                       c(pRAIDS$specificSNV$snp.position,
                         snp.position,
                         pRAIDS$specificSNV$snp.position),
                       c(-1 * seq_len(nrow(pRAIDS$specificSNV)),
                         rep(0, length(snp.position)),
                         seq_len(nrow(pRAIDS$specificSNV))),
                       c(rep(0, nrow(pRAIDS$specificSNV)),
                         seq_len(length(snp.position)),
                         rep(0, nrow(pRAIDS$specificSNV))))
            z <- z[order(z[,1], z[,2], z[,3]),]
            pRAIDS$specificSNV$snvKeep <- rep(-1, nrow(pRAIDS$specificSNV))
            pRAIDS$specificSNV$snvKeep[-1*cumsum(z[,3])[z[,3] == 0]] <- z[cumsum(z[,3]) < 0 & z[,3] == 0, 4]
            closefn.gds(gdsReference)
        }
    }

    return(pRAIDS)
}


