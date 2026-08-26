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


#' @title An S4 class to represent the RAIDS parameters
#'
#' @slot studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor). Default: 
#' \code{data.frame(study.id="NotDef", study.desc="NotDef", 
#' study.platform="NotDef", stringsAsFactors=FALSE)}.
#' 
#' @slot studyDFSyn a \code{data.frame} containing the information about the
#' synthetic data to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor). Default: \code{NULL}.
#' 
#' @slot pedStudy a \code{data.frame} containing the information TODO 
#' with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", and "Source". All columns must be in
#' \code{character} strings (no factor). All row names should correspond to the
#' values in the "Name.ID" column.
#' Default: \code{data.frame(Name.ID=c("ProfileId"), Case.ID=c("ProfileId"), 
#' Sample.Type=c("type"), Diagnosis="NotDef", Source=c("NotDef"), 
#' stringsAsFactors=FALSE, row.names = c("ProfileId"))}.
#' 
#' @slot studyType a \code{character} string representing the type of study.
#' The possible choices are: "LD" and "GeneAware". The type of study affects 
#' the way the estimation of the allelic fraction is done. 
#' Default: \code{"LD"}.
#' 
#' @slot genoSource TODO
#' 
#' @slot blockTypeId  a single \code{character} string corresponding to 
#' the block type used to extract the block identifiers. The block type must 
#' be present in the GDS Reference Annotation file. 
#' Default: \code{"GeneS.Ensembl.Hsapiens.v86"}.
#' 
#' @slot reference a \code{character} string with two possible values:
#' '1KGv1.0', '1k_hgdpV0.1'. It specifies the type of inference. 
#' Default: \code{"1KGc1.0"}.
#' 
#' @slot genome a \code{character} string with one possible value:
#' 'HG38'. It specifies the genome used. Default: \code{"HG38"}.
#' 
#' @slot chrInfo a \code{vector} of positive \code{integer} values
#' representing the length of the chromosomes. See 'details' section. 
#' Default: \code{seqlengths(Hsapiens)[seq_len(25)]}.
#' 
#' @slot paramAncestry a \code{list} of parameters related to ancestry. 
#' The \code{list} should contain those three entries: \code{ScanBamParam}, 
#' \code{PileupParam}, and \code{yieldSize}.  
#' Default: \code{list(ScanBamParam=NULL, PileupParam=NULL, 
#' yieldSize=10000000)}.
#' 
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
#' @slot prefix a single \code{character} string representing TODO.
#' Default: \code{"1"}.
#' 
#' @slot nbSim a positive \code{integer} representing TODO. 
#' Default: \code{1L}.
#' 
#' @slot offset a single \code{integer} that is added to the SNP position to
#' switch from 0-based to 1-based coordinate when needed (or reverse).
#' Default: \code{-1L}.
#' 
#' @slot minCov a single positive \code{integer} representing the minimum
#' required coverage. Default: \code{10L}.
#' 
#' @slot minProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability that the calculated genotype call is correct.
#' Default: \code{0.999}.
#' 
#' @slot seqError a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error. 
#' Default: \code{0.001}.
#' 
#' @slot seqErrorSyn a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error for synthetic.
#' Default: \code{0.001}.
#' 
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
#' 
#' @slot phase a \code{logical} indicating TODO. Default: \code{FALSE}.
#' 
#' @slot PCAmissingRate TODO
#' @slot PCAalgorithm TODO
#' @slot eigenCount TODO
#' @slot eigenCountSyn TODO
#' 
#' @slot kList a \code{vector} of positive \code{integer} representing the 
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis. 
#' Default: \code{seq(2L, 15L, 1L)}.
#' 
#' @slot pcaList a \code{vector} of positive \code{integer} representing 
#' the values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis. 
#' Default: \code{seq(2L, 15L, 1L)}.
#' 
#' @slot fieldPopInRef a \code{character} string representing the name of the
#' column that contains the known ancestry for the reference profiles in
#' the Population Reference GDS file (corresponding to the 
#' \code{fileReferenceGDS} parameter). Default: \code{"superPop"}.
#' 
#' @slot fieldPopInfAnc a \code{character} string representing the name of 
#' the column in the data frame that contains the inferred 
#' super-population ancestry for the samples. The column should be 
#' present in the data frame. Default: \code{"superPop"}.
#' 
#' @slot fieldSubPop a \code{character} string representing the name of 
#' the column in the Population Reference GDS file (corresponding to the 
#' \code{fileReferenceGDS} parameter) that contains the sub-population 
#' information for the samples. The column should be present in the Population 
#' Reference GDS file. Default: \code{"pop.group"}.
#' 
#' @slot verbose a \code{logical} indicating if message information should be
#' printed. Default: \code{FALSE}.
#' 
#' 
#' @details
#'
#' The `chrInfo` slot contains the length of the chromosomes. The
#' length of the chromosomes can be obtain through the
#' \code{\link[Seqinfo]{seqlengths}} library.
#'
#' As example, for hg38 genome:
#'
#' ```
#'
#' if (requireNamespace("Seqinfo", quietly=TRUE) &&
#'      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#'      chrInfo <- Seqinfo::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]
#' }
#'
#' ```
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @aliases RAIDSparam-class
#' @name RAIDSparam-class
#' @rdname RAIDSparam-class
#' 
#' @keywords classes
#' @exportClass RAIDSparam
#' @importFrom Seqinfo seqlengths
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 Hsapiens
#' @export
setClass("RAIDSparam",
  slots = c(
    studyDF = "data.frame",
    studyDFSyn = "data.frame",
    pedStudy="data.frame",
    studyType="character",
    genoSource="character",
    blockTypeId="character",
    reference="character",
    genome="character",
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
    studyDF=data.frame(study.id="NotDef",
                                study.desc="NotDef",
                                study.platform="NotDef",
                                stringsAsFactors=FALSE),
    studyDFSyn=data.frame(study.id="NotDef.Synthetic",
                        study.desc="NotDef synthetic data",
                        study.platform="Synthetic", stringsAsFactors=FALSE),
    pedStudy=data.frame(Name.ID=c("ProfileId"),
                                Case.ID=c("ProfileId"),
                                Sample.Type=c("type"),
                                Diagnosis="NotDef",
                                Source=c("NotDef"),
                                stringsAsFactors=FALSE, 
                                row.names=c("ProfileId")),
    studyType="LD",
    genoSource=NULL,
    blockTypeId="GeneS.Ensembl.Hsapiens.v86",
    reference="1KGc1.0",
    genome="HG38",
    chrInfo=seqlengths(Hsapiens)[seq_len(25)],
    paramAncestry=list(ScanBamParam=NULL, PileupParam=NULL,
                                yieldSize=10000000),
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
        if (!(is.data.frame(object@studyDF) && 
                all(c("study.id", "study.desc", "study.platform") %in% 
                    colnames(object@studyDF)))) {
            return(paste0("'studyDF' slot must be a data.frame with ", 
                "those 3 columns: \"study.id\", \"study.desc\", ", 
                "\"study.platform\"."))
        }
      
        ## Validate the studyDFSyn parameter
        if (!(is.data.frame(object@studyDFSyn) && 
                all(c("study.id", "study.desc", "study.platform") %in% 
                    colnames(object@studyDFSyn)))) {
            return(paste0("'studyDFSyn' slot must be a data.frame ", 
                "with those 3 columns: \"study.id\", \"study.desc\", ", 
                "\"study.platform\"."))
        }

        ## Validate the pedStudy parameter
        if (!(is.data.frame(object@pedStudy) && 
                all(c("Name.ID", "Case.ID", "Sample.Type", "Diagnosis", 
                        "Source") %in% colnames(object@pedStudy)) && 
                all(rownames(object@pedStudy) == object@pedStudy$Name.ID))) {
            return(paste0("'pedStudy' slot must be a data.frame ", 
                "with those 3 columns: \"Name.ID\", \"Case.ID\", ", 
                "\"Sample.Type\", \"Diagnosis\", and \"Source\". All row ", 
                "names should correspond to the Name.ID values."))
        }
      
        ## Validate the studyType parameter
        if (!(is.character(object@studyType) && 
                length(object@studyType) != 1 || 
                !object@studyType %in% c("LD", "GeneAware"))) {
            return(paste0("'studyType' slot must have one character", 
                " string within those 2 choices: \"LD\" and \"GeneAware\"."))
        }
      
        ## Validate the genoSource parameter TODO
      
        ## Validate thte blockTypeId parameter
        if (length(object@blockTypeId != 1)) {
            return("'blockTypeId' slot must have one character string.")
        }

        ## Validate the reference parameter
        if (length(object@reference) != 1 || 
            !object@reference %in% c("1KGv1.0", "1k_hgdpV0.1")) {
            return(paste0("'reference' slot must have one character string", 
                " within those 2 choices: \"1KGv1.0\", \"1k_hgdpV0.1\"."))
        }

        ## Validate the genome parameter
        if (length(object@genome != 1) || object@genome != "HG38") {
            return("'genome' slot must be the character string \"HG38\".")
        }

        ## Validate the chrInfo parameter TODO

        ## Validate the paramAncestry parameter 
        if (length(object@paramAncestry) != 3 || 
                    !all(c("ScanBamParam", "PileupParam", "yieldSize") %in% 
                      object@paramAncestry)) {
            return("'paramAncestry' slot must be a list with those three ", 
                        "entries: \"ScanBamParam\", \"PileupParam\", ", 
                        "and \"yieldSize\".")
        }
      
        ## Validate the profileFile parameter TODO
        ## Validate the profileFileGeno parameter TODO
        ## Validate the pathProfileGDS parameter TODO
        ## Validate the fileReferenceGDS parameter TODO
        ## Validate the fileReferenceAnnotGDS parameter TODO
        ## Validate the inferenceType parameter TODO
        ## Validate the sampleRef parameter TODO

        ## Validate the batch parameter
        if (length(object@batch != 1) || object@batch < 1) {
            return("'batch' slot must have one positive integer.")
        }

        ## Validate the prefix parameter
        if (length(object@prefix != 1)) {
            return("'prefix' slot must have one character string.")
        }

        ## Validate the nbSim parameter
        if (length(object@nbSim != 1) || object@nbSim < 1L) {
            return("'offset' slot must have one positive integer.")
        }

        ## Validate the offset parameter
        if (length(object@offset != 1)) {
            return("'offset' slot must have one integer.")
        }
        
        ## Validate the minCov parameter
        if (length(object@minCov != 1) || object@minCov < 1L) {
            return("'minCov' slot must have one positive integer.")
        }

        ## Validate the minProb parameter
        if (length(object@minProb != 1) || 
                !(object@minProb > 0 && object@minProb < 1)) {
            return(paste0("'minProb' slot must have one positive numeric", 
                        " between 0 and 1."))
        }

        ## Validate the seqError parameter
        if (length(object@seqError != 1) || 
                !(object@seqError > 0 && object@seqError < 1)) {
            return(paste0("'seqError' slot must have one positive numeric", 
                        " between 0 and 1."))
        }

        ## Validate the seqErrorSyn parameter
        if (length(object@seqErrorSyn != 1) || 
                !(object@seqErrorSyn > 0 && object@seqErrorSyn < 1)) {
            return(paste0("'seqErrorSyn' slot must have one positive numeric", 
                        " between 0 and 1."))
        }

        ## Validate the pRecomb parameter TODO
        ## Validate the np parameter TODO
        ## Validate the listPos parameter TODO
        ## Validate the syntheticRefDF parameter TODO

        ## Validate the pruningMethod parameter
        if (length(object@pruningMethod) != 1 || 
              !object@pruningMethod %in% c("corr", "r", "dprime",
                  "composite")) {
          return(paste0("'pruningMethod' slot must have one character string ", 
              "within those 4 choices: \"corr\", \"r\", \"dprime\", ", 
              "\"composite\"."))
        }

        ## Validate the slideWindowMaxBP parameter TODO
        ## Validate the thresholdLD parameter TODO
        ## Validate the specificSNV parameter TODO
        ## Validate the genoType parameter TODO
        ## Validate the phaseType parameter TODO
        
        ## Validate the phase parameter
        if (length(object@phase) != 1) {
          return("'verbose' slot must have one logical value.")
        }

        ## Validate the PCAmissingRate parameter TODO
        ## Validate the PCAalgorithm parameter TODO
        ## Validate the eigenCount parameter TODO
        ## Validate the eigenCountSyn parameter TODO

        ## Validate the kList parameter
        if (!all(object@kList > 0)) {
          return(paste0("'kList' slot must have one or more positive ", 
                    "integer values."))
        }

        ## Validate the pcaList parameter
        if (!all(object@pcaList > 0)) {
          return(paste0("'pcaList' slot must have one or more positive ", 
                    "integer values."))
        }

        ## Validate the fieldPopInRef parameter
        if (length(object@fieldPopInRef) != 1) {
          return("'fieldPopInRef' slot must have one character string.")
        }

        ## Validate the fieldPopInfAnc parameter
        if (length(object@fieldPopInfAnc) != 1) {
          return("'fieldPopInfAnc' slot must have one character string.")
        }

        ## Validate the fieldSubPop parameter
        if (length(object@fieldSubPop) != 1) {
          return("'fieldSubPop' slot must have one character string.")
        }

        ## Validate the verbose parameter
        if (length(object@verbose) != 1) {
          return("'verbose' slot must have one logical value.")
        }

        TRUE
    }
)

#' @title Create a RAIDSparam object 
#'
#' @description The function uses all parameters to create a 
#' \code{RAIDSparam} class object that contains all the parameters used 
#' in the RAIDS workflow.
#' 
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor). If \code{NULL}, the 
#' following will be assigned: 
#' \code{data.frame(study.id="NotDef", study.desc="NotDef", 
#' study.platform="NotDef", stringsAsFactors=FALSE)}.
#' Default: \code{NULL}.
#' 
#' @param studyDFSyn a \code{data.frame} containing the information about the
#' synthetic data to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor). If \code{NULL}, the 
#' following will be assigned: \code{data.frame(study.id="NotDef.Synthetic", 
#' study.desc="NotDef synthetic data", study.platform="Synthetic", 
#' stringsAsFactors=FALSE)}.
#' Default: \code{NULL}.
#' 
#' @param pedStudy a \code{data.frame} containing the information TODO 
#' with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", and "Source". All columns must be in
#' \code{character} strings (no factor). All row names should correspond to the
#' values in the "Name.ID" column. If \code{NULL}, the 
#' following will be assigned:  \code{data.frame(Name.ID=c("ProfileId"), 
#' Case.ID=c("ProfileId"), Sample.Type=c("type"), Diagnosis="NotDef", 
#' Source=c("NotDef"), stringsAsFactors=FALSE, row.names = c("ProfileId"))}.
#' Default: \code{NULL}.
#' 
#' @param studyType a \code{character} string representing the type of study.
#' The possible choices are: "LD" and "GeneAware". The type of study affects 
#' the way the estimation of the allelic fraction is done. 
#' Default: \code{"LD"}.
#' 
#' @param genoSource TODO. The valid options are: "VCF", "generic", 
#' "snp-pileup", and "bam". Default: \code{NULL}.
#' 
#' @param blockTypeId a single \code{character} string corresponding to 
#' the block type used to extract the block identifiers. The block type must 
#' be present in the GDS Reference Annotation file. 
#' Default: \code{"GeneS.Ensembl.Hsapiens.v86"}.
#' 
#' @param reference a \code{character} string with two possible values:
#' '1KGv1.0', '1k_hgdpV0.1'. It specifies the type of inference. 
#' Default: \code{"1KGc1.0"}.
#' 
#' @param genome a \code{character} string with one possible value:
#' 'HG38'. It specifies the genome used. Default: \code{"HG38"}.
#' 
#' @param chrInfo TODO
#' 
#' @param paramAncestry a \code{list} of parameters related to ancestry. 
#' The \code{list} should contain those three entries: \code{ScanBamParam}, 
#' \code{PileupParam}, and \code{yieldSize}. If \code{NULL}, the default value 
#' \code{list(ScanBamParam=NULL, PileupParam=NULL, yieldSize=10000000)} will 
#' be assigned to the parameter. Default: \code{NULL}.
#' 
#' @param profileFile a \code{character} string representing the path to the
#' file with genotype and the allele information of the profile. The format 
#' of the accepted file is determined by the value in the 'genoSource' slot.
#' A profile would have a file with extension "*vcf.gz" when 'genoSource' slot 
#' is "VCF". If 'genoSource' is "generic" or "snp-pileup", then "*.txt.gz". 
#' If 'genoSource' is "bam", then "*.bam" (the file need to be indexed with a
#' existing corresponding "*.bai" file). Default: \code{NULL}.
#' 
#' @param profileFileGeno TODO
#' @param pathProfileGDS TODO
#' @param fileReferenceGDS TODO
#' @param fileReferenceAnnotGDS TODO
#' @param inferenceType TODO
#' @param sampleRef TODO
#' 
#' @param batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore. 
#' Default: \code{1L}.
#' 
#' @param prefix a single \code{character} string representing TODO.
#' Default: \code{"1"}.
#' 
#' @param nbSim a positive \code{integer} representing TODO. 
#' Default: \code{1L}.
#' 
#' @param offset a single \code{integer} that is added to the SNP position to
#' switch from 0-based to 1-based coordinate when needed (or reverse).
#' Default: \code{-1L}.
#' 
#' @param minCov a single positive \code{integer} representing the minimum
#' required coverage. Default: \code{10L}.
#' 
#' @param minProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability that the calculated genotype call is correct.
#' Default: \code{0.999}.
#' 
#' @param seqError a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error. 
#' Default: \code{0.001}.
#' 
#' @param seqErrorSyn a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error for synthetic.
#' Default: \code{0.001}.
#' 
#' @param pRecomb TODO
#' @param np TODO
#' @param listPos TODO
#' @param syntheticRefDF TODO
#' 
#' @param pruningMethod a \code{character} string representing the method that 
#' will be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The 4 possible values
#' are: "corr", "r", "dprime", and "composite". Default: \code{"corr"}.
#' 
#' @param slideWindowMaxBP TODO
#' @param thresholdLD TODO
#' @param specificSNV TODO
#' @param genoType TODO
#' @param phaseType TODO
#' 
#' @param phase a \code{logical} indicating TODO. Default: \code{FALSE}.
#' 
#' @param PCAmissingRate TODO
#' @param PCAalgorithm TODO
#' @param eigenCount TODO
#' @param eigenCountSyn TODO
#' 
#' @param kList a \code{vector} of positive \code{integer} representing the 
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis. 
#' Default: \code{seq(2L, 15L, 1L)}.
#' 
#' @param pcaList a \code{vector} of positive \code{integer} representing 
#' the values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis. 
#' Default: \code{seq(2L, 15L, 1L)}.
#' 
#' @param fieldPopInRef a \code{character} string representing the name of the
#' column that contains the known ancestry for the reference profiles in
#' the Population Reference GDS file (corresponding to the 
#' \code{fileReferenceGDS} parameter). Default: \code{"superPop"}.
#' 
#' @param fieldPopInfAnc a \code{character} string representing the name of 
#' the column in the data frame that contains the inferred 
#' super-population ancestry for the samples. The column should be 
#' present in the data frame. Default: \code{"superPop"}.
#' 
#' @param fieldSubPop a \code{character} string representing the name of 
#' the column in the Population Reference GDS file (corresponding to the 
#' \code{fileReferenceGDS} parameter) that contains the sub-population 
#' information for the samples. The column should be present in the Population 
#' Reference GDS file. Default: \code{"pop.group"}.
#' 
#' @param verbose a \code{logical} indicating if message information should be
#' printed. Default: \code{FALSE}
#' 
#' @return an object of class \code{RAIDSparam} that contains all the 
#' required parameters needed by the RAIDS workflow.
#' 
#' @examples
#'
#'
#'  ## TODO
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @importFrom rlang arg_match
#' @importFrom methods new
#' @export
RAIDSparam <- function(studyDF=NULL, studyDFSyn=NULL, pedStudy=NULL, 
    studyType=c("LD", "GeneAware"), genoSource=NULL, 
    blockTypeId="GeneS.Ensembl.Hsapiens.v86", reference="1KGc1.0", 
    genome="HG38", chrInfo=NULL, paramAncestry=NULL, profileFile=NULL,
    profileFileGeno=NULL, pathProfileGDS=NULL, fileReferenceGDS=NULL,
    fileReferenceAnnotGDS=NULL, inferenceType="PCAknn", sampleRef=NULL,
    batch=1L, prefix="1", nbSim=1L, offset=-1L, minCov=10L, minProb=0.999,
    seqError=0.001, seqErrorSyn=0.001, pRecomb=0.01, np=1L, listPos=NULL, 
    syntheticRefDF=NULL, pruningMethod=c("corr", "r", "dprime", "composite"), 
    slideWindowMaxBP=500000L, thresholdLD=sqrt(0.1), 
    specificSNV=NULL, genoType="geno.ref",
    phaseType="phase.ref", phase=FALSE, PCAmissingRate=0.025, 
    PCAalgorithm="exact", eigenCount=32L, eigenCountSyn=15L, 
    kList=seq(2L, 15L, 1L), pcaList=seq(2L, 15L, 1L), fieldPopInRef="superPop",
    fieldPopInfAnc="superPop", fieldSubPop="pop.group", verbose=FALSE) {
    
    batch <- as.integer(batch)
    nbSim <- as.integer(nbSim)
    offset <- as.integer(offset)
    minCov <- as.integer(minCov)
    np <- as.integer(np)
    slideWindowMaxBP <- as.integer(slideWindowMaxBP)
    eigenCount <- as.integer(eigenCount)
    eigenCountSyn <- as.integer(eigenCountSyn)
    kList <- unique(as.integer(kList))
    pcaList <- unique(as.integer(pcaList))
    verbose <- as.logical(verbose)

    studyType <- arg_match(studyType)
    pruningMethod <- arg_match(pruningMethod)
    
    if (is.null(studyDF)) {
        studyDF <- data.frame(study.id="NotDef",
                                study.desc="NotDef",
                                study.platform="NotDef",
                                stringsAsFactors=FALSE)
    }

    if (is.null(studyDFSyn)) {
        studyDFSyn <- data.frame(study.id=paste0(studyDF$study.id, 
            ".Synthetic"),
            study.desc=paste0(studyDF$study.id, " synthetic data"),
            study.platform="Synthetic", stringsAsFactors=FALSE)
    }

    if(is.null(pedStudy)){
        pedStudy <- data.frame(Name.ID=c("ProfileId"),
                                Case.ID=c("ProfileId"),
                                Sample.Type=c("type"),
                                Diagnosis="NotDef",
                                Source=c("NotDef"),
                                stringsAsFactors=FALSE,
                                row.names=c("ProfileId"))
    }
  
    ## Assign default ancestry parameters when not assigned by user
    if(is.null(paramAncestry)) {
      paramAncestry <- list(ScanBamParam=NULL, PileupParam=NULL, 
                                yieldSize=10000000)
    }

    new("RAIDSparam", studyDF=studyDF, studyDFSyn=studyDFSyn,
    studyType=studyType, genoSource=genoSource, blockTypeId=blockTypeId,
    reference=reference, chrInfo=chrInfo, paramAncestry=paramAncestry,
    profileFile=profileFile, profileFileGeno=profileFileGeno,
    pathProfileGDS=pathProfileGDS, fileReferenceGDS=fileReferenceGDS,
    fileReferenceAnnotGDS=fileReferenceAnnotGDS,
    inferenceType=inferenceType, sampleRef=sampleRef, batch=batch, 
    prefix=prefix, nbSim=nbSim, offset=offset, minCov=minCov, minProb=minProb,
    seqError=seqError, seqErrorSyn=seqErrorSyn, pRecomb=pRecomb, np=np,
    listPos=listPos, syntheticRefDF=syntheticRefDF, pruningMethod=pruningMethod,
    slideWindowMaxBP=slideWindowMaxBP, thresholdLD=thresholdLD,
    specificSNV=specificSNV, genoType=genoType, phaseType=phaseType, 
    phase=phase, PCAmissingRate=PCAmissingRate, PCAalgorithm=PCAalgorithm,
    eigenCount=eigenCount, eigenCountSyn=eigenCountSyn, kList=kList,
    pcaList=pcaList, fieldPopInRef=fieldPopInRef, fieldPopInfAnc=fieldPopInfAnc,
    fieldSubPop=fieldSubPop, verbose=verbose)
}