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
#' @slot studyType a single \code{character} string representing the type of 
#' study. The possible choices are: "LD" and "GeneAware". The type of study 
#' affects the way the estimation of the allelic fraction is done. 
#' Default: \code{"LD"}.
#' 
#' @slot genoSource a single \code{character} string corresponding to the type 
#' of file with the genotype and the allele information of the profile that 
#' will be provided in the 'profileFile' slot. The valid options are: "VCF", 
#' "generic", "snp-pileup", and "bam". The "generic" format CSV file
#' must have at least these columns: 'Chromosome', 'Position', 'Ref', 'Alt', 
#' 'Count', 'File1R', and 'File1A'. The 'Count' is the depth at the 
#' specified position; 'FileR' is the depth of the reference allele, and
#' 'File1A' is the depth of the specific alternative allele.
#' Finally, in the case of a "VCF" file, the file must have at least those 
#' genotype fields: GT, AD, and DP. Default: \code{NULL}.
#' 
#' @slot blockTypeId  a single \code{character} string corresponding to 
#' the block type used to extract the block identifiers. The block type must 
#' be present in the GDS Reference Annotation file. 
#' Default: \code{"GeneS.Ensembl.Hsapiens.v86"}.
#' 
#' @slot reference a \code{character} string with two possible values:
#' '1KGv1.0', '1k_hgdpV0.1'. It specifies the type of inference. 
#' Default: \code{"1KGv1.0"}.
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
#' @slot profileFile a \code{character} string representing the path to the
#' file with genotype and the allele information of the profile. The format 
#' of the accepted file is determined by the value in the 'genoSource' slot.
#' A profile would have a file with extension "*vcf.gz" when 'genoSource' slot 
#' is "VCF". If 'genoSource' is "generic" or "snp-pileup", then "*.txt.gz". 
#' If 'genoSource' is "bam", then "*.bam" (the file needs to be indexed with a
#' existing corresponding "*.bai" file). Default: \code{NULL}.
#' 
#' @slot profileFileGeno TODO
#' 
#' @slot pathProfileGDS a single \code{character} string representing the 
#' path to the directory where the GDS Profile files will be created. If 
#' specified, the directory must exist. 
#' Default: \code{NULL}.
#' 
#' @slot fileReferenceGDS a single \code{character} string representing 
#' the file name of the Reference GDS file. If specified, the file  
#' must exist. Default: \code{NULL}.
#' 
#' @slot fileReferenceAnnotGDS a \code{character} string representing the
#' file name of the Population Reference GDS Annotation file. If specified,
#' the file must exist. Default: \code{NULL}.
#' 
#' @slot inferenceType a single \code{character} string representing the 
#' genotype ancestry inference method. The two possible values:
#' 'PCAknn' and 'haploAdmixture'. Default: \code{"PCAknn"}.
#' 
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
#' @slot pRecomb a single positive \code{numeric} between \code{0} and 
#' \code{1} that represents the frequency of phase switching in the 
#' synthetic profiles. Default: \code{0.01}.
#' 
#' @slot np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#' 
#' @slot listPos a \code{data.frame} containing 2 columns named: 
#' "snp.chromosome" and "snp.position". The first column,
#' called "snp.chromosome", contains the name of the chromosome where the
#' SNV is located. The second column, called "snp.position", contains the
#' position of the SNV on the chromosome. Default: \code{NULL}.
#' 
#' @slot syntheticRefDF a \code{data.frame} containing a subset of
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
#' @slot pruningMethod a \code{character} string representing the method that 
#' will be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The 4 possible values
#' are: "corr", "r", "dprime", and "composite". Default: \code{"corr"}.
#' 
#' @slot slideWindowMaxBP a single positive \code{integer} that represents
#' the maximum basepairs (bp) in the sliding window. This parameter is used
#' for the LD pruning done by the \code{\link[SNPRelate]{snpgdsLDpruning}} 
#' function. Default: \code{500000L}.
#' 
#' @slot thresholdLD a single positive \code{numeric} value that represents 
#' the LD threshold used in the \code{\link[SNPRelate]{snpgdsLDpruning}} 
#' function. Default: \code{sqrt(0.1)}.
#' 
#' @slot specificSNV a \code{data.frame} containing 2 columns. The first 
#' column, called "snp.chromosome" contains the name of the chromosome where 
#' the SNV is located. The second column, called "snp.position", contains the
#' position of the SNV on the chromosome. Optionally, the column "snvKeep" 
#' contains the SNV index in the reference or -1 if not in the reference. It 
#' is used during the pruning step. Default: \code{NULL}.
#' 
#' @slot genoType a TODO. Default: \code{"geno.ref"}.
#' 
#' @slot phaseType a TODO. Default: \code{"phase.ref"}.
#' 
#' @slot phase a \code{logical} indicating TODO. Default: \code{FALSE}.
#' 
#' @slot PCAmissingRate a positive \code{numeric} representing the maximum 
#' missing rate retained accepted to use SNPs in the PCA analysis done 
#' with the the \link[SNPRelate]{snpgdsPCA} function. If \code{NaN}, no 
#' missing threshold. Default: \code{0.025}.
#' 
#' @slot PCAalgorithm a single \code{character} string representing the 
#' algorithm to use with the \link[SNPRelate]{snpgdsPCA} function. The 
#' algorithm must be implemented and available to the 
#' \link[SNPRelate]{snpgdsPCA} function. The current options are:
#' 'exact' or 'randomized'.
#' Default: \code{"exact"}.
#' 
#' @slot eigenCount a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigenCount' <= 0, then all eigenvectors are returned. 
#' Default: \code{32L}.
#' 
#' @slot eigenCountSyn a TODO. Default: \code{15L}.
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
#' @importFrom stringr str_detect
#' @export
setClass("RAIDSparam",
  slots = c(
    studyDF = "data.frame",
    studyDFSyn = "data.frame",
    pedStudy="data.frame",
    studyType="character",
    genoSource="CharacterOrNULL",
    blockTypeId="character",
    reference="character",
    genome="character",
    chrInfo="integer",
    paramAncestry="list",
    profileFile="CharacterOrNULL",
    profileFileGeno="CharacterOrNULL",
    pathProfileGDS="CharacterOrNULL",
    fileReferenceGDS="CharacterOrNULL",
    fileReferenceAnnotGDS="CharacterOrNULL",
    inferenceType="character",
    sampleRef="CharacterOrNULL",
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
    listPos="DataFrameOrNULL",
    syntheticRefDF="DataFrameOrNULL",
    pruningMethod="character",
    slideWindowMaxBP="integer",
    thresholdLD="numeric",
    specificSNV="DataFrameOrNULL",
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
    reference="1KGv1.0",
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
        if (length(object@studyType) != 1 || 
                !object@studyType %in% c("LD", "GeneAware")) {
            return(paste0("'studyType' slot must have one character", 
                " string within those 2 choices: \"LD\" and \"GeneAware\"."))
        }
      
        ## Validate the genoSource parameter 
        if (!is.null(object@genoSource) && !(length(ojbect@genoSource) == 1 && 
            object@genSource %in% c("VCF", "bam", "generic", "snp-pileup"))) {
            return(paste0("'genoSource' slot must have one character ", 
                "string or NULL. The valid options are: \"VCF\", \"bam\", ", 
                "\"generic\", or \"snp-pileup\"."))
        }
      
        ## Validate thte blockTypeId parameter
        if (length(object@blockTypeId) != 1) {
            return("'blockTypeId' slot must have one character string.")
        }

        ## Validate the reference parameter
        if (length(object@reference) != 1 || 
            !object@reference %in% c("1KGv1.0", "1k_hgdpV0.1")) {
            return(paste0("'reference' slot must have one character string", 
                " within those 2 choices: \"1KGv1.0\", \"1k_hgdpV0.1\"."))
        }

        ## Validate the genome parameter
        if (length(object@genome) != 1 || object@genome != "HG38") {
            return("'genome' slot must be the character string \"HG38\".")
        }

        ## Validate the chrInfo parameter TODO

        ## Validate the paramAncestry parameter 
        if (length(object@paramAncestry) != 3 || 
                    !all(c("ScanBamParam", "PileupParam", "yieldSize") %in% 
                        names(object@paramAncestry))) {
            return(paste0("'paramAncestry' slot must be a list with those ", 
                        "three entries: \"ScanBamParam\", \"PileupParam\", ", 
                        "and \"yieldSize\"."))
        }
      
        ## Validate the profileFile parameter
        if (!is.null(object@profileFile) && length(object@profileFile) != 1) {
            return("'profileFile' slot must have one character string.")
        } 
      
        ## Validate the profileFile parameter using genoSource information
        if (!is.null(object@profileFile) && length(object@profileFile) == 1) {
            genoSource <- object@genoSource
            if (genoSource == "bam" && 
                    !stringr::str_detect(object@profileFile, ".bam$")) {
                return(paste0("'profileFile' slot must have one character ", 
                    "string representing a file with extension '.bam' ", 
                    "according to 'genoSource' slot."))
            } else if (genoSource == "VCF" && 
                    !stringr::str_detect(object@profileFile, ".vcf.gz$")) {
                return(paste0("'profileFile' slot must have one character ", 
                    "string representing a file with extension '.vcf.gz' ", 
                    "according to 'genoSource' slot."))
            } else if (genoSource %in% c("snp-pileup", "generic") && 
                    !stringr::str_detect(object@profileFile, ".txt.gz$")) {
                return(paste0("'profileFile' slot must have one character ", 
                    "string representing a file with extension '.txt.gz' ", 
                    "according to 'genoSource' slot."))
            } 
        } 
        
        ## Validate the profileFile file exists when not null
        if (!is.null(object@profileFile) && !dir.exists(object@profileFile)) {
            return(paste0("'profileFile' slot must have one character ", 
                    "string representing an existing file."))
        }
      
        ## Validate the profileFileGeno parameter TODO
      
        ## Validate the pathProfileGDS parameter
        if ((!is.null(object@pathProfileGDS) && 
                    length(object@pathProfileGDS) != 1) || 
            (!is.null(object@pathProfileGDS) && 
                    !dir.exists(object@pathProfileGDS))) {
            return(paste0("'pathProfileGDS' slot must have one character ", 
                    "string representing an existing directory."))
        }
      
        ## Validate the fileReferenceGDS parameter 
        if ((!is.null(object@fileReferenceGDS) && 
                    length(object@fileReferenceGDS) != 1) || 
            (!is.null(object@fileReferenceGDS) && 
                    !file.exists(object@fileReferenceGDS))) {
            return(paste0("'fileReferenceGDS' slot must have one character ", 
                    "string representing an existing file."))
        }
      
        ## Validate the fileReferenceAnnotGDS parameter
        if ((!is.null(object@fileReferenceAnnotGDS) && 
                    length(object@fileReferenceAnnotGDS) != 1) || 
            (!is.null(object@fileReferenceAnnotGDS) && 
                    !file.exists(object@fileReferenceAnnotGDS))) {
            return(paste0("'fileReferenceAnnotGDS' slot must have one ", 
                    "character string representing an existing file."))
        }
      
        ## Validate the inferenceType parameter 
        if (length(object@inferenceType) != 1 || 
                !object@inferenceType %in% c("PCAknn", "haploAdmixture")) {
            return(paste0("'inferenceType' slot must be a single character" , 
                " string. The valid options are: 'PCAknn' and ", 
                "'haploAdmixture'."))
        }

        ## Validate the sampleRef parameter TODO

        ## Validate the batch parameter
        if (length(object@batch) != 1 || object@batch < 1) {
            return("'batch' slot must have one positive integer.")
        }

        ## Validate the prefix parameter
        if (length(object@prefix) != 1) {
            return("'prefix' slot must have one character string.")
        }

        ## Validate the nbSim parameter
        if (length(object@nbSim) != 1 || object@nbSim < 1L) {
            return("'nbSim' slot must have one positive integer.")
        }

        ## Validate the offset parameter
        if (length(object@offset) != 1) {
            return("'offset' slot must have one integer.")
        }
        
        ## Validate the minCov parameter
        if (length(object@minCov) != 1 || object@minCov < 1L) {
            return("'minCov' slot must have one positive integer.")
        }

        ## Validate the minProb parameter
        if (length(object@minProb) != 1 || 
                !(object@minProb > 0 && object@minProb < 1)) {
            return(paste0("'minProb' slot must have one positive numeric", 
                        " between 0 and 1."))
        }

        ## Validate the seqError parameter
        if (length(object@seqError) != 1 || 
                !(object@seqError > 0 && object@seqError < 1)) {
            return(paste0("'seqError' slot must have one positive numeric", 
                        " between 0 and 1."))
        }

        ## Validate the seqErrorSyn parameter
        if (length(object@seqErrorSyn) != 1 || 
                !(object@seqErrorSyn > 0 && object@seqErrorSyn < 1)) {
            return(paste0("'seqErrorSyn' slot must have one positive numeric", 
                        " between 0 and 1."))
        }

        ## Validate the pRecomb parameter
        if (length(object@pRecomb) != 1 || 
                !(object@pRecomb > 0 && object@pRecomb < 1)) {
            return(paste0("'pRecomb' slot must have one positive numeric", 
                        " between 0 and 1."))
        }
      
        ## Validate the np parameter
        if (length(object@np) != 1 || 
                !(object@np > 0)) {
            return("'np' slot must have one positive integer.")
        }

        ## Validate the listPos parameter
        if (!(is.null(object@listPos) || (ncol(object@listPos) == 2 && 
                all(colnames(object@listPos) %in% c("snp.chromosome", 
                "snp.position"))))) {
            return(paste0("'listPos' slot must be NULL or a data.frame with ", 
                "2 columns named \"snp.chromosome\" and \"snp.position\"."))
        }

        ## Validate the syntheticRefDF parameter
        if (!(is.null(object@syntheticRefDF) || 
                (ncol(object@syntheticRefDF) >= 3 && 
                    all(c("sample.id", "pop.group", "superPop") %in% 
                            colnames(object@listPos) )))) {
            return(paste0("'syntheticRefDF' slot must be NULL or a ", 
                "data.frame with 3 columns named \"sample.id\", ", 
                "\"pop.group\", and \"superPop\"."))
        }

        ## Validate the pruningMethod parameter
        if (length(object@pruningMethod) != 1 || 
              !object@pruningMethod %in% c("corr", "r", "dprime",
                  "composite")) {
          return(paste0("'pruningMethod' slot must have one character string ", 
              "within those 4 choices: \"corr\", \"r\", \"dprime\", ", 
              "\"composite\"."))
        }

        ## Validate the slideWindowMaxBP parameter TODO
        if (length(object@slideWindowMaxBP) != 1 || 
                            object@slideWindowMaxBP < 0) {
          return(paste0("'slideWindowMaxBP' slot must have one positive ", 
                "integer value."))
        }
      
        ## Validate the thresholdLD parameter
        if (length(object@thresholdLD) != 1 || object@thresholdLD < 0) {
          return("'thresholdLD' slot must have one positive numeric value.")
        }
      
        ## Validate the specificSNV parameter
        if (!(is.null(object@specificSNV) || 
                (ncol(object@specificSNV) >= 2 && 
                    all(c("snp.chromosome", "snp.position") %in% 
                            colnames(object@specificSNV) )))) {
            return(paste0("'specificSNV' slot must be NULL or a ", 
                "data.frame with 2 columns named \"snp.chromosome\"", 
                " and \"snp.position\"."))
        }

        ## Validate the genoType parameter TODO

        ## Validate the phaseType parameter
        if (length(object@phaseType) != 1) {
            return("'phaseType' slot must be one character string.")
        }
        
        ## Validate the phase parameter
        if (length(object@phase) != 1) {
          return("'verbose' slot must have one logical value.")
        }

        ## Validate the PCAmissingRate parameter
        if (length(object@PCAmissingRate) != 1 || object@PCAmissingRate < 0) {
          return("'PCAmissingRate' slot must have one positive numeric value.")
        }
      
        ## Validate the PCAalgorithm parameter
        if (length(object@PCAalgorithm) != 1 || 
                !object@PCAalgorithm %in% c("exact", "randomized")) {
          return("'PCAalgorithm' slot must have one character string. ", 
                "The valid options are: \"exact\" or \"randomized\".")
        }
      
        ## Validate the eigenCount parameter
        if (length(object@eigenCount) != 1) {
          return("'eigenCount' slot must have one integer value.")
        }
      
        ## Validate the eigenCountSyn parameter 
        if (length(object@eigenCountSyn) != 1) {
          return("'eigenCountSyn' slot must have one integer value.")
        }

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
#' @param studyType a single \code{character} string representing the type 
#' of study. The possible choices are: "LD" and "GeneAware". The type of 
#' study affects the way the estimation of the allelic fraction is done. 
#' Default: \code{"LD"}.
#' 
#' @param genoSource a single \code{character} string corresponding to the type 
#' of file with the genotype and the allele information of the profile that 
#' will be provided in the 'profileFile' slot. The valid options are: "VCF", 
#' "generic", "snp-pileup", and "bam". The "generic" format CSV file
#' must have at least these columns: 'Chromosome', 'Position', 'Ref', 'Alt', 
#' 'Count', 'File1R', and 'File1A'. The 'Count' is the depth at the 
#' specified position; 'FileR' is the depth of the reference allele, and
#' 'File1A' is the depth of the specific alternative allele.
#' Finally, in the case of a "VCF" file, the file must have at least those 
#' genotype fields: GT, AD, and DP. Default: \code{NULL}.
#' 
#' @param blockTypeId a single \code{character} string corresponding to 
#' the block type used to extract the block identifiers. The block type must 
#' be present in the GDS Reference Annotation file. 
#' Default: \code{"GeneS.Ensembl.Hsapiens.v86"}.
#' 
#' @param reference a \code{character} string with two possible values:
#' '1KGv1.0', '1k_hgdpV0.1'. It specifies the type of inference. 
#' Default: \code{"1KGv1.0"}.
#' 
#' @param genome a \code{character} string with one possible value:
#' 'HG38'. It specifies the genome used. Default: \code{"HG38"}.
#' 
#' @param chrInfo a \code{vector} of positive \code{integer} values
#' representing the length of the chromosomes. See 'details' section. 
#' If \code{NULL}, the following value is 
#' assigned: \code{seqlengths(Hsapiens)[seq_len(25)]}. Default: \code{NULL}.
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
#' 
#' @param pathProfileGDS  a single \code{character} string 
#' representing the path to the directory where the GDS Profile files will 
#' be created. If specified, the directory must exist. 
#' Default: \code{NULL}.
#' 
#' @param fileReferenceGDS a single \code{character} string representing 
#' the file name of the Reference GDS file. If specified, the file  
#' must exist. Default: \code{NULL}.
#' 
#' @param fileReferenceAnnotGDS a \code{character} string representing the
#' file name of the Population Reference GDS Annotation file. If specified,
#' the file must exist. Default: \code{NULL}.
#' 
#' @param inferenceType a single \code{character} string representing the 
#' genotype ancestry inference method. The two possible values:
#' 'PCAknn' and 'haploAdmixture'. Default: \code{"PCAknn"}.
#' 
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
#' @param pRecomb a single positive \code{numeric} between \code{0} and 
#' \code{1} that represents the frequency of phase switching in the 
#' synthetic profiles. Default: \code{0.01}.
#' 
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#' 
#' @param listPos a \code{data.frame} containing 2 columns named: 
#' "snp.chromosome" and "snp.position". The first column,
#' called "snp.chromosome", contains the name of the chromosome where the
#' SNV is located. The second column, called "snp.position", contains the
#' position of the SNV on the chromosome. Default: \code{NULL}.
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
#' @param pruningMethod a \code{character} string representing the method that 
#' will be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The 4 possible values
#' are: "corr", "r", "dprime", and "composite". Default: \code{"corr"}.
#' 
#' @param slideWindowMaxBP a single positive \code{integer} that represents
#' the maximum basepairs (bp) in the sliding window. This parameter is used
#' for the LD pruning done by the \code{\link[SNPRelate]{snpgdsLDpruning}} 
#' function. Default: \code{500000L}.
#' 
#' @param thresholdLD a single positive \code{numeric} value that represents 
#' the LD threshold used in the \code{\link[SNPRelate]{snpgdsLDpruning}} 
#' function. Default: \code{sqrt(0.1)}.
#' 
#' @param specificSNV a \code{data.frame} containing 2 columns. The first 
#' column, called "snp.chromosome" contains the name of the chromosome where 
#' the SNV is located. The second column, called "snp.position", contains the
#' position of the SNV on the chromosome. Optionally, the column "snvKeep" 
#' contains the SNV index in the reference or -1 if not in the reference. It 
#' is used during the pruning step. Default: \code{NULL}.
#' 
#' @param genoType TODO
#' 
#' @param phaseType a TODO. Default: \code{"phase.ref"}.
#' 
#' @param phase a \code{logical} indicating TODO. Default: \code{FALSE}.
#' 
#' @param PCAmissingRate a positive \code{numeric} representing the maximum 
#' missing rate retained accepted to use SNPs in the PCA analysis done 
#' with the the \link[SNPRelate]{snpgdsPCA} function. If \code{NaN}, no 
#' missing threshold. Default: \code{0.025}.
#' 
#' @param PCAalgorithm a \code{character} string representing the algorithm 
#' to use with the \link[SNPRelate]{snpgdsPCA} function. The algorithm must 
#' be implemented and available to the \link[SNPRelate]{snpgdsPCA} function.
#' Default: \code{"exact"}.
#' 
#' @param eigenCount a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigenCount' <= 0, then all eigenvectors are returned. 
#' Default: \code{32L}.
#' 
#' @param eigenCountSyn a TODO. Default: \code{15L}.
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
#' @details
#'
#' The `chrInfo` parameter contains the length of the chromosomes. The
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
#' @examples
#'
#' ## New object of class "RAIDSparam" with default parameters
#' newParam <- RAIDSparam()
#' 
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @importFrom rlang arg_match
#' @importFrom methods new
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 Hsapiens
#' @importFrom Seqinfo seqlengths
#' @export
RAIDSparam <- function(studyDF=NULL, studyDFSyn=NULL, pedStudy=NULL, 
    studyType=c("LD", "GeneAware"), genoSource=NULL, 
    blockTypeId="GeneS.Ensembl.Hsapiens.v86", reference="1KGv1.0", 
    genome="HG38", chrInfo=NULL, paramAncestry=NULL, profileFile=NULL,
    profileFileGeno=NULL, pathProfileGDS=NULL, fileReferenceGDS=NULL,
    fileReferenceAnnotGDS=NULL, inferenceType="PCAknn", sampleRef=NULL,
    batch=1L, prefix="1", nbSim=1L, offset=-1L, minCov=10L, minProb=0.999,
    seqError=0.001, seqErrorSyn=0.001, pRecomb=0.01, np=1L, listPos=NULL, 
    syntheticRefDF=NULL, pruningMethod=c("corr", "r", "dprime", "composite"), 
    slideWindowMaxBP=500000L, thresholdLD=sqrt(0.1), specificSNV=NULL, 
    genoType="geno.ref", phaseType="phase.ref", phase=FALSE, 
    PCAmissingRate=0.025, PCAalgorithm="exact", eigenCount=32L, 
    eigenCountSyn=15L, kList=seq(2L, 15L, 1L), pcaList=seq(2L, 15L, 1L), 
    fieldPopInRef="superPop", fieldPopInfAnc="superPop", 
    fieldSubPop="pop.group", verbose=FALSE) {
    
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

    ## Assigne default chromosome information whe not assigned by user
    if (is.null(chrInfo)) {
        chrInfo <- seqlengths(Hsapiens)[seq_len(25)]
    }

    new("RAIDSparam", studyDF=studyDF, studyDFSyn=studyDFSyn,
    studyType=studyType, genoSource=genoSource, blockTypeId=blockTypeId,
    reference=reference, genome=genome, chrInfo=chrInfo, 
    paramAncestry=paramAncestry, profileFile=profileFile, 
    profileFileGeno=profileFileGeno,
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