#' @title Generate a VCF with the information from the SNPs that pass
#' a cut-off threshold
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the 1KG GDS file.
#'
#' @param fileOut a \code{character} string representing the path and file
#' name of the VCF file that will be created wit the retained SNP information.
#' The file should have the ".vcf" extension.
#'
#' @param offset a single \code{integer} that is added to the SNP position to
#' switch from 0-based to 1-based coordinate when needed (or reverse).
#' Default: \code{0L}.
#'
#' @param freqCutoff a single positive \code{numeric} specifying the cut-off to
#' keep a SNP. If \code{NULL}, all SNPs are retained. Default: \code{NULL}.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Demo 1KG Reference GDS file
#' fileGDS <- openfn.gds(file.path(dataDir,
#'                     "PopulationReferenceDemo.gds"))
#'
#' ## Output VCF file that will be created (temporary)
#' vcfFile <- file.path(tempdir(), "Demo_TMP_01.vcf")
#'
#' ## Create a VCF file with the SNV dataset present in the GDS file
#' ## No cutoff on frequency, so all SNVs are saved
#' snvListVCF(gdsReference=fileGDS, fileOut=vcfFile, offset=0L,
#'                     freqCutoff=NULL)
#'
#' ## Close GDS file (IMPORTANT)
#' closefn.gds(fileGDS)
#'
#' ## Remove temporary VCF file
#' unlink(vcfFile, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn ls.gdsn
#' @importFrom methods is
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
snvListVCF <- function(gdsReference, fileOut, offset=0L, freqCutoff=NULL) {

    ## Validate that gdsReference is an object of class gds.class
    if (!inherits(gdsReference, "gds.class")) {
        stop("The \'gdsReference\' must be an object of class \'gds.class\'.")
    }

    ## Validate that offset is a single integer
    if (! isSingleNumber(offset)) {
        stop("The \'offset\' must be a single integer.")
    }

    ## Validate that freqCutoff is a single numeric or NULL
    if (! isSingleNumber(freqCutoff) && ! is.null(freqCutoff)) {
        stop("The \'freqCutoff\' must be a single numeric or NULL.")
    }

    snpChromosome <- read.gdsn(index.gdsn(gdsReference, "snp.chromosome"))
    snpPosition <- read.gdsn(index.gdsn(gdsReference, "snp.position"))
    snpAllele <- read.gdsn(index.gdsn(gdsReference, "snp.allele"))

    allele <- matrix(unlist(strsplit(snpAllele, "\\/")), nrow=2)

    df <- NULL

    if(is.null(freqCutoff)){
        snpAF <- read.gdsn(index.gdsn(gdsReference, "snp.AF"))
        df <- data.frame(CHROM=snpChromosome,
                            POS=as.integer(snpPosition + offset),
                            ID=rep(".", length(snpChromosome)),
                            REF=allele[1,],
                            ALT=allele[2,],
                            QUAL=rep(".", length(snpChromosome)),
                            FILTER=rep(".", length(snpChromosome)),
                            INFO=paste0("AF=", snpAF),
                            stringsAsFactors=FALSE)
    } else {
        if( length(which(paste0("snp.",
                                c("EAS", "EUR",
                                  "AFR", "AMR",
                                  "SAS"),
                                "_AF") %in% ls.gdsn(gdsReference))) == 5){
            freqDF <- data.frame(
                    snp.AF=read.gdsn(index.gdsn(gdsReference, "snp.AF")),
                    snp.EAS_AF=read.gdsn(index.gdsn(gdsReference, "snp.EAS_AF")),
                    snp.EUR_AF=read.gdsn(index.gdsn(gdsReference, "snp.EUR_AF")),
                    snp.AFR_AF=read.gdsn(index.gdsn(gdsReference, "snp.AFR_AF")),
                    snp.AMR_AF=read.gdsn(index.gdsn(gdsReference, "snp.AMR_AF")),
                    snp.SAS_AF=read.gdsn(index.gdsn(gdsReference, "snp.SAS_AF")))

            listKeep <- which(rowSums(freqDF[,2:6] >= freqCutoff &
                                            freqDF[,2:6] <= 1 - freqCutoff) > 0)
            df <- data.frame(CHROM=snpChromosome[listKeep],
                                POS=as.integer(snpPosition[listKeep] + offset),
                                ID=rep(".", length(listKeep)),
                                REF=allele[1,listKeep],
                                ALT=allele[2,listKeep],
                                QUAL=rep(".", length(listKeep)),
                                FILTER=rep(".", length(listKeep)),
                                INFO=paste0("AF=", freqDF$snp.AF[listKeep]),
                                stringsAsFactors=FALSE)
        }else if("AF.superPop" %in% ls.gdsn(gdsReference)){

            matAF <- index.gdsn(gdsReference, "AF.superPop")
            snpAF <- read.gdsn(matAF)
            freqDF <- cbind(read.gdsn(index.gdsn(gdsReference, "snp.AF")),snpAF)
            rm(snpAF)
            listKeep <- which(rowSums(freqDF[,2:ncol(freqDF)] >= freqCutoff &
                                          freqDF[,2:ncol(freqDF)] <= 1 - freqCutoff) > 0)
            df <- data.frame(CHROM=snpChromosome[listKeep],
                             POS=as.integer(snpPosition[listKeep] + offset),
                             ID=rep(".", length(listKeep)),
                             REF=allele[1,listKeep],
                             ALT=allele[2,listKeep],
                             QUAL=rep(".", length(listKeep)),
                             FILTER=rep(".", length(listKeep)),
                             INFO=paste0("AF=", freqDF[listKeep,1]),
                             stringsAsFactors=FALSE)
        }
    }
    df <- df[df$REF %in% c("A", "C", "G", "T") & df$ALT %in% c("A", "C", "G", "T"),]
    df$CHROM <- paste0("chr", df$CHROM)

    ## Add the header
    ##fileformat=VCFv4.3
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency
    ##                            in the range (0,1)">
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

    cat(paste0('##fileformat=VCFv4.3', "\n"), file=fileOut)
    cat(paste0('##FILTER=<ID=PASS,Description="All filters passed">',
                "\n"), file=fileOut, append=TRUE)
    cat(paste0('##INFO=<ID=AF,Number=A,Type=Float,',
                'Description="Estimated allele frequency in the range (0,1)">',
                "\n"), file=fileOut, append=TRUE)
    cat('#', file=fileOut, append=TRUE)

    write.table(df, file=fileOut, sep="\t", append=TRUE, row.names=FALSE,
                    col.names=TRUE, quote=FALSE)

    ## Successful
    return(0L)
}

#' @title Generate a VCF with the information from the SNPs that pass
#' a cut-off threshold
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param pRAIDS a TODO.
#'
#' @param fileOut a \code{character} string representing the path and file
#' name of the VCF file that will be created wit the retained SNP information.
#' The file should have the ".vcf" extension.
#' 
#' @param freqCutoff a single positive \code{numeric} specifying the cut-off to
#' keep a SNP. If \code{NULL}, all SNPs are retained. Default: \code{NULL}.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn ls.gdsn openfn.gds closefn.gds
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
snvListVCFRef <- function(pRAIDS, fileOut, freqCutoff=NULL) {
    # gdsReference, fileOut, offset=0L, freqCutoff=NULL
    ## Validate that gdsReference is an object of class gds.class
    # if (!inherits(gdsReference, "gds.class")) {
    #     stop("The \'gdsReference\' must be an object of class \'gds.class\'.")
    # }
    #
    # ## Validate that offset is a single integer
    # if (! isSingleNumber(offset)) {
    #     stop("The \'offset\' must be a single integer.")
    # }
    #
    # ## Validate that freqCutoff is a single numeric or NULL
    # if (! isSingleNumber(freqCutoff) && ! is.null(freqCutoff)) {
    #     stop("The \'freqCutoff\' must be a single numeric or NULL.")
    # }
    # if reference HGDP1kg no offset

    gdsReference <- openfn.gds(pRAIDS$fileReferenceGDS)

    snpChromosome <- read.gdsn(index.gdsn(gdsReference, "snp.chromosome"))
    snpPosition <- read.gdsn(index.gdsn(gdsReference, "snp.position"))
    snpAllele <- read.gdsn(index.gdsn(gdsReference, "snp.allele"))
    listSNPKepp <- NULL
    allele <- matrix(unlist(strsplit(snpAllele, "\\/")), nrow=2)
    if("snp.KeepDefault" %in% ls.gdsn(gdsReference)){
        listSNPKepp <- read.gdsn(index.gdsn(gdsReference, "snp.KeepDefault"))
        snpChromosome <- snpChromosome[listSNPKepp]
        snpPosition <- snpPosition[listSNPKepp]
        allele <- allele[,listSNPKepp]
    }

    df <- NULL

    if(is.null(freqCutoff)){
        snpAF <- read.gdsn(index.gdsn(gdsReference, "snp.AF"))
        if(! is.null(listSNPKepp)){
            snpAF <- snpAF[listSNPKepp]
        }
        df <- data.frame(CHROM=snpChromosome,
                         POS=as.integer(snpPosition + pRAIDS$offset),
                         ID=rep(".", length(snpChromosome)),
                         REF=allele[1,],
                         ALT=allele[2,],
                         QUAL=rep(".", length(snpChromosome)),
                         FILTER=rep(".", length(snpChromosome)),
                         INFO=paste0("AF=", snpAF),
                         stringsAsFactors=FALSE)
    } else {
        if( length(which(paste0("snp.",
                                c("EAS", "EUR",
                                  "AFR", "AMR",
                                  "SAS"),
                                "_AF") %in% ls.gdsn(gdsReference))) == 5){
            freqDF <- data.frame(
                snp.AF=read.gdsn(index.gdsn(gdsReference, "snp.AF")),
                snp.EAS_AF=read.gdsn(index.gdsn(gdsReference, "snp.EAS_AF")),
                snp.EUR_AF=read.gdsn(index.gdsn(gdsReference, "snp.EUR_AF")),
                snp.AFR_AF=read.gdsn(index.gdsn(gdsReference, "snp.AFR_AF")),
                snp.AMR_AF=read.gdsn(index.gdsn(gdsReference, "snp.AMR_AF")),
                snp.SAS_AF=read.gdsn(index.gdsn(gdsReference, "snp.SAS_AF")))

            listKeep <- which(rowSums(freqDF[,2:6] >= freqCutoff &
                                          freqDF[,2:6] <= 1 - freqCutoff) > 0)
            df <- data.frame(CHROM=snpChromosome[listKeep],
                             POS=as.integer(snpPosition[listKeep] + pRAIDS$offset),
                             ID=rep(".", length(listKeep)),
                             REF=allele[1,listKeep],
                             ALT=allele[2,listKeep],
                             QUAL=rep(".", length(listKeep)),
                             FILTER=rep(".", length(listKeep)),
                             INFO=paste0("AF=", freqDF$snp.AF[listKeep]),
                             stringsAsFactors=FALSE)
        }else if("AF.superPop" %in% ls.gdsn(gdsReference)){

            matAF <- index.gdsn(gdsReference, "AF.superPop")
            snpAF <- read.gdsn(matAF)
            freqDF <- cbind(read.gdsn(index.gdsn(gdsReference, "snp.AF")),snpAF)
            rm(snpAF)
            listKeep <- which(rowSums(freqDF[,2:ncol(freqDF)] >= freqCutoff &
                                          freqDF[,2:ncol(freqDF)] <= 1 - freqCutoff) > 0)
            df <- data.frame(CHROM=snpChromosome[listKeep],
                             POS=as.integer(snpPosition[listKeep] + pRAIDS$offset),
                             ID=rep(".", length(listKeep)),
                             REF=allele[1,listKeep],
                             ALT=allele[2,listKeep],
                             QUAL=rep(".", length(listKeep)),
                             FILTER=rep(".", length(listKeep)),
                             INFO=paste0("AF=", freqDF[listKeep,1]),
                             stringsAsFactors=FALSE)
        }
    }
    df <- df[df$REF %in% c("A", "C", "G", "T") & df$ALT %in% c("A", "C", "G", "T"),]
    df$CHROM <- paste0("chr", df$CHROM)

    closefn.gds(gdsReference)

    ## Add the header
    ##fileformat=VCFv4.3
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency
    ##                            in the range (0,1)">
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

    cat(paste0('##fileformat=VCFv4.3', "\n"), file=fileOut)
    cat(paste0('##FILTER=<ID=PASS,Description="All filters passed">',
               "\n"), file=fileOut, append=TRUE)
    cat(paste0('##INFO=<ID=AF,Number=A,Type=Float,',
               'Description="Estimated allele frequency in the range (0,1)">',
               "\n"), file=fileOut, append=TRUE)
    cat('#', file=fileOut, append=TRUE)

    write.table(df, file=fileOut, sep="\t", append=TRUE, row.names=FALSE,
                col.names=TRUE, quote=FALSE)

    ## Successful
    return(0L)
}

#' @title Merge the genotyping files per chromosome into one file
#'
#' @description This function merge all the genotyping files associated to one
#' specific sample into one file. That merged VCF file will be saved in a
#' specified directory and will have the name of the sample. It will also be
#' compressed (bzip). The function will merge the
#' files for all samples present in the input directory.
#'
#' @param pathGenoChr a \code{character} string representing the path where
#' the genotyping files for each sample and chromosome are located. The path
#' must contains sub-directories (one per chromosome) and the genotyping files
#' must be present in those sub-directories.
#' The path must exists.
#'
#' @param pathOut a \code{character} string representing the path where
#' the merged genotyping files for each sample will be created.
#' The path must exists.
#'
#' @return The integer \code{0L} when successful or \code{FALSE} if not.
#'
#' @examples
#'
#' ## Path to the demo vcf files in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#' pathGenoTar <- file.path(dataDir, "demoGenoChr", "demoGenoChr.tar")
#'
#' ## Path where the chromosomes files will be located
#' pathGeno <- file.path(tempdir(), "tempGeno")
#' dir.create(pathGeno, showWarnings=FALSE)
#'
#' ## Untar the file that contains the VCF files for 3 samples split by
#' ## chromosome (one directory per chromosome)
#' untar(tarfile=pathGenoTar, exdir=pathGeno)
#'
#' ## Path where the output VCF file will be created is
#' ## the same where the split VCF are (pathGeno)
#'
#' ## The files must not exist
#' if (!file.exists(file.path(pathGeno, "NA12003.csv.bz2")) &&
#'         !file.exists(file.path(pathGeno, "NA12004.csv.bz2")) &&
#'         !file.exists(file.path(pathGeno, "NA12005.csv.bz2"))) {
#'
#'         ## Return 0 when successful
#'         ## The files "NA12003.csv.bz2", "NA12004.csv.bz2" and
#'         ## "NA12005.csv.bz2" should not be present in the current directory
#'         groupChr1KGSNV(pathGenoChr=pathGeno, pathOut=pathGeno)
#'
#'         ## Validate that files have been created
#'         file.exists(file.path(pathGeno, "NA12003.csv.bz2"))
#'         file.exists(file.path(pathGeno, "NA12004.csv.bz2"))
#'         file.exists(file.path(pathGeno, "NA12005.csv.bz2"))
#'
#' }
#'
#' ## Remove temporary directory
#' unlink(pathGeno, recursive=TRUE, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom utils write.csv2 read.csv2
#' @encoding UTF-8
#' @export
groupChr1KGSNV <- function(pathGenoChr, pathOut) {

    ## Validate that the input path for the genotyping files exists
    if (! file.exists(pathGenoChr)) {
        stop("The path \'", pathGenoChr, "\' does not exist.")
    }

    ## Validate that the output path for the genotyping files exists
    if (! file.exists(pathOut)) {
        stop("The path \'", pathOut, "\' does not exist.")
    }

    ## Obtain the comprehensive list of samples
    listFiles <- dir(file.path(pathGenoChr, "chr1"), ".+\\.chr1\\.vcf\\.bz2")
    listSamples <- gsub("\\.chr1\\.vcf\\.bz2", "", listFiles)


    ## Merge files associated to each samples into one csv file
    results <- lapply(X=listSamples, FUN=function(sampleId, pathOut) {

        ## For each chromosome, read genotyping file and append the information
        listGeno <- lapply(seq_len(22), function(chr, sampleId) {
            geno <- read.csv2(file.path(pathGenoChr, paste0("chr", chr),
                        paste0(sampleId, ".chr", chr,".vcf.bz2")),
                            sep="\t", row.names=NULL)
            return(geno)}, sampleId=sampleId)

        genoAll <- do.call(rbind, listGeno)

        ## Save the genotyping information into one file
        write.csv2(genoAll, file=bzfile(file.path(pathOut,
            paste0(sampleId, ".csv.bz2"))), row.names=FALSE)

        return(TRUE)}, pathOut=pathOut)

    ## Successful or not
    return(ifelse(all(unlist(results)), 0L, FALSE))
}


#' @title Generate a VCF with the information from the SNPs that pass
#' a cut-off threshold
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param freqCutoff a single positive \code{numeric} specifying the cut-off to
#' keep a SNP. If \code{NULL}, all SNPs are retained. Default: \code{NULL}.
#'
#' @param chr a \code{integer} between 1 to 22
#' chr1:22 specify SNP retained are on chr.
#' If \code{NULL}, all SNPs are retained. Default: \code{NULL}.
#'
#' @param fileOut a \code{character} string representing the path and
#' the file name of the new VCF.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Demo 1KG Reference GDS file
#' fileGDS <- openfn.gds(file.path(dataDir,
#'                     "PopulationReferenceDemo.gds"))
#'
#' ## Output VCF file that will be created (temporary)
#' vcfFile <- file.path(tempdir(), "Demo_TMP_01.vcf")
#'
#' ## Create a VCF file with the SNV dataset present in the GDS file
#' ## No cutoff on frequency, so all SNVs are saved
#' snvListVCF(gdsReference=fileGDS, fileOut=vcfFile, offset=0L,
#'                     freqCutoff=NULL)
#'
#' ## Close GDS file (IMPORTANT)
#' closefn.gds(fileGDS)
#'
#' ## Remove temporary VCF file
#' unlink(vcfFile, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn
#' @importFrom methods is
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
generateVCF <- function(freqCutoff=NULL, chr=NULL , fileOut, pRAIDS) {
    # gdsReference, fileOut, offset=0L,
    gdsReference <- snpgdsOpen(filename=pRAIDS$fileReferenceGDS)



    ## Validate that freqCutoff is a single numeric or NULL
    if (! isSingleNumber(freqCutoff) && ! is.null(freqCutoff)) {
        stop("The \'freqCutoff\' must be a single numeric or NULL.")
    }
    if(! is.null(chr)){
        chr <- as.integer(gsub("chr", "", chr))
        if(is.na(chr)){
            stop(paste0("The \'chr\':", chr, " is not an integer or chr1-chr22\n"))
        }
    }

    snpChromosome <- read.gdsn(index.gdsn(gdsReference, "snp.chromosome"))
    snpPosition <- read.gdsn(index.gdsn(gdsReference, "snp.position"))
    snpAllele <- read.gdsn(index.gdsn(gdsReference, "snp.allele"))

    allele <- matrix(unlist(strsplit(snpAllele, "\\/")), nrow=2)
    listSNP <- seq_len(length(snpPosition))
    if(!is.null(chr)){
        listSNP <- listSNP[snpChromosome == chr]
    }
    if(! is.null(freqCutoff)){
        snpAF <- read.gdsn(index.gdsn(gdsReference, "snp.AF"))[listSNP]
        listSNP <- listSNP[snpAF > freqCutoff &  snpAF < 1-freqCutoff]
    }
    closefn.gds(gdsReference)
    snpVCF <- data.frame(CHROM=snpChromosome[listSNP],
                         # not - the offset because if you load the vcf is + offset
                         POS=as.integer(snpPosition[listSNP] - pRAIDS$offset),
                         ID = rep(".", length(listSNP)),
                         REF=allele[1,listSNP],
                         ALT=allele[2,listSNP],
                         QUAL = rep(".", length(listSNP)),
                         FILTER = rep("PASS", length(listSNP)),
                         INFO = rep(".", length(listSNP)),
                         FORMAT = rep("GT", length(listSNP)),
                         stringsAsFactors = FALSE)

    ## Profile GDS file name
    fileProfileGDS <-  validateProfileGDSExist(pathProfile=pRAIDS$pathProfileGDS,
                                               profile=pRAIDS$pedStudy$Name.ID[1])
    gdsProfile <- openfn.gds(fileProfileGDS)

    listGeno <- list()
    g <- read.gdsn(index.gdsn(gdsProfile, "geno.ref"))[listSNP]
    closefn.gds(gdsProfile)
    geno <- rep("./.", length(listSNP))
    geno[which(g==0)] <- "0/0"
    geno[which(g==1)] <- "0/1"
    geno[which(g==2)] <- "1/1"
    snpVCF$QUAL[which(g == 3)] <- "."

    ## Add the header
    ##fileformat=VCFv4.3
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency
    ##                            in the range (0,1)">
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

    cat(paste0('##fileformat=VCFv4.3', "\n"), file = fileOut)
    cat(paste0('##FILTER=<ID=PASS,Description="All filters passed">',
               "\n"), file = fileOut, append=TRUE)
    cat(paste0('##INFO=<ID=GT,Number=1,Type=String,',
               'Description="Phased Genotype">',
               "\n"), file = fileOut, append=TRUE)
    listColName <- c('#CHROM', "POS", "ID",
                     "REF", "ALT", "QUAL",
                     "FILTER", "INFO","FORMAT",
                     pRAIDS$pedStudy$Name.ID[1])
    cat(paste0(paste(listColName,
                     collapse = "\t"), "\n"),
        file = fileOut, append=TRUE)

    snpgVCF <- cbind(snpVCF, geno)

    snpgVCF$CHROM <- paste0("chr", snpgVCF$CHROM)
    write.table(x = snpgVCF, file = fileOut,
                quote = FALSE, sep = "\t",
                col.names = FALSE, row.names = FALSE,
                append = TRUE)


    ## Successful
    return(0L)
}

#' @title Generate a plink bim file
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param fileOut a \code{character} string representing the path and
#' the file name of the new bim file.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Demo 1KG Reference GDS file
#' fileGDS <- openfn.gds(file.path(dataDir,
#'                     "PopulationReferenceDemo.gds"))
#'
#' ## Output VCF file that will be created (temporary)
#' vcfFile <- file.path(tempdir(), "Demo_TMP_01.vcf")
#'
#' ## Create a VCF file with the SNV dataset present in the GDS file
#' ## No cutoff on frequency, so all SNVs are saved
#' snvListVCF(gdsReference=fileGDS, fileOut=vcfFile, offset=0L,
#'                     freqCutoff=NULL)
#'
#' ## Close GDS file (IMPORTANT)
#' closefn.gds(fileGDS)
#'
#' ## Remove temporary VCF file
#' unlink(vcfFile, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn
#' @importFrom genio write_bim
#' @importFrom methods is
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
writeBimPruned <- function(fileOut, pRAIDS){
    fileGDSProofile <- file.path(pRAIDS$pathProfileGDS, paste0(pRAIDS$pedStudy$Name.ID[1], ".gds"))
    gdsProfile <- openfn.gds(fileGDSProofile)

    snpInfo <- data.frame(chr = read.gdsn(index.gdsn(gdsProfile, "snp.chromosome")),
                          id = read.gdsn(index.gdsn(gdsProfile, "snp.id")),
                          posg = 0,
                          pos = read.gdsn(index.gdsn(gdsProfile, "snp.position")),
                          alt = "A",
                          ref = "B",
                          stringsAsFactors = FALSE)
    closefn.gds(gdsProfile)
    write_bim(fileOut, snpInfo)
    return(0)
}


#' @title Generate a plink bed file TODO
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param fileOut a \code{character} string representing the path and
#' the file name of the new bim file.
#'
#' @param listProfile a \code{vector} of \code{character} representing the list
#' of profile from gdsProfile to keep
#'
#' @param listRM a \code{list} TODO. Default: \code{NULL}.
#' 
#' @param profileOnly a \code{boolean} TODO. Default: \code{FALSE}.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn
#' @importFrom genio write_bed
#' @importFrom methods is
#' @importFrom SNPRelate snpgdsGetGeno snpgdsClose snpgdsOpen
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
writeBedProfile <- function(fileOut, listProfile, listRM=NULL, 
    profileOnly=FALSE, pRAIDS) {
    
    fileGDSProfile <- file.path(pRAIDS$pathProfileGDS, 
                                    paste0(pRAIDS$pedStudy$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(fileGDSProfile)
    snpId <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))
    g <- NULL
    if(! profileOnly){
        gdsReference <- snpgdsOpen(pRAIDS$fileReferenceGDS)

        sampleId <- read.gdsn(index.gdsn(gdsReference,
                    "sample.id"))[
                        read.gdsn(index.gdsn(gdsReference, "sample.ref")) == 1]
        if(! is.null(listRM)){
            sampleId <- sampleId[! sampleId %in% listRM]
        }
        g <- snpgdsGetGeno(gdsReference, sample.id=sampleId, snp.id=snpId, 
                                snpfirstdim=TRUE)
        snpgdsClose(gdsReference)
    }
    #studyAn <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))
    gS <- snpgdsGetGeno(gdsProfile, sample.id=listProfile, snp.id=snpId, 
                                snpfirstdim=TRUE)
    if( profileOnly){
        g<- gS
    }else{
        g <- cbind(g,gS)
    }

    write_bed(fileOut, g, verbose=TRUE, append=FALSE)
    snpgdsClose(gdsProfile)

    return(0)
}


#' @title Generate a plink bed file TODO
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param pathOut TODO
#' 
#' @param fileP a \code{character} string representing the path and
#' the file name of the new bim file.
#'
#' @param listProfile a \code{vector} of \code{character} representing the list
#' of profiles from gdsProfile to keep.
#'
#' @param indexS a \code{integer} TODO. Default: \code{1}.
#' 
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom genio write_bed write_bim
#' @importFrom methods is
#' @importFrom SNPRelate snpgdsGetGeno snpgdsOpen snpgdsClose
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
writeBedBimFilesProfileFilter <- function(pathOut, fileP, listProfile, 
    indexS=1, pRAIDS) {
    
    fileGDSProfile <- file.path(pRAIDS$pathProfileGDS, 
                                    paste0(pRAIDS$pedStudy$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(fileGDSProfile)
    snpId <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))
    g <- NULL

    #studyAn <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))
    gS <- snpgdsGetGeno(gdsProfile, sample.id=listProfile, snp.id=snpId, 
                            snpfirstdim=TRUE)
    nbS <- ncol(gS)
    list2RM <- which(rowSums(is.na(gS))/nbS > 0.25)
    gS <- gS[-1 * list2RM,]

    snpInfo <- data.frame(chr=read.gdsn(index.gdsn(gdsProfile,
                                                "snp.chromosome")),
                          id=read.gdsn(index.gdsn(gdsProfile, "snp.id")),
                          posg=0,
                          pos=read.gdsn(index.gdsn(gdsProfile, "snp.position")),
                          alt="A",
                          ref="B",
                          stringsAsFactors=FALSE)
    snpInfo <- snpInfo[-1 * list2RM,]
    write_bim(file.path(pathOut, paste0(pRAIDS$pedStudy$Name.ID[1], ".syn.", 
                    indexS,".sv.bim")), snpInfo)
    write_bed(file.path(pathOut, paste0(pRAIDS$pedStudy$Name.ID[1], ".syn.", 
                    indexS,".sv.bed")), gS, verbose=TRUE, append=FALSE)

    snpgdsClose(gdsProfile)

    resP <- read.delim(fileP, sep=" ", header=FALSE)
    resP <- resP[-1 * list2RM,]
    # 068ba2ae-288c-446d-8d17-72445ce4f788.syn.1.sv.5.P.in
    write.table(resP, file.path(pathOut, paste0(pRAIDS$pedStudy$Name.ID[1], 
                                        ".syn.", indexS,".sv.5.P.in")),
                col.names=FALSE, row.names=FALSE, quote=FALSE, sep=" ")
    return(0)
}


#' @title Generate a plink fam file
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param fileOut a \code{character} string representing the path and
#' the file name of the new fam file.
#'
#' @param listProfile a \code{vector} of \code{character} representing the list
#' of profile from gdsProfile to keep.
#'
#' @param listRM a TODO. Default: \code{NULL}.
#' 
#' @param profileOnly a \code{integer} TODO. Default: \code{FALSE}.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn
#' @importFrom genio write_fam
#' @importFrom methods is
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
writeFamProfile <- function(fileOut, listProfile, listRM=NULL, 
    profileOnly=FALSE, pRAIDS){
    
    fileGDSProfile <- file.path(pRAIDS$pathProfileGDS, 
                            paste0(pRAIDS$pedStudy$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(fileGDSProfile)
    sampleId <- NULL
    if(! profileOnly){
        gdsReference <- snpgdsOpen(pRAIDS$fileReferenceGDS)

        sampleId <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[read.gdsn(index.gdsn(gdsReference, "sample.ref")) == 1]
        if(! is.null(listRM)){
            sampleId <- sampleId[! sampleId %in% listRM ]
        }
        snpgdsClose(gdsReference)
    }

    sujetFam <- data.frame(fam=seq_len(length(sampleId) + length(listProfile)),
                           id=c(sampleId, listProfile),
                           pat=0,
                           mat=0,
                           sex=1,
                           pheno=0,
                           stringsAsFactors=FALSE)

    write_fam(fileOut, sujetFam)
    snpgdsClose(gdsProfile)

    return(0)
}

#' @title Generate a matrix of synthetic group per pop
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param pRAIDS a \code{parametersRAIDS} object with all the RAIDS
#' parameters
#'
#' @return a \code{matrix} containing the sample identifiers and where
#' each column is the name of a subcontinental population. The number of
#' row corresponds to the number of samples for each subcontinental population.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## TODO example
#' 
#' ## Path to the demo pedigree file is located in this package
#' ## <- system.file("extdata", package="RAIDS")
#'
#' ## Demo 1KG Reference GDS file
#' ##fileGDS <- openfn.gds(file.path(dataDir,
#' ##                    "PopulationReferenceDemo.gds"))
#'
#' ## Output VCF file that will be created (temporary)
#' ##vcfFile <- file.path(tempdir(), "Demo_TMP_01.vcf")
#'
#' ## Create a VCF file with the SNV dataset present in the GDS file
#' ## No cutoff on frequency, so all SNVs are saved
#' ##snvListVCF(gdsReference=fileGDS, fileOut=vcfFile, offset=0L,
#' ##                    freqCutoff=NULL)
#'
#' ## Close GDS file (IMPORTANT)
#' ##closefn.gds(fileGDS)
#'
#' ## Remove temporary VCF file
#' ##unlink(vcfFile, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn
#' @importFrom genio write_bed
#' @importFrom methods is
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
getMatrixPopSynthetic <- function(pRAIDS) {
    fileGDSProfile <- file.path(pRAIDS$pathProfileGDS, 
                                paste0(pRAIDS$pedStudy$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(fileGDSProfile)
    gdsReference <- snpgdsOpen(pRAIDS$fileReferenceGDS)
    tmp <- which(read.gdsn(index.gdsn(gdsReference, "sample.ref")) == 1)
    sampleId <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[tmp]

    sampleAnnot <- read.gdsn(index.gdsn(gdsReference, "sample.annot"))[tmp,]
    # studySyn <- read.gdsn(index.gdsn(gdsProfile,""))
    # pRAIDS$studyDFSyn$study.id[1]
    profileSyn <- read.gdsn(index.gdsn(gdsProfile,"study.annot"))
    profileSyn <- profileSyn[profileSyn$study.id == 
                                        pRAIDS$studyDFSyn$study.id[1],]

    tmp <- which(sampleId %in% profileSyn$case.id)
    listPop <- unique(sampleAnnot[tmp, pRAIDS$fieldSubPop])
    listGr <- lapply(listPop,
                    FUN=function(x,sampleSyn, sampleAnnotSyn,fieldSubPop){
                        return(sampleSyn[sampleAnnotSyn[,fieldSubPop] == x])
                },
                sampleSyn=sampleId[tmp],
                sampleAnnotSyn=sampleAnnot[tmp,],
                fieldSubPop=pRAIDS$fieldSubPop)
    matGr <- do.call(cbind, listGr)
    snpgdsClose(gdsProfile)
    snpgdsClose(gdsReference)

    return(matGr)
}

#' @title Generate a matrix of synthetic group per pop
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param matGr a \code{matrix} containing the sample identifiers and where
#' each column is the name of a subcontinental population. The number of
#' row corresponds to the number of samples for each subcontinental population.
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return a \code{matrix} containing the data.id for synthetic data and where
#' each column is the name of a subcontinental population. The number of
#' row corresponds to the number of samples for each subcontinental population.
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @encoding UTF-8
#' @export
getMatrixDataId <- function(matGr, pRAIDS){
    fileGDSProfile <- file.path(pRAIDS$pathProfileGDS, 
                            paste0(pRAIDS$pedStudy$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(fileGDSProfile)

    profileSyn <- read.gdsn(index.gdsn(gdsProfile,"study.annot"))
    profileSyn <- profileSyn[profileSyn$study.id == 
                                        pRAIDS$studyDFSyn$study.id[1],]


    listGr <- lapply(seq_len(nrow(matGr)),
                     FUN=function(x,matGr, profileSyn){
                         return(profileSyn$data.id[profileSyn$case.id %in% 
                                                        matGr[x,]])
                     },
                     matGr=matGr,
                     profileSyn=profileSyn)
    matGrDataId <- do.call(rbind, listGr)

    snpgdsClose(gdsProfile)

    return(matGrDataId)
}

#' @title Generate a plink bed file
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param fileOut a \code{character} string representing the path and
#' the file name of the new bim file.
#'
#' @param listRM a \code{vector} of \code{character} strings containing the
#' identifiers for the reference samples that need to be removed for the
#' PCA analysis.
#'
#' @param subset a \code{logical} extract the genotypes from the reference in
#' listRM (TRUE) or extract from the reference not in listRM (FALSE).
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom genio write_bed
#' @importFrom SNPRelate snpgdsOpen snpgdsClose snpgdsGetGeno
#' @encoding UTF-8
#' @export
writeBedRef <- function(fileOut, listRM=NULL, subset=FALSE, pRAIDS){
    fileGDSProfile <- file.path(pRAIDS$pathProfileGDS, 
                            paste0(pRAIDS$pedStudy$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(fileGDSProfile)
    gdsReference <- snpgdsOpen(pRAIDS$fileReferenceGDS)

    snpId <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))
    sampleId <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[
                        read.gdsn(index.gdsn(gdsReference, "sample.ref")) == 1]
    if(subset){
        if(! is.null(listRM)){
            sampleId <- sampleId[sampleId %in% listRM]
        }
    } else{
        if(! is.null(listRM)){
            sampleId <- sampleId[! sampleId %in% listRM]
        }
    }
    g <- snpgdsGetGeno(gdsReference, sample.id=sampleId, snp.id=snpId, 
            snpfirstdim=TRUE)

    write_bed(fileOut, g, verbose=TRUE, append=FALSE)
    snpgdsClose(gdsProfile)
    snpgdsClose(gdsReference)
    return(0)
}

#' @title Generate a plink fam file
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param fileOut a \code{character} string representing the path and
#' the file name of the new fam file.
#'
#' @param listRM a \code{vector} of \code{character} strings containing the
#' identifiers for the reference samples that need to be removed for the
#' PCA analysis.
#'
#' @param subset a \code{logical} extract the profiles from the reference in
#' listRM (TRUE) or extract from the reference not in listRM (FALSE).
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom genio write_fam
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @encoding UTF-8
#' @export
writeFamRef <- function(fileOut, listRM=NULL, subset=FALSE, pRAIDS){

    gdsReference <- snpgdsOpen(pRAIDS$fileReferenceGDS)

    sampleId <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[read.gdsn(index.gdsn(gdsReference, "sample.ref")) == 1]
    if(subset){
        if(! is.null(listRM)){
            sampleId <- sampleId[sampleId %in% listRM ]
        }
    } else{
        if(! is.null(listRM)){
            sampleId <- sampleId[! sampleId %in% listRM ]
        }
    }

    sujetFam <- data.frame(fam=seq_len(length(sampleId)),
                           id=c(sampleId),
                           pat=0,
                           mat=0,
                           sex=1,
                           pheno=0,
                           stringsAsFactors=FALSE)

    write_fam(fileOut, sujetFam)
    snpgdsClose(gdsReference)
    return(0)
}

