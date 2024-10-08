% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allelicFraction_internal.R
\encoding{UTF-8}
\name{calcAFMLRNA}
\alias{calcAFMLRNA}
\title{Compute the log likelihood ratio based on the coverage of
each allele in a specific block (gene in the case of RNA-seq)}
\usage{
calcAFMLRNA(snpPosHetero)
}
\arguments{
\item{snpPosHetero}{a \code{data.frame}
containing the SNV information for a specific block (gene if RNA-seq).
The \code{data.frame} must contain those columns:
\describe{
\item{cnt.ref}{ a single \code{integer} representing the coverage for
the reference allele.}
\item{cnt.alt}{ a single \code{integer} representing the coverage for
the alternative allele.}
\item{phase}{ a single \code{integer} indicating the phase of the variant
if known, \code{3} if not known}
}}
}
\value{
a \code{list} for the block with the information
relative to the heterozygotes.
The \code{list} contains:
\describe{
\item{lR}{ a single \code{numeric} representing the sum of the log of
read depth of the lowest depth divide by the total depth of the position
minus of likelihood of the allelic fraction of 0.5.}
\item{aFraction}{ a single \code{numeric} representing the allele
fraction estimation.}
\item{sumAlleleLow}{ a \code{integer} representing the
sum of the allele read depth of the lowest read allele depth}
\item{sumAlleleHigh}{ a \code{integer} representing the
sum of the allele read depth
of the highest read allele depth}
}
}
\description{
This function sums the log of read depth of the lowest depth
divide by the total depth of the position minus of likelihood of the allelic
fraction of 0.5 for a block. If the phase is known, the SNVs in the same
haplotype are grouped together.
}
\examples{

## Loading demo dataset containing SNV information
data(snpPositionDemo)

## Only use a subset of heterozygote SNVs related to one block
subset <- snpPositionDemo[which(snpPositionDemo$block.id == 2750 &
                  snpPositionDemo$hetero), c("cnt.ref", "cnt.alt", "phase")]

result <- RAIDS:::calcAFMLRNA(subset)

head(result)

}
\author{
Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
}
\keyword{internal}
