% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RAIDS.R
\docType{data}
\name{pedSynthetic}
\alias{pedSynthetic}
\title{A small \code{data.frame} containing the information related to
synthetic profiles. The ancestry of the profiles used to generate the
synthetic profiles must be present.}
\format{
The \code{data.frame} containing the information about the
synthetic profiles.  The row names of
the \code{data.frame} correspond to the profile unique identifiers.
The \code{data.frame} contains 7 columns:
\describe{
\item{\code{data.id}}{ a \code{character} string representing the unique
synthetic profile identifier.}
\item{\code{case.id}}{ a \code{character} string representing the unique
profile identifier that was used to generate the synthetic profile.}
\item{\code{sample.type}}{ a \code{character} string representing the type
of profile. }
\item{\code{diagnosis}}{ a \code{character} string representing the
diagnosis of profile that was used to generate the synthetic profile. }
\item{\code{source}}{ a \code{character} string representing the
source of the synthetic profile. }
\item{\code{study.id}}{ a \code{character} string representing the
name of the study to which the synthetic profile is associated. }
\item{\code{superPop}}{ a \code{character} string representing the
super population of the profile that was used to generate the synthetic
profile. }
}
}
\usage{
data(pedSynthetic)
}
\value{
The \code{data.frame} containing the information about the
synthetic profiles. The row names of
the \code{data.frame} correspond to the profile unique identifiers.
The \code{data.frame} contains 7 columns:
\describe{
\item{\code{data.id}}{ a \code{character} string representing the unique
synthetic profile identifier.}
\item{\code{case.id}}{ a \code{character} string representing the unique
profile identifier that was used to generate the synthetic profile.}
\item{\code{sample.type}}{ a \code{character} string representing the type
of profile.}
\item{\code{diagnosis}}{ a \code{character} string representing the
diagnosis of profile that was used to generate the synthetic profile. }
\item{\code{source}}{ a \code{character} string representing the
source of the synthetic profile. }
\item{\code{study.id}}{ a \code{character} string representing the
name of the study to which the synthetic profile is associated. }
\item{\code{superPop}}{ a \code{character} string representing the
super population of the profile that was used to generate the synthetic
profile. }
}
}
\description{
The object is a \code{data.frame} with 7 columns. The row names of
the \code{data.frame} must be the profile unique identifiers.
}
\details{
This dataset can be
used to test the \code{\link{computeSyntheticROC}} function.
}
\examples{

## Loading demo dataset containing pedigree information for synthetic
## profiles
data(pedSynthetic)

## Loading demo dataset containing the inferred ancestry results
## for the synthetic data
data(matKNNSynthetic)

## Retain one K and one D value
matKNN <- matKNNSynthetic[matKNNSynthetic$D == 5 & matKNNSynthetic$K == 4, ]

## Compile statistics from the
## synthetic profiles for fixed values of D and K
results <- RAIDS:::computeSyntheticROC(matKNN=matKNN,
    matKNNAncestryColumn="SuperPop",
    pedCall=pedSynthetic, pedCallAncestryColumn="superPop",
    listCall=c("EAS", "EUR", "AFR", "AMR", "SAS"))

results$matAUROC.All
results$matAUROC.Call
results$listROC.Call


}
\seealso{
\describe{
\item{\code{\link{computeSyntheticROC}}}{ for calculating the AUROC of
the inferences for specific values of D and K using the inferred
ancestry results from the synthetic profiles}
}
}
\keyword{datasets}
