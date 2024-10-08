% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processStudy.R
\encoding{UTF-8}
\name{computePCARefSample}
\alias{computePCARefSample}
\title{Project specified profile onto PCA axes generated using known
reference profiles}
\usage{
computePCARefSample(
  gdsProfile,
  currentProfile,
  studyIDRef = "Ref.1KG",
  np = 1L,
  algorithm = c("exact", "randomized"),
  eigenCount = 32L,
  missingRate = NaN,
  verbose = FALSE
)
}
\arguments{
\item{gdsProfile}{an object of class \link[gdsfmt]{gds.class},
an opened Profile GDS file.}

\item{currentProfile}{a single \code{character} string representing
the profile identifier.}

\item{studyIDRef}{a single \code{character} string representing the
study identifier.}

\item{np}{a single positive \code{integer} representing the number of CPU
that will be used. Default: \code{1L}.}

\item{algorithm}{a \code{character} string representing the algorithm used
to calculate the PCA. The 2 choices are "exact" (traditional exact
calculation) and "randomized" (fast PCA with randomized algorithm
introduced in Galinsky et al. 2016). Default: \code{"exact"}.}

\item{eigenCount}{a single \code{integer} indicating the number of
eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
Default: \code{32L}.}

\item{missingRate}{a \code{numeric} value representing the threshold
missing rate at with the SNVs are discarded; the SNVs are retained in the
\link[SNPRelate]{snpgdsPCA}
with "<= missingRate" only; if \code{NaN}, no missing threshold.
Default: \code{NaN}.}

\item{verbose}{a \code{logical} indicating if messages should be printed
to show how the different steps in the function. Default: \code{FALSE}.}
}
\value{
a \code{list} containing 3 entries:
\describe{
\item{\code{sample.id}}{ a \code{character} string representing the unique
identifier of the analyzed profile.}
\item{\code{eigenvector.ref}}{ a \code{matrix} of \code{numeric}
representing the eigenvectors of the reference profiles. }
\item{\code{eigenvector}}{ a \code{matrix} of \code{numeric} representing
the eigenvectors of the analyzed profile. }
}
}
\description{
This function generates a PCA using the know reference
profiles. Them, it projects the specified profile onto the PCA axes.
}
\examples{

## Required library
library(gdsfmt)

## Path to the demo Profile GDS file is located in this package
dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")

## Open the Profile GDS file
gdsProfile <- snpgdsOpen(file.path(dataDir, "ex1.gds"))

## Project a profile onto a PCA generated using reference profiles
## The reference profiles come from 1KG
resPCA <- computePCARefSample(gdsProfile=gdsProfile,
    currentProfile=c("ex1"), studyIDRef="Ref.1KG", np=1L, verbose=FALSE)
resPCA$sample.id
resPCA$eigenvector

## Close the GDS files (important)
closefn.gds(gdsProfile)


}
\references{
Galinsky KJ, Bhatia G, Loh PR, Georgiev S, Mukherjee S, Patterson NJ,
Price AL. Fast Principal-Component Analysis Reveals Convergent Evolution
of ADH1B in Europe and East Asia. Am J Hum Genet. 2016 Mar 3;98(3):456-72.
doi: 10.1016/j.ajhg.2015.12.022. Epub 2016 Feb 25.
}
\author{
Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
}
