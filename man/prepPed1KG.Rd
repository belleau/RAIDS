% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process1KG.R
\encoding{UTF-8}
\name{prepPed1KG}
\alias{prepPed1KG}
\title{Prepare the pedigree file using pedigree information from Reference}
\usage{
prepPed1KG(filePed, pathGeno = file.path("data", "sampleGeno"), batch = 0L)
}
\arguments{
\item{filePed}{a \code{character} string representing the path and
file name of the pedigree file (PED file) that contains the information
related to the profiles present in the Reference GDS file. The PED file must
exist.}

\item{pathGeno}{a \code{character} string representing the path where
the Reference genotyping files for each profile are located. Only the
profiles with associated genotyping files are retained in the creation of
the final \code{data.frame}. The name of the genotyping files must
correspond to the individual identification (Individual.ID) in the
pedigree file (PED file).
Default: \code{"./data/sampleGeno"}.}

\item{batch}{a\code{integer} that uniquely identifies the source of the
pedigree information. The Reference is usually \code{0L}.
Default: \code{0L}.}
}
\value{
a \code{data.frame} containing the needed pedigree information
from Reference. The \code{data.frame} contains those columns:
\describe{
\item{sample.id}{ a \code{character} string representing the profile unique
ID.}
\item{Name.ID}{ a \code{character} string representing the profile name.}
\item{sex}{a \code{character} string representing the sex of the profile.}
\item{pop.group}{ a \code{character} string representing the
sub-continental ancestry of the profile.}
\item{superPop }{ a \code{character} string representing the continental
ancestry of the profile.}
\item{superPop }{ a \code{integer} representing the batch of the profile.}
}
}
\description{
Using the pedigree file from Reference, this function extracts
needed information and formats it into a \code{data.frame} so in can
be used in following steps of the ancestry inference process. The
function also requires that the genotyping files associated to each
sample be available in a specified directory.
}
\examples{

## Path to the demo pedigree file is located in this package
dataDir <- system.file("extdata", package="RAIDS")

## Path where the demo genotype CSV files are located
pathGeno <- file.path(dataDir, "demoProfileGenotypes")

## Demo pedigree file
pedDemoFile <- file.path(dataDir, "PedigreeDemo.ped")

## Create a data.frame containing the information of the retained
## samples (samples with existing genotyping files)
prepPed1KG(filePed=pedDemoFile, pathGeno=pathGeno, batch=0L)


}
\author{
Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
}
