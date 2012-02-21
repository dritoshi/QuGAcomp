\name{loadBedFile}
\alias{loadBedFile}
\title{Import BED file}
\usage{
  loadBedFile(bed.file, genome.length.file)
}
\arguments{
  \item{bed.file}{bed file of peak}

  \item{genome.length.file}{genome length file. Tab
  delimited file. First column is chromosome name. Second
  column is choromosome size.}
}
\description{
  This function is importer of BED file to Genome Ranges
  object
}
\examples{
\dontrun{
bed.file <- file.path(
##  system.file(package="QuGAcomp"),
  "qugacomp",
  "inst",
  "data",
  "GSM288346_ES_Oct4.mm9.header.bed")
genome.length.file <- file.path(
##  system.file(package="QuGAcomp"),
  "qugacomp",
  "inst",
  "data",
  "mm9.info")
oct4.gr <- loadBedFile(bed.file, genome.length.file)
}
}

