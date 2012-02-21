#' Import BED file
#'
#' This function is importer of BED file to Genome Ranges object
#'
#' @usage loadBedFile(bed.file, genome.length.file)
#' @param bed.file bed file of peak
#' @param genome.length.file genome length file. Tab delimited file.
#' First column is chromosome name. Second column is choromosome size.
#' @export
#' @examples \dontrun{
#' bed.file <- file.path(
#' ##  system.file(package="QuGAcomp"),
#'   "qugacomp",
#'   "inst",
#'   "data",
#'   "GSM288346_ES_Oct4.mm9.header.bed")
#' genome.length.file <- file.path(
#' ##  system.file(package="QuGAcomp"),
#'   "qugacomp",
#'   "inst",
#'   "data",
#'   "mm9.info")
#' oct4.gr <- loadBedFile(bed.file, genome.length.file) 
#' }
loadBedFile <- function(bed.file, genome.length.file) {
  x <- read.table(bed.file, skip=1)

  # score
  if (class(try(x[,5])) == "try-error") {
    x <- cbind(x, rep(0, nrow(x)) )
  }
  # strand
  if (class(try(x[,6])) == "try-error") {
    x <- cbind(x, as.factor(rep("+", nrow(x))) )
  }

  print( dim(x) )

  gr <- GRanges(
    seqnames = Rle(x[,1]),
    ranges   = IRanges(x[,2], end = x[,3], names = x[,4]),
    strand   = x[,6],
    score    = x[,5]
  )
  
  # load genome length file
  genome.length <- loadGenomeLength(genome.length.file)
  seqlengths(gr) <- getGenomeLength(gr, genome.length)

  return(gr)
}

loadGenomeLength <- function(file) {
  mm9.length <- read.table(file, row.names=1)
  return(mm9.length) 
}

getGenomeLength <- function(x, genome.length)  {
  genome.length[names(seqlengths(x)),] 
}
