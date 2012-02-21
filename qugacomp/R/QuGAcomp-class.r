#' Quantitative Genome Annotation Comparison Tool
#' 
#' This R package can quantitative compare genome annotations such as
#' peaks of two ChIP-seq data.
#' 
#' @name QuGAcomp-package
#' @docType package
#' @aliases qugacomp-package package-qugacomp
#' @references Itoshi NIAIDO and Hiroki R. Ueda., QuGAcpmp: 
#' Quantitative genome annoation comparison Tool for high-throughput
#' sequencing data. (in prep.)
NULL

#' QuGAcomp object of Sox2 and Oct4 ChIP-seq data
#' 
#' This object is Oct4 and Sox2 ChIP-seq data in QuGAcomp object
#'
#' @docType data
#' @keywords datasets
#' @name quga
#' @usage data(quga)
#' @format QuGAcomp object
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11431}
#' @references Xi Chen, et. al. Integration of External Signaling
#'  Pathways with the Core Transcriptional Network in Embryonic Stem 
#'  Cells, Cell, Volume 133, Issue 6, 13 June 2008, Pages 1106-1117, ISSN
#'  0092-8674, 10.1016/j.cell.2008.04.043.
NULL

#' GenomicRanges object of Sox2 ChIP-seq data
#' 
#' This object is Sox2 ChIP-seq data in GenomicRanges object
#'
#' @docType data
#' @keywords datasets
#' @name sox2.gr
#' @usage data(sox2.gr)
#' @format GenomicRanges object
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11431}
#' @references Xi Chen, et. al. Integration of External Signaling
#'  Pathways with the Core Transcriptional Network in Embryonic Stem 
#'  Cells, Cell, Volume 133, Issue 6, 13 June 2008, Pages 1106-1117, ISSN
#'  0092-8674, 10.1016/j.cell.2008.04.043.
NULL

#' GenomicRanges object of Oct4 ChIP-seq data
#' 
#' This object is Oct4 ChIP-seq data in GenomicRanges object
#'
#' @docType data
#' @keywords datasets
#' @name oct4.gr
#' @usage data(oct4.gr)
#' @format GenomicRanges object
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11431}
#' @references Xi Chen, et. al. Integration of External Signaling
#'  Pathways with the Core Transcriptional Network in Embryonic Stem 
#'  Cells, Cell, Volume 133, Issue 6, 13 June 2008, Pages 1106-1117, ISSN
#'  0092-8674, 10.1016/j.cell.2008.04.043.
NULL

setClass(
  "QuGAcomp",
  representation(
    data = "RleList"
  )
)

#' Construct QuGAcomp object
#'
#' This function make QuGAcomp object from Rle objects
#'
#' @usage qugacomp(gr1, gr2)
#' @param gr1 GenomicRanges object
#' @param gr2 GenomicRanges object
#' @export
#' @examples
#' data(oct4.gr)
#' data(sox2.gr)
#'
#' oct4.fat <- fat(oct4.gr, 200)
#' sox2.fat <- fat(sox2.gr, 200)
#'
#' oct4.unistd <- unifyStrand(oct4.fat)
#' sox2.unistd <- unifyStrand(sox2.fat)
#' 
#' oct4.cov <- coverage(oct4.unistd)
#' sox2.cov <- coverage(sox2.unistd)
#' 
#' oct4.bin500 <- lapply( oct4.cov, function(x) rleBinning(x, 500) )
#' sox2.bin500 <- lapply( sox2.cov, function(x) rleBinning(x, 500) )
#' 
#' oct4.bin500 <- flatRleList(oct4.bin500)
#' sox2.bin500 <- flatRleList(sox2.bin500)
#' 
#' quga <- qugacomp(oct4.bin500, sox2.bin500)
qugacomp <- function(gr1, gr2) {
  new("QuGAcomp", data = RleList(gr1, gr2))
}
