#' Expand both end of peak
#'
#' This function is expansion of start and end of peak position
#'
#' @usage fat(x, width.size)
#' @param x Genomic Ranges
#' @param width.size numeric value of expanding length
#' @export
#' @examples
#' data(oct4.gr)
#' ## load(file.path("qugacomp", "data", "oct4.gr.rda"))
#' gr.fat <- fat(oct4.gr, 200)
fat <- function(x, width.size) {
  start(x) <- start(x) - width.size
  end(x)   <- end(x)   + width.size
  return(x)
}

#' Unify strand of Genomic Ranges object
#'
#' This function unifies genomic strand of peaks in Genomic Ranges object
#'
#' @usage unifyStrand(x, strand.name = "+")
#' @param x Genomic Ranges object
#' @param strand.name character of genomic strand ("+" or "-")
#' @export
#' @examples
#' data(oct4.gr)
#' ## load(file.path("qugacomp", "data", "oct4.gr.rda"))
#' gr.unifyStrand <- unifyStrand(oct4.gr, "+")
unifyStrand <- function(x, strand.name = "+") {
  strand(x) <- "+"
  return(x)
}
