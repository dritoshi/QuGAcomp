#' Unlist RleList to Rle object
#'
#' flat RleList to Rle object
#'
#' @usage flatRleList(x)
#' @param x RleList
#' @export
#' @examples
#' data(oct4.gr)
#' ## load(file.path("qugacomp", "data", "oct4.gr.rda"))
#'
#' oct4.fat <- fat(oct4.gr, 200)
#' oct4.unistd <- unifyStrand(oct4.fat)
#' oct4.cov <- coverage(oct4.unistd)
#'
#' bin.size <- 500
#' oct4.bin500 <- lapply( oct4.cov, function(x) rleBinning(x, bin.size) )
#' oct4.bin500.flt <- flatRleList(oct4.bin500)
flatRleList <- function(x) {
  x.names <- names(x)
  joinedRleList <- x[[1]]
  for (i in 2:length(x.names) ) {
    joinedRleList <- append(joinedRleList, x[[ x.names[i] ]])
  }
  return(joinedRleList)
}

#' Binning of Rle data
#'
#' Binning of Rle data
#'
#' @usage rleBinning(x, bin.size)
#' @param x GenomicRanges object
#' @param bin.size bin size (numeric)
#' @export
#' @examples
#' data(oct4.gr)
#' ## load(file.path("qugacomp", "data", "oct4.gr.rda"))
#'
#' oct4.fat <- fat(oct4.gr, 200)
#' oct4.unistd <- unifyStrand(oct4.fat)
#' oct4.cov <- coverage(oct4.unistd)
#'
#' bin.size <- 500
#' oct4.bin500 <- lapply( oct4.cov, function(x) rleBinning(x, bin.size) )
rleBinning <- function(x, bin.size) {
  
  run.length <- runLength(x)
  num.runs   <- length(run.length)
  
  run.ends   <- cumsum(run.length)
  run.starts <- c(1, run.ends+1)[1:num.runs]
  
  bin.num <- ceiling(length(x) / bin.size)
  x.bin <- rep(0, bin.num)

  run.value  <- runValue(x)
  
  for (i in 1:num.runs) {
    
    bin.start <- ceiling(run.starts[i] / bin.size)
    bin.end   <- ceiling(run.ends[i]   / bin.size)

    x.bin[bin.start:bin.end] <- x.bin[bin.start:bin.end] + run.value[i]
    
    #cat(i, run.starts[i], run.ends[i], run.value[i], bin.start, bin.end, "\n")
  }
  x.bin[x.bin > 0] <- 1
  return(Rle(x.bin))
}

## very slow
.rleBinning <- function(x, bin.size) {

    x.length <- length(x)
  
    starts <- seq.int(1,        x.length, bin.size)
    ends   <- seq.int(bin.size, x.length, bin.size)
    ends   <- c(ends, x.length)
    
    aggregate(as.vector(x), start = starts, end = ends, FUN = sum)
}
