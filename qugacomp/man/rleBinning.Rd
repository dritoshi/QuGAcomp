\name{rleBinning}
\alias{rleBinning}
\title{Binning of Rle data}
\usage{
  rleBinning(x, bin.size)
}
\arguments{
  \item{x}{GenomicRanges object}

  \item{bin.size}{bin size (numeric)}
}
\description{
  Binning of Rle data
}
\examples{
data(oct4.gr)
## load(file.path("qugacomp", "data", "oct4.gr.rda"))

oct4.fat <- fat(oct4.gr, 200)
oct4.unistd <- unifyStrand(oct4.fat)
oct4.cov <- coverage(oct4.unistd)

bin.size <- 500
oct4.bin500 <- lapply( oct4.cov, function(x) rleBinning(x, bin.size) )
}

