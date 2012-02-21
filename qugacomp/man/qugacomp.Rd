\name{qugacomp}
\alias{qugacomp}
\title{Construct QuGAcomp object}
\usage{
  qugacomp(gr1, gr2)
}
\arguments{
  \item{gr1}{GenomicRanges object}

  \item{gr2}{GenomicRanges object}
}
\description{
  This function make QuGAcomp object from Rle objects
}
\examples{
data(oct4.gr)
data(sox2.gr)

oct4.fat <- fat(oct4.gr, 200)
sox2.fat <- fat(sox2.gr, 200)

oct4.unistd <- unifyStrand(oct4.fat)
sox2.unistd <- unifyStrand(sox2.fat)

oct4.cov <- coverage(oct4.unistd)
sox2.cov <- coverage(sox2.unistd)

oct4.bin500 <- lapply( oct4.cov, function(x) rleBinning(x, 500) )
sox2.bin500 <- lapply( sox2.cov, function(x) rleBinning(x, 500) )

oct4.bin500 <- flatRleList(oct4.bin500)
sox2.bin500 <- flatRleList(sox2.bin500)

quga <- qugacomp(oct4.bin500, sox2.bin500)
}

