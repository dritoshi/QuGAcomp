\name{flatRleList}
\alias{flatRleList}
\title{Unlist RleList to Rle object}
\usage{
  flatRleList(x)
}
\arguments{
  \item{x}{RleList}
}
\description{
  flat RleList to Rle object
}
\examples{
data(oct4.gr)
## load(file.path("qugacomp", "data", "oct4.gr.rda"))

oct4.fat <- fat(oct4.gr, 200)
oct4.unistd <- unifyStrand(oct4.fat)
oct4.cov <- coverage(oct4.unistd)

bin.size <- 500
oct4.bin500 <- lapply( oct4.cov, function(x) rleBinning(x, bin.size) )
oct4.bin500.flt <- flatRleList(oct4.bin500)
}

