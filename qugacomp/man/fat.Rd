\name{fat}
\alias{fat}
\title{Expand both end of peak}
\usage{
  fat(x, width.size)
}
\arguments{
  \item{x}{Genomic Ranges}

  \item{width.size}{numeric value of expanding length}
}
\description{
  This function is expansion of start and end of peak
  position
}
\examples{
data(oct4.gr)
## load(file.path("qugacomp", "data", "oct4.gr.rda"))
gr.fat <- fat(oct4.gr, 200)
}

