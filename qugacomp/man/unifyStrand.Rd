\name{unifyStrand}
\alias{unifyStrand}
\title{Unify strand of Genomic Ranges object}
\usage{
  unifyStrand(x, strand.name = "+")
}
\arguments{
  \item{x}{Genomic Ranges object}

  \item{strand.name}{character of genomic strand ("+" or
  "-")}
}
\description{
  This function unifies genomic strand of peaks in Genomic
  Ranges object
}
\examples{
data(oct4.gr)
## load(file.path("qugacomp", "data", "oct4.gr.rda"))
gr.unifyStrand <- unifyStrand(oct4.gr, "+")
}

