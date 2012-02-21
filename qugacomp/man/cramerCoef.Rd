\name{cramerCoef}
\alias{cramerCoef}
\title{Caluculate cramer's coefficient V of cross table}
\arguments{
  \item{quga}{QuGAcomp object}
}
\description{
  This function caluculate cramer's coefficient V of cross
  table (contingency table).
}
\examples{
data(quga)

quga.cramer <- cramerCoef(quga)
}

