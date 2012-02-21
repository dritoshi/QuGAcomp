\name{phiCoef}
\alias{phiCoef}
\title{Caluculate Phi coefficient of cross table}
\arguments{
  \item{quga}{QuGAcomp object}
}
\description{
  This function caluculate Phi coefficient of cross table
  (contingency table). (= mean square contingency
  coefficient) (= Pearson correlation coefficient for two
  binary variables)
}
\examples{
data(quga)

quga.phi <- phiCoef(quga)
}

