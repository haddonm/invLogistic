% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funIL.R
\name{invlog}
\alias{invlog}
\title{invlog calculates the mean growth increment for the IL}
\usage{
invlog(p, Lt)
}
\arguments{
\item{p}{a vector of parameters with at least c(MaxDL,L50,L95)}

\item{Lt}{the initial lengths, Lt, for which to calculate the
predicted growth increment}
}
\value{
a vector of growth increments the same length as x
}
\description{
invlog is the simplest equation dealing with the
    inverse logistic equation
}
\examples{
 data(midg)
 pars <- initpars(midg$Lt,midg$DL)
 invlog(pars,sort(midg$Lt[1:10]))
}
