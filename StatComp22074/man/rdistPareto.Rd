% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Setting.R
\name{rdistPareto}
\alias{rdistPareto}
\title{Return random numbers of Pareto(a,b) distributes}
\usage{
rdistPareto(n = 1000, a = 2, b = 2, picture = FALSE)
}
\arguments{
\item{n}{random sample size}

\item{a}{parameters a}

\item{b}{parameters b}

\item{picture}{plot histogram}
}
\value{
random numbers of Pareto(a,b) distributes
}
\description{
Return random numbers of Pareto(a,b) distributes
}
\examples{
\dontrun{
x <- rdistPareto(5,2,2,FALSE)
print(x)
}

}
