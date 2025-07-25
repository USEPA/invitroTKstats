% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{initfunction_clint}
\alias{initfunction_clint}
\title{Set Initial Values for Intrinsic Hepatic Clearance (Clint) Bayesian Model}
\usage{
initfunction_clint(mydata, chain)
}
\arguments{
\item{mydata}{(List) Output of \code{build_mydata_clint}.}

\item{chain}{(Numeric) The number of Markov Chains to use.}
}
\value{
A list of initial values.
}
\description{
Sets the initial values of arguments required for JAGS such as assumed initial probability
distributions. The list is used as an argument to JAGS during level-4 processing.
}
