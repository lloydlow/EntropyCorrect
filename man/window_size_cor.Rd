% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/window_size_cor.R
\name{window_size_cor}
\alias{window_size_cor}
\title{window_size_cor function}
\usage{
window_size_cor(MSA3, n, randomise)
}
\arguments{
\item{MSA3}{The input is a multiple sequence alignment in FASTA format.}

\item{n}{By default, n = 1, which means only a single character per alignment
column per sequence. For example, users may change to n = 9 for nonamer 
scoring instead of monomer.}

\item{randomise}{This parameter allows users to try other corrected entropy
values. By default, set.seed is used in entropy correction and hence, same
values are produced, which made it easier for the empiral study to determine
the effects of correction in various sample and population sizes.}
}
\description{
This function is one of the functions that makes EntropyMSA works.
}
\examples{
#see EntropyMSA function
}
\keyword{entropy}

