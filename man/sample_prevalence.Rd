% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_prevalence.R
\name{sample_prevalence}
\alias{sample_prevalence}
\title{Sample Prevalence}
\usage{
sample_prevalence(prevalence, reporting_prob)
}
\arguments{
\item{prevalence}{data frame of prevalence per time.}

\item{reporting_prob}{proportion of cases observed.}
}
\value{
data frame of observed prevalence per day.
}
\description{
Generates a partially observed prevalence.
}
\examples{
sample_prevalence(prevalence = prev, reporting_prob = 0.2)
}
