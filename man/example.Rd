% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example.R
\name{example}
\alias{example}
\title{Get an example}
\usage{
example(
  name = c("blocked", "CPP_sim", "CPP", "nolapse", "recovery", "regression", "simple"),
  type = c(NULL, "jags", "stan"),
  filename = NULL
)
}
\arguments{
\item{name}{The name of the example}

\item{type}{The name of sampler template}

\item{filename}{What filename you want to save (defaults to the name of the exercise in the working directory)}
}
\value{
Saves a file to the working directory (or path from filename)
}
\description{
Get an example
}
\examples{
\dontrun{
example("regression", "jags") # get example for the regression jags
example("nolapse", "jags") # get example for the nolapse stan
example("nolapse", "stan", "exercises/simulation_study.Rmd") # save into exercises directory
}
}
