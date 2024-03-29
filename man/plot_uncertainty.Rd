% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_uncertainty.R
\name{plot_uncertainty}
\alias{plot_uncertainty}
\title{Plot uncertainty of tissue labels from fit_banyan()}
\usage{
plot_uncertainty(fit)
}
\arguments{
\item{fit}{A list returned by fit_banyan(). Must have slots named "coords" and "U_scores".}
}
\value{
A ggplot object
}
\description{
This function allows you to visualize the uncertainty levels inferred cell sub-populations after running fit_banyan() and get_scores()
}
\keyword{Bayesian}
\keyword{Gibbs}
\keyword{MLSBM}
\keyword{SBM}
\keyword{expression}
\keyword{gene}
\keyword{networks}
\keyword{spatial}
