% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/service.check_plausibility.R
\name{service.check_plausibility.creatinine}
\alias{service.check_plausibility.creatinine}
\title{Check and modify if necessary the creatinine values.}
\usage{
service.check_plausibility.creatinine(creatinine)
}
\arguments{
\item{creatinine}{Numeric. The value to be checked.}
}
\value{
numeric Vector with controlled values.
}
\description{
Check and modify if necessary the creatinine values.
}
\details{
Service function which check whether creatinine is in biologically plausible boundaries, shows to user warnings if any, and substitute unplausable values.
}
