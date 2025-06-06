% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix.R
\name{matrix.cross_table}
\alias{matrix.cross_table}
\title{Creates a named matrix from two variables.}
\usage{
matrix.cross_table(var1, var2, predefined_levels = NA)
}
\arguments{
\item{var1}{Character vector. Values representing first variable.}

\item{var2}{Character vector. Values representing second variable.}

\item{predefined_levels}{Character vector. Levels for var1 and var2. If omitted, the variables just coded according to the levels they have. If contains vector, the variables are coded according to predefined_levels values, the latter could be useful if var1 and var2 contain not all levels of interest}
}
\value{
matrix with cross-tabulation of var1 and var2.
}
\description{
Creates a named matrix from two variables.
}
\details{
Creates a named matrix from two variables.
}
