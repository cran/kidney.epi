% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/service.general.R
\name{service.strict_to_numeric_threshhold_lower}
\alias{service.strict_to_numeric_threshhold_lower}
\title{Select only numeric values greater than defined threshhold.}
\usage{
service.strict_to_numeric_threshhold_lower(x, threshhold)
}
\arguments{
\item{x}{the vector to be checked.}

\item{threshhold}{numeric the threshhold to compare with.}
}
\value{
numeric returns only numeric values greater than threshhold.
}
\description{
Select only numeric values greater than defined threshhold.
}
\details{
Select only numeric values greater than defined threshhold, and substitute other values with NA. 

Programming: Boris Bikbov \email{boris@bikbov.ru}.
}
\examples{
myvals <- c(1, 8, -5, "oggi", NA)
# ruturn to myvals2 only numeric values greater than defined threshhold (0 in this case)
#    and susbstitute non-numeric or negative values with NA
myvals2 <- service.strict_to_numeric_threshhold_lower(myvals, 0)
myvals2 # 1, 8, NA, NA, NA

}
