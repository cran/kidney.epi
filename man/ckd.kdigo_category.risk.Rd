% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ckd_classify.R
\name{ckd.kdigo_category.risk}
\alias{ckd.kdigo_category.risk}
\title{Calculate KDIGO risk categories}
\usage{
ckd.kdigo_category.risk(gfr_cat, alb_cat)
}
\arguments{
\item{gfr_cat}{Character vector. eGFR categories coded as G1, G2, G3a, G3b, G4, G5.}

\item{alb_cat}{Character vector. Albuminuria categories coded as A1, A2, A3.}
}
\value{
string risk category.
}
\description{
Calculate KDIGO risk categories
}
\details{
Calculate KDIGO risk of complications categories (Low, Moderate, High, Very high) based on eGFR and albuminuria grades.

Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
}
\examples{
# for a single patient
ckd.kdigo_category.risk(gfr_cat = "G2", alb_cat = "A3")
# for a dataset - see vignettes for details
# ckd.kdigo_category.risk(gfr_cat = dta$gfr_cat, alb_cat = dta$alb_cat)
}
\author{
Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
}
