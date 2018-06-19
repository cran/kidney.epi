#' Sample dataset with kidney transplant patients.
#'
#' A dataset containing 10 records for kidney transplant patients, including information for deceased donors.
#'
#' @format A data frame with 10 rows and 12 variables:
#' \describe{
#'   \item{ptid}{patient identifier}
#'   \item{rec.age}{age of the recipient, in years}
#'   \item{don.age}{age of the donor, in years}
#'   \item{don.height}{height of the donor, in cm}
#'   \item{don.weight}{weight of the donor, in kg}
#'   \item{don.ethnicity}{ethnicity of the donor}
#'   \item{don.hypertension}{history of hypertension for the donor}
#'   \item{don.diabetes}{history of diabetes for the donor}
#'   \item{don.causeofdeath}{cause of death for the donor}
#'   \item{don.creatinine}{serum creatinine of the donor, in mg/dL}
#'   \item{don.hcv}{hepatitis c virus status of the donor}
#'   \item{don.dcdstatus}{donation after circulatory death status of the donor}
#' }
#' @source Generation from different patients' records
"ktx"