#' Calculate eGFR by the revised Lund-Malmö creatinine-based equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by the revised Lund-Malmö creatinine-based equation.
#'
#' Reference to the equation: Björk J, Grubb A, Sterner G, Nyman U. Revised equations for estimating glomerular filtration rate based on the Lund-Malmö Study cohort. Scand J Clin Lab Invest. 2011;71: 232-239.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.lm.cr
#' @examples
#' # for a single patient
#' egfr.lm.cr (creatinine = 1.4, age = 60, sex = "Male", 
#'   creatinine_units = "mg/dl")
#' # for a dataset - see vignettes for details
#' # egfr.lm.cr (creatinine = dta$scr, age = dta$age, sex = dta$sex, 
#' #  creatinine_units = "mg/dl")

egfr.lm.cr <- function(
  # variables for calculation of eGFR
  creatinine, age, sex,
  # creatinine measurement units
  creatinine_units = "micromol/l",
  # custom labels for factor parameters and their unknown values - more explanations are available in the vignette
    # label for definition male sex in data set
    label_sex_male = c("Male", 1),
    # label for definition female sex in data set
    label_sex_female = c("Female", 0),
  max_age = 100
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  if(!service.check_equal_length(creatinine, age, sex)) stop("The length of variables should be equal.")

  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("creatinine", "age", "sex") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)

  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)

  # check the range of creatinine_units
  service.check_param_arguments(creatinine_units, c("mg/dl", "micromol/l", "mmol/l"))

  # Check whether numerical arguments inputed by user are fine (weight, height etc have to be positive numbers)
  service.check_params_numeric(age, creatinine)
  
  
  # check plausible biologic boundaries
  age <- service.check_and_correct_numeric(age, "age",
    rules = c(
      non_negative = TRUE,
      lower_than = max_age,
      greater_than = 18,
      adult_equation = 18
    )
  )
  creatinine <- service.check_and_correct_numeric(creatinine, "creatinine")

  # CHECK FUNCTION INPUT: END
  ##################################################################


  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  creatinine_units <- rep(creatinine_units, length(creatinine))
  
  #
  # convert creatinine units
  # since the equation is developed using micromol/L, use them as a reference
  creatinine <- service.convert_creatinine(creatinine, creatinine_units, creatinine_reference_units = "micromol/l")

  # x coefficient
  x <- ifelse(sex %in% label_sex_female, 
    ifelse(creatinine < 150, 2.5 + 0.0121 * (150 - creatinine), 2.5 - 0.926 * log(creatinine / 150)),
  ifelse(sex %in% label_sex_male, 
    ifelse(creatinine < 180, 2.56 + 0.00968 * (180 - creatinine), 2.56 - 0.926 * log(creatinine / 180)), NA
  )
  )

  eGFR <- exp(x - 0.0158 * age + 0.438 * log(age))

return (round(eGFR, 2))

}


# FUNCTION: END
##################################################################





#' Calculate creatinine clearance by the Cockcroft-Gault equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate Cockcroft-Gault by the Cockcroft-Gault equation.
#'
#' Reference to the equation: Cockcroft, DW, Gault MH. Prediction of creatinine clearance from serum creatinine. Nephron. 1976. 16(1):31-41.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param weight Numeric vector. Weight, kg.
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric creatinine clearance expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.cg.cr
#' @examples
#' # for a single patient
#' egfr.cg.cr (creatinine = 1.4, age = 60, sex = "Male",
#'   weight = 80, creatinine_units = "mg/dl")
#' # for a dataset - see vignettes for details

egfr.cg.cr <- function(
  # variables for calculation of eGFR
  creatinine, age, sex, weight,
  # creatinine measurement units
  creatinine_units = "micromol/l",
  # custom labels for factor parameters and their unknown values - more explanations are available in the vignette
    # label for definition male sex in data set
    label_sex_male = c ("Male", 1),
    # label for definition female sex in data set
    label_sex_female = c ("Female", 0),
  max_age = 100
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN

  if(!service.check_equal_length(creatinine, age, sex, weight)) stop("The length of variables should be equal.")
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("creatinine", "age", "sex", "weight") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)

  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)

  # check the range of creatinine_units
  service.check_param_arguments(creatinine_units, c("mg/dl", "micromol/l", "mmol/l"))

  
  # check plausible biologic boundaries
  age <- service.check_and_correct_numeric(age, "age",
    rules = c(
      non_negative = TRUE,
      lower_than = max_age,
      greater_than = 18,
      adult_equation = 18
    )
  )
  creatinine <- service.check_and_correct_numeric(creatinine, "creatinine")
  weight <- service.check_and_correct_numeric(weight, "weight")

  # CHECK FUNCTION INPUT: END
  ##################################################################


  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  creatinine_units <- rep(creatinine_units, length(creatinine))
  
  #
  # convert creatinine units
  # since the equation is developed using mg/dL, use them as a reference
  creatinine <- service.convert_creatinine(creatinine, creatinine_units, creatinine_reference_units = "mg/dl")

  eGFR <- (140 - age) * weight / (72 * creatinine)
  eGFR <- ifelse(sex %in% label_sex_female, 
    eGFR * 0.85,
    ifelse(sex %in% label_sex_male, 
      eGFR, NA
    )
  )

return (round(eGFR, 2))

}


# FUNCTION: END
##################################################################








#' Calculate eGFR by the Berlin Initiative Study (BIS1) creatinine-based equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by the Berlin Initiative Study (BIS1) creatinine-based equation.
#'
#' Reference to the equation: Schaeffner ES, Ebert N, Delanaye P et al. Two novel equations to estimate kidney function in persons aged 70 years or older. Ann Intern Med 2012; 157: 471–481  doi: 10.7326/0003-4819-157-7-201210020-00003.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.bis.cr
#' @examples
#' # for a single patient
#' egfr.bis.cr (creatinine = 1.4, age = 80, sex = "Male", 
#'   creatinine_units = "mg/dl")
#' # for a dataset - see vignettes for details
#' # egfr.bis.cr (creatinine = dta$scr, age = dta$age, sex = dta$sex, 
#' #  creatinine_units = "mg/dl")

egfr.bis.cr <- function(
  # variables for calculation of eGFR
  creatinine, age, sex,
  # creatinine measurement units
  creatinine_units = "micromol/l",
  # custom labels for factor parameters and their unknown values - more explanations are available in the vignette
    # label for definition male sex in data set
    label_sex_male = c ("Male", 1),
    # label for definition female sex in data set
    label_sex_female = c ("Female", 0),
  max_age = 100
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN

  if(!service.check_equal_length(creatinine, age, sex)) stop("The length of variables should be equal.")
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("creatinine", "age", "sex") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)

  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)

  # check the range of creatinine_units
  service.check_param_arguments(creatinine_units, c("mg/dl", "micromol/l", "mmol/l"))

  # Check whether numerical arguments inputed by user are fine (weight, height etc have to be positive numbers)
  service.check_params_numeric(age, creatinine)
  
  
  # check plausible biologic boundaries
  age <- service.check_and_correct_numeric(age, "age",
    rules = c(
      non_negative = TRUE,
      lower_than = max_age,
      greater_or_equal_than = 70
    )
  )
  creatinine <- service.check_and_correct_numeric(creatinine, "creatinine")

  # CHECK FUNCTION INPUT: END
  ##################################################################


  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  creatinine_units <- rep(creatinine_units, length(creatinine))
  
  #
  # convert creatinine units
  # since the equation is developed using micromol/L, use them as a reference
  creatinine <- service.convert_creatinine(creatinine, creatinine_units, creatinine_reference_units = "mg/dl")

  eGFR <- 3736 * creatinine^-0.87 * age^-0.95
  eGFR <- ifelse(sex %in% label_sex_female, eGFR * 0.82, eGFR)

return (round(eGFR, 2))
}


# FUNCTION: END
##################################################################








#' Calculate eGFR by the Berlin Initiative Study (BIS2) creatinine-cystatin C-based equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by the Berlin Initiative Study (BIS2) creatinine-cystatin C-based equation.
#'
#' Reference to the equation: Schaeffner ES, Ebert N, Delanaye P et al. Two novel equations to estimate kidney function in persons aged 70 years or older. Ann Intern Med 2012; 157: 471–481  doi: 10.7326/0003-4819-157-7-201210020-00003.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param cystatin Numeric vector. Serum cystatin, could be expressed in "mg/L" or "nanomol/L". Units of measurement should be defined in variable cystatin_units (if not defined explicitly by user, the default value is "mg/L").
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param cystatin_units Character string. Units in which serum cystatin is expressed. Could be one of the following: "mg/L" or "nanomol/L"
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.bis.cr_cys
#' @examples
#' # for a single patient
#' egfr.bis.cr_cys (creatinine = 1.4, cystatin = 1.0, age = 80,
#'    sex = "Male", creatinine_units = "mg/dl",
#'    cystatin_units = "mg/L")
#' # for a dataset - see vignettes for details
#' # egfr.bis.cr_cys (creatinine = dta$scr, cystatin = dta$cys, 
#' #   age = dta$age, sex = dta$sex, 
#' #   creatinine_units = "mg/dl", cystatin_units = "mg/l")

egfr.bis.cr_cys <- function(
  # variables for calculation of eGFR
  creatinine, cystatin, age, sex,
  # creatinine measurement units
  creatinine_units = "micromol/l", cystatin_units = "mg/l",
  # custom labels for factor parameters and their unknown values - more explanations are available in the vignette
    # label for definition male sex in data set
    label_sex_male = c ("Male", 1),
    # label for definition female sex in data set
    label_sex_female = c ("Female", 0),
  max_age = 100
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN

  if(!service.check_equal_length(creatinine, cystatin, age, sex)) stop("The length of variables should be equal.")
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("creatinine", "cystatin", "age", "sex") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)

  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)

  # check the range of creatinine_units
  service.check_param_arguments(creatinine_units, c("mg/dl", "micromol/l", "mmol/l"))

  # check that user defined a single cystatin_units
  service.check_param_number(cystatin_units)

  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  cystatin_units <- tolower(cystatin_units)

  # check the range of cystatin_units
  service.check_param_arguments(cystatin_units, c("mg/l", "nanomol/l") )
  
  # check plausible biologic boundaries
  age <- service.check_and_correct_numeric(age, "age",
    rules = c(
      non_negative = TRUE,
      lower_than = max_age,
      greater_or_equal_than = 70
    )
  )
  creatinine <- service.check_and_correct_numeric(creatinine, "creatinine")
  cystatin <- service.check_and_correct_numeric(cystatin, "cystatin C")

  # CHECK FUNCTION INPUT: END
  ##################################################################


  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  creatinine_units <- rep(creatinine_units, length(creatinine))
  cystatin_units <- rep(cystatin_units, length(cystatin))
 
  #
  # convert creatinine units
  # since the equation is developed using micromol/L, use them as a reference
  creatinine <- service.convert_creatinine(creatinine, creatinine_units, creatinine_reference_units = "mg/dl")
  cystatin <- service.convert_cystatin(cystatin, cystatin_units)

  eGFR <- 767 * cystatin^-0.61 * creatinine^-0.40 * age^-0.57
  eGFR <- ifelse(sex %in% label_sex_female, eGFR * 0.87, eGFR)

return (round(eGFR, 2))
}


# FUNCTION: END
##################################################################



