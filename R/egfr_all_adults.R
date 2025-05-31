#' Calculate eGFR by all creatinine-based equations for adults
##################################################################
# FUNCTION: BEGIN
#' @details Calculate eGFR by all creatinine-based equations for adults available in the kidney.epi package.
#'
#' References to the equations are available in single functions of the kidney.epi package.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param ethnicity Vector. Ethnicity. If no ethnicity will be defined, the calculation will use coefficients for Caucasian subjects. Specify ethnicity labels in the function parameter label_african.
#' @param creatinine_units Character string. Units in which serum creatinine is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @param label_african List. Label(s) for African ethnicity. Required only by race-specific equations.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.all_adults.cr
#' @examples
#' # for a single patient
#' egfr.all_adults.cr (creatinine = 1.4, age = 60, sex = "Male", 
#'   creatinine_units = "mg/dl")
#' # for a dataset - see vignettes for details
#' # egfr.all_adults.cr (creatinine = dta$scr, age = dta$age, sex = dta$sex, 
#' #  creatinine_units = "mg/dl")

egfr.all_adults.cr <- function(
  # variables for calculation of eGFR
  creatinine, age, sex, ethnicity = NA,
  # creatinine measurement units
  creatinine_units = "micromol/l",
  # custom labels for factor parameters and their unknown values - more eexplanations are available in the vignette
    # label for definition male sex in data set
    label_sex_male = c ("Male", 1),
    # label for definition female sex in data set
    label_sex_female = c ("Female", 0),
    # label for Afroamerican ethnicity
    label_african = c ("Afroamerican"),
  max_age = 100
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
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


	egfr.ckdepi.cr.2009 <- egfr.ckdepi.cr.2009(creatinine = creatinine, age = age, sex = sex, ethnicity = ethnicity, creatinine_units = creatinine_units, label_sex_male = label_sex_male, label_sex_female = label_sex_female, label_afroamerican = label_african, max_age = max_age)

	egfr.ckdepi.cr.2021 <- egfr.ckdepi.cr.2021(creatinine = creatinine, age = age, sex = sex, creatinine_units = creatinine_units, label_sex_male = label_sex_male, label_sex_female = label_sex_female, max_age = max_age)

	egfr.ekfc.cr <- egfr.ekfc.cr(creatinine = creatinine, age = age, sex = sex, ethnicity = ethnicity, creatinine_units = creatinine_units, label_sex_male = label_sex_male, label_sex_female = label_sex_female, max_age = max_age)

	egfr.fas.cr <- egfr.fas.cr(creatinine = creatinine, age = age, sex = sex, creatinine_units = creatinine_units, label_sex_male = label_sex_male, label_sex_female = label_sex_female, max_age = max_age)

	egfr.mdrd4 <- egfr.mdrd4(creatinine = creatinine, age = age, sex = sex, ethnicity = ethnicity, creatinine_units = creatinine_units, label_sex_male = label_sex_male, label_sex_female = label_sex_female, label_afroamerican = label_african, max_age = max_age)

	egfr.lm.cr <- egfr.lm.cr(creatinine = creatinine, age = age, sex = sex, creatinine_units = creatinine_units, label_sex_male = label_sex_male, label_sex_female = label_sex_female, max_age = max_age)

	egfr.bis.cr <- egfr.bis.cr(creatinine = creatinine, age = age, sex = sex, creatinine_units = creatinine_units, label_sex_male = label_sex_male, label_sex_female = label_sex_female, max_age = max_age)

return (data.frame(
	egfr.ckdepi.cr.2009,
	egfr.ckdepi.cr.2021,
	egfr.ekfc.cr,
	egfr.fas.cr,
	egfr.mdrd4,
	egfr.lm.cr,
	egfr.bis.cr
	))

}


# FUNCTION: END
##################################################################







#' Calculate eGFR by all cystatin-based equations for adults
##################################################################
# FUNCTION: BEGIN
#' @details Calculate eGFR by all cystatin-based equations for adults available in the kidney.epi package.
#'
#' References to the equations are available in single functions of the kidney.epi package.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param cystatin Numeric vector. Serum cystatin, could be expressed in "mg/L" or "nanomol/L". Units of measurement should be defined in variable cystatin_units (if not defined explicitly by user, the default value is "mg/L").
#' @param age Numeric vector. Age, in years.
#' @param cystatin_units Character string. Units in which serum cystatin is expressed. Could be one of the following: "mg/L" or "nanomol/L"
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.all_adults.cys
#' @examples
#' # for a single patient
#' egfr.all_adults.cys (cystatin = 1.4, age = 60,
#'   cystatin_units = "mg/L")
#' # for a dataset - see vignettes for details
#' # egfr.all_adults.cys (cystatin = dta$cys, age = dta$age, 
#' #  cystatin_units = "mg/L")

egfr.all_adults.cys <- function(
  cystatin, age,
  cystatin_units = "mg/L",
  max_age = 100
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("cystatin", "age") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single cystatin_units
  service.check_param_number(cystatin_units)

  cystatin_units <- tolower(cystatin_units)

  # check the range of cystatin_units
  service.check_param_arguments(cystatin_units, c("mg/l", "nanomol/l"))

  # check plausible biologic boundaries
  age <- service.check_and_correct_numeric(age, "age",
    rules = c(
      non_negative = TRUE,
      lower_than = max_age,
      greater_than = 18,
      adult_equation = 18
    )
  )
  cystatin <- service.check_and_correct_numeric(cystatin, "cystatin C")


  # CHECK FUNCTION INPUT: END
  ##################################################################

	egfr.ekfc.cys <- egfr.ekfc.cys(cystatin = cystatin, age = age, cystatin_units = cystatin_units, max_age = max_age)
	egfr.fas.cys <- egfr.fas.cys(cystatin = cystatin, age = age, cystatin_units = cystatin_units, max_age = max_age)

return (data.frame(
	egfr.ekfc.cys,
	egfr.fas.cys
	))
}


# FUNCTION: END
##################################################################





#' Calculate eGFR by all creatinine-cystatin-based equations for adults
##################################################################
# FUNCTION: BEGIN
#' @details Calculate eGFR by all creatinine-cystatin-based equations for adults available in the kidney.epi package.
#'
#' References to the equations are available in single functions of the kidney.epi package.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param creatinine_units Character string. Units in which serum creatinine is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param cystatin Numeric vector. Serum cystatin, could be expressed in "mg/L" or "nanomol/L". Units of measurement should be defined in variable cystatin_units (if not defined explicitly by user, the default value is "mg/L").
#' @param cystatin_units Character string. Units in which serum cystatin is expressed. Could be one of the following: "mg/L" or "nanomol/L"
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.all_adults.cr_cys

egfr.all_adults.cr_cys <- function(
  # variables for calculation of eGFR
  creatinine, cystatin, age, sex,
  # creatinine measurement units
  creatinine_units = "micromol/l", cystatin_units = "mg/L",
  # custom labels for factor parameters and their unknown values - more explanations are available in the vignette
    # label for definition male sex in data set
    label_sex_male = c ("Male", 1),
    # label for definition female sex in data set
    label_sex_female = c ("Female", 0),
  max_age = 100
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("creatinine", "age", "sex") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)
  service.check_param_number(cystatin_units)
  
  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)
  cystatin_units <- tolower(cystatin_units)
  
  # check the range of creatinine_units
  service.check_param_arguments(creatinine_units, c("mg/dl", "micromol/l", "mmol/l"))
  service.check_param_arguments(cystatin_units, c("mg/l", "nanomol/l"))

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
  cystatin <- service.check_and_correct_numeric(cystatin, "cystatin C")

  # CHECK FUNCTION INPUT: END
  ##################################################################

  egfr.ckdepi.cr_cys.2021 <- egfr.ckdepi.cr_cys.2021 (creatinine = creatinine, cystatin = cystatin, age = age, sex = sex, creatinine_units = creatinine_units, cystatin_units = cystatin_units, label_sex_male = label_sex_male, label_sex_female = label_sex_female, max_age = max_age)

  egfr.fas.cr_cys <- egfr.fas.cr_cys (creatinine = creatinine, cystatin = cystatin, age = age, sex = sex, creatinine_units = creatinine_units, cystatin_units = cystatin_units, label_sex_male = label_sex_male, label_sex_female = label_sex_female, max_age = max_age)

  egfr.bis.cr_cys <- egfr.bis.cr_cys(creatinine = creatinine, cystatin = cystatin, age = age, sex = sex, creatinine_units = creatinine_units, cystatin_units = cystatin_units, label_sex_male = label_sex_male, label_sex_female = label_sex_female, max_age = max_age)

return (data.frame(
	egfr.ckdepi.cr_cys.2021,
	egfr.fas.cr_cys,
	egfr.bis.cr_cys
	))
}


# FUNCTION: END
##################################################################

