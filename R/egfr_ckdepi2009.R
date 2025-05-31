#' Calculate eGFR based on CKD-EPI 2009 creatinine-based equation
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by CKD-EPI 2009 creatinine-based equation.
#'
#' Reference to the equation: Levey AS, Stevens LA, Schmid CH et al. A New Equation to Estimate Glomerular Filtration Rate. Ann Intern Med 2009;150:604–12.
#'
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param ethnicity Vector. Ethnicity. If no ethnicity will be defined, the calculation will use coefficients for White subjects. Specify ethnicity if a study includes African-American subjects, and define the the values of variable in the parameter label_afroamerican.
#' @param creatinine_units Character string. Units in which serum creatinine is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param label_afroamerican List. Label(s) for Afroamerican ethnicity.
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @param max_age Numeric. Maximal age suitable for the equation application, in years. By default is 100 years, but change this value in case you would like to apply equation to older persons.
#' @return numeric eGFR expressed in ml/min/1.73m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}}.
#' @export
#' @name egfr.ckdepi.cr.2009
#' @examples
#' # for a single patient
#' egfr.ckdepi.cr.2009 (creatinine = 1.4, age = 60, sex = "Male", ethnicity = "White",
#'   creatinine_units = "mg/dl")
#' # for a dataset - see vignettes for details
#' # egfr.ckdepi.cr.2009 (creatinine = dta$scr, age = dta$age, sex = dta$sex,
#' #   ethnicity = dta$race, creatinine_units = "mg/dl")

egfr.ckdepi.cr.2009 <- function(
  # variables for eGFR calculation
  creatinine, age, sex, ethnicity = NA, 
  # creatinine measurement units
  creatinine_units = "micromol/l",
  # labels used in the data set - more explanations are available in the vignette
    # label used to define Afroamerican ethnicity in the dataset
    label_afroamerican = c ("Afroamerican"),
    # label(s) used to define male sex in the dataset
    label_sex_male = c ("Male", 1),
    # label(s) used to define female sex in the dataset
    label_sex_female = c ("Female", 0),
  max_age = 100
) {
  
  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN

  # if ethnicity column is not defined
  if(length(ethnicity) == 1) ethnicity <- rep("none", length(creatinine))  

  if(!service.check_equal_length(creatinine, age, ethnicity, sex)) stop("The length of variables should be equal.")
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("creatinine", "age", "sex") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)

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


  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  creatinine_units <- rep(creatinine_units, length(creatinine))
  
  #
  # convert creatinine units if necessary
  #
  creatinine <- service.convert_creatinine(creatinine, creatinine_units)

  

  # apply first set of coefficients for sex and creatinine
  eGFR <- ifelse( sex %in% label_sex_female,
                ifelse( creatinine <= 0.7,
                       ( (creatinine / 0.7)^(-0.329) ) * (0.993^age),
                       ( (creatinine / 0.7)^(-1.209) ) * (0.993^age)
                ),
                  ifelse( sex %in% label_sex_male,
                    ifelse( creatinine <= 0.9,
                      ( (creatinine / 0.9)^(-0.411) ) * (0.993^age),
                      ( (creatinine / 0.9)^(-1.209) ) * (0.993^age)
                    ),
                    NA # if sex value is not corresponding neither to male nor female labels)
                  )
                )



  # apply second set of coefficients for sex and race
  eGFR <- ifelse( ethnicity %in% label_afroamerican,
                ifelse( sex %in% label_sex_female,
                        eGFR * 166,
                        ifelse( sex %in% label_sex_male,
                                eGFR * 163,
                                NA
                                )
            ),
                ifelse( sex %in% label_sex_female,
                        eGFR * 144,
                        ifelse( sex %in% label_sex_male,
                                eGFR * 141,
                      NA # if sex value is not corresponding neither to male nor female labels)
                              )
                       )
                )

return (round(eGFR, 2))
}


# FUNCTION: END
##################################################################


