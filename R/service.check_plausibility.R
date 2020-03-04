#' Service functions for data check on biological plausibility and biochemistry conversion which could be applied in any function of the package or externally
#'
##################################################################
# FUNCTION: BEGIN
#' Check and modify if necessary the age values. 
#' @details Service function which check whether age is in biologically plausible boundaries, shows to user warnings if any, and substitute unplausable values.
#' 
#' @param age Numeric. The value to be checked.
#' @return numeric Vector with controlled values.
#' @export
#' @name service.check_plausibility.age
service.check_plausibility.age <- function (age){
  # logic for this and other similiar functions to check plausible biologic boundaries is the following:
  #   check and inform user whether any values out of boundaries were substituted by NA
  #   after the check change the value to boundaries in the possible range (i.e. age > 0 and < 100)

  suspiciosly_low <- service.count_lowerequal_threshhold(age, 0)
  if(suspiciosly_low > 0) warning(service.output_message(suspiciosly_low, "negative values for age", "NA"))
  suspiciosly_high <- service.count_greater_threshhold(age, 100)
  if(suspiciosly_high > 0) warning(service.output_message(suspiciosly_high, "age >100 years", "NA"))
  age <- service.strict_to_numeric_threshhold_lower(age, 0)
  age <- service.strict_to_numeric_threshhold_greater(age, 100)

  return (age)
}
# FUNCTION: END
##################################################################



##################################################################
# FUNCTION: BEGIN
#' Check and modify if necessary the creatinine values. 
#' @details Service function which check whether creatinine is in biologically plausible boundaries, shows to user warnings if any, and substitute unplausable values.
#' 
#' @param creatinine Numeric. The value to be checked.
#' @return numeric Vector with controlled values.
#' @export
#' @name service.check_plausibility.creatinine
service.check_plausibility.creatinine <- function (creatinine){

  suspiciosly_low <- service.count_lowerequal_threshhold(creatinine, 0)
  if(suspiciosly_low > 0) warning(service.output_message(suspiciosly_low, "creatinine <=0", "NA"))
  creatinine <- service.strict_to_numeric_threshhold_lower(creatinine, 0)

  return (creatinine)
}
# FUNCTION: END
##################################################################





##################################################################
# FUNCTION: BEGIN
#' Convert creatinine values if necessary (depending on the mesurement units). 
#' @details Service function which check mesurement units and convert creatinine values if necessary.
#' 
#' @param creatinine Numeric. The creatinine value from data set.
#' @param creatinine_units Character. Creatinine mesurement units defined by user.
#' @return numeric Vector with converted values.
#' @export
#' @name service.convert_creatinine
service.convert_creatinine <- function (creatinine, creatinine_units){

  creatinine <- ifelse (  creatinine_units == "mg/dl",
              creatinine, # do nothing 
              ifelse(  creatinine_units == "micromol/l",
                  creatinine / 88.4,
                  ifelse(  creatinine_units == "mmol/l",
                    1000 * creatinine / 88.4,
                    NA # if any other undefined units is used, assume NA
                  )
              )
            )

  return (creatinine)
}
# FUNCTION: END
##################################################################

