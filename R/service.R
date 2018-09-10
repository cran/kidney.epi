#' Service functions for data check which could be applied in any function of the package or externally
#'
#' @details Verifies whether the single value is among the values of the vector. 
#' Function is useful to check whether the argument of the function defined by the user is among the possible arguments recognized inside the function.
#' 
#' Programming: Boris Bikbov \email{boris@@bikbov.ru}.
#'
#' @param param2check Numeric value or character string. The single value to be verified.
#' @param possible_params Vector. The vector of values which contains all possible values.
#' @return logic returns TRUE if argument param2check is foundin possible values possible_params, and FALSE if it is not.
#' @export
#' @name service.is.param_possible
#' @examples
#' possible_params = c("KDPI", " KDRI_Rao", "KDRI_median")
#' service.is.param_possible("KDZO", possible_params) # return FALSE
#' service.is.param_possible("KDPI", possible_params) # return TRUE

##################################################################
# FUNCTION: BEGIN
service.is.param_possible <- function(param2check, possible_params){
  if(param2check %in% possible_params){
    # OK, given value is among possible values
	return (TRUE)
  }else{
    # given value is NOT among possible values
    return (FALSE)
  }
}
# FUNCTION: END
##################################################################



#' Select only numeric values greater than defined threshhold.
#' @details Select only numeric values greater than defined threshhold, and substitute other values with NA. 
#' 
#' Programming: Boris Bikbov \email{boris@@bikbov.ru}.
#'
#' @param x the vector to be checked.
#' @param threshhold numeric the threshhold to compare with.
#' @return numeric returns only numeric values greater than threshhold.
#' @export
#' @name service.strict_to_numeric_threshhold_lower
#' @examples
#' myvals <- c(1, 8, -5, "oggi", NA)
#' # ruturn to myvals2 only numeric values greater than defined threshhold (0 in this case)
#' #    and susbstitute non-numeric or negative values with NA
#' myvals2 <- service.strict_to_numeric_threshhold_lower(myvals, 0)
#' myvals2 # 1, 8, NA, NA, NA
#' 
##################################################################
# FUNCTION: BEGIN
service.strict_to_numeric_threshhold_lower <- function(x, threshhold){
  y <- ifelse( is.numeric(x) & x > threshhold,
              x,
			  NA)
  return (y)
}
# FUNCTION: END
##################################################################



#' Select only numeric values lower than defined threshhold
#' @details Select only numeric values lower than defined threshhold, and substitute other values with NA. 
#' 
#' Programming: Boris Bikbov \email{boris@@bikbov.ru}.
#'
#' @param x the vector to be checked.
#' @param threshhold numeric the threshhold to compare with.
#' @return numeric returns only numeric values lower than threshhold.
#' @export
#' @name service.strict_to_numeric_threshhold_greater
#' @examples
#' myvals <- c(1, 8, -5, "oggi", NA)
#' # ruturn to myvals2 only numeric values lower than threshhold  (3 in this case)
#' #   susbstitute non-numeric or negative values with NA
#' myvals2 <- service.strict_to_numeric_threshhold_greater(myvals, 3)
#' myvals2 # 1, NA, -5, NA, NA
#' 
##################################################################
# FUNCTION: BEGIN
service.strict_to_numeric_threshhold_greater <- function(x, threshhold){
  y <- ifelse( is.numeric(x) & x < threshhold,
              x,
			  NA)
  return (y)
}
# FUNCTION: END
##################################################################




#' Count how many values are less or equal than the defined threshhold.
#' @details Count how many values are less or equal than the defined threshhold. 
#' 
#' Programming: Boris Bikbov \email{boris@@bikbov.ru}.
#'
#' @param x the vector to be checked.
#' @param threshhold numeric the threshhold to compare with.
#' @return numeric returns number of numeric values less or equal to the threshhold.
#' @export
#' @name service.count_lowerequal_threshhold
#' @examples
#' myvals <- c(1, 8, -5, "oggi", NA)
#' myvals2 <- service.count_lowerequal_threshhold(myvals, 0)
#' myvals2 # 1
##################################################################
# FUNCTION: BEGIN
service.count_lowerequal_threshhold <- function(x, threshhold){
  mycounter <- sum (is.numeric(x) & x <= threshhold & !is.na(x))
  return (mycounter)
}
# FUNCTION: END
##################################################################


#' Count how many values are greater than the defined threshhold. 
#' @details Count how many values are greater than the defined threshhold. 
#' 
#' Programming: Boris Bikbov \email{boris@@bikbov.ru}.
#'
#' @param x the vector to be checked.
#' @param threshhold numeric the threshhold to compare with.
#' @return numeric returns number of numeric values greater or equal to the threshhold.
#' @export
#' @name service.count_greater_threshhold
#' @examples
#' myvals <- c(1, 8, -5, "oggi", NA)
#' myvals2 <- service.count_greater_threshhold(myvals, 0)
#' myvals2 # 2
##################################################################
# FUNCTION: BEGIN
service.count_greater_threshhold <- function(x, threshhold){
  mycounter <- sum (is.numeric(x) & x > threshhold & !is.na(x))
  return (mycounter)
}
# FUNCTION: END
##################################################################




#' Check whether a vector is numeric.
#' @details Check whether a vector is numeric. 
#' 
#' Programming: Boris Bikbov \email{boris@@bikbov.ru}.
#'
#' @param x the vector to be checked.
#' @return logic whether vector x is numeric or not.
#' @name service.is_numeric
# @examples
# myvals <- c(1, 8, -5, "oggi", NA)
# service.is_numeric(myvals) # FALSE
##################################################################
# FUNCTION: BEGIN
service.is_numeric <- function(x){
  return (is.numeric(x))
}
# FUNCTION: END
##################################################################




#' Form output message in singular or plural.
#' @details Provide different output for constructing messages in singular or plural. 
#' 
#' Programming: Boris Bikbov \email{boris@@bikbov.ru}.
#'
#' @param x Numeric. The value to be checked (usualy a counter of some variable).
#' @param singular Character string. The value to be returned in case of singular form (usualy a string, but could be any type).
#' @param plural Character string. The value to be returned in case of plural form (usualy a string, but could be any type).
#' @return Character string. Returns a value for constructing messages in singular or plural form.
#' @export
#' @name service.singular_or_plural
#' @examples
#' service.singular_or_plural(1, "This value was", "These values were") # "This value was"
#' service.singular_or_plural(99, "This value was", "These values were") # "These values were"
##################################################################
# FUNCTION: BEGIN
service.singular_or_plural <- function(x, singular, plural){
  if(x == 1){
    return (singular)
  }else{
    return (plural)
  }
}
# FUNCTION: END
##################################################################





#' Produce message for warning or cat
#' @details Produce message that is used by warning or cat in the ktx.kdpi.optn function. 
#' Service function that will not be exported to user, and used only in the ktx.kdpi.optn function.
#' 
#' Programming: Boris Bikbov \email{boris@@bikbov.ru}.
#'
#' @param x Numeric. The value to be checked (usualy a counter of some variable).
#' @param custom_phrase Character string. Custom message to be inserted in the middle of standard message.
#' @param warning_type Character string. The type of message: warning (with substitution to NA) or cat (with leave as is).
#' @return Character string. Returns a phrase.
#' @name service.kdpi.optn.output_message
# @examples
# service.kdpi.optn.output_message(suspiciosly_high, "age >100 years", "NA")
# 
##################################################################
# FUNCTION: BEGIN
service.kdpi.optn.output_message <- function(x, custom_phrase, warning_type){
  if(warning_type == "NA"){
    last_sentence = paste(" ", service.singular_or_plural(x, "This value was", "These values were"), " substituted to NA.", sep = "")
  }else if(warning_type == "as is"){
    last_sentence = paste(" ", service.singular_or_plural(x, "This value was", "These values were"), " kept as is.", sep = "")
  }
  
  whole_phrase = paste("There", service.singular_or_plural(x, " is ", " are "), x, " donor", service.singular_or_plural(x, "", "s"), " with ", custom_phrase, ". ", last_sentence, "\n", sep = "")
  
  return(whole_phrase)
}
# FUNCTION: END
##################################################################

