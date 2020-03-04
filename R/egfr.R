#' Calculate estimated glomerular filtration rate (eGFR) by different equations
#'
#' MDRD
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by MDRD equation.
#' Reference to the equation: Levey AS, Coresh J, Greene T, et al. Using standardized serum creatinine values in the modification of diet in renal disease study equation for estimating glomerular filtration rate. Annals of Internal Medicine 2006;145:247–54.
#' 
#' Programming: Boris Bikbov \email{boris@@bikbov.ru}.
#' Citation: Bikbov B. R open source programming code for calculation of the Kidney Donor Profile Index and Kidney Donor Risk Index. Kidney Diseases, 2018. DOI: 10.1159/000492427 (citation for the whole package)
#'
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param ethnicity Vector. Ethnicity, specify in case of African-American patients. The value of variable refers to the parameter label_afroamerican.
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param creatinine_method Character string. Creatinine standartisation method in a laboratory. Could be either "IDMS" or "non-IDMS". If not explicitly defined by user, the default assumption is "non-IDMS".
#' @param label_afroamerican List. Label(s) for Afroamerican ethnicity.
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @return numeric eGFR expressed in ml/min/1.73m<sup>2</sup>.
#' @export
#' @name egfr.mdrd4
#' @examples
#' egfr.mdrd4 (creatinine = 1.4, age = 60, sex = "Male", ethnicity = "White",
#'   creatinine_units = "mg/dl")

egfr.mdrd4 <- function (
              # variables for calculation of eGFR
              creatinine, age, sex, ethnicity, 
              # creatinine measurement units
              creatinine_units = "micromol/l",
              # creatinine standartisation method in a laboratory: "IDMS" in case of isotope dilution mass spectrometry standartisation, "non-IDMS" in other cases
              creatinine_method = "non-IDMS",
              # custom labels for factor parameters and their unknown values
                  # introduce variables' labels which any user can adapt to their own labeling in the data file
                  # all labels could be a character, numeric, vector, or list - whatever you prefer
                  # you have to change the labels according your data file
                  #  for example:
                  #    if males in your data file defined as "Male", you have to change below the variable label_sex_male = "Male";
                  #    if you use 1 for definition of male sex, you have to change below the variable label_sex_male = 1;
                  #    if males in your data file defined either as "Male" or as "MALE" or as "Males" or as "M", you have to change below the variable label_sex_male = c ("Male", "MALE", "Males", "M"), but in this case it is worth to normalize your data dictionary;
                  #
                  # label for Afroamerican ethnicity
                  label_afroamerican = c ("Afroamerican"),
                  # label for definition male sex in data set
                  label_sex_male = c ("Male", 1),
                  # label for definition female sex in data set
                  label_sex_female = c ("Female", 0)
              # end of function parameters definition
              ){
  

  
  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("creatinine", "age", "ethnicity", "sex") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)
  # check that user defined a single creatinine_method
  service.check_param_number(creatinine_method)

  # check the range of creatinine_units
  service.check_param_arguments(creatinine_units, c("mg/dl", "micromol/l", "mmol/l") )
  # check the range of creatinine_method
  # it is completely out of logic and has no explanations, but check for creatinine_method produce error in any case, while the check for creatinine_units works fine
  service.check_param_arguments(creatinine_method, c("IDMS", "non-IDMS") )

  # Check whether numerical arguments inputed by user are fine (weight, height etc have to be positive numbers)
  dta <- data.frame(cbind(age, creatinine))
  service.check_params_numeric(age, creatinine)
  
  
  # check plausible biologic boundaries (by functions in the service.check_plausibility.R):
  #   check and inform user whether any values out of boundaries were substituted by NA
  #   after the check change the value to boundaries in the possible range (i.e. age > 0 and < 100)
  # age
  # first: general check and tidy: age <0 OR age >100
  service.check_plausibility.age(age)
  # second: since this eGFR equation was developed and validated for adults only, notify user if any children were found, and exclude them from calculation
  suspiciosly_low <- service.count_lowerequal_threshhold(age, 17)
  if(suspiciosly_low > 0) cat(service.output_message(suspiciosly_low, "age <=17 years", "NA"))
  age <- service.strict_to_numeric_threshhold_lower(age, 17)
  # creatinine
  service.check_plausibility.creatinine(creatinine)

  # CHECK FUNCTION INPUT: END
  ##################################################################


  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)
  creatinine_units <- rep(creatinine_units, length(creatinine))  
  # repeat creatinine_method
  creatinine_method <- tolower(creatinine_method)
  creatinine_method <- rep(creatinine_method, length(creatinine))  
  
  #
  # convert creatinine units if necessary
  #
  creatinine <- service.convert_creatinine(creatinine, creatinine_units)

  # start eGFR calculation
  eGFR <- creatinine^(-1.154)*age^(-0.203)
  

  # apply sex coefficients
  eGFR <- ifelse( sex %in% label_sex_male,
                eGFR, # males
                ifelse( sex %in% label_sex_female,
                        eGFR * 0.742, # females
                        NA # if sex value is not corresponding neither to male nor female labels)
                       )
                )

  # apply ethnicity coefficient
  eGFR <- ifelse(  ethnicity %in% label_afroamerican,
                eGFR * 1.212,
                eGFR
                )

  # apply creatinine standartisation method coefficient
  eGFR <- ifelse(  creatinine_method == "idms",
                eGFR * 175,
                eGFR * 186
                )


return (round(eGFR, 2))

}


# FUNCTION: END
##################################################################










#' CKD-EPI
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by CKD-EPI equation
#' Reference to the equation: Levey AS, Stevens LA, Schmid CH et al. A New Equation to Estimate Glomerular Filtration Rate. Ann Intern Med 2009;150:604–12.
#' 
#' Programming: Boris Bikbov \email{boris@@bikbov.ru}.
#' Citation: Bikbov B. R open source programming code for calculation of the Kidney Donor Profile Index and Kidney Donor Risk Index. Kidney Diseases, 2018. DOI: 10.1159/000492427 (citation for the whole package)
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param age Numeric vector. Age, in years.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female.
#' @param ethnicity Vector. Ethnicity, specify in case of African-American patients. The value of variable refers to the parameter label_afroamerican.
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param label_afroamerican List. Label(s) for Afroamerican ethnicity.
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @return numeric eGFR expressed in ml/min/1.73m<sup>2</sup>.
#' @export
#' @name egfr.ckdepi
#' @examples
#' egfr.ckdepi (creatinine = 1.4, age = 60, sex = "Male", ethnicity = "White",
#'   creatinine_units = "mg/dl")

egfr.ckdepi <- function (
              # variables for calculation of eGFR
              creatinine, age, sex, ethnicity, 
              # creatinine measurement units
              creatinine_units = "micromol/l",
              # custom labels for factor parameters and their unknown values - more eexplanations are available at the MDRD function
                  #
                  # label for Afroamerican ethnicity
                  label_afroamerican = c ("Afroamerican"),
                  # label for definition male sex in data set
                  label_sex_male = c ("Male", 1),
                  # label for definition female sex in data set
                  label_sex_female = c ("Female", 0)
              # end of function parameters definition
              ){
  

  
  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("creatinine", "age", "ethnicity", "sex") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)

  # check the range of creatinine_units
  service.check_param_arguments(creatinine_units, c("mg/dl", "micromol/l", "mmol/l") )

  # Check whether numerical arguments inputed by user are fine (weight, height etc have to be positive numbers)
  dta <- data.frame(cbind(age, creatinine))
  service.check_params_numeric(age, creatinine)
  
  
  # check plausible biologic boundaries (by functions in the service.check_plausibility.R):
  #   check and inform user whether any values out of boundaries were substituted by NA
  #   after the check change the value to boundaries in the possible range (i.e. age > 0 and < 100)
  # age
  # first: general check and tidy: age <0 OR age >100
  service.check_plausibility.age(age)
  # second: since this eGFR equation was developed and validated for adults only, notify user if any children were found, and exclude them from calculation
  suspiciosly_low <- service.count_lowerequal_threshhold(age, 17)
  if(suspiciosly_low > 0) cat(service.output_message(suspiciosly_low, "age <=17 years", "NA"))
  age <- service.strict_to_numeric_threshhold_lower(age, 17)
  # creatinine
  service.check_plausibility.creatinine(creatinine)

  # CHECK FUNCTION INPUT: END
  ##################################################################


  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)
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














#' Schwartz (for children only)
##################################################################
# FUNCTION: BEGIN
#' @details Calculate estimated glomerular filtration rate (eGFR) by Schwartz equation
#' Reference to the equation: Gao A, Cachat F, Faouzi M et al. Comparison of the glomerular filtration rate in children by the new revised Schwartz formula and a new generalized formula. Kidney International 2013;83:524–30.
#' 
#' Programming: Boris Bikbov \email{boris@@bikbov.ru}.
#' Citation: Bikbov B. R open source programming code for calculation of the Kidney Donor Profile Index and Kidney Donor Risk Index. Kidney Diseases, 2018. DOI: 10.1159/000492427 (citation for the whole package)
#' @param creatinine Numeric vector. Serum creatinine, could be expressed in "micromol/L", "mmol/L" or "mg/dL". Units of measurement should be defined in variable creatinine_units (if not defined explicitly by user, the default value is "micromol/L").
#' @param age Numeric vector. Age, in years. Age does not accounted in Schwartz equation, but used in the function to check whether Schwartz equation could be applied to a given patient.
#' @param sex Vector. The value of variable refers to the parameters label_sex_male and label_sex_female. Required only in case of quadratic Schwartz equation. 
#' @param height_cm Numeric vector. Could be defined either as height_cm if is measured in cm, or as height_ft and height_inch if is measured in feet and inches.
#'    If the parameter height_cm is greater than 0, the function uses cm, otherwise - feet and inches.
#' @param height_ft see height_cm
#' @param height_inch see height_cm
#' @param creatinine_units Character string. Units in which serum creatinne is expressed. Could be one of the following: "micromol/L", "mmol/L" or "mg/dL".
#' @param equation_type Character string. Define whether to calculate eGFR either by classic Schwartz or quadratic Schwartz equation. Could be one of the following: "classic", "quadratic". If not explicitly defined by user, the default assumption is "classic".
#' @param label_sex_male List. Label(s) for definition(s) of male sex.
#' @param label_sex_female List. Label(s) for definition(s) of female sex.
#' @return numeric eGFR expressed in ml/min/1.73m<sup>2</sup>.
#' @export
#' @name egfr.schwartz
#' @examples
#' egfr.schwartz (creatinine = 1.4, age = 10, height_cm = 90, sex = "Male",
#'   creatinine_units = "mg/dl")
#' egfr.schwartz (creatinine = 1.4, age = 10, height_cm = 90, sex = "Male",
#'   creatinine_units = "mg/dl", equation_type = "quadratic")

egfr.schwartz <- function (
              # variables for calculation of eGFR
              creatinine, age, sex, height_cm = 0, height_ft = 0, height_inch = 0, 
              # creatinine measurement units
              creatinine_units = "micromol/l",
			  # Schwartz equation type
			  equation_type = "classic",
              # custom labels for factor parameters and their unknown values - more eexplanations are available at the MDRD function
                  # label for definition male sex in data set
                  label_sex_male = c ("Male", 1),
                  # label for definition female sex in data set
                  label_sex_female = c ("Female", 0)
              # end of function parameters definition
              ){
  

  
  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # check whether all obligatory argument(s) is(are) defined by user
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  # at the beginning assume that all function arguments are present, and in the following code change it to FALSE if any of the obligatory argument(s) is(are) absent
  fx_params_resulting <- TRUE
  # introduce variable for understanding how user defined height
  local_height_type <- "unknown"
  # Check of height is complex because they could be expressed in different variables by different dimentions
  if ("height_cm" %in% args){
    # do nothing
	local_height_type <- "cm"
  }else{
    # if height_cm not presented check whether height_ft and height_inch are presented
    if( ("height_ft" %in% args) && ("height_inch" %in% args)){
      # do nothing
	  local_height_type <- "ft"
    }else{
      # if neither height_cm nor height_ft and height_inch are presented - there is no data for height at all
      warning("Obligatory argument for height is not defined by user in the function arguments", "\n")
      err_num <- err_num + 1
      fx_params_resulting <- FALSE
    }
  }
  
  # List of obligatory function params which have to be defined by user at the function call
  if(equation_type == "classic"){fx_params <- c("creatinine")} 
  if(equation_type == "quadratic"){fx_params <- c("creatinine", "sex")} 
 
  service.check_obligatory_params(fx_params, args, fx_params_resulting)

  # Check the type of arguments inputed by user
  # check that user defined a single creatinine_units
  service.check_param_number(creatinine_units)

  # check the range of creatinine_units
  service.check_param_arguments(creatinine_units, c("mg/dl", "micromol/l", "mmol/l") )

  # Check whether numerical arguments inputed by user are fine (weight, height etc have to be positive numbers)
  service.check_params_numeric(age, creatinine)
  if(equation_type == "quadratic"){
    if(local_height_type == "cm"){service.check_params_numeric(height_cm)}
    if(local_height_type == "ft"){service.check_params_numeric(height_ft, height_inch)}
  } 
  
  
  # check plausible biologic boundaries (by functions in the service.check_plausibility.R):
  #   check and inform user whether any values out of boundaries were substituted by NA
  #   after the check change the value to boundaries in the possible range (i.e. age > 0 and < 100)
  # age
  # first: general check and tidy: age <0 OR age >100
  service.check_plausibility.age(age)
  # second: since this eGFR equation was developed and validated only for children older 2 years, notify user if younger children were found, and exclude them from calculation
  suspiciosly_low <- service.count_lowerequal_threshhold(age, 2)
  if(suspiciosly_low > 0){
    cat(service.output_message(suspiciosly_low, "age <=2 years", "NA"))
    cat("eGFR for children younger 2 years has to be defined according nomograms\n")
  }
  age <- service.strict_to_numeric_threshhold_lower(age, 2)
  # third: since this eGFR equation was developed and validated only for children, notify user if adults were found, and exclude them from calculation
  suspiciosly_high <- service.count_greater_threshhold(age, 18)
  if(suspiciosly_high > 0) warning(service.output_message(suspiciosly_high, "age >18 years", "NA"))
  age <- service.strict_to_numeric_threshhold_greater(age, 18)

  # creatinine
  service.check_plausibility.creatinine(creatinine)
  
  # I don't check height, since it should be done by user
  

  # CHECK FUNCTION INPUT: END
  ##################################################################


  #
  # repeat creatinine_units according to the length of the file, since the creatinine_units is defined either by user or by default value as a single value for the whole function
  #
  # convert to lower case to avoid any troubles with possible definitions as "mg/dl" or "mg/dL"
  creatinine_units <- tolower(creatinine_units)
  creatinine_units <- rep(creatinine_units, length(creatinine))  
  
  #
  # convert creatinine units if necessary
  #
  creatinine <- service.convert_creatinine(creatinine, creatinine_units)

  #
  # convert height units if necessary
  #
  height <-   ifelse(  height_cm > 0,
            height_cm, # if height is defined in cm, just take this value
            ifelse(height_ft > 0 | height_inch > 0, # if height is defined NOT in cm, check whether height in feets or inches is defined
            height_ft * 30.48 + height_inch * 2.54, # if they contain number than convert to cm
            NA) # if height in feets or inches is not defined, assume NA
        )
  

  # if user required classic equation
  if(equation_type == "classic"){
    eGFR <- 41.3 * (0.01 * height / creatinine )
  }

  # if user required quadratic equation
  if(equation_type == "quadratic"){
  eGFR <- 0.68 * ( height / creatinine ) -0.0008 * (height / creatinine ) * (height / creatinine ) + 0.48 * age
  eGFR <- ifelse( sex %in% label_sex_female,
                        eGFR -25.68,
                        ifelse( sex %in% label_sex_male,
                                eGFR  -21.53,
                                NA
                                )
                )
  
  }


return (round(eGFR, 2))

}


# FUNCTION: END
##################################################################




