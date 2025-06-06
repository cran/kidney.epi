# Calculate categories of albuminuria and eGFR based on laboratory values
#' Calculate KDIGO albuminuria categories
##################################################################
# FUNCTION: BEGIN
#' @details Calculate albuminuria categories (A1, A2, A3) based on albuminuria values.
#'
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param albuminuria Numeric vector. Urine albumin, could be expressed in "mg/day" (24-hour urine collection), "mg/mmol" (UACR) or "mg/g" (UACR). Units of measurement should be defined in variable albuminuria_units (if not defined explicitly by user, the default value is "mg/g").
#' @param albuminuria_units Character string. Units in which urine albumin is measured. Could be one of the following: "mg/day", "mg/mmol" or "mg/g".
#' @param semiquantitative_values Character string. Defines whether semiquantitative values are allowed in the data. If "any", all semiquantitative values ('<30', '30-300', '>300') and any numeric values (29, 30, 35, etc) will be classified into A categories (NB! both '30-300' and '30-299' will be classified as A2). If "only_limits", only limiting semiquantitative values ('<30', '>300') will be classified into A categories, but middle semiquantitative values ('30-300') will be omitted; but numeric values (29, 30, 35, etc) will be classified into A categories. If "forbidden", only numeric values will be classified into A categories.
#' @return string albuminuria category.
#' @export
#' @name ckd.kdigo_category.albuminuria
#' @examples
#' # for a single patient
#' ckd.kdigo_category.albuminuria(albuminuria = 25, albuminuria_units = "mg/g")
#' # for a dataset - see vignettes for details
#' # ckd.kdigo_category.albuminuria(albuminuria = dta$alb, albuminuria_units = "mg/g")

ckd.kdigo_category.albuminuria <- function(
  albuminuria, albuminuria_units = "mg/g",
  semiquantitative_values = "forbidden"
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("albuminuria") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single albuminuria_units
  service.check_param_number(albuminuria_units)

  # check the range of albuminuria_units
  albuminuria_units <- tolower(albuminuria_units)
  service.check_param_arguments(albuminuria_units, c("mg/day", "mg/mmol", "mg/g"))

  # check the range of semiquantitative_values
  semiquantitative_values <- tolower(semiquantitative_values)
  service.check_param_arguments(semiquantitative_values, c("any", "only_limits", "forbidden"))

  # check plausible biologic boundaries
  #albuminuria <- service.check_and_correct_numeric(albuminuria, "albuminuria")

  # CHECK FUNCTION INPUT: END
  ##################################################################

  # define albuminuria categories (A1, A2, A3) tresholds
  if (albuminuria_units == "mg/day") {
    A1_lower_threshold <- 0
    A2_lower_threshold <- 30
    A3_lower_threshold <- 300
  }

  if (albuminuria_units == "mg/mmol") {
    A1_lower_threshold <- 0
    A2_lower_threshold <- 3
    A3_lower_threshold <- 30
  }

  if (albuminuria_units == "mg/g") {
    A1_lower_threshold <- 0
    A2_lower_threshold <- 30
    A3_lower_threshold <- 300
  }



  if (semiquantitative_values == "forbidden") {
    # remove internally values with signs or with letters
    albuminuria[ grepl("[A-Za-z]|[-<>=]", albuminuria ) ] <- NA
    albuminuria <- as.numeric(albuminuria)
    # classify cleaned values
    albuminuria_category <- NA
    albuminuria_category[albuminuria >= A1_lower_threshold & albuminuria < A2_lower_threshold] <- "A1"
    albuminuria_category[albuminuria >= A2_lower_threshold & albuminuria < A3_lower_threshold] <- "A2"
    albuminuria_category[albuminuria >= A3_lower_threshold] <- "A3"
  }

  if (semiquantitative_values == "any" || semiquantitative_values == "only_limits") {

    if (semiquantitative_values == "any") {
      # if '30-300' or '30-299' - just assign internally a value corresponding to A2 
      albuminuria [ albuminuria == paste0(A2_lower_threshold, "-", A3_lower_threshold)] <- A3_lower_threshold - 1
      albuminuria [ albuminuria == paste0(A2_lower_threshold, "-", A3_lower_threshold - 1)] <- A3_lower_threshold - 1
    }
	
	if (semiquantitative_values == "only_limits") {
	  albuminuria <- ifelse(grepl("-", albuminuria), NA, albuminuria)
	}
    
    # proceed both recoded "any" and "only_limits" in the same way
    albuminuria_category <- rep(NA_character_, length(albuminuria)) 
    albuminuria_value <- as.numeric(gsub("[<>=]", "", albuminuria))
    albuminuria_operator <- gsub("[0-9.]+", "", albuminuria)
    
    albuminuria_category[albuminuria_value < A1_lower_threshold] <- NA
    # code A3 first
    albuminuria_category[(albuminuria_operator == ">" & albuminuria_value > A3_lower_threshold & albuminuria_value >= A1_lower_threshold) | (albuminuria_value >= A3_lower_threshold & albuminuria_value >= A1_lower_threshold)] <- "A3"
    # code A2 first to replace it later if <30
    albuminuria_category[(albuminuria_value >= A2_lower_threshold & albuminuria_value < A3_lower_threshold & albuminuria_value >= A1_lower_threshold)] <- "A2"
    albuminuria_category[(albuminuria_operator == "<" & albuminuria_value >= A1_lower_threshold) | (albuminuria_value < A2_lower_threshold & albuminuria_value >= A1_lower_threshold)] <- "A1"
    albuminuria_category[(albuminuria_operator == "<=" & albuminuria_value <= A2_lower_threshold & albuminuria_value >= A1_lower_threshold)] <- "A1"
  }
  
  
return(albuminuria_category)
}

# FUNCTION: END
##################################################################

#' Legacy function to calculate albuminuria categories
##################################################################
#' @details Legacy function to calculate albuminuria categories (A1, A2, A3) based on albuminuria values. Alias to the ckd.kdigo_category.albuminuria() function. This function is for legacy use only.
#' @param ...  all arguments for the ckd.kdigo_category.albuminuria() function.
#' @export
#' @name nephro.albuminuria_category
nephro.albuminuria_category <- function(...) {
	ckd.kdigo_category.albuminuria(...)
}




#'
#' Calculate proteinuria categories
##################################################################
# FUNCTION: BEGIN
#' @details Calculate albuminuria categories (A1, A2, A3) based on proteinuria values.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#'
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param proteinuria Numeric vector. Urine protein, could be expressed in "mg/day" (24-hour urine collection), "mg/mmol" (UPCR) or "mg/g" (UPCR). Units of measurement should be defined in variable proteinuria_units (if not defined explicitly by user, the default value is "mg/g").
#' @param proteinuria_units Character string. Units in which urine protein is measured. Could be one of the following: "mg/day", "mg/mmol" or "mg/g".
#' @param semiquantitative_values Character string. Defines whether semiquantitative values are allowed in the data. If "any", all semiquantitative values ('<30', '30-300', '>300') and any numeric values (29, 30, 35, etc) will be classified into A categories (NB! both '30-300' and '30-299' will be classified as A2). If "only_limits", only limiting semiquantitative values ('<30', '>300') will be classified into A categories, but middle semiquantitative values ('30-300') will be omitted; but numeric values (29, 30, 35, etc) will be classified into A categories. If "forbidden", only numeric values will be classified into A categories. 
#' @return string albuminuria category.
#' @export
#' @name ckd.kdigo_category.proteinuria
#' @examples
#' # for a single patient
#' ckd.kdigo_category.proteinuria(proteinuria = 25, proteinuria_units = "mg/g")
#' # for a dataset - see vignettes for details
#' # ckd.kdigo_category.proteinuria(proteinuria = dta$alb, proteinuria_units = "mg/g")

ckd.kdigo_category.proteinuria <- function(
  proteinuria, proteinuria_units = "mg/g",
  semiquantitative_values = "forbidden"
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("proteinuria") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # Check the type of arguments inputed by user
  # check that user defined a single albuminuria_units
  service.check_param_number(proteinuria_units)

  # check the range of albuminuria_units
  proteinuria_units <- tolower(proteinuria_units)
  service.check_param_arguments(proteinuria_units, c("mg/day", "mg/mmol", "mg/g"))

  # check the range of semiquantitative_values
  semiquantitative_values <- tolower(semiquantitative_values)
  service.check_param_arguments(semiquantitative_values, c("any", "only_limits", "forbidden"))

  # check plausible biologic boundaries
  #proteinuria <- service.check_and_correct_numeric(proteinuria, "proteinuria")

  # CHECK FUNCTION INPUT: END
  ##################################################################

  # define albuminuria categories (A1, A2, A3) tresholds
  if (proteinuria_units == "mg/day") {
    A1_lower_threshold <- 0
    A2_lower_threshold <- 150
    A3_lower_threshold <- 500
  }

  if (proteinuria_units == "mg/mmol") {
    A1_lower_threshold <- 0
    A2_lower_threshold <- 15
    A3_lower_threshold <- 50
  }

  if (proteinuria_units == "mg/g") {
    A1_lower_threshold <- 0
    A2_lower_threshold <- 150
    A3_lower_threshold <- 500
  }

  if (semiquantitative_values == "forbidden") {
    # remove internally values with signs or with letters
    proteinuria[ grepl("[A-Za-z]|[-<>=]", proteinuria ) ] <- NA
    proteinuria <- as.numeric(proteinuria)
    # classify cleaned values
    albuminuria_category <- NA
    albuminuria_category[proteinuria >= A1_lower_threshold & proteinuria < A2_lower_threshold] <- "A1"
    albuminuria_category[proteinuria >= A2_lower_threshold & proteinuria < A3_lower_threshold] <- "A2"
    albuminuria_category[proteinuria >= A3_lower_threshold] <- "A3"
  }

  if (semiquantitative_values == "any" || semiquantitative_values == "only_limits") {

    if (semiquantitative_values == "any") {
      # if '30-300' or '30-299' - just assign internally a value corresponding to A2 
      proteinuria [ proteinuria == paste0(A2_lower_threshold, "-", A3_lower_threshold)] <- A3_lower_threshold - 1
      proteinuria [ proteinuria == paste0(A2_lower_threshold, "-", A3_lower_threshold - 1)] <- A3_lower_threshold - 1
    }

	if (semiquantitative_values == "only_limits") {
	  proteinuria <- ifelse(grepl("-", proteinuria), NA, proteinuria)
	}
    
    # proceed both recoded "any" and "only_limits" in the same way
    albuminuria_category <- rep(NA_character_, length(proteinuria)) 
    proteinuria_value <- as.numeric(gsub("[<>=]", "", proteinuria))
    proteinuria_operator <- gsub("[0-9.]+", "", proteinuria)
    
    albuminuria_category[proteinuria_value < A1_lower_threshold] <- NA
    # code A3
    albuminuria_category[(proteinuria_operator == ">" & proteinuria_value > A3_lower_threshold & proteinuria_value >= A1_lower_threshold) | (proteinuria_value >= A3_lower_threshold & proteinuria_value >= A1_lower_threshold)] <- "A3"
    # code A2 first to replace it later if <150
    albuminuria_category[(proteinuria_value >= A2_lower_threshold & proteinuria_value < A3_lower_threshold & proteinuria_value >= A1_lower_threshold)] <- "A2"
    albuminuria_category[(proteinuria_operator == "<" & proteinuria_value >= A1_lower_threshold) | (proteinuria_value < A2_lower_threshold & proteinuria_value >= A1_lower_threshold)] <- "A1"
    albuminuria_category[(proteinuria_operator == "<=" & proteinuria_value <= A2_lower_threshold & proteinuria_value >= A1_lower_threshold)] <- "A1"
  }

return(albuminuria_category)
}

# FUNCTION: END
##################################################################

#' Legacy function to calculate proteinuria categories
##################################################################
#' @details Legacy function to calculate albuminuria categories (A1, A2, A3) based on proteinuria values. Alias to the ckd.kdigo_category.proteinuria() function. This function is for legacy use only.
#' @param ...  all arguments for the ckd.kdigo_category.proteinuria() function.
#' @export
#' @name nephro.proteinuria_category
nephro.proteinuria_category <- function(...) {
	ckd.kdigo_category.proteinuria(...)
}



#'
#' Calculate eGFR categories
##################################################################
# FUNCTION: BEGIN
#' @details Calculate eGFR categories (G1, G2, G3a, G3b, G4, G5) based on eGFR values.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#'
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param gfr Numeric vector. eGFR, expressed in "ml/min/1.73m2".
#' @return string gfr category.
#' @export
#' @name ckd.kdigo_category.gfr
#' @examples
#' # for a single patient
#' ckd.kdigo_category.gfr(gfr = 25)
#' # for a dataset - see vignettes for details
#' # ckd.kdigo_category.gfr(gfr = dta$egfr)

ckd.kdigo_category.gfr <- function(
  gfr
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("gfr") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  # check plausible biologic boundaries
  gfr <- service.check_and_correct_numeric(gfr, "GFR")

  # CHECK FUNCTION INPUT: END
  ##################################################################

  # define gfr categories (G1, G2, G3a, G3b, G4, G5) thresholds
  G1_lower_threshold <- 90
  G2_lower_threshold <- 60
  G3a_lower_threshold <- 45
  G3b_lower_threshold <- 30
  G4_lower_threshold <- 15
  G5_lower_threshold <- 0

  gfr_category <- 
    ifelse( gfr >= G5_lower_threshold & gfr < G4_lower_threshold, "G5",
      ifelse( gfr >= G4_lower_threshold & gfr < G3b_lower_threshold, "G4", 
      ifelse( gfr >= G3b_lower_threshold & gfr < G3a_lower_threshold, "G3b", 
      ifelse( gfr >= G3a_lower_threshold & gfr < G2_lower_threshold, "G3a", 
      ifelse( gfr >= G2_lower_threshold & gfr < G1_lower_threshold, "G2", 
      ifelse( gfr >= G1_lower_threshold, "G1", "" 
      ))))))

return(gfr_category)
}

# FUNCTION: END
##################################################################


#' Legacy function to calculate KDIGO GFR categories
##################################################################
#' @details Legacy function to calculate KDIGO GFR categories. Alias to the ckd.kdigo_category.gfr() function. This function is for legacy use only.
#' @param ...  all arguments for the ckd.kdigo_category.gfr() function.
#' @export
#' @name nephro.gfr_category
nephro.gfr_category <- function(...) {
	ckd.kdigo_category.gfr(...)
}




#'
#' Calculate KDIGO risk categories
##################################################################
# FUNCTION: BEGIN
#' @details Calculate KDIGO risk of complications categories (Low, Moderate, High, Very high) based on eGFR and albuminuria grades.
#' 
#' Citation: Bikbov B. kidney.epi: Kidney-Related Functions for Clinical and Epidemiological Research. Scientific-Tools.Org, https://Scientific-Tools.Org. DOI: 10.32614/CRAN.package.kidney.epi
#'
#' @author Programming: Boris Bikbov https://www.linkedin.com/in/boris-bikbov.
#' @param gfr_cat Character vector. eGFR categories coded as G1, G2, G3a, G3b, G4, G5.
#' @param alb_cat Character vector. Albuminuria categories coded as A1, A2, A3.
#' @return string risk category.
#' @export
#' @name ckd.kdigo_category.risk
#' @examples
#' # for a single patient
#' ckd.kdigo_category.risk(gfr_cat = "G2", alb_cat = "A3")
#' # for a dataset - see vignettes for details
#' # ckd.kdigo_category.risk(gfr_cat = dta$gfr_cat, alb_cat = dta$alb_cat)

ckd.kdigo_category.risk <- function(
  gfr_cat, alb_cat
) {

  ##################################################################
  # CHECK FUNCTION INPUT: BEGIN

  if(!service.check_equal_length(gfr_cat, alb_cat)) stop("The length of variables should be equal.")
  
  # check whether all obligatory argument(s) is(are) defined by user
  fx_params <- c("alb_cat", "gfr_cat") # List of obligatory function params which have to be defined by user at the function call
  args <- names(as.list(match.call())[-1]) # take all params defined by user from the function
  service.check_obligatory_params(fx_params, args)
  
  service.check_string_values(unique(gfr_cat), c("G1", "G2", "G3a", "G3b", "G4", "G5"))
  service.check_string_values(unique(alb_cat), c("A1", "A2", "A3"))

  risk <- character(length(gfr_cat))
  
  risk [
    (gfr_cat == "G1" & alb_cat == "A1") | 
    (gfr_cat == "G2" & alb_cat == "A1")
  ] <- "No CKD"
  
  risk [
    (gfr_cat == "G1" & alb_cat == "A2") | 
    (gfr_cat == "G2" & alb_cat == "A2") | 
    (gfr_cat == "G3a" & alb_cat == "A1") 
  ] <- "Moderate"
  
  risk [
    (gfr_cat == "G1" & alb_cat == "A3") | 
    (gfr_cat == "G2" & alb_cat == "A3") | 
    (gfr_cat == "G3a" & alb_cat == "A2") | 
    (gfr_cat == "G3b" & alb_cat == "A1")
  ] <- "High"
  
  risk [
    (gfr_cat == "G4" & alb_cat == "A1") | 
    (gfr_cat == "G5" & alb_cat == "A1") | 
    (gfr_cat == "G3b" & alb_cat == "A2") | 
    (gfr_cat == "G4" & alb_cat == "A2") | 
    (gfr_cat == "G5" & alb_cat == "A2") | 
    (gfr_cat == "G3a" & alb_cat == "A3") | 
    (gfr_cat == "G3b" & alb_cat == "A3") | 
    (gfr_cat == "G4" & alb_cat == "A3") | 
    (gfr_cat == "G5" & alb_cat == "A3") 
  ] <- "Very high"

return(risk)
}

# FUNCTION: END
##################################################################



#' Legacy function to calculate KDIGO risk categories
##################################################################
#' @details Legacy function to calculate KDIGO risk categories. Alias to the ckd.kdigo_category.risk() function. This function is for legacy use only.
#' @param ...  all arguments for the ckd.kdigo_category.risk() function.
#' @export
#' @name nephro.gfr_category
nephro.kdigo_risk_category <- function(...) {
	ckd.kdigo_category.risk(...)
}
