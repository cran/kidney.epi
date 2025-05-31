## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(kidney.epi)
head(ckd.data)

## -----------------------------------------------------------------------------
# call egfr.ckdepi.cr.2009 function, and directly set parameters values
egfr.ckdepi.cr.2009(
  creatinine = 1.4,  
  age = 60,  
  sex = "Male", 
  ethnicity = "White", 
  creatinine_units = "mg/dl", 
  label_afroamerican = c("Afroamerican"), 
  label_sex_male = c("Male"), 
  label_sex_female = c("Female")
)

# Definitions of the labels for sex and race are optional if you use the same labels defined as default in the function. The following also works well:
egfr.ckdepi.cr.2009(
  creatinine = 1.4,  
  age = 60,  
  sex = "Male", 
  ethnicity = "White", 
  creatinine_units = "mg/dl"
)

# If you measure creatinine in micromol/l, it is possible to omit also 'creatinine_units' since the default value is "micromol/l":
egfr.ckdepi.cr.2009(
  creatinine = 103, # creatinine is in micromol/l
  age = 60,  
  sex = "Male", 
  ethnicity = "White"
)

## -----------------------------------------------------------------------------
# copy as an example the internal dataframe ckd.data from R package to your dataframe
mydata <- ckd.data

# calculate eGFR by CKD-EPI equation
mydata$ckdepi <- egfr.ckdepi.cr.2009(
  creatinine = mydata$cr, age = mydata$age,
  sex = mydata$sex, ethnicity = mydata$ethnicity,
  creatinine_units = "micromol/L",
  # customize all labels for those used in the data frame
  # label(s) used to define male sex in the dataset
  label_sex_male = c("Male"), 
  # label(s) used to define female sex in the dataset
  label_sex_female = c("Female"),
  # label used to define Afroamerican ethnicity in the dataset
  label_afroamerican = c("Black")
) 

# show descriptive stat for the calculated values
# note that synthetic data set ckd.data contains input parameters for both adults and children, and since the CKD-EPI equation was developed and validated for adults only, the resulting eGFR values for children will be NA. Use children-specific eGFR equations when necessary.
summary(mydata$ckdepi)

