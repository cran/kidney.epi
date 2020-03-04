## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(kidney.epi)
head(ktx)

## ------------------------------------------------------------------------
# call egfr.ckdepi function, and directly set parameters values
egfr.ckdepi (creatinine = 1.4, age = 60, sex = "Male", ethnicity = "White", creatinine_units = "mg/dl", label_afroamerican = c ("Afroamerican"), label_sex_male = c ("Male"), label_sex_female = c ("Female"))

## ------------------------------------------------------------------------
# copy as an example the internal dataframe ktx from R package to your dataframe
mydata <- ktx

# calculate eGFR by CKD-EPI equation
mydata$ckdepi <- egfr.ckdepi ( creatinine = mydata$don.creatinine, age = mydata$don.age,
  sex = mydata$don.sex, ethnicity = mydata$don.ethnicity, creatinine_units = "mg/dl",
  # customize all labels used in the dataframe
  label_afroamerican = c ("Afroamerican"),
  label_sex_male = c ("Male"), label_sex_female = c ("Female")) 

# show descriptive stat for the calculated values
summary(mydata$ckdepi)

