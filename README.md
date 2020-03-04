# Kidney.epi package for R

This is the source code for the "kidney.epi" package in R. It gets posted to the comprehensive R archive (CRAN) at intervals, each such posting preceded a throrough test. In general, each new push to CRAN with new function(s) will update the second term of the version number, e.g. 1.2.0 to 1.3.0. Updates only to the code of existing functions increment the third term of the version number, e.g. 1.2.0 to 1.2.1.

Home page of the project: http://www.kidneyepidemiology.org/r  
Twitter of the project: https://twitter.com/KidneyEpiOrg  
Maintainer web-site: http://boris.bikbov.ru/english/  
                     https://www.linkedin.com/in/boris-bikbov-56339a1b/  


## Agreement on prefix of function names and internal datafarames
	nephro. - related to nephrology in general
	ktx. - related to kidney transplantation
	hd. - related to hemodialysis
	pd. - related to peritoneal dialysis
	egfr. - related to equations for calculation of estimated glomerular filtration rate
	epi. - related to epidemiology in general
	service. - service functions or internal purposes of the package (data check, etc)
	
## Internal datafarames
Internal datafarames are contained in the R/sysdata.rda of the source code, used in the R package functions *but not accessible to the user*.  
	ktx.kdpi_mapping_table - contains data with mapping KDPI and KDRI reported by OPTN for the years 2015-2016  
	ktx.kdpi_coefficients - contains data with coefficients used by OPTN for the calculation (KDRI scaling factor, chances of hypertension and diabetes in case if they were unknown for donor, etc)  

## External datafarames 
	ktx - contains data for 10 kidney transplant patients.
	
## Vignettes
	List of vignettes are available via browseVignettes(package = "kidney.epi").

