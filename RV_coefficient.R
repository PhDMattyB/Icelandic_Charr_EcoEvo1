##############################
## RV coefficients
##
## Matt Brachmann (PhDMattyB)
##
## 2021-03-09
##
##############################

# packages to load --------------------------------------------------------


library(tidyverse)
library(FactoMineR)


# data --------------------------------------------------------------------


data = read_csv('~/PhD/SNP Demographic modelling/RV_coeff/Rv_coefficient_dataframe.csv')


# lapply madness ----------------------------------------------------------


q1_lapply = lapply(data_morphpair, function(x){
  lm(cbind(PW1X, PW1Y, PW2X, PW2Y,PW3X, PW3Y, PW4X, PW4Y, PW5X, PW5Y, PW6X, 
           PW6Y,PW7X, PW7Y, PW8X, PW8Y,PW9X, PW9Y, PW10X, PW10Y,PW11X, PW11Y, 
           PW12X, PW12Y,PW13X, PW13Y, PW14X, PW14Y,
           PW15X, PW15Y, PW16X, PW16Y,
           PW17X, PW17Y, PW18X, PW18Y, 
           PW19X, PW19Y, UNIX, UNIY) ~ q1, data = x)
})

fitted_q1 = q1_lapply %>% 
  lapply(function(x){
    fitted.values(x)
  })

data_lapply_carbon = lapply(data_morphpair, function(x){
  lm(cbind(PW1X, PW1Y, PW2X, PW2Y,PW3X, PW3Y, PW4X, PW4Y, PW5X, PW5Y, PW6X, 
           PW6Y,PW7X, PW7Y, PW8X, PW8Y,PW9X, PW9Y, PW10X, PW10Y,PW11X, PW11Y, 
           PW12X, PW12Y,PW13X, PW13Y, PW14X, PW14Y,
           PW15X, PW15Y, PW16X, PW16Y,
           PW17X, PW17Y, PW18X, PW18Y, 
           PW19X, PW19Y, UNIX, UNIY) ~ carbon, data = x)
})
## grab the fitted values for the matrices
fitted_carbon = data_lapply_carbon %>% 
  lapply(function(x){
    fitted.values(x)
  })
# dim(fitted_carbon)
## nitrogen
data_lapply_nitrogen = lapply(data_morphpair, function(x){
  lm(cbind(PW1X, PW1Y, PW2X, PW2Y,PW3X, PW3Y, PW4X, PW4Y, PW5X, PW5Y, PW6X, 
           PW6Y,PW7X, PW7Y, PW8X, PW8Y,PW9X, PW9Y, PW10X, PW10Y,PW11X, PW11Y, 
           PW12X, PW12Y,PW13X, PW13Y, PW14X, PW14Y,
           PW15X, PW15Y, PW16X, PW16Y,
           PW17X, PW17Y, PW18X, PW18Y, 
           PW19X, PW19Y, UNIX, UNIY) ~ nitrogen, data = x)
})
## grab the fitted values for the matrices
fitted_nitrogen = data_lapply_nitrogen%>% 
  lapply(function(x){
    fitted.values(x)
  })

## isotope RV coeffcients
coeffRV(fitted_carbon$GSBPI, fitted_nitrogen$GSBPI)
coeffRV(fitted_carbon$SLGBPEL, fitted_nitrogen$SLGBPEL)
coeffRV(fitted_carbon$TLGBPL, fitted_nitrogen$TLGBPL)
coeffRV(fitted_carbon$TSBPL, fitted_nitrogen$TSBPL)
coeffRV(fitted_carbon$VBRSIL, fitted_nitrogen$VBRSIL)
coeffRV(fitted_carbon$FLJ, fitted_nitrogen$FLJ)
coeffRV(fitted_carbon$MJOA, fitted_nitrogen$MJOA)


## carbon - gene flow rv coefficients
coeffRV(fitted_carbon$GSBPI, fitted_q1$GSBPI)
coeffRV(fitted_carbon$SLGBPEL, fitted_q1$SLGBPEL)
coeffRV(fitted_carbon$TLGBPL, fitted_q1$TLGBPL)
coeffRV(fitted_carbon$TSBPL, fitted_q1$TSBPL)
coeffRV(fitted_carbon$VBRSIL, fitted_q1$VBRSIL)
coeffRV(fitted_carbon$FLJ, fitted_q1$FLJ)
coeffRV(fitted_carbon$MJOA, fitted_q1$MJOA)


## nitrogen - gene flow rv coefficients
coeffRV(fitted_nitrogen$GSBPI, fitted_q1$GSBPI)
coeffRV(fitted_nitrogen$SLGBPEL, fitted_q1$SLGBPEL)
coeffRV(fitted_nitrogen$TLGBPL, fitted_q1$TLGBPL)
coeffRV(fitted_nitrogen$TSBPL, fitted_q1$TSBPL)
coeffRV(fitted_nitrogen$VBRSIL, fitted_q1$VBRSIL)
coeffRV(fitted_nitrogen$FLJ, fitted_q1$FLJ)
coeffRV(fitted_nitrogen$MJOA, fitted_q1$MJOA)

