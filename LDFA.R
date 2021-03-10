##############################
## Linear discriminant function analyses
##
## Matt Brachmann (PhDMattyB)
##
## 09-03-2021
##
##############################

# load in packages --------------------------------------------------------
library(MASS)
library(tidyverse)

# read in data ------------------------------------------------------------
## This data set has body size and all affine and non-affine components
## of shape variables. The LDFA will calculate an axis of divergence
## related to ALL of the body shape variables

data = read_csv('/~/PhD/Morphometrics/Working Directory/SIASamples_PWS_StandSexAllo_GSTV.csv')

## filter out the lake or morph pair that we want to perform the 
## analysis on. We want to do this on each benthic and pelagic 
## morph pair
pws = data %>% 
  filter(Lake == 'Svinavatn', 
         Morph %in% c('Piscivorous', 
                      'Planktivorous'))
  # filter(Lake == 'Thingvallavatn',
  #        Morph %in% c('Small Benthic',
  #                     'Planktivorous'))
  # filter(Lake == 'Thingvallavatn',
  #        Morph %in% c('Large Benthic',
  #                     'Planktivorous'))
  # filter(Lake == 'Galtabol')
  # filter(Lake == 'Svinavatn')
  # filter(Lake == 'Vatnshlidarvatn')

phenotype = pws %>% 
  dplyr::select(contains('PW'), 
                UNIX, 
                UNIY)

## SVINAVATN ONLY!!!
## Need to combine the two pelagic morphs from svinavatn
## this is based on the population genetics of the population
pws = mutate(.data = pws, 
             Morph = as.factor(case_when(
               Morph == 'Large Benthic' ~ 'Large benthic', 
               Morph == 'Piscivorous' ~ 'Pelagic', 
               Morph == 'Planktivorous' ~ 'Pelagic'
             )))
# LDFA analysis -----------------------------------------------------------

ldfa = lda(pws$Morph ~ pws$PW1X+ pws$PW1Y+ pws$PW2X+ pws$PW2Y+ 
                        pws$PW3X+ pws$PW3Y+ pws$PW4X+ pws$PW4Y+
                        pws$PW5X+ pws$PW5Y+ pws$PW6X+ pws$PW6Y+
                        pws$PW7X+ pws$PW7Y+ pws$PW8X+ pws$PW8Y+ 
                        pws$PW9X+ pws$PW9Y+ pws$PW10X+ pws$PW10Y+
                        pws$PW11X+ pws$PW11Y+ pws$PW12X+ pws$PW12Y+
                        pws$PW13X+ pws$PW13Y+ pws$PW14X+ pws$PW14Y+
                        pws$PW15X+ pws$PW15Y+ pws$PW16X+ pws$PW16Y+
                        pws$PW17X+ pws$PW17Y+ pws$PW18X+ pws$PW18Y+ 
                        pws$PW19X+ pws$PW19Y+ pws$UNIX+ pws$UNIY, 
           CV = T)

##Predicted values based on the ldfa
## CV CANNOT BE TRUE FOR THIS!!!!

ldfa_predict = predict(ldfa)
apply(ldfa_predict$posterior, MARGIN = 1, FUN = max)

##Plot the first linear discriminant axis!
## these plots show how the analysis is splitting up the data
## it also shows you how well the analysis splits up the data
## based on the grouping variables. 

plot(ldfa)
plot(ldfa, dimen = 1, type = 'both')

## tells us how well the analysis did without a cross validation
table = table(pws$Morph, ldfa_predict$class)
sum(table[row(table) == col(table)])/sum(table)

# Cross Validation --------------------------------------------------------
#leave one out cross validation
## NEED TO ADD CV = TRUE IN THE LDFA
CV = table(pws$Morph, ldfa$class)

##Calculate the re-substitution error for the cross validation
sum(CV[row(CV) == col(CV)])/sum(CV)