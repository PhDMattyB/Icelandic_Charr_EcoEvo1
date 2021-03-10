##############################
## Linear and non linear regression
##
## Matt Brachmann (PhDMattyB)
##
## 09-03-2021
##
##############################
## This script is to compare linear and non-linear regression models
## for each population of Icelandic Arctic charr

library(tidyverse)

# data --------------------------------------------------------------------
le_data = read_csv('~/PhD/SNP Demographic Modelling/Working Directory/GSTVMF_Morph_Eco_Geno.csv')%>%
  filter(Lake %in% c('Galtabol', 'Svinavatn', 'Thingvallavatn',
                     'Vatnshlidarvatn')) 

le_data = mutate(.data = le_data,
                         Vector2 = as.factor(case_when(
                           Vector == "GSBPI" ~ "GSBPI",
                           Vector == 'SLGBPI' ~ 'SLGBPEL',
                           Vector == 'SLGBPL' ~ 'SLGBPEL',
                           Vector == 'TLGBPL'~ 'TLGBPL',
                           Vector == 'TSBPL' ~ 'TSBPL',
                           Vector == 'VBRSIL' ~ 'VBRSIL')))

## LM NLM Compare #####
sub_data = le_data %>% group_by(Vector) %>%
  dplyr::select(Specimen.ID, Lake, Morph, BP, Vector, carbon, nitrogen,
                Sex, AGE, Fork.length, eco_vec, BP_ld1_noallo) %>%
  ## filter the vector name for the morph comparisons to check out
  filter(Vector == 'GSBPI')

shape_models = function(data){
  shape.lm1 = lm(BP_ld1_noallo ~ carbon + nitrogen, data = data)
  summary(shape.lm1)
  AIC1 = AIC(shape.lm1)
  shape.nlm3 = lm(BP_ld1_noallo ~ poly(carbon, 2) + poly(nitrogen, 2), 
                  data = data)
  summary(shape.nlm3)
  AIC2 = AIC(shape.nlm3) 
  shape.nlm4 = lm(BP_ld1_noallo ~ poly(carbon, 3) + poly(nitrogen, 3), 
                  data = data)
  summary(shape.nlm4)
  AIC3 = AIC(shape.nlm4)
  # output = list(AIC1, 
  #               AIC2, 
  #               AIC3)
  # return(output)
  print(paste(AIC1, 'AIC value for linear model'))
  print(paste(AIC2, 'AIC value for quadratic model'))
  print(paste(AIC3, 'AIC value for polynomial model'))
  
}

shape_models(sub_data)
## Here are the model pieces

# ## Shape models
# ##Fit a linear relationship
# shape.lm1 = lm(BP_ld1_noallo ~ carbon + nitrogen, data = sub_data)
# summary(shape.lm1)
# AIC(shape.lm1)
# 
# ##Quadratic fit to the data. 
# shape.nlm3 = lm(BP_ld1_noallo ~ poly(carbon, 2) + poly(nitrogen, 2), 
#                 data = sub_data)
# summary(shape.nlm3)
# AIC(shape.nlm3) 
# 
# ## Polynomial fit to the data
# shape.nlm4 = lm(BP_ld1_noallo ~ poly(carbon, 3) + poly(nitrogen, 3), 
#                 data = sub_data)
# summary(shape.nlm4)
# AIC(shape.nlm4)


size_models = function(data){
  size.lm1 = lm(Fork.length ~ carbon + nitrogen, data = sub_data)
  summary(size.lm1)
  AIC1 = AIC(size.lm1)
  
  size.nlm3 = lm(Fork.length ~ poly(carbon, 2) + poly(nitrogen, 2), data = sub_data)
  summary(size.nlm3)
  AIC2 = AIC(size.nlm3)
  
  size.nlm4 = lm(Fork.length ~ poly(carbon, 3) + poly(nitrogen, 3), data = sub_data)
  summary(size.nlm4)
  AIC3 = AIC(size.nlm4)
  print(paste(AIC1, 'AIC value for linear model'))
  print(paste(AIC2, 'AIC value for quadratic model'))
  print(paste(AIC3, 'AIC value for polynomial model'))
  
}

size_models(sub_data)

## Here are the model pieces
# ## size models
# size.lm1 = lm(Fork.length ~ carbon + nitrogen, data = sub_data)
# summary(size.lm1)
# AIC(size.lm1)
#  
# size.nlm3 = lm(Fork.length ~ poly(carbon, 2) + poly(nitrogen, 2), data = sub_data)
# summary(size.nlm3)
# AIC(size.nlm3)
# 
# size.nlm4 = lm(Fork.length ~ poly(carbon, 3) + poly(nitrogen, 3), data = sub_data)
# summary(size.nlm4)
# AIC(size.nlm4)
## Take the predicted values for plottinn using predict(modelname)
