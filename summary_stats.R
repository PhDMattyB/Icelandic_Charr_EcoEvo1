##############################
## Morpho-eco-geno relationships
##
## Matt Brachmann (PhDMattyB)
##
## 2021-03-09
##
##############################

# install packages --------------------------------------------------------
library(patchwork)
library(lme4)
library(tidyverse)

# Read in the data --------------------------------------------------------
le_data = read_csv('~/PhD/SNP Demographic Modelling/WOrking Directory/GSTVMF_Morph_Eco_Geno.csv ')%>% 
  filter(Lake %in% c('Galtabol', 'Svinavatn', 
                     'Thingvallavatn', 'Vatnshlidarvatn')) %>% 
  filter(LaMorph != 'S.LGB2',
         LaMorph != 'T.PL2') %>% 
  group_by(LaMorph)

## Anovas and ttests ####

data = mutate(.data = data, 
              Full_name = as.factor(case_when(
                Vector == 'GSBPI' ~ 'GSBPI', 
                Vector == 'SLGBPI' ~ 'SLGBPEL', 
                Vector == 'SLGBPL' ~ 'SLGBPEL', 
                Vector == 'TLGBPL' ~ 'TLGBPL', 
                Vector == 'TSBPL' ~ 'TSBPL', 
                Vector == 'VBRSIL' ~ 'VBRSIL'
              )))


sub_data = le_data %>% 
  ## Specify which lake to perform anova on
  filter(Lake == 'Thingvallavatn')

aov_model = function(data){
  aov_model1 = aov(data = data, 
                  carbon ~ Morph)
  aov_model2 = aov(data = data, 
                   nitrogen ~ Morph)
  sum1 = summary(aov_model1)
  sum2 = summary(aov_model2)
  Tukey1 = TukeyHSD(aov_model1)
  Tukey2 = TukeyHSD(aov_model2)
  output = list(aov_model1, 
                sum1, 
                Tukey1, 
                aov_model2, 
                sum2, 
                Tukey2)
 
  return(output)
  
}

aov_model(sub_data)

## For lake with two populations
## perform the t-tests
sub_data = le_data %>% 
  filter(Lake == 'Galtabol')

t.test_output = function(data, Lake){
  require(dplyr)
  data %>%
    group_by(BP)
  
  test1 = t.test(data = data, 
                 Fork.length ~ BP)
  test2 = t.test(data = data, 
                 carbon ~ BP)
  test3 = t.test(data = data, 
                 nitrogen ~ BP)
  
  output = list(test1,
                test2, 
                test3)
  
  return(output)
    
}

t.test_output(sub_data)

## Hybrid percentage #####


hybrid_num = function(data){
  h1 = data %>% dplyr::select(Specimen.ID,
                         LaMorph, 
                         q1) %>% 
    na.omit() %>%
    group_by(LaMorph) %>% 
    count(q1<0.75)
  h2 = data %>% dplyr::select(Specimen.ID,
                              LaMorph, 
                              q1) %>% 
    na.omit() %>%
    group_by(LaMorph) %>% 
    count(q1>0.75)
  # h1 = (h1/312)*100
  # h2 = (h2/312)*100
 output = list(h1, 
               h2)
 return(output)
}

test = le_data %>% 
  filter(Vector == 'VBRSIL')

hybrid_num(test)

# (19/312)*100
## 6.1% of all individuals have hybrid genotypes
## G 0%
## S (2/90)*100 = 2.2%
## TLGBPL (2/63)*100 = 3.2%
## TSBPL (6/63) *100 = 9.5%
## T (8/96)*100 = 8.3%
## V (10/64)*100 = 15.6%
