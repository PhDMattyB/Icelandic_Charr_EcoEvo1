##############################
## Calculating Bhattacharyya distance
##
## Matt Brachmann (PhDMattyB)
##
## 2020-12-29
##
##############################

library(fpc)
library(Hotelling)
library(patchwork)
library(tidyverse)

# Data --------------------------------------------------------------------
isotopes = read_csv('~/PhD/SNP Demographic Modelling/Working Directory/GSTVMF_Morph_Eco_Geno.csv') %>% 
  filter(Lake %in% c('Galtabol',
                     'Svinavatn', 
                     'Thingvallavatn', 
                     'Vatnshlidarvatn')) %>% 
  select(Specimen.ID, 
         Lake, 
         Morph, 
         BP, 
         Vector, 
         LaMorph, 
         carbon, 
         nitrogen)

# Analysis ----------------------------------------------------------------

# Galtabol ----------------------------------------------------------------

G.benthic = isotopes %>% 
  filter(Lake == 'Galtabol', 
         BP == 'Benthic') %>%
  summarise(mean_carbon = mean(carbon),
            mean_nitrogen = mean(nitrogen),
            covariance = cov(x = carbon,
                             y = nitrogen,
                             method = 'pearson'))

G.pelagic = isotopes %>% 
  filter(Lake == 'Galtabol', 
         BP == 'Pelagic') %>%
  summarise(mean_carbon = mean(carbon),
            mean_nitrogen = mean(nitrogen),
            covariance = cov(x = carbon,
                             y = nitrogen,
                             method = 'pearson'))

cov_b = diag(2)
cov_b[1,2] = -0.905
cov_b[2,1] = -0.905

cov_p = diag(2)
cov_p[1,2] = -1.62
cov_p[2,1] = -1.62


Galtabol_dist = bhattacharyya.dist(c(-18.0, 6.80), 
                   c(-20.1, 7.54), 
                   Sigma1 = cov_b, 
                   Sigma2 = cov_p) %>% 
  as_tibble()

### 1000 permutation Hotelling T-test
### This is a non-parametric multivariate t-test
G.benthic_stat = isotopes %>% 
  filter(Lake == 'Galtabol', 
         BP == 'Benthic') %>% 
  select(carbon, nitrogen)

G.pelagic_stat = isotopes %>% 
  filter(Lake == 'Galtabol', 
         BP == 'Pelagic') %>%
  select(carbon, nitrogen)

g.hoteltest = hotelling.test(G.benthic_stat, 
                             G.pelagic_stat, 
                             perm = TRUE, 
                             B = 1000)
# Svinavatn ---------------------------------------------------------------

S.benthic = isotopes %>% 
  filter(Lake == 'Svinavatn', 
         BP == 'Benthic') %>%
  summarise(mean_carbon = mean(carbon),
            mean_nitrogen = mean(nitrogen),
            covariance = cov(x = carbon,
                             y = nitrogen,
                             method = 'pearson'))

S.pelagic = isotopes %>% 
  filter(Lake == 'Svinavatn', 
         BP == 'Pelagic') %>%
  summarise(mean_carbon = mean(carbon),
            mean_nitrogen = mean(nitrogen),
            covariance = cov(x = carbon,
                             y = nitrogen,
                             method = 'pearson'))

cov_b = diag(2)
cov_b[1,2] = -0.523
cov_b[2,1] = -0.523

cov_p = diag(2)
cov_p[1,2] = 0.176
cov_p[2,1] = 0.176


Svinavatn_dist = bhattacharyya.dist(c(-22.0, 9.07), 
                   c(-28.4, 9.60), 
                   Sigma1 = cov_b, 
                   Sigma2 = cov_p)%>% 
  as_tibble()

### 1000 permutation Hotelling T-test
### This is a non-parametric multivariate t-test
s.benthic_stat = isotopes %>% 
  filter(Lake == 'Svinavatn', 
         BP == 'Benthic') %>% 
  select(carbon, nitrogen)

s.pelagic_stat = isotopes %>% 
  filter(Lake == 'Svinavatn', 
         BP == 'Pelagic') %>%
  select(carbon, nitrogen)

s.hoteltest = hotelling.test(s.benthic_stat, 
                             s.pelagic_stat, 
                             perm = TRUE, 
                             B = 1000)


# Thingvallavatn 1 --------------------------------------------------------
T.benthic1 = isotopes %>% 
  filter(Lake == 'Thingvallavatn', 
         Vector == 'TLGBPL', 
         BP == 'Benthic') %>%
  summarise(mean_carbon = mean(carbon),
            mean_nitrogen = mean(nitrogen),
            covariance = cov(x = carbon,
                             y = nitrogen,
                             method = 'pearson'))

T.pelagic = isotopes %>% 
  filter(Lake == 'Thingvallavatn',
         Vector == 'TLGBPL',
         BP == 'Pelagic') %>%
  summarise(mean_carbon = mean(carbon),
            mean_nitrogen = mean(nitrogen),
            covariance = cov(x = carbon,
                             y = nitrogen,
                             method = 'pearson'))

cov_b = diag(2)
cov_b[1,2] = -0.295
cov_b[2,1] = -0.295

cov_p = diag(2)
cov_p[1,2] = 0.678
cov_p[2,1] = 0.678


Thingvallavatn1_dist = bhattacharyya.dist(c(-11.8, 5.68), 
                                    c(-28.7, 6.15), 
                                    Sigma1 = cov_b, 
                                    Sigma2 = cov_p)%>% 
  as_tibble()


### 1000 permutation Hotelling T-test
### This is a non-parametric multivariate t-test
t1.benthic_stat = isotopes %>% 
  filter(Lake == 'Thingvallavatn',
         Vector == 'TLGBPL',
         BP == 'Benthic') %>%
  select(carbon, nitrogen)

t1.pelagic_stat = isotopes %>% 
  filter(Lake == 'Thingvallavatn',
         Vector == 'TLGBPL',
         BP == 'Pelagic') %>%
  select(carbon, nitrogen)

t1.hoteltest = hotelling.test(t1.benthic_stat, 
                             t1.pelagic_stat, 
                             perm = TRUE, 
                             B = 1000)

# Thingvallavatn 2 --------------------------------------------------------
T.benthic2 = isotopes %>% 
  filter(Lake == 'Thingvallavatn', 
         Vector == 'TSBPL', 
         BP == 'Benthic') %>%
  summarise(mean_carbon = mean(carbon),
            mean_nitrogen = mean(nitrogen),
            covariance = cov(x = carbon,
                             y = nitrogen,
                             method = 'pearson'))

T.pelagic2 = isotopes %>% 
  filter(Lake == 'Thingvallavatn',
         Vector == 'TSBPL',
         BP == 'Pelagic') %>%
  summarise(mean_carbon = mean(carbon),
            mean_nitrogen = mean(nitrogen),
            covariance = cov(x = carbon,
                             y = nitrogen,
                             method = 'pearson'))

cov_b = diag(2)
cov_b[1,2] = -4.34
cov_b[2,1] = -4.34

cov_p = diag(2)
cov_p[1,2] = 0.678
cov_p[2,1] = 0.678


Thingvallavatn2_dist = bhattacharyya.dist(c(-12.8, 4.98), 
                                          c(-28.7, 6.15), 
                                          Sigma1 = cov_b, 
                                          Sigma2 = cov_p)%>% 
  as_tibble()


### 1000 permutation Hotelling T-test
### This is a non-parametric multivariate t-test
t2.benthic_stat = isotopes %>% 
  filter(Lake == 'Thingvallavatn',
         Vector == 'TSBPL',
         BP == 'Benthic') %>%
  select(carbon, nitrogen)

t2.pelagic_stat = isotopes %>% 
  filter(Lake == 'Thingvallavatn',
         Vector == 'TSBPL',
         BP == 'Pelagic') %>%
  select(carbon, nitrogen)

t2.hoteltest = hotelling.test(t2.benthic_stat, 
                              t2.pelagic_stat, 
                              perm = TRUE, 
                              B = 1000)

# Vatnshlidarvatn ---------------------------------------------------------

v.benthic = isotopes %>% 
  filter(Lake == 'Vatnshlidarvatn', 
         BP == 'Benthic') %>%
  summarise(mean_carbon = mean(carbon),
            mean_nitrogen = mean(nitrogen),
            covariance = cov(x = carbon,
                             y = nitrogen,
                             method = 'pearson'))

v.pelagic = isotopes %>% 
  filter(Lake == 'Vatnshlidarvatn', 
         BP == 'Pelagic') %>%
  summarise(mean_carbon = mean(carbon),
            mean_nitrogen = mean(nitrogen),
            covariance = cov(x = carbon,
                             y = nitrogen,
                             method = 'pearson'))

cov_b = diag(2)
cov_b[1,2] = -0.189
cov_b[2,1] = -0.189

cov_p = diag(2)
cov_p[1,2] = 0.0644
cov_p[2,1] = 0.0644


Vatnshlidarvatn_dist = bhattacharyya.dist(c(-23.2, 7.27), 
                                    c(-24.5, 7.45), 
                                    Sigma1 = cov_b, 
                                    Sigma2 = cov_p)%>% 
  as_tibble()


### 1000 permutation Hotelling T-test
### This is a non-parametric multivariate t-test
v.benthic_stat = isotopes %>% 
  filter(Lake == 'Vatnshlidarvatn',
         BP == 'Benthic') %>%
  select(carbon, nitrogen)

v.pelagic_stat = isotopes %>% 
  filter(Lake == 'Vatnshlidarvatn',
         BP == 'Pelagic') %>%
  select(carbon, nitrogen)

v.hoteltest = hotelling.test(v.benthic_stat, 
                              v.pelagic_stat, 
                              perm = TRUE, 
                              B = 1000)


# Combine the data into a single df ---------------------------------------

label = rep('Galtabol', length(Galtabol_dist$value)) %>% 
  as_tibble()
Galtabol_dist = bind_cols(Galtabol_dist, label)

label = rep('Svinavatn', length(Svinavatn_dist$value)) %>% 
  as_tibble()
Svinavatn_dist = bind_cols(Svinavatn_dist, label)

label = rep('Thingvallavatn1', length(Thingvallavatn1_dist$value)) %>% 
  as_tibble()
Thingvallavatn1_dist = bind_cols(Thingvallavatn1_dist, label)

label = rep('Thingvallavatn2', length(Thingvallavatn2_dist$value)) %>% 
  as_tibble()
Thingvallavatn2_dist = bind_cols(Thingvallavatn2_dist, label)

label = rep('Vatnshlidarvatn', length(Vatnshlidarvatn_dist$value)) %>% 
  as_tibble()
Vatnshlidarvatn_dist = bind_cols(Vatnshlidarvatn_dist, label)

dist_df = bind_rows(Galtabol_dist, 
                    Svinavatn_dist)
dist_df = bind_rows(dist_df, 
                    Thingvallavatn1_dist)
dist_df = bind_rows(dist_df, 
                     Thingvallavatn2_dist)
dist_df = bind_rows(dist_df, 
                    Vatnshlidarvatn_dist)

dist_df = dist_df %>% 
  rename(Bhattacharryya_distance = value, 
         Lake = value1) %>% 
  select(Lake, 
         Bhattacharryya_distance)


# write the data ----------------------------------------------------------

write_csv(dist_df, 
          'Bhattacharryya_distances_morph_pairs_isotopes.csv')

read_csv('Bhattacharryya_distances_morph_pairs_isotopes.csv')

# DISTANCE COMPARISON -----------------------------------------------------

df_iso = read_csv('Bhattacharryya_distances_morph_pairs_isotopes.csv')

# df_genetic = read_csv('Neis_genetic_distance_allpops.csv')
df_genetic = read_csv('Clean_Neis_distance_polypops.csv') %>% 
  arrange(Lake)

df_genetic = read_csv('Clean_Provesti_genetic_distance.csv') %>% 
  arrange(Lake)

df_genetic = read_csv('Pairwise_Fst_polypops.csv')

dist_df = left_join(df_iso, 
                    df_genetic, 
                    by = 'Lake')

dist_df[4,2] = 10.0

summary(lm(Neis_distance ~ Bhattacharryya_distance, data = dist_df))
summary(lm(Provesti_distance ~ Bhattacharryya_distance, data = dist_df))
summary(lm(Fst ~ Bhattacharryya_distance, data = dist_df))

summary(lm(Bhattacharryya_distance ~ Neis_distance, 
           data = dist_df))
summary(lm(Bhattacharryya_distance ~ Provesti_distance, 
           data = dist_df))
summary(lm(Bhattacharryya_distance ~ Fst, 
           data = dist_df))

colours = c('#B72939',
            '#372673',
            '#FAD421',
            '#0171A6',
            '#DB520D')

ggplot(data = dist_df, 
       aes(x = Bhattacharryya_distance, 
           y = Fst, 
           col = Lake))+
  geom_point(size = 4)+
  geom_smooth(col = 'black', 
              method = 'lm', 
              se = F, 
              size = 2)+
  scale_color_manual(values = colours)+
  labs(x = 'Ecological distance', 
       y = 'Fst')+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        # axis.line = element_line(size = 1), 
        axis.ticks = element_line(size = 1))

ggsave('Population_Distance_regression.tiff', 
       plot = last_plot(),
       dpi = 'retina',
       units = 'cm')

# Malhanobis distance -----------------------------------------------------

isotopes = read_csv('GSTVMF_Morph_Eco_Geno.csv') %>% 
  filter(Lake %in% c('Galtabol',
                     'Svinavatn', 
                     'Thingvallavatn', 
                     'Vatnshlidarvatn')) %>% 
  select(Specimen.ID, 
         Lake, 
         Morph, 
         BP, 
         Vector, 
         LaMorph, 
         carbon, 
         nitrogen)

base_df = isotopes %>% 
  select(carbon, nitrogen)

distances = mahalanobis(base_df, 
            colMeans(base_df), 
            cov(base_df)) %>% 
  as_tibble() %>% 
  rename(Mahalanobis_dist = value)

clean_df = bind_cols(isotopes, 
                     distances)


# Eco distance vs hybrids -------------------------------------------------

df = read_csv('Bhattacharryya_distances_morph_pairs_isotopes.csv')
# df[4,2] = 10.0

summary(lm(Percentage_hybrids ~ Bhattacharryya_distance, data = df))
## Relationship is not significant at all. 
## The whole model isn't significant and neither is the predictor
## variable. 

summary(lm(Hybrid_percent_90 ~ Bhattacharryya_distance, data = df))