##############################
## Levenes test for variation in stable isotope signatures
##
## Matt Brachmann (PhDMattyB)
##
## 2020-08-25
##
##############################

# working directory -------------------------------------------------------
setwd('~/PhD/SNP Demographic Modelling/Working Directory')

# install packages --------------------------------------------------------
library(patchwork)
library(car)
library(tidyverse)

# read in the data --------------------------------------------------------

le_data = read_csv('GSTVMF_Morph_Eco_Geno.csv')


# data cleaning -----------------------------------------------------------
poly_mono = le_data %>% 
  # filter(LaMorph != 'S.LGB2',
  #        LaMorph != 'T.PL2',
  #        Lake %in% c('Svinavatn', 'Mjoavatn')) %>% 
  # arrange(desc(Lake))
  filter(LaMorph != 'S.LGB2',
         LaMorph != 'T.PL2')

poly_mono = mutate(.data = poly_mono, 
                   Morph_type = as.factor(case_when(
                     LaMorph == "G.SB" ~ "G: Benthic",
                     LaMorph == 'G.PI' ~ 'G: Pelagic',
                     LaMorph == 'V.BR' ~ 'V: Pelagic',
                     LaMorph == 'V.SIL'~ 'V: Benthic',
                     LaMorph == 'T.LGB' ~ 'T: Benthic (1)',
                     LaMorph == 'T.SB' ~ 'T: Benthic (2)', 
                     LaMorph == 'T.PL' ~ 'T: Pelagic',
                     LaMorph == 'S.LGB' ~ 'S: Benthic',
                     LaMorph == 'S.PI' ~ 'S: Pelagic',
                     LaMorph == 'S.PL' ~ 'S: Pelagic', 
                     LaMorph == 'F.ana' ~ 'F: Pelagic',
                     LaMorph == 'M.single' ~ 'M: Benthic')))


# calculate st.dev --------------------------------------------------------
st.dev = poly_mono %>% 
  group_by(Lake) %>% 
  summarise(sd_carbon = sd(carbon),
            sd_nitrogen = sd(nitrogen))


# perfom levenetest -------------------------------------------------------
leveneTest(poly_mono$carbon, 
           poly_mono$Lake, 
           data = poly_mono)

leveneTest(poly_mono$nitrogen, 
           poly_mono$Lake, 
           data = poly_mono)


# leventest svinavatn pelagic ---------------------------------------------

Svin_pelagic = poly_mono %>% 
  filter(Lake == 'Svinavatn', 
         Morph %in% c('Piscivorous', 
                      'Planktivorous'))

leveneTest(Svin_pelagic$carbon, 
           Svin_pelagic$Morph, 
           data = Svin_pelagic)

# Graph that shit! --------------------------------------------------------
theme_set(theme_bw())

# colour palette
box_palette = c('#ED1F17',
                '#009EFA',
                '#1A1D45',
                '#FA4613',
                '#B0B0EB',
                '#02B04D')


## Need to manipulate the box plot levels
box_plot1 = poly_mono
box_plot1$Lake = factor(box_plot1$Lake, 
                        levels = c('Galtabol',
                                   'Svinavatn', 
                                   'Thingvallavatn', 
                                   'Vatnshlidarvatn', 
                                   'Fljotaa', 
                                   'Mjoavatn'))

b1 = ggplot(data = box_plot1, 
            aes(x = Lake, 
                y = carbon))+
  geom_boxplot(aes(fill = Lake), 
               col = 'black')+
  scale_fill_manual(values = box_palette)+
  labs(x = 'Population', 
       y = 'Carbon 13 (???)', 
       title = 'A)')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        legend.position = 'none', 
        plot.title = element_text(size = 12, 
                                  hjust = 0), 
        axis.ticks.y = element_line(size = 1)) 

## Need to manipulate the levels for the second boxplot
box_plot2 = poly_mono
box_plot2$Lake = factor(box_plot2$Lake, 
                        levels = c('Galtabol',
                                   'Svinavatn', 
                                   'Thingvallavatn', 
                                   'Vatnshlidarvatn', 
                                   'Fljotaa', 
                                   'Mjoavatn'))

b2 = ggplot(data = box_plot2, 
            aes(x = Lake, 
                y = nitrogen))+
  geom_boxplot(aes(fill = Lake), 
               col = 'black')+
  scale_fill_manual(values = box_palette)+
  labs(x = 'Population', 
       y = 'Nitrogen 15(???)', 
       title = 'B)')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        legend.position = 'none',
        plot.title = element_text(size = 12,
                                  hjust = 0), 
        axis.ticks.y = element_line(size = 1)) 


combo_plot = b1+b2

ggsave('Fig.3_Isovar_24.08.2020.tiff', 
       units = 'cm', 
       # width = 50, 
       # height = 20, 
       plot = combo_plot)