##############################
## Phylogenetic tree
##
## Matt Brachmann (PhDMattyB)
##
## 09-03-2021
##
##############################

setwd('~/PhD/SNP Demographic modelling/Working Directory/')

library(tidyverse)
library(poppr)
library(adegenet)
library(ape)
library(devtools)
library(ggtree)

SNP_iceland = read_tsv('ALLPOPS_Plink_input.ped', 
                       col_names = T) %>% 
  select(-PaternalID,
         -MaternalID, 
         -Sex, 
         -Phenotype) %>% 
  arrange(`#FamilyID`)

# Combine S pelagic -------------------------------------------------------

# S.PI = SNP_iceland %>% 
#   filter(`#FamilyID` == '5') %>% 
#   dplyr::select(`#FamilyID`, 
#                 IndividualID)
# 
# value = rep(4, length(S.PI$IndividualID)) %>% 
#   as.numeric() %>% 
#   as_tibble()
# 
# S.PI_id = bind_cols(S.PI, value) %>% 
#   select(-`#FamilyID`) %>% 
#   rename(`#FamilyID` = value)
# 
# S.PI_Geno = SNP_iceland %>% 
#   filter(`#FamilyID` == '5') %>% 
#   select(contains('Affx-'))
# 
# S.PI_data = bind_cols(S.PI_id, 
#                       S.PI_Geno) %>% 
#   select(`#FamilyID`, 
#          IndividualID, 
#          contains('Affx-'))
# 
# SNP_iceland = SNP_iceland %>% 
#   filter(`#FamilyID` != '5')
# 
# SNP_iceland = bind_rows(SNP_iceland, 
#                         S.PI_data)
# back to analysis --------------------------------------------------------

SNP_iceland = mutate(.data = SNP_iceland,
              POP_name = as.factor(case_when(
                `#FamilyID` == "1" ~ "T.LGB",
                `#FamilyID` == '2' ~ 'V.BR',
                `#FamilyID` == '3' ~ 'V.SIL',
                `#FamilyID` == '4' ~ 'S.PL', 
                `#FamilyID` == '5' ~ 'S.PI',
                # `#FamilyID` == '4'~ 'S.PEL',
                # `#FamilyID` == '5' ~ 'S.PEL',
                `#FamilyID` == '6' ~ 'T.PL',
                `#FamilyID` == '7' ~ 'T.SB',
                `#FamilyID` == '8' ~ 'S.LGB',
                `#FamilyID` == '9' ~ 'G.SB',
                `#FamilyID` == '10' ~ 'G.PI',
                `#FamilyID` == '11' ~ 'Mjoavatn',
                `#FamilyID` == '12' ~ 'Fljotaa')))

SNP_iceland = SNP_iceland %>% 
  select(IndividualID, POP_name, contains("Affx-"))

Classifiers = SNP_iceland %>% select(POP_name)

SNP_iceland = as.data.frame(SNP_iceland)

## MAKE DA TREE #####
iceland_genind2 = df2genind(SNP_iceland[,3:length(SNP_iceland)], 
                            ploidy = 2, 
                           ind.names = SNP_iceland[,1], 
                           sep = '\t',
                           strata = Classifiers)

nei_dist = nei.dist(iceland_genind2)

melt(as.matrix(nei_dist), 
     varnames = c('row', 'column')) %>% 
  as_tibble()

## testing out the poppr package to genind2genalex
# test = genind2genalex(iceland_genind2, 
#                       filename = 'AllPops_Genalex', 
#                       allstrata = TRUE)


set.seed(666)
# set.seed(1738)

nameStrata(iceland_genind2) = ~LaMorph

CHARR_TREE = iceland_genind2 %>% 
  genind2genpop(pop = ~LaMorph) %>% 
  aboot(cutoff = 50,
        quiet = FALSE,
        sample = 1000,
        distance = nei.dist,
        tree = 'upgma',
        showtree = F)

nei.dist(CHARR_TREE)

dist.gene(CHARR_TREE)

# str(CHARR_TREE)
write.tree(CHARR_TREE, file = "ALLPOPS_CHARR_TREE_SPELNOCOMBO")

## READ DA TREE ####
FISHTREE = read.tree('ALLPOPS_CHARR_TREE_SPELNOCOMBO')

## GGTREE #####
## read in the data
phylo_data = read_csv('ALLPOPS_NJtree_SPELNOCOMBO.csv') %>% 
  select(parent, 
         node, 
         branch.length, 
         large_label, 
         label3)
## data needs to be class phlyo
phylo_data = as.phylo(Joined_data)
## check the structure to make sure
str(phylo_data)

ggtree(phylo_data, ladderize = T) +
  #geom_text(aes(label= boot), hjust=-.3)+
  geom_text2(aes(label = label,
                 subset = !is.na(as.numeric(label)) &
                   as.numeric(label) > 50),
             vjust = 1.2,
             hjust = 1.2)+
  theme_tree2()