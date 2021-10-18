library(tidyverse)

admix_res <- read_csv('../splitSpecies/ADMIXTURE/MiSeq_lightSpecies2.10.1.2.results.csv') 


admix_res %>%
  select(-admix_assignment) %>%
  rowwise %>%
  mutate(max_prob = max(c_across(starts_with('Cluster'))),
         which_clust = which.max(c_across(starts_with('Cluster')))) %>%
  ungroup %>%
  mutate(admix_id = if_else(round(max_prob, 4) >= 1, 
                            as.character(which_clust - 1), 'U')) %>%
  # filter(admix_id == 'U') %>%
  arrange(-max_prob) %>%
  select(ID, admix_id) %>%
  rename(Indiv = ID) %>%
  count(admix_id)

