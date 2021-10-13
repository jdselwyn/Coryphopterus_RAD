library(tidyverse)

mtdna <- read_csv('../Mitochondrial_Mapping/blast_speciesID.csv') 
admix_data <- read_csv('../splitSpecies/ADMIXTURE/MiSeq_lightSpecies2.10.1.2.results.csv')
newhybrids_subsetLoci <- read_csv('../splitSpecies/newHybrids/newHybrids_fullResults.csv')


newhybrids_subsetLoci %>%
  filter(is.na(Z),
         prior == 'Jeffreys') %>%
  filter(Indiv == 'COPE_1023.fp2.repr')

newhybrids_subsetLoci %>%
  filter(prior == 'Jeffreys') %>%
  select(-prior) %>%
  group_by(Indiv) %>%
  mutate(hybrid_type = hybrid_class[mean_posterior_probability == max(mean_posterior_probability)],
         post_prob = max(mean_posterior_probability)) %>%
  ungroup %>%
  select(-starts_with('C'), -sd_posterior_probability, -cv_posterior_probability) %>%
  pivot_wider(names_from = hybrid_class, values_from = mean_posterior_probability) %>%
  arrange(post_prob) %>%
  filter(str_detect(hybrid_type, 'pure', negate = TRUE))


tmp <- newhybrids_subsetLoci %>%
  filter(prior == 'Jeffreys') %>%
  select(-prior) %>%
  mutate(Indiv = str_remove(Indiv, '.fp2.repr')) %>%
  group_by(Indiv) %>%
  mutate(hybrid_type = hybrid_class[mean_posterior_probability == max(mean_posterior_probability)],
         post_prob = max(mean_posterior_probability)) %>%
  ungroup %>%
  select(-starts_with('C'), -sd_posterior_probability, -cv_posterior_probability) %>%
  pivot_wider(names_from = hybrid_class, values_from = mean_posterior_probability) %>%
  left_join(mtdna, by = c('Indiv' = 'ID')) 


tmp %>%
  filter(evalue < -65) %>%
  filter(!is.na(species),
         str_detect(hybrid_type, 'pure')) %>%
  # count(species, hybrid_type)
  mutate(match = case_when(hybrid_type == 'pure_b' & species == 'Coryphopterus hyalinus' ~ 'match_chya',
                           hybrid_type == 'pure_a' & species == 'Coryphopterus personatus' ~ 'match_cper',
                           hybrid_type == 'pure_a' & species == 'Coryphopterus hyalinus' ~ 'mtDNA_chya-nuc_cper',
                           hybrid_type == 'pure_b' & species == 'Coryphopterus personatus' ~ 'mtDNA_cper-nuc_chya',
                           TRUE ~ 'error')) %>%
  
  ggplot(aes(x = match, y = bitscore)) +
  geom_boxplot()

  
  