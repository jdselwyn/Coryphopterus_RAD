suppressMessages(library(tidyverse))


hybrid_match <- read_csv('splitSpecies/newHybrids_best/newHybrids_fullResults.csv',
                         col_types = cols(.default = col_double(),
                                          prior = col_character(),
                                          Indiv = col_character(),
                                          Z = col_character(),
                                          hybrid_class = col_character())) %>%
  filter(prior == 'Jeffreys') %>%
  select(-prior) %>%
  mutate(Indiv = str_remove(Indiv, '.fp2.repr')) %>%
  group_by(Indiv) %>%
  mutate(hybrid_type = hybrid_class[mean_posterior_probability == max(mean_posterior_probability)],
         post_prob = max(mean_posterior_probability)) %>%
  ungroup %>%
  select(-starts_with('C'), -sd_posterior_probability, -cv_posterior_probability) %>%
  pivot_wider(names_from = hybrid_class, values_from = mean_posterior_probability) %>%
  select(Indiv, hybrid_type) %>%
  mutate(hybrid_type = case_when(hybrid_type == 'pure_a' ~ 'CPER',
                                 hybrid_type == 'pure_b' ~ 'CHYA',
                                 
                                 hybrid_type == 'a_bx' ~ 'CPERbx',
                                 hybrid_type == 'b_bx' ~ 'CHYAbx',
                                 
                                 hybrid_type == 'F1' ~ 'F1',
                                 hybrid_type == 'F2' ~ 'F2',
                                 TRUE ~ NA_character_))

# species_match <- read_csv('splitSpecies/DAPC/MiSeq_lightSpecies2_dapc_all_cluster_pca.csv') %>%
#   select(ID, dapc_species) %>%
#   mutate(ID = str_remove(ID, '.fp2.repr'),
#          dapc_species = if_else(dapc_species == 'hyalinus', 'CHYA', 'CPER'))

all_files <- list.files(path = "NCBI_upload", pattern = '*gz', full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(ID = str_extract(file, '[CF][A-Z0-9bx]+_[0-9]+'))

rename_scheme %>%
  filter(str_detect(file, '0489'))

rename_scheme <- full_join(all_files, hybrid_match, by = c('ID' = 'Indiv')) %>% 
  mutate(hybrid = if_else(is.na(hybrid_type), 'CSP', hybrid_type),
         new_file = str_replace(file, 'COPE', hybrid)) %>%
  select(-ID, -contains('hybrid'))

with(rename_scheme, file.rename(from = file, to = new_file))

list.files(path = "NCBI_upload", pattern = '*gz', full.names = FALSE) %>%
  tibble(file_name = .) %>%
  write_delim('sample_list.txt', delim = '\t', col_names = FALSE)