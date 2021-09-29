library(tidyverse)

species_match <- read_csv('splitSpecies/MiSeq_lightSpecies_mini_dapc_all_cluster_pca.csv') %>%
  select(ID, dapc_species) %>%
  mutate(ID = str_remove(ID, '.fp2.repr'),
         dapc_species = if_else(dapc_species == 'hyalinus', 'CHYA', 'CPER'))

all_files <- list.files(path = "NCBI_upload", pattern = '*gz', full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(ID = str_extract(file, 'COPE_[0-9]+'))

rename_scheme <- full_join(all_files, species_match, by = 'ID') %>%
  mutate(species = if_else(is.na(dapc_species), 'CSP', dapc_species),
         new_file = str_replace(file, 'COPE', species)) %>%
  select(-ID, -contains('species'))

with(rename_scheme, file.rename(from = file, to = new_file))