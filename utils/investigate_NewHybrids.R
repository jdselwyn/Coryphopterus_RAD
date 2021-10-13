library(tidyverse)


#### Read in Data ####
preprocess_data <- list.files('..', pattern = 'sample_preprocess.csv$', 
                              recursive = TRUE, full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(sequencing_type = str_extract(file, 'MiSeq|NovaSeq'),
         preprocess_step = str_extract(file, '(/[a-z0-9_]+/){1}') %>% str_remove_all('/')) %>%
  filter(preprocess_step != 'fq_fp1_clmp_fp2_fqscrn') %>% #not needed since repair is the same
  rowwise(sequencing_type, preprocess_step) %>%
  summarise(read_csv(file, col_types = cols(.default = col_character(), 
                                            number_reads = col_integer())), 
            .groups = 'drop') %>%
  mutate(preprocess_step = factor(preprocess_step, levels = c("demultiplexed_seqs", "fq_fp1", "fq_fp1_fp2", "fq_fp1_fp2_fqscrn_repaired"), ordered = TRUE)) %>%
  filter(preprocess_step == 'fq_fp1_fp2_fqscrn_repaired') %>%
  select(sequencing_type, sample, number_reads) %>%
  mutate(number_reads = if_else(number_reads < 10000, NA_integer_, number_reads)) %>%
  na.omit() %>%
  pivot_wider(names_from = 'sequencing_type',
              values_from = 'number_reads')


mtdna <- read_csv('../Mitochondrial_Mapping/blast_speciesID.csv') %>%
  rename(mt_species = species)

admix_data <- read_csv('../splitSpecies/ADMIXTURE/MiSeq_lightSpecies2.10.1.2.results.csv') %>%
  mutate(ID = str_remove(ID, '.fp2.repr')) %>%
  rowwise %>%
  mutate(admix_post_prob = max(c_across(starts_with('Cluster')))) %>%
  ungroup

newhybrids_subsetLoci <- read_csv('../splitSpecies/newHybrids/newHybrids_fullResults.csv') %>%
  filter(prior == 'Jeffreys') %>%
  select(-prior) %>%
  mutate(Indiv = str_remove(Indiv, '.fp2.repr')) %>%
  group_by(Indiv) %>%
  mutate(hybrid_type = hybrid_class[mean_posterior_probability == max(mean_posterior_probability)],
         post_prob = max(mean_posterior_probability)) %>%
  ungroup %>%
  select(-starts_with('C'), -sd_posterior_probability, -cv_posterior_probability) %>%
  pivot_wider(names_from = hybrid_class, values_from = mean_posterior_probability) 


newhybrids_allLoci <- full_join(read_csv('../splitSpecies/newHybrids/All Loci/newhybrid_data.csv') %>%
            # select(-starts_with('L')) %>%
            mutate(id = row_number()), 
          read_delim('../splitSpecies/newHybrids/All Loci/aa-PofZ.txt', delim = '\t', 
                           col_types = cols(.default = col_double(),
                                            tmp = col_character()),
                           skip = 1,
                           col_names = c('id', 'tmp', 'pure_a', 'pure_b', 'F1', 'F2', 'a_bx', 'b_bx')) %>%
            select(-tmp), 
          by = 'id') %>% 
  mutate(Z = str_extract(ID, 'z[0-1]')) %>%
  select(Indiv, Z, everything()) %>%
  select(-id, -ID) %>%
  mutate(Indiv = str_remove(Indiv, '.fp2.repr')) %>%
  pivot_longer(cols = -c(Indiv, Z, starts_with('L'))) %>%
  group_by(Indiv) %>%
  mutate(post_prob = max(value),
         hybrid_type = name[value == max(value)]) %>%
  ungroup %>%
  pivot_wider()

#### Compare 200 loci vs All loci ####
mixup_ind <- newhybrids_subsetLoci %>%
  pivot_longer(cols = -c(Indiv, Z, hybrid_type, post_prob)) %>%
  
  full_join(newhybrids_allLoci %>% 
              select(-starts_with('L')) %>%
              pivot_longer(cols = -c(Indiv, Z, hybrid_type, post_prob)), 
            by = c('Indiv', 'Z', 'name'),
            suffix = c('.subset', '.all')) %>%
  mutate(diff = abs(value.subset - value.all)) %>%
  filter(hybrid_type.all != hybrid_type.subset) %>%
  
  group_by(Indiv) %>%
  filter(value.subset != post_prob.subset) %>%
  filter(value.subset == max(value.subset)) %>%
  pull(Indiv)


full_join(
  select(newhybrids_subsetLoci, Indiv, hybrid_type, post_prob),
  select(newhybrids_allLoci, Indiv, hybrid_type, post_prob),
  by = c('Indiv'),
  suffix = c('.subset', '.all')
) %>%
  filter(hybrid_type.all != hybrid_type.subset)


read_csv('../splitSpecies/newHybrids/newHybrids_fullResults.csv') %>%
  filter(prior == 'Jeffreys') %>%
  select(-prior) %>%
  mutate(Indiv = str_remove(Indiv, '.fp2.repr')) %>%
  filter(Indiv %in% c(tmp))
  

newhybrids_allLoci %>%
  filter(Indiv %in% tmp) %>%
  rowwise %>%
  mutate(missing = rowSums(across(starts_with('L'), ~. == '0'))) %>%
  ungroup %>%
  select(Indiv, missing) %>%
  mutate(missing = missing / length(colnames(select(newhybrids_allLoci, starts_with('L')))) * 100) %>%
  arrange(-missing)


  
  
#### Join Data ####
full_data <- newhybrids_allLoci %>%
  left_join(mtdna, by = c('Indiv' = 'ID')) %>%
  left_join(preprocess_data, by = c('Indiv' = 'sample')) %>%
  left_join(admix_data, by = c('Indiv' = 'ID')) %>%
  mutate(match = case_when(is.na(mt_species) ~ NA_character_,
                           hybrid_type == 'pure_b' & mt_species == 'Coryphopterus hyalinus' ~ 'match_chya',
                           hybrid_type == 'pure_a' & mt_species == 'Coryphopterus personatus' ~ 'match_cper',
                           hybrid_type == 'pure_a' & mt_species == 'Coryphopterus hyalinus' ~ 'mtDNA_chya-nuc_cper',
                           hybrid_type == 'pure_b' & mt_species == 'Coryphopterus personatus' ~ 'mtDNA_cper-nuc_chya',
                           TRUE ~ 'hybrid')) %>%
  mutate(hybrid = case_when(hybrid_type == 'pure_b' ~ 'chya',
                            hybrid_type == 'pure_a' ~ 'cper',
                            hybrid_type == 'a_bx' ~ 'cper_bx',
                            hybrid_type == 'b_bx' ~ 'chya_bx',
                            TRUE ~ hybrid_type)) 
  
full_data %>%
  ggplot(aes(x = hybrid_type, y = NovaSeq)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(y = 'Number of Reads')
  
full_data %>%
  filter(!is.na(mt_species)) %>%
  ggplot(aes(x = match, y = NovaSeq)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(y = 'Number of Reads')

full_data %>%
  filter(!is.na(mt_species)) %>%
  ggplot(aes(x = match, y = evalue)) +
  geom_boxplot() +
  labs(y = 'log10(E-Value)')

full_data %>%
  count(match)

full_data %>%
  filter(is.na(hybrid))
  count(hybrid)
