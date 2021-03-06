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

newhybrids_bestLoci <- read_csv('../splitSpecies/newHybrids_best/newHybrids_fullResults.csv') %>%
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
mixup_ind_all <- newhybrids_subsetLoci %>%
  
  full_join(newhybrids_allLoci %>% 
              select(-starts_with('L')), 
            by = c('Indiv', 'Z'),
            suffix = c('.best', '.all')) %>%
  filter(hybrid_type.all != hybrid_type.best) %>%
  pull(Indiv)

all_nan <- newhybrids_allLoci %>%
  filter(is.nan(post_prob)) %>%
  pull(Indiv)

full_join(
  select(newhybrids_subsetLoci, Indiv, hybrid_type, post_prob),
  select(newhybrids_allLoci, Indiv, hybrid_type, post_prob),
  by = c('Indiv'),
  suffix = c('.subset', '.all')
) %>%
  filter(hybrid_type.all != hybrid_type.subset)

#### Compare 200 loci vs Best loci ####
mixup_ind_best <- newhybrids_subsetLoci %>%
  
  full_join(newhybrids_bestLoci, 
            by = c('Indiv', 'Z'),
            suffix = c('.subset', '.best')) %>%
  select(-Z) %>%
  filter(hybrid_type.best != hybrid_type.subset) %>%
  pull(Indiv)

mixup_best_all <- newhybrids_bestLoci %>%
  
  full_join(newhybrids_allLoci %>% 
              select(-starts_with('L')), 
            by = c('Indiv', 'Z'),
            suffix = c('.best', '.all')) %>%
  filter(hybrid_type.all != hybrid_type.best) %>%
  pull(Indiv)

#### PCA ####
library(adegenet)
library(vcfR)
vcf <- read.vcfR('../tmp_dir/MiSeq_lightSpecies2.10.1.Fltr20.8.randSNPperLoc.vcf', verbose = TRUE)

vcf_genind <- vcfR2genind(vcf)

full_pca <- dudi.pca(tab(vcf_genind, NA.method = "mean"), scannf=FALSE, 
                     center = TRUE, scale=TRUE, nf = nLoc(vcf_genind) - 1)

#### Join Data ####
full_data <- newhybrids_subsetLoci %>%
  full_join(newhybrids_bestLoci, 
            by = c('Indiv'),
            suffix = c('.subset', '.best')) %>%
  full_join(select(newhybrids_allLoci, -starts_with('L')) %>%
              rename_with(~str_c(., '.all'), -c(Indiv)),
            by = c('Indiv')) %>%
  left_join(mtdna, by = c('Indiv' = 'ID')) %>%
  left_join(preprocess_data, by = c('Indiv' = 'sample')) %>%
  left_join(admix_data, by = c('Indiv' = 'ID')) %>%
  mutate(mixup_mt_best = case_when(is.na(mt_species) ~ NA,
                           hybrid_type.best == 'pure_b' & mt_species == 'Coryphopterus hyalinus' ~ FALSE,
                           hybrid_type.best == 'pure_a' & mt_species == 'Coryphopterus personatus' ~ FALSE,
                           hybrid_type.best == 'pure_a' & mt_species == 'Coryphopterus hyalinus' ~ TRUE,
                           hybrid_type.best == 'pure_b' & mt_species == 'Coryphopterus personatus' ~ TRUE,
                           TRUE ~ NA),
         mixup_mt_sub = case_when(is.na(mt_species) ~ NA,
                                   hybrid_type.subset == 'pure_b' & mt_species == 'Coryphopterus hyalinus' ~ FALSE,
                                  hybrid_type.subset == 'pure_a' & mt_species == 'Coryphopterus personatus' ~ FALSE,
                                  hybrid_type.subset == 'pure_a' & mt_species == 'Coryphopterus hyalinus' ~ TRUE,
                                  hybrid_type.subset == 'pure_b' & mt_species == 'Coryphopterus personatus' ~ TRUE,
                                   TRUE ~ NA),
         mixup_mt_all = case_when(is.na(mt_species) ~ NA,
                                  hybrid_type.all == 'pure_b' & mt_species == 'Coryphopterus hyalinus' ~ FALSE,
                                  hybrid_type.all == 'pure_a' & mt_species == 'Coryphopterus personatus' ~ FALSE,
                                  hybrid_type.all == 'pure_a' & mt_species == 'Coryphopterus hyalinus' ~ TRUE,
                                  hybrid_type.all == 'pure_b' & mt_species == 'Coryphopterus personatus' ~ TRUE,
                                  TRUE ~ NA),
         nan_all = is.nan(post_prob.all),
         mixup_all_sub = hybrid_type.subset != hybrid_type.all,
         mixup_best_sub = hybrid_type.best != hybrid_type.subset,
         mixup_best_all = hybrid_type.best != hybrid_type.all) %>%
  inner_join(as_tibble(full_pca$li, rownames = 'Indiv') %>%
               mutate(Indiv = str_remove(Indiv, '.fp2.repr')) %>%
               select(Indiv, Axis1, Axis2), 
             by = 'Indiv')

#### Plot of interesting/disagreements ####
full_data %>%
  select(Indiv, Axis1, Axis2, contains('hybrid_type')) %>%
  mutate(across(starts_with('hybrid_type'), is.na)) %>%
  pivot_longer(cols = starts_with('hybrid_type')) %>%
  mutate(order_var = if_else(!value | is.na(value), 1, 2),
         name = str_remove(name, 'hybrid_type.') %>%
           str_replace('best', 'Top Fst - 200') %>%
           str_replace('subset', 'random - 200') %>%
           str_replace('all', 'all - 1,726') %>%
           str_to_title()) %>%
  arrange(order_var) %>%
  filter(Axis2 > -20) %>%
  ggplot(aes(x = Axis1, y = Axis2, colour = value, size = value)) +
  geom_point(show.legend = FALSE) +
  scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'gray25'), na.value = 'grey75') +
  scale_size_manual(values = c('TRUE' = 3, 'FALSE' = 1), na.value = 1) +
  facet_wrap(~name) +
  theme_classic() +
  theme(text = element_text(size = 32),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = 'transparent', colour = 'black'),
        strip.background = element_blank())


full_data %>%
  select(Indiv, Axis1, Axis2, nan_all, starts_with('mixup')) %>%
  pivot_longer(cols = c(nan_all, starts_with('mixup'))) %>%
  mutate(order_var = if_else(!value | is.na(value), 1, 2)) %>%
  arrange(order_var) %>%
  filter(Axis2 > -20) %>%
  
  ggplot(aes(x = Axis1, y = Axis2, colour = value, size = value)) +
  geom_point(show.legend = FALSE) +
  scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'gray25'), na.value = 'grey75') +
  scale_size_manual(values = c('TRUE' = 3, 'FALSE' = 1), na.value = 1) +
  facet_wrap(~name) +
  theme_classic()


full_data %>%
  filter(mixup_mt_best | mixup_mt_sub | mixup_mt_all) %>%
  select(Indiv, starts_with('hybrid_type'), starts_with('mixup'))

#### Investigate Mismatches ####
inspect <- full_data %>%
  select(Indiv, Axis1, Axis2, nan_all, starts_with('mixup')) %>%
  pivot_longer(cols = c(nan_all, starts_with('mixup'))) %>%
  mutate(order_var = if_else(!value | is.na(value), 1, 2)) %>%
  arrange(order_var) %>%
  filter(Axis2 > -20,
         str_detect(name, 'mixup'),
         str_detect(name, 'mt', negate = TRUE)) %>%
  mutate(name = str_replace_all(name, c('mixup_all_sub' = 'All v Random',
                                        'mixup_best_all' = 'All v Top Fst',
                                        'mixup_best_sub' = 'Top Fst v Random'))) %>%
  filter(value) %>%
  
  filter(str_detect(name, 'All', negate = TRUE)) %>%
  
  select(Indiv, name) %>%
  mutate(exisit = 1) %>%
  pivot_wider(names_from = 'name',
              values_from = 'exisit') %>%
  pull(Indiv)


full_data %>%
  filter(Indiv %in% inspect) %>%
  select(Indiv, ends_with('subset'), ends_with('all'), ends_with('best')) %>%
  select(-starts_with('Z'), -starts_with('mixup'), -nan_all) %>%
  pivot_longer(cols = -c(Indiv),
               names_to = c('.value', "loci"),
               names_pattern = '(.*)\\.(.*)') %>%
  pivot_longer(cols = pure_a:b_bx) %>%
  
  # filter(Indiv %in% c('COPE_0489', "COPE_0502")) %>%
  
  ggplot(aes(x = loci, y = value, fill = name)) +
  geom_col() +
  facet_wrap(~ Indiv)


#### Plot Hybrid types ####
full_data %>% 
  select(Indiv, Axis1, Axis2, starts_with('hybrid_type')) %>%
  pivot_longer(cols = starts_with('hybrid_type')) %>%
  mutate(name = str_remove(name, 'hybrid_type.')) %>%
  filter(Axis2 > -20) %>%
  filter(name == 'best') %>%
  ggplot(aes(x = Axis1, y = Axis2, colour = value, size = value)) +
  geom_point(show.legend = TRUE) +
  scale_colour_manual(values = c('a_bx' = 'red', 'b_bx' = 'green', 
                                 'F1' = 'blue', 'F2' = 'purple',
                                 'pure_a' = 'gray10', 'pure_b' = 'gray70'), 
                      na.value = 'pink') +
  scale_size_manual(values = c('a_bx' = 2, 'b_bx' = 2, 
                               'F1' = 2, 'F2' = 2,
                               'pure_a' = 1, 'pure_b' = 1), 
                    na.value = 2) +
  
  # facet_wrap(~name) +
  theme_classic()

full_data %>%
  select(Indiv, ends_with('best')) %>%
  count(hybrid_type.best)

full_data %>%
  select(Indiv, mt_species, ends_with('best')) %>%
  filter(!is.na(mt_species),
         str_detect(hybrid_type.best, 'pure')) %>%
  count(hybrid_type.best, mt_species)

full_data %>%
  select(Indiv, ends_with('best')) %>%
  filter(mixup_mt_best)


full_data %>%
  select(Indiv, ends_with('best')) %>%
  filter(is.na(hybrid_type.best))


full_data %>%
  select(Indiv, ends_with('best')) %>%
  select(-Z.best)

full_data %>%
  select(Indiv, ends_with('best')) %>%
  arrange(post_prob.best)

(4+4+3+7) / (4+4+3+7+133+625)



full_data %>%
  filter(!is.na(MiSeq)) %>%
  select(Indiv, MiSeq, ends_with('best'))

#### Output dataframe with ID/Hybrid type ####
full_data %>%
  select(Indiv, hybrid_type.best) %>%
  rename(ID = Indiv,
         hybrid = hybrid_type.best) %>%
  mutate(hybrid = case_when(hybrid == 'pure_a' ~ 'CPER',
                            hybrid == 'pure_b' ~ 'CHYA',
                            
                            hybrid == 'a_bx' ~ 'CPERbx',
                            hybrid == 'b_bx' ~ 'CHYAbx',
                            
                            hybrid == 'F1' ~ 'F1',
                            hybrid == 'F2' ~ 'F2',
                            TRUE ~ NA_character_),
         ID = str_c(ID, '.fp2.repr')) %>%
  write_csv('../splitSpecies/hybrid_types.csv')
