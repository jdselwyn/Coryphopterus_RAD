library(tidyverse)
library(magrittr)
library(broom)
library(broom.mixed)
library(emmeans)
library(lme4)
library(afex)

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
  mutate(sample = str_replace(sample, '_', '-')) %>%
  filter(sequencing_type != 'MiSeq') %>%
  filter(preprocess_step == 'fq_fp1_fp2_fqscrn_repaired') %>%
  select(sample, number_reads)


mapping_data <- read_csv('../fltrBAM_MiSeq/sample_bammapping.csv') %>%
  mutate(sample = str_extract(file, 'COPE_[0-9]+') %>% str_replace('_', '-')) %>%
  select(-file)

species_id = read_csv('../splitSpecies/hybrid_types.csv') %>%
  rename(sample = ID) %>%
  mutate(sample = str_extract(sample, 'COPE_[0-9]+') %>% str_replace('_', '-')) 


inner_join(preprocess_data, mapping_data, by = 'sample') %>%
  summarise(across(where(is.numeric), list(mean = mean, sd = sd))) %>%
  mutate(across(ends_with('sd'), ~./sqrt(778)))


inner_join(preprocess_data, mapping_data, by = 'sample') %>%
  inner_join(species_id, by = 'sample') %>%
  ggplot(aes(x = number_reads, y = mapped_reads, colour = hybrid)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)
