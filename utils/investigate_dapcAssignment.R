library(tidyverse)
library(magrittr)
library(ggdist)
library(broom)
library(ggbeeswarm)
library(emmeans)
count <- dplyr::count

#### Get Data ####
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
  select(sequencing_type, sample, number_reads)

mitochondrial_id <- read_csv('../Mitochondrial_Mapping/blast_speciesID.csv') 

# 'MiSeq_lightSpecies_mini_dapc_all_cluster_pca.csv'
dapc_id <- read_csv('../splitSpecies/MiSeq_lightSpecies_mini_dapc_all_cluster_pca.csv') %>%
  mutate(ID = str_remove(ID, '.fp2.repr'))

full_data <- full_join(dapc_id, preprocess_data, by = c('ID' = 'sample')) %>%
  filter(number_reads >= 10000) %>%
  mutate(number_reads = log(number_reads, base = 10)) %>%
  
  group_by(ID) %>%
  mutate(has_miseq = any(sequencing_type == 'MiSeq')) %>%
  ungroup %>%
  filter(sequencing_type != 'MiSeq') %>%
  left_join(select(mitochondrial_id, -species), by = 'ID')

#### Where are MiSeq specimens in PCA space ####
ggplot(full_data, aes(x = Axis1, y = Axis2, colour = has_miseq)) +
  geom_point() +
  ylim(c(-10,NA))

#### Do Reads relate to species? ####
full_data %>%
  filter(str_detect(id_match, 'pred', negate = TRUE)) %>%
  mutate(id_match = if_else(species == dapc_species, 'match', 'fail')) %>%
  select(-starts_with('Axis')) %>%
  mutate(success = as.integer(id_match == 'match')) %>%
  
  glm(success ~ number_reads, data = ., family = 'binomial') %>% 
  summary

full_data %>%
  count(dapc_species, species) %>%
  filter(!is.na(species))


full_data %>%
  filter(!is.na(species), 
         dapc_species != species) %>%
  select(-starts_with('Axis')) %>%
  ggplot(aes(x = certainty, y = evalue)) +
  geom_point()


#### Make List of C. hyalinus to keep in VCF ####
full_data %>%
  filter(dapc_species == 'hyalinus', 
         species != 'personatus' | is.na(species))  %>%
  select(ID) %>%
  mutate(ID = str_c(ID, '.fp2.repr')) %>%
  pull(ID) %>% 
  write_lines('../splitSpecies/CHYA.list')


#### Get most extreme individuals which agree on mito & nuc and find diagnostic SNPs ####
full_data %>%
  select(-starts_with('Axis')) %>%
  filter(dapc_species == species)
