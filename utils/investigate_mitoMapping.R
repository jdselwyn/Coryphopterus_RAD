library(tidyverse)


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
  filter(sequencing_type == 'NovaSeq',
         preprocess_step == 'fq_fp1_fp2_fqscrn_repaired') %>%
  select(sample, number_reads)



species_id <- read_csv('../Mitochondrial_Mapping/summarized_blast_results.csv') %>%
  filter(ScientificName != 'Coryphopterus urospilus') %>%
  group_by(ID) %>%
  filter(n == max(n)) %>%
  filter(n() == 1) %>%
  ungroup %>%
  select(ScientificName, ID, n)

full_join(species_id, preprocess_data, by = c('ID' = 'sample')) %>%
  filter(number_reads >= 10000) %>%
  mutate(has_id = !is.na(ScientificName)) %>%
  
  glm(has_id ~ log(number_reads, base = 10), data = ., family = 'binomial') %>%
  summary

prop_test <- function(x){
  a <- 1;b <- 1
  x <- stats::na.omit(x)
  
  tmp_test <- binom.test(x = sum(x), n = length(x))
  
  ggplot2:::new_data_frame(list(y = tmp_test$estimate, ymin = tmp_test$conf.int[1], 
                      ymax = tmp_test$conf.int[2]), n = 1)
}


full_join(species_id, preprocess_data, by = c('ID' = 'sample')) %>%
  filter(number_reads >= 10000) %>%
  mutate(has_id = !is.na(ScientificName)) %>%
  
  ggplot(aes(x = log(number_reads, base = 10), y = as.numeric(has_id))) +
  geom_point(colour = 'grey50', size = 0.7) +
  stat_summary_bin(fun.data = prop_test, bins = 30) +
  geom_smooth(method = 'glm', formula = y ~ x, method.args = list(family = 'binomial')) +
  theme_classic()

