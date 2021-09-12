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
  mutate(sample = str_replace(sample, '_', '-') %>% str_c(sequencing_type, sep = '_')) 

contamination_data <- list.files(path = '..', pattern = "screen.txt$", recursive = TRUE, full.names = TRUE) %>%
  str_subset('multiqc', negate = TRUE) %>%
  tibble(file = .) %>%
  rowwise %>%
  summarise(sequencing_type = str_extract(file, 'MiSeq|NovaSeq'),
            individual = str_extract(file, 'COPE-[0-9]+|Blank-[0-9]'),
            direction = str_extract(file, 'r[12]'),
            percent_no_genome = read_tsv(file, skip = 14, n_max = 1, col_names = 'no_hit', 
                                         col_types = cols(no_hit = col_character())) %>%
              pull(no_hit) %>%
              str_extract('[0-9\\.]+') %>%
              as.numeric() %>%
              divide_by(100),
            read_tsv(file, skip = 1, col_types = cols(.default = col_double(), Genome = col_character()), n_max = 12),
            .groups = 'drop')

#### Both Sequencing Types alone ####
model_out <- preprocess_data %>%
  select(-file) %>%
  group_by(sequencing_type) %>%
  summarise(model = list(glmer(number_reads ~ preprocess_step + (1 | sample), 
                               family = poisson(link = 'log'))),
            .groups = 'rowwise') %>%
  mutate(#model_out = tidy(car::Anova(model)),
         plotting = list(emmeans(model, ~preprocess_step, type = 'response') %>% 
                           multcomp::cld(alpha = 0.10, Letters = LETTERS) %>% tidy(conf.int = TRUE) %>% 
                           mutate(.group = str_trim(.group))),
         pairs = list(emmeans(model, ~preprocess_step) %>% pairs %>% tidy)) %>%
  ungroup 

model_out$model[[1]] %>% tidy
model_out$model[[2]] %>% tidy

model_out %>%
  select(sequencing_type, plotting) %>%
  unnest(plotting) %>%
  mutate(letter_shift = if_else(sequencing_type == 'MiSeq', 1.1 * asymp.UCL, 1.1 * asymp.UCL)) %>%

  ggplot(aes(x = preprocess_step, y = rate)) +
  
  geom_path(data = preprocess_data,
            aes(y = number_reads, group = sample), 
            position = position_jitter(width = 0.25, seed = 123), 
            colour = 'grey70') +
  geom_point(data = preprocess_data,
             aes(y = number_reads),
             position = position_jitter(width = 0.25, seed = 123), 
             size = 0.1, colour = 'grey40') +
  
  geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL), size = 1) +
  
  geom_text(aes(y = letter_shift, label = .group), vjust = 0,
            position = position_dodge(0.2), show.legend = FALSE) +
  scale_y_log10(labels = scales::comma_format()) +
  facet_wrap(~sequencing_type, scales = 'fixed') +
  theme_classic()


preprocess_data %>%
  count(sequencing_type, preprocess_step, sample) %>%
  filter(n > 1)

preprocess_data %>%
  filter(sequencing_type == 'NovaSeq') %>%
  count(sample) %>%
  filter(n < 3)


#### Complete Model ####
full_model <- mixed(number_reads ~ sequencing_type * preprocess_step + (1 | sample),
                    data = preprocess_data, 
                    family = poisson(link = 'log'),
                    method = 'LRT', 
                    all_fit = FALSE)
full_model


emmeans(full_model, ~preprocess_step | sequencing_type, type = 'link') %>% contrast("poly")

emmeans(full_model, ~sequencing_type, type = 'link') %>% contrast("pairwise")

emmeans(full_model, ~preprocess_step | sequencing_type, type = 'response') %>% 
  multcomp::cld(alpha = 0.1, Letters = LETTERS) %>%
  tidy(conf.int = TRUE) %>%
  mutate(.group = str_trim(.group),
         letter_shift = if_else(sequencing_type == 'MiSeq', 1.1 * asymp.UCL, 1.01 * asymp.UCL)) %>%
  
  ggplot(aes(x = preprocess_step, y = rate)) +
  
  geom_path(data = preprocess_data,
            aes(y = number_reads, group = sample), 
            position = position_jitter(width = 0.25, seed = 123), 
            colour = 'grey70') +
  geom_point(data = preprocess_data,
             aes(y = number_reads),
             position = position_jitter(width = 0.25, seed = 123), 
             size = 0.1, colour = 'grey40') +
  
  geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge(0.2)) +
  geom_text(aes(y = letter_shift, label = .group), vjust = 0,
            position = position_dodge(0.2), show.legend = FALSE) +
  scale_y_log10(labels = scales::comma_format()) +
  facet_wrap(~sequencing_type) +
  theme_classic()


#### Decide which Samples to toss ####
## - samples with too little sequencing effort or too much contamination

## MiSeq ##
library(outliers)

preprocess_data %>%
  select(-file) %>%
  filter(preprocess_step == 'fq_fp1_fp2_fqscrn_repaired') %>%
  
  filter(!number_reads %in% c(6018)) %>%
  
  summarise(tidy(dixon.test(number_reads)))
  


#NovaSeq
preprocess_data %>%
  filter(sequencing_type == 'NovaSeq',
         preprocess_step == 'fq_fp1_fp2_fqscrn_repaired') %>%
  filter(str_detect(sample, 'Blank', negate = TRUE)) %>%
  arrange(number_reads) %>%
  
  filter(!(number_reads >= 10000))



samples_to_exclude <- preprocess_data %>%
  select(-file) %>%
  mutate(sample = str_remove(sample, str_c('_', sequencing_type))) %>%
  filter(sequencing_type == 'NovaSeq') %>%
  filter(preprocess_step == 'fq_fp1_fp2_fqscrn_repaired') %>%
  arrange(number_reads) %>%
  filter(number_reads < 10000 | str_detect(sample, 'Blank')) 




samples_to_exclude %>%
  mutate(sample = str_replace(sample, '-', '_')) %>%
  mutate(command = str_c('mv mkREF_NovaSeq/', sample, '* NovaSeq/fq_fp1_fp2_fqscrn_repaired')) %>%
  pull(command) %>%
  str_c(collapse = '; ')











