library(tidyverse)
library(magrittr)
library(seqinr)
library(ggdist)
library(broom)
library(ggbeeswarm)
library(emmeans)
count <- dplyr::count

#### Get Data ####
reference_sequences <- read.fasta('../Reference_Sequence/Coryphopterus_mtDNA.fasta') %>%
  tibble(sseqid = names(.),
         sequence = .) %>%
  rowwise %>%
  mutate(ref_length = length(sequence),
         extra = attr(sequence,"Annot"),
         sequence = str_c(sequence, collapse = '')) %>%
  mutate(species = str_extract(extra, 'hyalinus|personatus')) %>%
  select(sseqid, species, ref_length, sequence)

reference_sequences %>%
  dplyr::count(species)

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

species_id <- read_csv('../Mitochondrial_Mapping/blast_speciesID.csv') 

full_data <- full_join(species_id, preprocess_data, by = c('ID' = 'sample')) %>%
  filter(number_reads >= 10000) %>%
  mutate(has_id = as.integer(!is.na(species)),
         number_reads = log(number_reads, base = 10),
         binary = factor(has_id, levels = c("0", "1")))

filter(full_data, sequencing_type == 'MiSeq')

#### See how ID relates to number of reads ####
fit_model <- glm(has_id ~ number_reads, family = 'binomial',
                 data = filter(full_data, sequencing_type == 'NovaSeq'))
summary(fit_model)

# compute the fitted lines and SE's
prediction_out <- predict(fit_model,
                          newdata = tibble(number_reads = modelr::seq_range(full_data$number_reads, 100)),
                          type = "link",
                          se.fit = TRUE) %>% 
  # wrangle
  data.frame() %>% 
  as_tibble %>%
  mutate(.lower = fit - 1.96 * se.fit,
         .upper = fit + 1.96 * se.fit) %>% 
  select(-residual.scale, -se.fit) %>% 
  mutate_all(plogis) %>%
  bind_cols(tibble(number_reads = modelr::seq_range(full_data$number_reads, 100)))

prediction_out %>%
  ggplot(aes(x = 10^number_reads, y = fit)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 1/2) +
  geom_line() +
  stat_slab(data = filter(full_data, sequencing_type == 'NovaSeq'),
            aes(y = has_id, side = ifelse(has_id == 0, "top", "bottom"),
                colour = binary, fill = binary),
            slab_type = "histogram",
            scale = 0.2, breaks = 50, show.legend = FALSE) +
  
  # stat_dots(data = filter(full_data, sequencing_type == 'NovaSeq'),
  #           aes(y = has_id, 
  #               side = ifelse(has_id == 0, "top", "bottom"),
  #               color = binary),
  #           scale = 0.2, shape = 19) +
  
  geom_jitter(data = filter(full_data, sequencing_type == 'MiSeq'), 
             aes(y = if_else(has_id == 1, 0.98, 0.02)),
             shape = 22, size = 2, 
             colour = 'black', fill = 'white',
             height = 0.02, width = 0.05) +
  
  geom_pointrange(data = filter(full_data, sequencing_type == 'NovaSeq') %>%
                    mutate(bins = case_when(number_reads < 4 ~ 1,
                                            number_reads >= 4 & number_reads < 4.5 ~ 2,
                                            number_reads >= 4.5 & number_reads < 5 ~ 3,
                                            number_reads >= 5 & number_reads < 5.5 ~ 4,
                                            number_reads >= 5.5 & number_reads < 6 ~ 5,
                                            number_reads >= 6 & number_reads < 6.5 ~ 6,
                                            number_reads >= 6.5 & number_reads < 7 ~ 7,
                                            number_reads >= 7 ~ 8)) %>%
                    group_by(bins) %>%
                    summarise(number_reads = (min(number_reads) + max(number_reads)) / 2,
                              n_id = sum(has_id),
                              n = n()) %>%
                    rowwise %>%
                    mutate(tidy(prop.test(x = n_id, n = n))),
                  
                  aes(y = estimate, ymin = conf.low, ymax = conf.high)) +
  
  
  scale_y_continuous(labels = scales::percent_format(1), expand = c(0, 0)) +
  scale_x_continuous(trans = scales::log10_trans(),
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides = 'b', outside = FALSE) +
  coord_cartesian(clip = "on") +
  labs(x = 'Number of Reads',
       y = 'Percent Individuals with Mitochondrial Species ID') +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black'),
        panel.border = element_rect(colour = 'black', fill = 'transparent'))

#### Difference in number of reads by species ID ###
t_test_res <- full_data %>%
  filter(sequencing_type == 'NovaSeq') %>%
  filter(!is.na(species)) %>%
  t.test(number_reads ~ species, data = .) %$%
  c(statistic, parameter, p.value) %>%
  round(2)

intervals <- full_data %>%
  filter(sequencing_type == 'NovaSeq') %>%
  filter(!is.na(species)) %>%
  group_by(species) %>%
  summarise(tidy(t.test(number_reads)))

t_res <- c(str_c('t', '[(', t_test_res[2], ')]', ' == ', t_test_res[1]),
         str_c('p == ', t_test_res[3]))

full_data %>%
  filter(sequencing_type == 'NovaSeq') %>%
  filter(!is.na(species)) %>%
  ggplot(aes(x = species, y = 10^number_reads)) +
  geom_beeswarm(alpha = 0.5, size = 0.2) +
  geom_pointrange(data = intervals, aes(y = 10^estimate, ymin = 10^conf.low, ymax = 10^conf.high,
                                        colour = species), 
                  show.legend = FALSE, size = 0.7) +
  annotate('text', x = 1.5, y = c(10^max(full_data$number_reads), 0.75 * 10^max(full_data$number_reads)),
           label = t_res, parse = TRUE, hjust = 0.5) +
  scale_y_continuous(trans = scales::log10_trans(),
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides = 'l', outside = FALSE) +
  coord_cartesian(clip = "on") +
  labs(x = NULL,
       y = 'Number of Reads') +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(face = 'italic'),
        panel.border = element_rect(colour = 'black', fill = 'transparent'))


#### The three groups ####
trio_anova <- full_data %>%
  filter(sequencing_type == 'NovaSeq') %>%
  mutate(species = if_else(is.na(species), 'unknown', species),
         species = str_remove(species, 'oryphopterus'),
         species = factor(species)) %>%
  aov(number_reads ~ species, data = .) 

aov_res <- c(str_c("F", "[\"(", anova(trio_anova)$Df[1], ", ", anova(trio_anova)$Df[2], ")\"]", 
        ' == ', round(anova(trio_anova)$`F value`[1], 2)),
  str_c("p < 0.001"))

source('~/R/R Functions/diagPlots.R')
diag.plots(trio_anova)
anova(trio_anova)

emmeans(trio_anova, ~species) %>% contrast('pairwise')
  
emmeans(trio_anova, ~species) %>%
  tidy(conf.int = TRUE) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~10^.)) %>%
  
  ggplot(aes(x = species, y = estimate)) +
  geom_beeswarm(data = filter(full_data, sequencing_type == 'NovaSeq') %>%
                  mutate(species = if_else(is.na(species), 'unknown', species),
                         species = str_remove(species, 'oryphopterus'),
                         species = factor(species)),
                aes(y = 10^number_reads),
                alpha = 0.5, size = 0.2) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high, colour = species), show.legend = FALSE) +

  annotate('text', x = 2, y = c(10^max(full_data$number_reads), 0.75 * 10^max(full_data$number_reads)),
           label = aov_res, parse = TRUE, hjust = 0.5) +
  scale_y_continuous(trans = scales::log10_trans(),
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides = 'l', outside = FALSE) +
  coord_cartesian(clip = "on") +
  labs(x = NULL,
       y = 'Number of Reads') +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(face = 'italic'),
        panel.border = element_rect(colour = 'black', fill = 'transparent'))

