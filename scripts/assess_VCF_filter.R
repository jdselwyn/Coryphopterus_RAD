##TODO find config files and read in the various settings

library(tidyverse)
library(broom)
library(emmeans)
library(effectsize)

tukey.nonadditivity.test <- function(the.aov) {
  the.model <- the.aov$model
  y <- the.model[,1]; f1 <- the.model[,2]; f2 <- the.model[,3]
  lm.1 <- lm(y~f1+f2)
  interact.term <- fitted(lm.1)^2
  lm.2 <- lm(y~f1+f2+interact.term)
  return(anova(lm.2)[3,])
}

vcf_summary_out <- list.files(path = '..', pattern = '.summary.csv$', recursive = TRUE, full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(stage = str_extract(file, 'mk|flt') %>% str_replace_all(c('mk' = 'preFilter', 'flt' = 'postFilter')),
         species = str_extract(file, 'h_[a-z]+'),
         param_set = str_extract(file, str_c(species, '_[A-Za-z]+')) %>% str_remove(str_c(species, '_'))) %>%
  group_by(stage, species, param_set) %>%
  summarise(read_csv(file, col_types = cols(.default = col_double())),
            .groups = 'drop') %>%
  mutate(param_set = factor(param_set),
         param_set = fct_relevel(param_set, 'default', after = 0L))


#### Very Basic Analysis ####
anova_results <- vcf_summary_out %>%
  filter(stage != 'preFilter') %>%
  select(-min_snp_per_contig, -min_missing_snp) %>% #its always 1
  pivot_longer(cols = -c(stage, species, param_set)) %>%
  group_by(stage, name) %>%
  summarise(anova = list(aov(value ~ param_set + species)),
            .groups = 'rowwise') %>%
  mutate(tukey.nonadditivity.test(anova) %>% tidy,
         param_set_emmean = list(emmeans(anova, ~ param_set) %>% tidy),
         species_emmean = list(emmeans(anova, ~ species) %>% tidy),
         interaction_emmean = list(emmeans(anova, ~ species + param_set) %>% tidy),
         effect_size = list(effectsize(anova)),
         anova = list(tidy(anova) %>% full_join(effect_size, by = c('term' = 'Parameter')))) %>%
  select(-effect_size) %>%
  ungroup() 


anova_results %>%
  select(name, param_set_emmean) %>%
  unnest(param_set_emmean) %>%
  ggplot(aes(x = param_set, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error)) +
  geom_pointrange() +
  facet_wrap(~name, scales = 'free_y')

anova_results %>%
  select(name, species_emmean) %>%
  unnest(species_emmean) %>%
  ggplot(aes(x = species, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error)) +
  geom_pointrange() +
  facet_wrap(~name, scales = 'free_y')


anova_results %>%
  # filter(p.value < 0.05) %>%
  select(name, interaction_emmean) %>%
  unnest(interaction_emmean) %>%
  ggplot(aes(x = param_set, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error, colour = species)) +
  geom_pointrange(position = position_dodge(0.5)) +
  facet_wrap(~name, scales = 'free_y')

#### Slightly better mixed model variant ####
library(afex)
library(broom.mixed)

mixed_results <- vcf_summary_out %>%
  filter(stage != 'preFilter') %>%
  select(-min_snp_per_contig, -min_missing_snp) %>% #its always 1
  pivot_longer(cols = -c(stage, species, param_set)) %>%
  nest_by(stage, name) %>%
  summarise(anova = list(mixed(value ~ param_set + (1 | species), data = data, progress = FALSE)),
            .groups = 'rowwise') %>%
  mutate(as_tibble(anova$anova_table),
         param_set_emmean = list(emmeans(anova, ~ param_set) %>% tidy)) %>%
  ungroup %>%
  janitor::clean_names() %>%
  mutate(p.adjusted = p.adjust(pr_f, 'holm'))


mixed_results %>%
  slice(19)

mixed_results %>%
  select(name, param_set_emmean) %>%
  unnest(param_set_emmean) %>%
  ggplot(aes(x = param_set, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error)) +
  geom_pointrange() +
  facet_wrap(~name, scales = 'free_y')



#### BRMS with SD incorporated ####

