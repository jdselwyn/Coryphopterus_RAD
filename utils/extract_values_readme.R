library(tidyverse)
library(magrittr)

the_readme <- read_lines('README.md')

miseq_ref_stats_line <- str_which(the_readme, 'MiSeq Reference Stats')


the_data <- the_readme[(miseq_ref_stats_line + 1):(miseq_ref_stats_line + 7)][-2]



map(the_data, ~tibble(text = unlist(str_split(.x, pattern = " \\| "))) %>%
    rowid_to_column(var = "line")) %>%
  bind_rows(.id = "page") %>%
  select(page, line, text) %>%
  pivot_wider(names_from = 'line',
              values_from = 'text') %>%
  select(-page) %>%
  as.data.frame() %>%
  mutate(`8` = str_remove(`8`, ' \\|')) %>%
  column_to_rownames(1) %>%
  t %>%
  as_tibble() %>%
  rename_with(~str_remove(., '\\| ')) %>%
  separate(Metric, into = c('cut1', 'cut2')) %>%
  separate('Mean Length', into = c('mean_length', 'sd_length', NA)) %>%
  separate('Range Length', into = c('min_length', 'max_length')) %>%
  rename(contigs = `Number Contigs`,
         total_length = `Total Length`) %>%
  select(-`Contigs with Central Ns`) %>%
  mutate(across(everything(), ~str_remove_all(., ',')),
         across(everything(), parse_integer)) %>%
  
  ggplot(aes(x = cut1, y = contigs, colour = as.character(cut2))) +
  geom_point()
