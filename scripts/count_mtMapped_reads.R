library(tidyverse)

all_files <- list.files('Mitochondrial_Mapping/fastqfiles', pattern = 'fasta$',
           full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(ID = str_extract(file, 'COPE_[0-9]+')) %>%
  rowwise() %>%
  mutate(number_reads = read_lines(file) %>% length)


all_files %>%
  ungroup %>%
  filter(number_reads > 0) %>%
  summarise(n = n(),
            total = sum(number_reads),
            mean = mean(number_reads),
            sd = sd(number_reads))
