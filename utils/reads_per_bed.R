library(tidyverse)

coverage <- read_delim('cov.10.1.stats', delim = '\t', 
           col_names = c('contig', 'start', 
                         'end', 'coverage'),
           col_type = cols(
             contig = col_character(),
             start = col_integer(),
             end = col_integer(),
             coverage = col_integer()
           ))

splits <- list.files(pattern = 'mapped.*.10.1.bed') %>%
  str_subset("mapped.10.1.bed", negate = TRUE) %>%
  tibble(file = .) %>%
  rowwise(file) %>%
  summarise(read_delim(file, delim = '\t',
                       col_names = c('contig', 'start', 'end'),
                       col_types = cols(
                         contig = col_character(),
                         start = col_integer(),
                         end = col_integer()
                       )),
            .groups = 'drop')

splits %>%
  filter(contig == 'dDocent_Contig_1',
         start == 0,
         end == 462)


anti_join(splits, coverage, by = c('contig', 'start', 'end'))
anti_join(coverage, splits, by = c('contig', 'start', 'end'))

inner_join(splits, coverage, by = c('contig', 'start', 'end')) %>%
  group_by(file) %>%
  summarise(coverage = sum(coverage)) %>%
  arrange(coverage)


%>%
  filter(is.na(file))
