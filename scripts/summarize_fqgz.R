print('Start')
args <- commandArgs(trailingOnly = TRUE)
print(args)

DIRECTORY <- args[1]
# DIRECTORY <- 'NovaSeq/fq_fp1_clmp_fp2_fqscrn_repaired'

print('libraries')
suppressWarnings(suppressMessages(library(tidyverse)))
print('A')
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(furrr)))

plan('multisession')

setwd(DIRECTORY)

chunk_reader <- function(x, pos){
  sum(str_detect(x, '^\\+$'))
}

chunk_reader2 <- function(x, pos){
  str_subset(x, '^@') %>%
    str_extract('copies=[0-9]+') %>%
    str_replace_na(replacement = '1') %>%
    parse_number() %>%
    sum
}



summarize_preprocessing <- list.files(full.names = FALSE, pattern = 'fq\\.gz$') %>%
  tibble(file = .) %>%
  mutate(sample = str_extract(file, '^[A-Za-z]+[_-][A-Za-z0-9]+')) %>%
  filter(str_detect(file, '[Rr]1|F')) %>%
  mutate(file = str_c(getwd(), file, sep = '/')) %>%
  mutate(number_reads = future_map_dbl(file, ~sum(read_lines_chunked(.x, DataFrameCallback$new(chunk_reader), chunk_size = 100000)),
                                       .options = furrr_options(seed = TRUE))) %>%
  write_csv('sample_preprocess.csv') %>%
  summarise(n = n(),
            mean = mean(number_reads),
            sd = sd(number_reads),
            min = min(number_reads),
            max = max(number_reads),
            .groups = 'drop') %T>%
  write_csv('summary_preprocess.csv')

print(summarize_preprocessing)