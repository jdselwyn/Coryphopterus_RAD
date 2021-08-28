args <- commandArgs(trailingOnly=TRUE)
target_dir <- args[1]
job_id <- args[2]

suppressMessages(library(tidyverse))

setwd(target_dir)

failed_out_files <- list.files(pattern = str_c('clmp.*', job_id)) %>%
  tibble(out_files = .) %>%
  mutate(successful = map_lgl(out_files, ~read_lines(.x) %>%
                                str_detect("Done!") %>%
                                any)) %>%
  filter(!successful) %>%
  pull(out_files)

if(length(failed_out_files) == 0){
  message('Clumpify Successfully worked on all samples')
} else {
  message('Clumpify failed on ', length(failed_out_files), ' samples. Inspect the following outfiles:\n', str_c(failed_out_files, collapse = '\n'))
}
