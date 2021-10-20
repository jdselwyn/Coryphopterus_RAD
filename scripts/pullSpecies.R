if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  
  DIRECTORY <- args[1]
  speciesList <- args[2]
  HYBRID_TYPE <- args[3]
} else {
  setwd("~/Coryphopterus/Bioinformatics/Coryphopterus_RAD")
  DIRECTORY <- 'mkBAM_MiSeq'
  speciesList <- 'splitSpecies/hybrid_types.csv'
  HYBRID_TYPE <- 'CHYA'
}

suppressWarnings(suppressMessages(library(tidyverse)))

dir.create(str_c(DIRECTORY, HYBRID_TYPE, sep = '_'),
           recursive = TRUE, showWarnings = FALSE)

right_hybrid_type <- read_csv(speciesList, col_types = cols(.default = col_character())) %>%
  filter(hybrid == HYBRID_TYPE) %>%
  mutate(ID = str_extract(ID, 'COPE_[0-9]{4}'))

out_dir <- list.files(path = DIRECTORY, full.names = TRUE, 
           pattern = ".fq.gz$") %>%
  tibble(file = .) %>%
  mutate(ID = str_extract(file, 'COPE_[0-9]{4}')) %>%
  inner_join(right_hybrid_type,
             by = 'ID') %>%
  select(-hybrid, -ID) %>%
  mutate(new_file = str_replace(file, DIRECTORY, 
                                str_c(DIRECTORY, HYBRID_TYPE, sep = '_')))

with(out_dir, walk2(file, new_file, ~file.copy(from = .x, to = .y)))