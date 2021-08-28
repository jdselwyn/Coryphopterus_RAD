args <- commandArgs(trailingOnly=TRUE)
ref_genome <- args[1]

# ref_genome <- 'mkREF/reference.2.2.fasta'

suppressMessages(library(tidyverse))

message('Contig summary stats:')
read_lines(ref_genome) %>%
  str_subset('^>', negate = TRUE) %>%
  tibble(contig = .) %>%
  mutate(length = str_length(contig),
         has_central_n = str_detect(contig, 
                                    'N[ATGC]+N+[ATGC]+N')) %>%
  summarise(n = n(),
            mean_length = mean(length),
            sd_length = sd(length),
            min_length = min(length),
            max_length = max(length),
            total_length = sum(length),
            number_with_central_n = sum(has_central_n)) %>%
  print(width = Inf)
