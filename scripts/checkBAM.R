args <- commandArgs(trailingOnly=TRUE)
folder_to_check <- args[1]
subset_chunk <- args[2]
# folder_to_check <- 'mkbam_testing'

suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(furrr))

if(Sys.info()['nodename'] == 'hpcm1'){
  plan('sequential')
} else {
  plan('multisession')
}


setwd(folder_to_check)

check_mapped_reads <- function(file){
  command <- str_c('module load samtools; samtools view -c', file, sep = ' ')
  
  out <- system(command, intern = TRUE) #, ignore.stdout = TRUE
  as.integer(out)
}

message('Mapped Reads summary stats')
bam_stats <- list.files(pattern = '.bam$') %>%
  str_subset(subset_chunk) %>%
  tibble(file = .) %>%
  mutate(file = str_c(getwd(), file, sep = '/'),
         mapped_reads = future_map_int(file, check_mapped_reads, 
                                         .options = furrr_options(seed = TRUE))) %>%
  mutate(file = str_remove(file, str_c(getwd(), '/'))) %T>%
  write_csv('sample_bammapping.csv') %>%
  
  summarise(mean_number = mean(mapped_reads), 
            sd_number = sd(mapped_reads),
            min_number = min(mapped_reads),
            max_number = max(mapped_reads)) %T>%
  write_csv('summary_bammapping.csv') %T>%
  print

