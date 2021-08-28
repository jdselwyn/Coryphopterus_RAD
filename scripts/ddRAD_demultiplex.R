args <- commandArgs(trailingOnly=T)

DIRECTORY<-args[1]
n_jobs <- as.integer(args[2])
# DIRECTORY <- 'MiSeq/fastq'

suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))

#### Make job indicies ####
indicies <- list.files(DIRECTORY, pattern = "*.demultiplex.txt$") %T>%
  write_lines(str_c(DIRECTORY, 'demulti_array', sep = '/')) %>%
  tibble(demulti_file = .) %>%
  mutate(index = str_extract(demulti_file, "Plate[0-9]*Pool[0-9]*"),
         index = str_remove_all(index, "late|ool")) %$%
  index %T>%
  write_lines(str_c(DIRECTORY, 'index_array', sep = '/'))
#Last job is a fake one that is meant to fail - write_lines adds a blank line at the end

#### Send Sbatch array ####
array_call <- str_c("sbatch --array=0-", length(indicies), '%', n_jobs, 
                    " --output=./SLURM_out/demultiplex.%A_%a.out ",
                    "./scripts/ddRAD_demultiplex.slurm ",
                    str_remove(DIRECTORY, '/fastq'))

print(array_call)
system(array_call)