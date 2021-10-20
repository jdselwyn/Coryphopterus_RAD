if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  
  DIRECTORY <- args[1]
  VCF <- args[2]
  DECODER <- args[3]
  HYBRID_TYPE <- args[4]
} else {
  setwd("~/Coryphopterus/Bioinformatics/Coryphopterus_RAD")
  DIRECTORY <- 'tmp_dir'
  VCF <- 'tmp_dir/MiSeq_lightSpecies2.10.1.Fltr20.8.randSNPperLoc.vcf'
  DECODER <- 'splitSpecies/hybrid_types.csv'
  HYBRID_TYPE <- 'CHYA'
}

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(magrittr)))


right_hybrid_type <- read_csv(DECODER, col_types = cols(.default = col_character())) %>%
  filter(hybrid == HYBRID_TYPE) %>%
  pull(ID) %T>%
  write_lines(str_c('splitSpecies/', HYBRID_TYPE, '.list', sep = ''))



str_c()

c('module load bcftools')

str_c('bgzip -c ', VCF, ' > ', VCF, '.gz')
str_c('tabix -p vcf ', VCF, '.gz')

str_c('bcftools view -S ', 
'splitSpecies/', HYBRID_TYPE, '.list ',
VCF, '.gz > ')
fltrVCF_MiSeq/MiSeq_CHYA_Initial.10.1.Fltr02.2.recode.vcf