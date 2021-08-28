args <- commandArgs(trailingOnly=TRUE)
vcf_file <- args[1]
# vcf_file <- '../mkVCF_test/TotalRawSNPs.5.1.vcf'


suppressMessages(library(tidyverse))
suppressMessages(library(vcfR))
suppressMessages(library(magrittr))

out_name <- str_replace(vcf_file, '.vcf', '.summary.csv')

vcf <- read.vcfR(vcf_file, verbose = TRUE) %>%
  vcfR2tidy()


for_missing_data <- vcf$gt %>%
  mutate(Indiv = str_remove(Indiv, ".clmp.fp2.repr")) %>%
  select(ChromKey, POS, Indiv, gt_GT) %>%
  pivot_wider(names_from = Indiv, values_from = gt_GT) 


missing_locus <- for_missing_data %>%
  rowwise %>%
  mutate(count_missing_indiv = mean(is.na(c_across(starts_with('COPE'))))) %>%
  ungroup %>%
  # filter(count_missing_indiv == max(count_missing_indiv))
  summarise(mean_missing_snp = mean(count_missing_indiv),
            sd_missing_snp = sd(count_missing_indiv),
            min_missing_snp = min(count_missing_indiv),
            max_missing_snp = max(count_missing_indiv),
            .groups = 'drop') %T>%
  print(width = Inf)
#For each locus they are missing on average from $MEAN% of individuals ranginging from $MIN% to $MAX%


missing_indiv <- for_missing_data %>%
  summarise(across(starts_with('COPE'), ~mean(is.na(.)))) %>%
  pivot_longer(cols = everything()) %>%
  summarise(n_ind = n(),
            mean_missing_ind = mean(value),
            sd_missing_ind = sd(value),
            min_missing_ind = min(value),
            max_missing_ind = max(value),
            .groups = 'drop') %T>%
  print(width = Inf)
#Each Individual is missing on average $MEAN% of the loci ranging from $MIN% to $MAX%

snp_per_contig <- vcf$fix %>%
  group_by(CHROM) %>%
  summarise(n = n()) %>%
  summarise(mean_snp_per_contig = mean(n),
            sd_snp_per_contig = sd(n),
            min_snp_per_contig = min(n),
            max_snp_per_contig = max(n),
            .groups = 'drop') %T>%
  print(width = Inf)

output <- vcf$fix %>%
  summarise(n_snp = n(),
            n_contig = n_distinct(CHROM),
            
            mean_depth = mean(DP),
            sd_depth = sd(DP),
            min_depth = min(DP),
            max_depth = max(DP),
            
            mean_phred = mean(QUAL),
            sd_phred = sd(QUAL),
            min_phred = min(QUAL),
            max_phred = max(QUAL)) %>%
  select(starts_with('n'), everything()) %T>%
  print(width = Inf) %>%
  bind_cols(missing_locus, missing_indiv, snp_per_contig) %T>%
  write_csv(out_name)

