library(tidyverse)
library(seqinr)
count <- dplyr::count

#### Get Data ####
reference_sequences <- read.fasta('../Reference_Sequence/Coryphopterus_mtDNA.fasta') %>%
  tibble(sseqid = names(.),
         sequence = .) %>%
  rowwise %>%
  mutate(ref_length = length(sequence),
         extra = attr(sequence,"Annot"),
         sequence = str_c(sequence, collapse = '')) %>%
  mutate(species = str_extract(extra, 'hyalinus|personatus')) %>%
  select(sseqid, species, ref_length, sequence)

reference_sequences %>%
  dplyr::count(species)

preprocess_data <- list.files('..', pattern = 'sample_preprocess.csv$', 
                              recursive = TRUE, full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(sequencing_type = str_extract(file, 'MiSeq|NovaSeq'),
         preprocess_step = str_extract(file, '(/[a-z0-9_]+/){1}') %>% str_remove_all('/')) %>%
  filter(preprocess_step != 'fq_fp1_clmp_fp2_fqscrn') %>% #not needed since repair is the same
  rowwise(sequencing_type, preprocess_step) %>%
  summarise(read_csv(file, col_types = cols(.default = col_character(), 
                                            number_reads = col_integer())), 
            .groups = 'drop') %>%
  mutate(preprocess_step = factor(preprocess_step, levels = c("demultiplexed_seqs", "fq_fp1", "fq_fp1_fp2", "fq_fp1_fp2_fqscrn_repaired"), ordered = TRUE)) %>%
  filter(sequencing_type == 'NovaSeq',
         preprocess_step == 'fq_fp1_fp2_fqscrn_repaired') %>%
  select(sample, number_reads)


all_results <- read_csv('../Mitochondrial_Mapping/raw_blast_results.csv') %>%
  arrange(ID) %>%
  group_by(ID) %>%
  mutate(certainty = if_else(n_distinct(ScientificName) == 1, 'certain', 'uncertain')) %>%
  ungroup %>%
  mutate(evalue = log(evalue, base = 10)) %>%
  select(ID, ScientificName, sseqid, certainty, 
                                                bitscore, evalue, pident, length) %>%
  group_by(ID, ScientificName, certainty, sseqid) %>%
  summarise(across(where(is.numeric), mean), .groups = 'drop') 


#### Identify Species of each ID which successfully mapped to mitochondrion & has a blast hit ####
species_id <- all_results %>%
  select(-sseqid) %>%
  distinct %>%
  group_by(ID) %>%
  filter(evalue == min(evalue)) %>%
  arrange(ID) %>%
  ungroup() %>%
  rename(species = ScientificName)

select(species_id, ID, species) %>%
  write_csv('../')

#### See how ID relates to number of reads ####
tmp <- full_join(species_id, preprocess_data, by = c('ID' = 'sample')) %>%
  filter(number_reads >= 10000) %>%
  mutate(has_id = !is.na(species),
         number_reads = log(number_reads, base = 10)) 

glm(has_id ~ number_reads, data = tmp, family = 'binomial') %>%
  summary

prop_test <- function(x){
  a <- 1;b <- 1
  x <- stats::na.omit(x)
  
  tmp_test <- binom.test(x = sum(x), n = length(x))
  
  ggplot2:::new_data_frame(list(y = tmp_test$estimate, ymin = tmp_test$conf.int[1], 
                      ymax = tmp_test$conf.int[2]), n = 1)
}

tst <- as.character(sort(unique(cut_number(tmp$number_reads, n = 50)))) %>%
  str_extract(',[0-9\\.]+') %>% parse_number()

tmp %>%
  ggplot(aes(x = number_reads, y = as.numeric(has_id))) +
  geom_point(colour = 'grey50', size = 0.7) +
  stat_summary_bin(fun.data = prop_test,
                   breaks = c(seq(4, 5, by = 0.25),
                              seq(5.1, 6.5, by = 0.1),
                              seq(6.75, 8, by = 0.5)),
                   colour = 'black') +
  # stat_summary_bin(fun.data = prop_test, 
  #                  breaks = tst,
  #                  colour = 'black') +
  geom_smooth(method = 'glm', formula = y ~ x, 
              method.args = list(family = 'binomial')) +
  theme_classic()

