suppressMessages(library(tidyverse))
suppressMessages(library(rentrez))
suppressMessages(library(XML))
suppressMessages(library(magrittr))

# Read in BLAST Results 
blast_out <- list.files("Mitochondrial_Mapping/blast_out", 
                        pattern = 'tsv$', full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(ID = str_extract(file, 'COPE_[0-9]+')) %>%
  rowwise(ID) %>%
  summarise(read_tsv(file, 
                     col_names = c('qseqid', 'sseqid', 'pident', 'length', 
                                   'mismatch', 'gapopen', 'qstart', 'qend', 
                                   'sstart', 'send', 'evalue', 'bitscore'), 
                     col_types = cols(.default = col_double(), qseqid = col_character(), 
                                      sseqid = col_character())),
            .groups = 'drop') %>%
  select(-qseqid)


#Associate Sequence IDs with Species IDs 
#Runs slower having in pipe and doing each ID individually rather than together
#keeps everything together though so makes life easier

#https://www.metagenomics.wiki/tools/blast/blastn-output-format-6

blast_results <- blast_out %>%
  nest_by(sseqid) %>%
  mutate(taxomic_id = entrez_link(dbfrom = "nucleotide", id = sseqid, db = 'taxonomy')$links$nuccore_taxonomy) %>%
  nest_by(taxomic_id) %>%
  mutate(entrez_fetch(db = 'taxonomy', id = taxomic_id, rettype = 'xml', parsed = TRUE) %>%
           xmlToDataFrame() %>%
           as_tibble) %>%
  unnest(data) %>%
  unnest(data) %>%
  ungroup %T>%
  write_csv('Mitochondrial_Mapping/raw_blast_results.csv')

# blast_results <- read_csv('Mitochondrial_Mapping/raw_blast_results.csv')

#### Identify Species ####
species_id <- blast_results %>%
  arrange(ID) %>%
  mutate(evalue = log(evalue, base = 10)) %>%
  filter(evalue < -65) %>%
  group_by(ID) %>%
  mutate(certainty = if_else(n_distinct(ScientificName) == 1, 'certain', 'uncertain')) %>%
  ungroup %>%
  
  select(ID, ScientificName, sseqid, certainty, 
         bitscore, evalue, pident, length) %>%
  group_by(ID, ScientificName, certainty, sseqid) %>%
  summarise(across(where(is.numeric), mean), .groups = 'drop') %>%
  
  select(-sseqid) %>%
  distinct %>%
  group_by(ID) %>%
  filter(evalue == min(evalue)) %>%
  arrange(ID) %>%
  ungroup() %>%
  rename(species = ScientificName) %T>%
  
  write_csv('Mitochondrial_Mapping/blast_speciesID.csv')
