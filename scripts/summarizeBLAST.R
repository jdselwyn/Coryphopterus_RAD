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


# Identify Species
summarized_blast_results <- blast_results %>%
  select(-taxomic_id, -qstart, -qend, -sstart, -send) %>%
  select(sseqid, ID, pident, length, mismatch, gapopen, evalue, bitscore, ScientificName) %>%
  group_by(ScientificName, ID, sseqid) %>%
  summarise(across(where(is.numeric), mean),
            n = n(),
            .groups = 'drop_last') %>%
  summarise(n = sum(n),
            across(where(is.numeric), mean),
            .groups = 'drop') %>%
  write_csv('Mitochondrial_Mapping/summarized_blast_results.csv')
