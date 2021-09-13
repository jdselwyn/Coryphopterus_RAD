## Download and write fasta all mitochondrial sequences of a coryphopterus fron Genbank

library(tidyverse)
library(rentrez)

all_cory_search <- 'Coryphopterus[Organism] AND mitochondrial'

NCBI_sequence_ids <- entrez_search(db="nucleotide", term = all_cory_search, retmax=9999, use_history = TRUE)
NCBI_sequences <- entrez_fetch(db = 'nucleotide', rettype = 'fasta', web_history = NCBI_sequence_ids$web_history) 
write_lines(NCBI_sequences, 'Reference_Sequence/Coryphopterus_mtDNA.fasta')
