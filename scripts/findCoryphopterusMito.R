## Download and write fasta all mitochondrial sequences of a coryphopterus fron Genbank

suppressMessages(library(tidyverse))
suppressMessages(library(rentrez))

all_cory_search <- 'Coryphopterus[Organism] AND mitochondrial'
just_relevant_cory <- '(Coryphopterus personatus[Organism] OR Coryphopterus hyalinus[Organism]) AND mitochondrial'

NCBI_sequence_ids <- entrez_search(db="nucleotide", term = just_relevant_cory, retmax=9999, use_history = TRUE)
NCBI_sequences <- entrez_fetch(db = 'nucleotide', rettype = 'fasta', web_history = NCBI_sequence_ids$web_history) 
write_lines(NCBI_sequences, 'Reference_Sequence/Coryphopterus_mtDNA.fasta')
