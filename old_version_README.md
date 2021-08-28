# Coryphopterus_RADseq - Prior to discovering bad interaction btw clumpify & dDocent
 Process Coryphopterus NovaSeq samples with dDocent

## Preprocessing Tracking
- MiSeq - Finished
- NovaSeq - Finished

## dDocent Tracking
- MiSeq only - VCF Filtering
- MiSeq Assembly, NovaSeq Mapping - Genotyping
- NovaSeq only - Assembly - paused for question for Chris

## ToDo
- Assemble novaseq based reference genomes
- Map novaseq reads to both miseq & novaseq reference genomes
- Map preprocessed reads to mtGenome & blast for species ID

## Step 1.  Demultiplex Sequences
```
cd /work/hobi/jselwyn/Coryphopterus_RAD

module load R/gcc/64/3.5.1

#Demultiplex MiSeq sequences as SBATCH array
Rscript scripts/ddRAD_demultiplex.R MiSeq/fastq 5

#Demultiplex NovaSeq sequences as SBATCH array
Rscript scripts/ddRAD_demultiplex.R NovaSeq/fastq 5
```

## Step 2. 1st fastp
```
#Arguments are inDir, outDir, minimum length
sbatch -p gpu -t 1-00:00:00 scripts/runFASTP_1st_trim.sbatch MiSeq/demultiplexed_seqs MiSeq/fq_fp1 280
sbatch --dependency=afterany:37752 scripts/runFASTP_1st_trim.sbatch NovaSeq/demultiplexed_seqs NovaSeq/fq_fp1 140
```

## Step 3. Clumpify
```
#Arguments are inDir, outDir, simultanious array jobs
bash scripts/runCLUMPIFY_r1r2_array.bash MiSeq/fq_fp1 MiSeq/fq_fp1_clmp 5
bash scripts/runCLUMPIFY_r1r2_array.bash NovaSeq/fq_fp1 NovaSeq/fq_fp1_clmp 10

module load R/gcc/64/3.5.1
Rscript scripts/checkClumpify.R SLURM_out 37765 #number is job ID of clumpify array${MiSeqJOBID}
Rscript scripts/checkClumpify.R SLURM_out 37866
```

## Step 4. Run fastp2
```
#Arguments are inDir, outDir, minimum length
sbatch -p gpu -t 1-00:00:00 scripts/runFASTP_2nd_trim.sbatch MiSeq/fq_fp1_clmp MiSeq/fq_fp1_clmp_fp2 280
sbatch scripts/runFASTP_2nd_trim.sbatch NovaSeq/fq_fp1_clmp NovaSeq/fq_fp1_clmp_fp2 140
```

## Step 5. Run fastq_screen
```
#Arguments are inDir, outDir, simultanious array jobs, node type, time limit
bash scripts/runFQSCRN_array.bash MiSeq/fq_fp1_clmp_fp2 MiSeq/fq_fp1_clmp_fp2_fqscrn 5
bash scripts/runFQSCRN_array.bash NovaSeq/fq_fp1_clmp_fp2 NovaSeq/fq_fp1_clmp_fp2_fqscrn 10
```

## Step 6. repair fastq_screen paired end files
```
sbatch scripts/runREPAIR.sbatch MiSeq/fq_fp1_clmp_fp2_fqscrn MiSeq/fq_fp1_clmp_fp2_fqscrn_repaired
sbatch scripts/runREPAIR.sbatch NovaSeq/fq_fp1_clmp_fp2_fqscrn NovaSeq/fq_fp1_clmp_fp2_fqscrn_repaired
```

## Step 7. Make naming convention work for dDocent
Make sure the files follow the dDocent naming convention
* only 1 underscore, and it delineates group from individiual
* `.r1.fq.gz` `.r2.fq.gz` suffixes requried
```
cd MiSeq/fq_fp1_clmp_fp2_fqscrn_repaired/
rename _ . ./*gz
rename - _ ./*gz
cd ../../

cd NovaSeq/fq_fp1_clmp_fp2_fqscrn_repaired/
rename _ . ./*gz
rename - _ ./*gz
cd ../../
```

## Step 8. Summarize post-preprocessing contigs
```
#Repeat this with each preprocessing folder as desired
sbatch -o SLURM_out/novaseq_preprocess_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/summarize_fqgz.R \
  NovaSeq/fq_fp1_clmp_fp2_fqscrn_repaired

```
| Metric | MiSeq | NovaSeq |
| --- | ----- | ----- |
| Number Samples | 8 | 799 |
| Mean Number Unique Reads | 212,872 ± 132,883 SD | 248,482 ± 214,491 |
| Range Number Unique Reads | 4,694 - 348,118 | 80 - 2,079,565 |

## Step 9. Get dDocent
I copied [dDocentHPC](https://github.com/cbirdlab/dDocentHPC) to `/work/hobi/jselwyn/Coryphopterus_RAD/scripts`, and added it to `.gitignore`.

## Step 10. Assemble *de novo* reference genomes
### MiSeq
```
mkdir mkREF_MiSeq
mv MiSeq/fq_fp1_clmp_fp2_fqscrn_repaired/*gz mkREF_MiSeq
```
Move poorly sequenced samples back to the preprocessing area so they aren't used moving forward
- COPE_1033 - 4,694 unique reads
- COPE_0773 - 63,259 unique reads
Keep all others since they have >100k unique reads

```
mv mkREF_MiSeq/COPE_1033* MiSeq/fq_fp1_clmp_fp2_fqscrn_repaired
mv mkREF_MiSeq/COPE_0773* MiSeq/fq_fp1_clmp_fp2_fqscrn_repaired
```

| Metric | Remaining Samples |
| --- | ----- |
| Number Samples | 6 |
| Mean Number Unique Reads | 272,504 ± 85,495 SD |
| Range Number Unique Reads | 115,481 - 348,118 |

Run multiple times with different cutoffs
Cutoff1 is the minimum coverage required to keep a contig
Cutoff2 is the minimum number of individuals a contig must be present in to keep

1. Cutoff1 = 1, Cutoff2 = 1
2. Cutoff1 = 2, Cutoff2 = 2
3. Cutoff1 = 2, Cutoff2 = 1
4. Cutoff1 = 10, Cutoff2 = 1
5. Cutoff1 = 5, Cutoff2 = 1

Must edit config file each time and then run below line
```
sbatch scripts/mkRef.sbatch mkREF_MiSeq config_files/mkREF_MiSeq.config
```

Check Reference Genome Stats
```
module load R/gcc/64/3.5.1
Rscript scripts/checkContigs.R  mkREF_MiSeq/reference.1.1.fasta
Rscript scripts/checkContigs.R  mkREF_MiSeq/reference.2.2.fasta
Rscript scripts/checkContigs.R  mkREF_MiSeq/reference.2.1.fasta
Rscript scripts/checkContigs.R  mkREF_MiSeq/reference.10.1.fasta
Rscript scripts/checkContigs.R  mkREF_MiSeq/reference.5.1.fasta
```
MiSeq Reference Stats
| Metric | 1.1 | 2.2 | 2.1 | 10.1 | 5.1 |
| --- | ----- | ----- | ----- | ----- | ----- |
| Number Contigs | 41,684 | 4,615 | 41,684 | 1,465 | 9,937 |
| Mean Length | 481 ± 50 SD | 491 ± 31 SD | 481 ± 50 SD | 501 ± 18 SD | 497 ± 21 SD |
| Range Length | 282 - 580 | 284 - 576 | 282 - 580 | 358 - 545 | 295 - 558 |
| Total Length | 20,033,447 | 2,265,735 | 20,033,447 | 734,461 | 4,943,468 |
| Contigs with Central Ns | 0 | 0 | 0 | 0 | 0 |

### NovaSeq
```
mkdir mkREF_NovaSeq
mv NovaSeq/fq_fp1_clmp_fp2_fqscrn_repaired/*gz mkREF_NovaSeq

mkdir mkREF_NovaSeq_OL
mv mkREF_NovaSeq/*gz mkREF_NovaSeq_OL
```
Move poorly sequenced samples back to the preprocessing area so they aren't used moving forward
- COPE-1057 - 80 unique_reads
- COPE-0533 - 101 unique_reads
- COPE-0881 - 148 unique_reads
- COPE-0513 - 293 unique_reads
- COPE-0516 - 404 unique_reads
- COPE-0581 - 434 unique_reads
- COPE-0564 - 439 unique_reads
- COPE-0959 - 687 unique_reads
- COPE-1133 - 762 unique_reads
- COPE-0574 - 858 unique_reads
- COPE-1083 - 868 unique_reads
- COPE-0609 - 962 unique_reads
- COPE-0691 - 2,210 unique_reads
- COPE-0590 - 2,320 unique_reads
- COPE-1022 - 2,655 unique_reads
- COPE-1192 - 2,696 unique_reads
- COPE-1262 - 2,703 unique_reads
- COPE-0844 - 2,862 unique_reads
- COPE-0714 - 3,211 unique_reads
- COPE-0689 - 3,269 unique_reads
- COPE-0733 - 3,633 unique_reads
- COPE-0502 - 3,705 unique_reads
- COPE-0777 - 5,437 unique_reads
- COPE-0598 - 5,992 unique_reads
- COPE-0641 - 6,384 unique_reads
- COPE-1066 - 7,950 unique_reads
- COPE-0779 - 8,650 unique_reads
- COPE-1247 - 9,298 unique_reads
- COPE-0827 - 9,446 unique_reads
- COPE-0562 - 9,562 unique_reads
- Blank-1 - 46,244 unique_reads
- Blank-2 - 66,109 unique_reads
Keep all others since they have >10k unique reads. Also exclude the two blanks. Command to move all samples generated by `utils/investigate_preprocessing.R`

| Metric | Remaining Samples |
| --- | ----- |
| Number Samples | 767 |
| Mean Number Unique Reads | 258,575 ± 213,012 |
| Range Number Unique Reads | 10,934 - 2,079,565 |

Run multiple times with different cutoffs
Cutoff1 is the minimum coverage required to keep a contig
Cutoff2 is the minimum number of individuals a contig must be present in to keep

1. Cutoff1 = 2, Cutoff2 = 2, PE
2. Cutoff1 = 1, Cutoff2 = 1, PE
3. Cutoff1 = 2, Cutoff2 = 1, PE
4. Cutoff1 = 5, Cutoff2 = 1, PE
5. Cutoff1 = 2, Cutoff2 = 2, OL
6. Cutoff1 = 1, Cutoff2 = 1, OL
7. Cutoff1 = 5, Cutoff2 = 1, OL
7. Cutoff1 = 1, Cutoff2 = 5, OL


Must edit config file each time and then run below line
```
sbatch scripts/mkRef.sbatch mkREF_NovaSeq config_files/mkREF_NovaSeq.config #PE
sbatch scripts/mkRef.sbatch mkREF_NovaSeq_OL config_files/mkREF_NovaSeq.config #OL
```

Check Reference Genome Stats
```
module load R/gcc/64/3.5.1
Rscript scripts/checkContigs.R  mkREF_NovaSeq/reference.2.2.fasta
Rscript scripts/checkContigs.R  mkREF_NovaSeq/reference.1.1.fasta
Rscript scripts/checkContigs.R  mkREF_NovaSeq/reference.2.1.fasta
Rscript scripts/checkContigs.R  mkREF_NovaSeq/reference.5.1.fasta
Rscript scripts/checkContigs.R  mkREF_NovaSeq_OL/reference.2.2.fasta
Rscript scripts/checkContigs.R  mkREF_NovaSeq_OL/reference.1.1.fasta
Rscript scripts/checkContigs.R  mkREF_NovaSeq_OL/reference.5.1.fasta
Rscript scripts/checkContigs.R  mkREF_NovaSeq_OL/reference.1.5.fasta
```
MiSeq Reference Stats
| Metric | 2.2 PE | 1.1 PE | 2.1 PE | 5.1 PE | 2.2 OL | 1.1 OL | 5.1 OL | 1.5 OL |
| --- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
| Number Contigs | 3,279 | 21,600 | 21,600 | 274 | 7,103 | 40,124 | 1,150 | 1,163 |
| Mean Length | 207 ± 33 SD | 209 ± 36 SD | 209 ± 36 SD | 218 ± 43 SD | 208 ± 39 SD | 214 ± 39 SD | 206 ± 51 SD | 207 ± 41 SD |
| Range Length | 142 - 347 | 142 - 494 | 142 - 494 | 142 - 445 | 52 - 285 | 52 - 285 | 52 - 282 | 52 - 285 |
| Total Length | 679,783 | 4,507,011 | 4,507,011 | 60,087 | 1,475,171 | 8,583,958 | 236,633 | 241,311 |
| Contigs with Central Ns | 12 | 145 | 145 | 6 | 0 | 0 | 0 | 0 |

## Step 11. Map reads to *de novo* reference genomes
### Choose Reference Genomes to use
1. MiSeq - 5.1
  - 9,937 contigs
  - Less variable length contigs but still many
  - Higher confidence that contigs are real because of higher necessary number of reads
2. NovaSeq

### Test scripts by mapping MiSeq reads to MiSeq reference Genome
```
mkdir mkBAM_test
cp mkREF_MiSeq/*gz mkBAM_test

#Rename from r1/r2 to R1/R2
rename .r1.fq.gz .R1.fq.gz ./mkBAM_test/*gz
rename .r2.fq.gz .R2.fq.gz ./mkBAM_test/*gz

sbatch scripts/mkBAM.sbatch \
  mkBAM_test \
  config_files/mkBAM_MiSeq.config \
  mkREF_MiSeq/reference.5.1.fasta

sbatch -o SLURM_out/bam_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/checkBAM.R \
  mkBAM_test RAW
```

Mapping Stats
| Metric | # Reads Mapped |
| --- | ----- |
| Mean | 863,946 ± 245,138 SD |
| Range | 397,755 - 1,068,099 |

### Map NovaSeq reads to MiSeq reference Genome
```
mkdir mkBAM_MiSeq
cp mkREF_NovaSeq/*gz mkBAM_MiSeq

#Rename from r1/r2 to R1/R2
rename .r1.fq.gz .R1.fq.gz ./mkBAM_MiSeq/*gz
rename .r2.fq.gz .R2.fq.gz ./mkBAM_MiSeq/*gz

sbatch scripts/mkBAM.sbatch \
  mkBAM_MiSeq \
  config_files/mkBAM_MiSeq.config \
  mkREF_MiSeq/reference.5.1.fasta

sbatch -o SLURM_out/bam_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/checkBAM.R \
  mkBAM_MiSeq RAW
```

Mapping Stats
| Metric | # Reads Mapped / individual |
| --- | ----- |
| Mean | 745,308 ± 616,456 SD |
| Range | 31,096 - 6,208,250 |

### Map NovaSeq reads to NovaSeq reference Genome
```
mkdir mkBAM_NovaSeq
cp mkBAM_MiSeq/*gz mkBAM_NovaSeq
```
## Step 12. Filter Mapped Reads
### Test scripts by filtering MiSeq reads mapped to MiSeq reference Genome
```
mkdir fltrBAM_test
mv mkBAM_test/*RAW* fltrBAM_test
mv mkBAM_test/reference.*.fasta fltrBAM_test

sbatch scripts/fltrBAM.sbatch \
  fltrBAM_test \
  config_files/fltrBAM_MiSeq.config

sbatch -o SLURM_out/bam_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/checkBAM.R \
  fltrBAM_test RG
```

Mapping Stats
| Metric | Value |
| --- | ----- |
| Mean Reads Mapped | 215,772 ± 90,136 SD |
| Range Reads Mapped | 59,855 - 301,911 |

### Filter reads mapped to MiSeq reference
```
mkdir fltrBAM_MiSeq
mv mkBAM_MiSeq/*RAW* fltrBAM_MiSeq
mv mkBAM_MiSeq/reference.*.fasta fltrBAM_MiSeq

sbatch scripts/fltrBAM.sbatch \
  fltrBAM_MiSeq \
  config_files/fltrBAM_MiSeq.config

sbatch -o SLURM_out/bam_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/checkBAM.R \
  fltrBAM_MiSeq RG
```

Mapping Stats
| Metric | Value |
| --- | ----- |
| Mean Reads Mapped | 90,130 ± 90,437 SD |
| Range Reads Mapped | 463 - 749,713 |

### Filter reads mapped to NovaSeq reference

## Step 13. Genotyping
### Test scripts by genotyping MiSeq reads mapped to MiSeq reference Genome
```
mkdir mkVCF_test
mv fltrBAM_test/*RG* mkVCF_test

sbatch scripts/mkVCF.sbatch \
  mkVCF_test \
  config_files/mkVCF_MiSeq.config \
  mkREF_MiSeq/reference.5.1.fasta

#Run on Head Node
module load R/gcc/64/3.5.1
Rscript scripts/summarizeVCF.R  mkVCF_test/TotalRawSNPs.5.1.vcf

#Run on Node
sbatch -o SLURM_out/vcf_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/summarizeVCF.R \
  mkVCF_test/TotalRawSNPs.5.1.vcf
```

Genotyping Stats
| Metric | Unfiltered VCF |
| --- | ----- |
| Number Individuals | 6 |
| Number SNPs | 18,859 |
| Number Contigs | 4,069 |
| Mean SNPs/Contig | 4.63 ± 4.56 SD |
| Range SNPs/Contig | 1 - 47 |
| Mean Coverage | 152 ± 133 SD |
| Range Coverage | 29 - 2,416 |
| Mean PHRED | 1,698 ± 2,168 SD |
| Range PHRED | 160 - 65,849 |
| Mean Missing (Ind) | 6.4% ± 3.9% |
| Range Missing (Ind) | 2.4% - 13.4% |
| Mean Missing (Loci) | 6.4% ± 13.2% |
| Range Missing (Loci) | 0% - 83.3% |

### Genotype reads mapped to MiSeq reference
```
mkdir mkVCF_MiSeq
mv fltrBAM_MiSeq/*RG* mkVCF_MiSeq

sbatch scripts/mkVCF.sbatch \
  mkVCF_MiSeq \
  config_files/mkVCF_MiSeq.config \
  mkREF_MiSeq/reference.5.1.fasta

#Run on Head Node
module load R/gcc/64/3.5.1
Rscript scripts/summarizeVCF.R  mkVCF_MiSeq/TotalRawSNPs.5.1.vcf

#Run on Node
sbatch -o SLURM_out/vcf_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/summarizeVCF.R \
  mkVCF_MiSeq/TotalRawSNPs.5.1.vcf
```

Genotyping Stats
| Metric | Unfiltered VCF |
| --- | ----- |
| Number Individuals | 767 |
| Number SNPs | 91,985 |
| Number Contigs | 5,075 |
| Mean SNPs/Contig |18.1  ± 14 SD |
| Range SNPs/Contig | 1 - 120 |
| Mean Coverage | 5,588 ± 5,120 SD |
| Range Coverage | 21 - 151,721 |
| Mean PHRED | 13,503 ± 37,921 SD |
| Range PHRED | 0 - 2,375,280 |
| Mean Missing (Ind) | 30.4% ± 24.4% |
| Range Missing (Ind) | 2.9% - 97.9% |
| Mean Missing (Loci) | 30.4% ± 15.8% |
| Range Missing (Loci) | 0% - 99.9% |


### Genotype reads mapped to NovaSeq reference
