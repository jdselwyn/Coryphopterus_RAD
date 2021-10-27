# Coryphopterus_RADseq
 Process Coryphopterus NovaSeq samples with dDocent

## Preprocessing Tracking
- MiSeq - Finished
- NovaSeq - Finished

## dDocent Tracking - unsplit species
- MiSeq only - Finished
- MiSeq Assembly, NovaSeq Mapping - Finished
- NovaSeq only - Finished

## ToDo
- Delete NCBI data - need to change species IDs post ADMIXTURE & NewHybrids
- Upload data to NCBI
- Rerun SNP filtering

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
skip clumpify step here

## Step 3. Run fastp2
```
#Arguments are inDir, outDir, minimum length
sbatch scripts/runFASTP_2nd_trim.sbatch MiSeq/fq_fp1 MiSeq/fq_fp1_fp2 280
sbatch scripts/runFASTP_2nd_trim.sbatch NovaSeq/fq_fp1 NovaSeq/fq_fp1_fp2 140
```

## Step 5. Run fastq_screen
```
#Arguments are inDir, outDir, simultanious array jobs, node type, time limit
bash scripts/runFQSCRN_array.bash MiSeq/fq_fp1_fp2 MiSeq/fq_fp1_fp2_fqscrn 10
bash scripts/runFQSCRN_array.bash NovaSeq/fq_fp1_fp2 NovaSeq/fq_fp1_fp2_fqscrn 15
```

## Step 6. repair fastq_screen paired end files
```
sbatch scripts/runREPAIR.sbatch MiSeq/fq_fp1_fp2_fqscrn MiSeq/fq_fp1_fp2_fqscrn_repaired
sbatch -t 4-00:00:00 -p normal,cbirdq scripts/runREPAIR.sbatch NovaSeq/fq_fp1_fp2_fqscrn NovaSeq/fq_fp1_fp2_fqscrn_repaired
```

## Step 7. Make naming convention work for dDocent
Make sure the files follow the dDocent naming convention
* only 1 underscore, and it delineates group from individiual
* `.r1.fq.gz` `.r2.fq.gz` suffixes requried
```
cd MiSeq/fq_fp1_fp2_fqscrn_repaired/
rename _ . ./*gz
rename - _ ./*gz
cd ../../

cd NovaSeq/fq_fp1_fp2_fqscrn_repaired/
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
  NovaSeq/fq_fp1_fp2_fqscrn_repaired
```
| Metric | MiSeq | NovaSeq |
| --- | ----- | ----- |
| Number Samples | 8 | 800 |
| Mean Number Reads | 835,009 ± 613,265 SD | 2,193,480 ± 2,653,565 SD |
| Range Number Reads | 6018 - 1,636,205 | 202 - 31,430,472 |

## Step 9. Species Identification
Index reference genome to be used - *Bathygobius cocosensis*
```
sbatch scripts/index_ref.slurm \
  Reference_Sequence/bathygobius_cocosensis_complete_mitochondrion.fasta
```
Make BLAST Database with all the Coryphopterus Mitochondrial DNA on GenBank
```
sbatch scripts/buildBlast.sbatch
```

Map to Mitochondrial Genome & BLAST Results
```
mkdir -p Mitochondrial_Mapping
all_prefix=$(ls mkREF_NovaSeq/*.r1.fq.gz | sed 's/.*\///' | sed "s/\..*//")
IFS=' ' read -ra all_prefix <<< $all_prefix
printf "%s\n" "${all_prefix[@]}" > Mitochondrial_Mapping/tmp_prefix_file

sbatch --array=0-$((${#all_prefix[@]}-1))%20 \
  --output=SLURM_out/mitoBLAST_%A_%a.out \
  scripts/mitoBLAST.sbatch \
  Reference_Sequence/bathygobius_cocosensis_complete_mitochondrion.fasta \
  Reference_Sequence/CoryphopterusBlast \
  mkREF_NovaSeq \
  Mitochondrial_Mapping

module load R/gcc/64/3.5.1
Rscript scripts/summarizeBLAST.R
rm Mitochondrial_Mapping/tmp_prefix_file
```


## Step 10. Get dDocent
I copied [dDocentHPC](https://github.com/cbirdlab/dDocentHPC) to `/work/hobi/jselwyn/Coryphopterus_RAD/scripts`, and added it to `.gitignore`.

## Step 11. Assemble *de novo* reference genomes
### MiSeq
```
mkdir mkREF_MiSeq
mv MiSeq/fq_fp1_fp2_fqscrn_repaired/*gz mkREF_MiSeq
```
Move poorly sequenced samples back to the preprocessing area so they aren't used moving forward
- COPE_1033 - 6,018 reads
Keep all others since they have >100k reads and seem reasonably in the ballpark of each other

```
mv mkREF_MiSeq/COPE_1033* MiSeq/fq_fp1_fp2_fqscrn_repaired

sbatch -o SLURM_out/miseqseq_preprocess_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/summarize_fqgz.R \
  mkREF_MiSeq
```

| Metric | Remaining Samples |
| --- | ----- |
| Number Samples | 7 |
| Mean Number Unique Reads | 953,435 ± 554,866 SD |
| Range Number Unique Reads | 149,273 - 1,636,205 |

Run multiple times with different cutoffs
Cutoff1 is the minimum coverage required to keep a contig
Cutoff2 is the minimum number of individuals a contig must be present in to keep

1. Cutoff1 = 1, Cutoff2 = 1
2. Cutoff1 = 2, Cutoff2 = 2
3. Cutoff1 = 10, Cutoff2 = 1
4. Cutoff1 = 5, Cutoff2 = 1
5. Cutoff1 = 5, Cutoff2 = 2
6. Cutoff1 = 10, Cutoff2 = 2
7. Cutoff1 = 15, Cutoff2 = 1
8. Cutoff1 = 2, Cutoff2 = 1

Must edit config file each time and then run below line
```
sbatch scripts/mkRef.sbatch mkREF_MiSeq config_files/MiSeq.config
```

Check Reference Genome Stats
```
module load R/gcc/64/3.5.1
Rscript scripts/checkContigs.R  mkREF_MiSeq/reference.1.1.fasta
Rscript scripts/checkContigs.R  mkREF_MiSeq/reference.2.2.fasta
Rscript scripts/checkContigs.R  mkREF_MiSeq/reference.10.1.fasta
Rscript scripts/checkContigs.R  mkREF_MiSeq/reference.5.1.fasta
Rscript scripts/checkContigs.R  mkREF_MiSeq/reference.5.2.fasta
Rscript scripts/checkContigs.R  mkREF_MiSeq/reference.10.2.fasta
Rscript scripts/checkContigs.R  mkREF_MiSeq/reference.15.1.fasta
Rscript scripts/checkContigs.R  mkREF_MiSeq/reference.2.1.fasta
```
MiSeq Reference Stats
| Metric | 1.1 | 2.2 | 10.1 | 5.1 | 5.2 | 10.2 | 15.1 | 2.1 |
| --- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
| Number Contigs | 88,540 | 13,330 | 33,685 | 47,443 | 8,097 | 5,459 | 26,910 | 88,540 |
| Mean Length | 470 ± 65 SD | 476 ± 56 SD | 497 ± 27 SD | 493 ± 33 SD | 492 ± 29 SD | 496 ± 25 SD | 499 ± 25 SD | 470 ± 65 SD |
| Range Length | 282 - 586 | 282 - 585 | 285 - 575 | 282 - 581 | 284 - 575 | 285 - 575 | 285 - 562 | 282 - 586 |
| Total Length | 41,644,725 | 6,339,890 | 16,733,979 | 23,379,362 | 3,986,179 | 2,706,419 | 13,416,533 | 41,644,725 |
| Contigs with Central Ns | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### NovaSeq
```
mkdir mkREF_NovaSeq
mv NovaSeq/fq_fp1_fp2_fqscrn_repaired/*gz mkREF_NovaSeq

sbatch -o SLURM_out/novaseq_preprocess_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/summarize_fqgz.R \
  mkREF_NovaSeq
```
Move poorly sequenced samples back to the preprocessing area so they aren't used moving forward. Remove the 20 samples with < 10,000 reads

Keep all others since they have >10k reads. Also exclude the two blanks. Command to move all samples generated by `utils/investigate_preprocessing.R`

| Metric | Remaining Samples |
| --- | ----- |
| Number Samples | 778 |
| Mean Number Unique Reads | 2,254,302 ± 2,665,634 |
| Range Number Unique Reads | 10,799 - 31,430,472 |

Run multiple times with different cutoffs
Cutoff1 is the minimum coverage required to keep a contig
Cutoff2 is the minimum number of individuals a contig must be present in to keep

1. Cutoff1 = 2, Cutoff2 = 2
2. Cutoff1 = 10, Cutoff2 = 10
3. Cutoff1 = 20, Cutoff2 = 20
4. Cutoff1 = 10, Cutoff2 = 20
5. Cutoff1 = 20, Cutoff2 = 10


Must edit config file each time and then run below line
```
sbatch scripts/mkRef.sbatch mkREF_NovaSeq config_files/NovaSeq.config
#took ~24 hours for 2,2
```

Check Reference Genome Stats
```
module load R/gcc/64/3.5.1
Rscript scripts/checkContigs.R  mkREF_NovaSeq/reference.2.2.fasta
Rscript scripts/checkContigs.R  mkREF_NovaSeq/reference.10.10.fasta
Rscript scripts/checkContigs.R  mkREF_NovaSeq/reference.20.20.fasta
Rscript scripts/checkContigs.R  mkREF_NovaSeq/reference.10.20.fasta
Rscript scripts/checkContigs.R  mkREF_NovaSeq/reference.20.10.fasta
```
NovaSeq Reference Stats
| Metric | 2.2 | 10.10 | 20.20 | 10.20 | 20.10 |
| --- | ----- | ----- | ----- | ----- | ----- |
| Number Contigs | 361,609 | 45,531 | 20,763 | 27,839 | 33,407 |
| Mean Length | 330 ± 49 SD | 337 ± 23 SD | 338 ± 12 SD | 337 ± 17 SD | 338 ± 17 SD |
| Range Length | 142 - 914 | 142 - 563 | 148 - 534 | 142 - 547 | 142 - 547 |
| Total Length | 119,462,309 | 15,330,836 | 7,010,970 | 9,393,569 | 11,277,928 |
| Contigs with Central Ns | 316,160 | 44,167 | 20,650 | 27,440 | 32,947 |

## Step 12. Map reads to *de novo* reference genomes
### Choose Reference Genomes to use
1. MiSeq - mkREF_MiSeq/reference.10.1.fasta
  - Balance confidence in existence of locus with not removing too many.
  - Discuss with others to see about choice
  - 33,685 contigs
2. NovaSeq - mkREF_NovaSeq/reference.20.10.fasta
  - Same reasons as above
  - 33,407 contigs

### Test scripts by mapping MiSeq reads to MiSeq reference Genome
```
mkdir mkBAM_test
cp mkREF_MiSeq/*gz mkBAM_test

#Rename from r1/r2 to R1/R2
rename .r1.fq.gz .R1.fq.gz ./mkBAM_test/*gz
rename .r2.fq.gz .R2.fq.gz ./mkBAM_test/*gz

sbatch scripts/mkBAM.sbatch \
  mkBAM_test \
  config_files/MiSeq.config \
  mkREF_MiSeq/reference.10.1.fasta

Rscript scripts/checkBAM.R mkBAM_test RAW
```
Mapping Stats
| Metric | # Reads Mapped / individual |
| --- | ----- |
| Mean | 2,485,075 ± 1,286,427 SD |
| Range | 484,360 - 3,960,753 |

### Map NovaSeq reads to MiSeq reference Genome
```
mkdir mkBAM_MiSeq
cp mkREF_NovaSeq/*gz mkBAM_MiSeq

#Rename from r1/r2 to R1/R2
rename .r1.fq.gz .R1.fq.gz ./mkBAM_MiSeq/*gz
rename .r2.fq.gz .R2.fq.gz ./mkBAM_MiSeq/*gz

sbatch scripts/mkBAM.sbatch \
  mkBAM_MiSeq \
  config_files/MiSeq.config \
  mkREF_MiSeq/reference.10.1.fasta

sbatch -o SLURM_out/bam_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/checkBAM.R \
  mkBAM_MiSeq RAW
```
Mapping Stats
| Metric | # Reads Mapped |
| --- | ----- |
| Mean | 6,891,175 ± 8,686,655 SD |
| Range | 31,687 - 119,213,337 |


### Map NovaSeq reads to NovaSeq reference Genome
```
mkdir mkBAM_NovaSeq
cp mkBAM_MiSeq/*gz mkBAM_NovaSeq

sbatch scripts/mkBAM.sbatch \
  mkBAM_NovaSeq \
  config_files/NovaSeq.config \
  mkREF_NovaSeq/reference.20.10.fasta

sbatch -o SLURM_out/bam_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/checkBAM.R \
  mkBAM_NovaSeq RAW

```
Mapping Stats
| Metric | # Reads Mapped |
| --- | ----- |
| Mean | 5,675,411 ± 6,753,900 SD |
| Range | 25,320 - 80,286,143 |


## Step 13. Filter Mapped Reads
### Test scripts by filtering MiSeq reads mapped to MiSeq reference Genome
```
mkdir fltrBAM_test
mv mkBAM_test/*RAW* fltrBAM_test

sbatch scripts/fltrBAM.sbatch \
  fltrBAM_test \
  config_files/MiSeq.config \
  mkREF_MiSeq/reference.10.1.fasta

Rscript scripts/checkBAM.R fltrBAM_test RG
```

Mapping Stats
| Metric | Value |
| --- | ----- |
| Mean Reads Mapped | 1,479,872 ± 922,700 SD |
| Range Reads Mapped | 194,095 - 2,647,685 |

### Filter reads mapped to MiSeq reference
```
mkdir fltrBAM_MiSeq
mv mkBAM_MiSeq/*RAW* fltrBAM_MiSeq

sbatch scripts/fltrBAM.sbatch \
  fltrBAM_MiSeq \
  config_files/MiSeq.config \
  mkREF_MiSeq/reference.10.1.fasta

sbatch -o SLURM_out/bam_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/checkBAM.R \
  fltrBAM_MiSeq RG
```

Mapping Stats
| Metric | Value |
| --- | ----- |
| Mean Reads Mapped | 1,850,867 ± 2,265,530 SD |
| Range Reads Mapped | 1,990 - 20,326,650 |

### Filter reads mapped to NovaSeq reference
```
mkdir fltrBAM_NovaSeq
mv mkBAM_NovaSeq/*RAW* fltrBAM_NovaSeq

sbatch scripts/fltrBAM.sbatch \
  fltrBAM_NovaSeq \
  config_files/NovaSeq.config \
  mkREF_NovaSeq/reference.20.10.fasta

sbatch -o SLURM_out/bam_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/checkBAM.R \
  fltrBAM_NovaSeq RG
```

Mapping Stats
| Metric | Value |
| --- | ----- |
| Mean Reads Mapped | 1,609,163 ± 1,905,412 SD |
| Range Reads Mapped | 3,844 - 21,480,937 |


## Step 14. Genotyping

### Test scripts by genotyping MiSeq reads mapped to MiSeq reference Genome
```
mkdir mkVCF_test
mv fltrBAM_test/*RG* mkVCF_test

sbatch scripts/mkVCF.sbatch \
  mkVCF_test \
  config_files/MiSeq.config \
  mkREF_MiSeq/reference.10.1.fasta

#Run on Head Node
Rscript scripts/summarizeVCF.R  mkVCF_test/TotalRawSNPs.10.1.vcf
```

Genotyping Stats
| Metric | Unfiltered VCF |
| --- | ----- |
| Number Individuals | 7 |
| Number SNPs | 185,312 |
| Number Contigs | 15,690 |
| Mean SNPs/Contig | 11.8 ± 9.4 SD |
| Range SNPs/Contig | 1 - 67 |
| Mean Coverage | 360 ± 316 SD |
| Range Coverage | 20 - 4,049 |
| Mean PHRED | 2,824 ± 3,960 SD |
| Range PHRED | 0 - 117,327 |
| Mean Missing (Ind) | 18.5% ± 8.7% |
| Range Missing (Ind) | 9.6% - 30.8% |
| Mean Missing (Loci) | 18.5% ± 22.3% |
| Range Missing (Loci) | 0% - 85.7% |

### Genotype reads mapped to MiSeq reference
```
mkdir mkVCF_MiSeq
mv fltrBAM_MiSeq/*RG* mkVCF_MiSeq

sbatch -p cbirdq,long -t 15-00:00:00 scripts/mkVCF.sbatch \
  mkVCF_MiSeq \
  config_files/MiSeq.config \
  mkREF_MiSeq/reference.10.1.fasta
#47809

#Run on Node
sbatch -o SLURM_out/vcf_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/summarizeVCF.R \
  mkVCF_MiSeq/TotalRawSNPs.10.1.vcf
48213
```

Genotyping Stats
| Metric | Unfiltered VCF |
| --- | ----- |
| Number Individuals | 778 |
| Number SNPs | 1,020,147 |
| Number Contigs | 22,217 |
| Mean SNPs/Contig | 46 ± 25 SD |
| Range SNPs/Contig | 1 - 162 |
| Mean Coverage | 30,681 ± 33,414 SD |
| Range Coverage | 20 - 708,357 |
| Mean PHRED | 37,422 ± 154,634 SD |
| Range PHRED | 0 - 16,328,900 |
| Mean Missing (Ind) | 49% ± 20% |
| Range Missing (Ind) | 18% - 99.8% |
| Mean Missing (Loci) | 49% ± 26% |
| Range Missing (Loci) | 1.1% - 99.9% |


### Genotype reads mapped to NovaSeq reference
```
mkdir mkVCF_NovaSeq
mv fltrBAM_NovaSeq/*RG* mkVCF_NovaSeq

sbatch -p long,cbirdq -t 15-00:00:00 scripts/mkVCF.sbatch \
  mkVCF_NovaSeq \
  config_files/NovaSeq.config \
  mkREF_NovaSeq/reference.20.10.fasta
#47810

#Run on Node
sbatch -o SLURM_out/vcf_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/summarizeVCF.R \
  mkVCF_NovaSeq/TotalRawSNPs.20.10.vcf
48214
```
Error - too big??
```
Error: Internal error in `dict_hash_with()`: Dictionary is full.
Backtrace:
    █
 1. ├─`%>%`(...)
 2. ├─tidyr::pivot_wider(., names_from = Indiv, values_from = gt_GT)
 3. └─tidyr:::pivot_wider.data.frame(., names_from = Indiv, values_from = gt_GT)
 4.   └─tidyr::build_wider_spec(...)
 5.     └─vctrs::vec_unique(data[names_from])
 6.       ├─vctrs::vec_slice(x, vec_unique_loc(x))
 7.       └─vctrs::vec_unique_loc(x)
 8.         └─(function () ...
```

Genotyping Stats
| Metric | Unfiltered VCF |
| --- | ----- |
| Number Individuals |  |
| Number SNPs |  |
| Number Contigs |  |
| Mean SNPs/Contig |  ±  SD |
| Range SNPs/Contig |  -  |
| Mean Coverage |  ±  SD |
| Range Coverage |  -  |
| Mean PHRED |  ±  SD |
| Range PHRED |  -  |
| Mean Missing (Ind) | % ± % |
| Range Missing (Ind) | % - % |
| Mean Missing (Loci) | % ± % |
| Range Missing (Loci) | % - % |


## Step 14. Filter Genotypes

### Test scripts by filtering vcf of MiSeq reads mapped to MiSeq reference
```
sbatch scripts/fltrVCF.sbatch \
	fltrVCF_test \
	mkVCF_test/TotalRawSNPs.10.1.vcf \
	config_files/fltrVCF_A.config \
	A

#Run on Head Node
module load R/gcc/64/3.5.1
Rscript scripts/summarizeVCF.R  fltrVCF_test/test_A.10.1.Fltr21.22.MostInformativeSNP.vcf

#Run on Node
sbatch -o SLURM_out/vcf_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/summarizeVCF.R \
  fltrVCF_test/test_A.10.1.Fltr21.22.MostInformativeSNP.vcf
```

Genotyping Stats
| Metric | [Filter Set A](config_files/fltrVCF_A.config) |
| --- | ----- |
| JobID | [`48267`](SLURM_out/fltrVCF-48267.out) |
| Summary Graph | [A](fltrVCF_test/test_A.fltrStats2.plots.pdf) |
| Number Individuals | 5 |
| Number SNPs | 6,451 |
| Number Contigs | 6,451 |
| Mean SNPs/Contig | 1 ± 0 SD |
| Range SNPs/Contig | 1 - 1 |
| Mean Coverage | 415 ± 234 SD |
| Range Coverage | 86 - 1,370 |
| Mean PHRED | 3,429 ± 2,874 SD |
| Range PHRED | 283 - 27,495 |
| Mean Missing (Ind) | 2.8% ± 1.2% |
| Range Missing (Ind) | 1.3% - 4.0% |
| Mean Missing (Loci) | 2.8% ± 7.5% |
| Range Missing (Loci) | 0% - 40% |

### Filter genotypes of reads mapped to MiSeq reference
```
sbatch scripts/fltrVCF.sbatch \
	fltrVCF_MiSeq \
	mkVCF_MiSeq/TotalRawSNPs.10.1.vcf \
	config_files/fltrVCF_A.config \
	A

sbatch scripts/fltrVCF.sbatch \
	fltrVCF_MiSeq \
	mkVCF_MiSeq/TotalRawSNPs.10.1.vcf \
	config_files/fltrVCF_lightSpecies.config \
	lightSpecies

sbatch scripts/fltrVCF.sbatch \
	fltrVCF_MiSeq \
	mkVCF_MiSeq/TotalRawSNPs.10.1.vcf \
	config_files/fltrVCF_lightSpecies2.config \
	lightSpecies2

sbatch scripts/fltrVCF.sbatch \
	fltrVCF_MiSeq \
	mkVCF_MiSeq/TotalRawSNPs.10.1.vcf \
	config_files/fltrVCF_initial.config \
	Initial

#Run on Head Node
module load R/gcc/64/3.5.1
Rscript scripts/summarizeVCF.R  fltrVCF_MiSeq/MiSeq_lightSpecies.10.1.Fltr20.7.randSNPperLoc.vcf

Rscript scripts/summarizeVCF.R  fltrVCF_MiSeq/MiSeq_lightSpecies2.10.1.Fltr20.8.randSNPperLoc.vcf

#Run on Node
sbatch -o SLURM_out/vcf_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/summarizeVCF.R \
  fltrVCF_MiSeq/MiSeq_Initial.10.1.Fltr02.2.recode.vcf
48368
```

Genotyping Stats
| Metric | [Filter Set A](config_files/fltrVCF_A.config) | [Species Filter](config_files/fltrVCF_lightSpecies.config) | [Species Filter 2](config_files/fltrVCF_lightSpecies2.config) |
| --- | ----- | ----- | ----- |
| JobID | [`48279`](SLURM_out/fltrVCF-48279.out) | [`48336`](SLURM_out/fltrVCF-48336.out) | [`49316`](SLURM_out/fltrVCF-49316.out) |
| Summary Graph | [A](fltrVCF_MiSeq/MiSeq_A.fltrStats2.plots.pdf) | [Species](fltrVCF_MiSeq/MiSeq_lightSpecies.fltrStats2.plots.pdf) | [Species Filter 2](fltrVCF_MiSeq/MiSeq_lightSpecies2.fltrStats2.plots.pdf) |
| Number Individuals | 516 | 778 | 778 |  |
| Number SNPs | 3,372 | 1,575 | 1,726 |  |
| Number Contigs | 3,372 | 1,575 | 1,726 |  |
| Mean SNPs/Contig | - | - | - |  |
| Range SNPs/Contig | - | - | - |  |
| Mean Coverage | 68,920 ± 27,365 SD | 76,556 ± 49,534 SD | 75,158 ± 48,132 SD |  |
| Range Coverage | 18,326 - 156,828 | 10,101 - 706,622 | 9,804 - 646,358 |  |
| Mean PHRED | 209,170 ± 376,042 SD | 485,037 ± 22,468 SD | 461,508 ± 52,904 SD |  |
| Range PHRED | 258 - 3,264,620 | 82,979 - 500,000 | 13,985 - 500,000 |  |
| Mean Missing (Ind) | 11% ± 10% | 16% ± 25% | 16.3% ± 24.9% |  |
| Range Missing (Ind) | 4.5% - 42% | 0% - 99.5% | 0% - 99.5% |  |
| Mean Missing (Loci) | 11% ± 7% | 16% ± 4.7% | 16.3% ± 4.8% |  |
| Range Missing (Loci) | 0% - 29% | 1.8% - 25% | 1.7% - 25% |  |

### Filter genotypes of reads mapped to NovaSeq reference
```
sbatch scripts/fltrVCF.sbatch \
	fltrVCF_NovaSeq \
	mkVCF_NovaSeq/TotalRawSNPs.20.10.vcf \
	config_files/fltrVCF_A.config \
	A

sbatch scripts/fltrVCF.sbatch \
	fltrVCF_NovaSeq \
	mkVCF_NovaSeq/TotalRawSNPs.20.10.vcf \
	config_files/fltrVCF_lightSpecies.config \
	lightSpecies

sbatch scripts/fltrVCF.sbatch \
	fltrVCF_NovaSeq \
	mkVCF_NovaSeq/TotalRawSNPs.20.10.vcf \
	config_files/fltrVCF_initial.config \
	Initial
48340

#Run on Head Node
module load R/gcc/64/3.5.1
Rscript scripts/summarizeVCF.R  fltrVCF_NovaSeq/NovaSeq_lightSpecies.20.10.Fltr20.7.randSNPperLoc.vcf

#Run on Node
sbatch -o SLURM_out/vcf_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/summarizeVCF.R \
  fltrVCF_NovaSeq/NovaSeq_Initial.20.10.Fltr02.2.recode.vcf
48370
```

Genotyping Stats
| Metric | [Filter Set A](config_files/fltrVCF_A.config) | [Species Filter](config_files/fltrVCF_lightSpecies.config) | [Initial Filter](config_files/fltrVCF_initial.config) |
| --- | ----- | ----- | ----- |
| JobID | [`48280`](SLURM_out/fltrVCF-48280.out) | [`48337`](SLURM_out/fltrVCF-48337.out) | [`48340`](SLURM_out/fltrVCF-48340.out) |
| Summary Graph | [A](fltrVCF_NovaSeq/NovaSeq_A.fltrStats2.plots.pdf) | [Species](fltrVCF_NovaSeq/NovaSeq_lightSpecies.fltrStats2.plots.pdf) | [Species](fltrVCF_NovaSeq/NovaSeq_Initial.fltrStats2.plots.pdf) |
| Number Individuals | 508 | 778 |  |
| Number SNPs | 2,261 | 2,861 |  |
| Number Contigs | 2,261 | 2,861 |  |
| Mean SNPs/Contig | - | - | - |
| Range SNPs/Contig | - | - | - |
| Mean Coverage | 61,620 ± 23,967 SD | 50,394 ± 28,253 SD |  ±  SD |
| Range Coverage | 18,311 - 152,024 | 9,551 - 175,262 |  -  |
| Mean PHRED | 214,921 ± 405,328 SD | 483,908 ± 21,055 SD |  ±  SD |
| Range PHRED | 300 - 3,222,750 | 219,499 - 499,994 |  -  |
| Mean Missing (Ind) | 12% ± 10% | 17% ± 22% | % ± % |
| Range Missing (Ind) | 0.5% - 40% | 0.3% - 98.6% | % - % |
| Mean Missing (Loci) | 12% ± 7% | 17% ± 4% | % ± % |
| Range Missing (Loci) | 0% - 31% | 3.5% - 25% | % - % |

## Step 15. DAPC Analysis
```
sbatch -o SLURM_out/dapc_miseq-%j.out \
  -p cbirdq \
  -t 15-00:00:00 \
  scripts/runRscript.sbatch \
  scripts/assignSpecies.R \
    Mitochondrial_Mapping/blast_speciesID.csv \
    fltrVCF_MiSeq/MiSeq_lightSpecies2.10.1.Fltr20.8.randSNPperLoc.vcf \
    500 \
    MiSeq_lightSpecies2
```
[Out File](SLURM_out/dapc_miseq-49364.out)

## Step 16. Admixture Analysis
Use admixture to find pure specimens to use in NewHybrids. From here on just use MiSeq assembly

Test k = 1 - 25 with 10-fold cross validation
```
sbatch scripts/runADMIXTURE.slurm \
  splitSpecies/ADMIXTURE \
  fltrVCF_MiSeq/MiSeq_lightSpecies2.10.1.Fltr20.8.randSNPperLoc.vcf \
  25 \
  10
```
[Out File](SLURM_out/admixture-49372.out)

## Step 17. NewHybrids Analysis
Use admixture results to identify pure specimens & run NewHybrids to determine recent hybrid status of all other individuals

```
sbatch -o SLURM_out/newHybrids_miseq-%j.out \
  --job-name=NewHybrids \
  -p cbirdq \
  -t 15-00:00:00 \
  scripts/runRscript.sbatch \
  scripts/runNewHybrids.R \
    splitSpecies/newHybrids_random \
    fltrVCF_MiSeq/MiSeq_lightSpecies2.10.1.Fltr20.8.randSNPperLoc.vcf \
    splitSpecies/ADMIXTURE/MiSeq_lightSpecies2.10.1.2.results.csv \
    100000  \
    1000000 \
    100 \
    100 \
    random

sbatch -o SLURM_out/newHybrids_miseq-%j.out \
  --job-name=NewHybrids \
  -p normal \
  -t 4-00:00:00 \
  scripts/runRscript.sbatch \
  scripts/runNewHybrids.R \
    splitSpecies/newHybrids_best \
    fltrVCF_MiSeq/MiSeq_lightSpecies2.10.1.Fltr20.8.randSNPperLoc.vcf \
    splitSpecies/ADMIXTURE/MiSeq_lightSpecies2.10.1.2.results.csv \
    100000  \
    1000000 \
    100 \
    5 \
    best

```
[Initial File](SLURM_out/newHybrids_miseq-49484.out)

[Random 100 Samples 200 Loci File](SLURM_out/newHybrids_miseq-49797.out)

[Largest Fst 200 Loci File](SLURM_out/newHybrids_miseq-49798.out)

## Step 18. Interpret NewHybrids
`utils/investigate_NewHybrids.R`


## Step 19. Upload to NCBI
Go to: https://submit.ncbi.nlm.nih.gov/subs/sra/ and click 'Aspera command line and FTP upload options' to request a preload folder
This is a helpful page: https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/SRA-Upload_Protocol/
  1. Make Folder with all files to upload
  ```
  mkdir NCBI_upload
  cp MiSeq/fq_fp1_fp2_fqscrn_repaired/COPE*gz NCBI_upload/
  cp mkREF_MiSeq/COPE*gz NCBI_upload/
  rename fp2.repr miseq NCBI_upload/*gz

  cp NovaSeq/fq_fp1_fp2_fqscrn_repaired/COPE*gz NCBI_upload/
  cp mkREF_NovaSeq/COPE*gz NCBI_upload/
  rename fp2.repr novaseq NCBI_upload/*gz
  ```
  2. Rename `COPE` portion to species ID from newHybrids
  ```
  module load R/gcc/64/3.5.1
  Rscript scripts/rename_ncbi.R
  ```
  3. Follow NCBI instructions
  ```
  cd NCBI_upload
  sshpass -p PWORD sftp UID@SFTP

  cd THEIR_FOLDER
  mkdir coryphopterus_upload
  cd coryphopterus_upload
  mput *fq.gz

  ```
  4. Make metadata using `Create NCBI metadata.R`

## Step 20. Remove CPERS and Filter Genotypes
```
module load vcftools
vcftools --vcf mkVCF_MiSeq/TotalRawSNPs.10.1.vcf --keep splitSpecies/CHYA.list --recode-INFO-all --recode --out fltrVCF_MiSeq/MiSeq_CHYA.10.1
```

```
sbatch scripts/fltrVCF.sbatch \
	fltrVCF_MiSeq \
	fltrVCF_MiSeq/MiSeq_CHYA.10.1.recode.vcf \
	config_files/fltrVCF_chya_A.config \
	chyaA

sbatch scripts/fltrVCF.sbatch \
	fltrVCF_MiSeq2 \
	fltrVCF_MiSeq/MiSeq_CHYA.10.1.recode.vcf  \
	config_files/fltrVCF_chya_B.config \
	chyaB

sbatch scripts/fltrVCF.sbatch \
	fltrVCF_MiSeq3 \
	fltrVCF_MiSeq/MiSeq_CHYA.10.1.recode.vcf  \
	config_files/fltrVCF_chya_C.config \
	chyaC

mv fltrVCF_MiSeq2/* fltrVCF_MiSeq/; rmdir fltrVCF_MiSeq2
mv fltrVCF_MiSeq3/* fltrVCF_MiSeq/; rmdir fltrVCF_MiSeq3

ls fltrVCF_MiSeq/MiSeq_chya*MostInformativeSNP.vcf

module load R/gcc/64/3.5.1
Rscript scripts/summarizeVCF.R fltrVCF_MiSeq/MiSeq_chyaA.10.1.Fltr21.26.MostInformativeSNP.vcf

Rscript scripts/summarizeVCF.R fltrVCF_MiSeq/MiSeq_chyaB.10.1.Fltr21.28.MostInformativeSNP.vcf

Rscript scripts/summarizeVCF.R fltrVCF_MiSeq/MiSeq_chyaC.10.1.Fltr21.30.MostInformativeSNP.vcf
```

Genotyping Stats
| Metric | [Initial](config_files/fltrVCF_initial.config) | [chyaA](config_files/fltrVCF_chya_A.config) | [chyaB](config_files/fltrVCF_chya_B.config) | [chyaC](config_files/fltrVCF_chya_C.config) |
| --- | ----- | ----- | ----- | ----- |
| JobID | [`48339`](SLURM_out/fltrVCF-48339.out) | [`50250`](SLURM_out/fltrVCF-50250.out) | [`50251`](SLURM_out/fltrVCF-50251.out) | [`50252`](SLURM_out/fltrVCF-50252.out) |
| Summary Graph | [Initial](fltrVCF_MiSeq/MiSeq_Initial.fltrStats2.plots.pdf) | [chyaA](fltrVCF_MiSeq/MiSeq_chyaA.fltrStats2.plots.pdf) | [chyaB](fltrVCF_MiSeq/MiSeq_chyaB.fltrStats2.plots.pdf) | [chyaC](fltrVCF_MiSeq/MiSeq_chyaC.fltrStats2.plots.pdf) |
| Number Individuals | 625 | 471 | 471 | 489 |
| Number SNPs | 802,220 | 2,366 | 2,315 | 2,476 |
| Number Contigs | 22,196 | 2,366 | 2,315 | 2,476 |
| Mean SNPs/Contig | 36.1 ± 18.9 SD | 1 ± 0 SD | 1 ± 0 SD  | 1 ± 0 SD  |
| Range SNPs/Contig | 1 - 113 | 1 - 1 | 1 - 1 | 1 - 1 |
| Mean Coverage | 30,097 ± 33,587 SD | 63,151 ± 20,366 SD | 62,872 ± 19,912 SD | 60,284 ± 19,318 SD |
| Range Coverage | 20 - 708,357 | 18,748 - 128,679 | 24,604 - 128,679 | 22,997 - 110,017 |
| Mean PHRED | 26,029 ± 123,271 SD | 352,606 ± 373,493 SD | 353,785 ± 374,820 SD | 337,504 ± 355,630 SD |
| Range PHRED | 0 - 16,289,400 | 6,608 - 2,463,370 | 6,608 - 2,463,370 | 5,909 - 2,463,370 |
| Mean Missing (Ind) | 48% ± 20% | 6% ± 8% | 6% ± 8% | 8% ± 11% |
| Range Missing (Ind) | 19% - 98% | 0% - 39% | 0% - 38% | 0% - 50% |
| Mean Missing (Loci) | 48% ± 27% | 6% ± 4% | 6% ± 4% | 8% ± 4% |
| Range Missing (Loci) | 1% - 100% | 0% - 15% | 0% - 15% | 0% - 15% |

Variant C seems like a marginally better "filter" in that it retains more individuals and more loci. However the starting point is low with only ~1/2 the contigs as we had before removing non CHYA. Remap to 2.1 MiSeq Reference

## Step 21. Redo Mapping -> genotyping with only CHYA
### Mapping
```
module load R/gcc/64/3.5.1
Rscript scripts/pullSpecies.R \
  mkBAM_MiSeq \
  splitSpecies/hybrid_types.csv \
  CHYA

sbatch scripts/mkBAM.sbatch \
  mkBAM_MiSeq_CHYA \
  config_files/MiSeq.config \
  mkREF_MiSeq/reference.2.1.fasta

sbatch -o SLURM_out/bam_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/checkBAM.R \
  mkBAM_MiSeq_CHYA RAW
```
Mapping Stats
| Metric | # Reads Mapped |
| --- | ----- |
| Mean | 7,072,244 ± 8,730,970 SD |
| Range | 33,464 - 111,946,159 |

### Filter BAM
```
mkdir fltrBAM_MiSeq_CHYA
mv mkBAM_MiSeq_CHYA/*RAW* fltrBAM_MiSeq_CHYA

sbatch scripts/fltrBAM.sbatch \
  fltrBAM_MiSeq_CHYA \
  config_files/MiSeq.config \
  mkREF_MiSeq/reference.2.1.fasta

sbatch -o SLURM_out/bam_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/checkBAM.R \
  fltrBAM_MiSeq_CHYA RG
```

Mapping Stats
| Metric | Value |
| --- | ----- |
| Mean Reads Mapped | 2,624,607 ± 2,907,942 SD |
| Range Reads Mapped | 10,216 - 24,371,260 |

### Genotyping
```
mkdir mkVCF_MiSeq_CHYA
mv fltrBAM_MiSeq_CHYA/*RG* mkVCF_MiSeq_CHYA

sbatch -p cbirdq,long -t 15-00:00:00 scripts/mkVCF.sbatch \
  mkVCF_MiSeq_CHYA \
  config_files/MiSeq.config \
  mkREF_MiSeq/reference.2.1.fasta
50369

#Run on Node
sbatch -o SLURM_out/vcf_summary-%j.out \
  --job-name=VCFsummary \
  scripts/runRscript.sbatch \
  scripts/summarizeVCF.R \
  mkVCF_MiSeq_CHYA/TotalRawSNPs.2.1.vcf
50580
```


Error: Internal error in `dict_hash_with()`: Dictionary is full.
```
Backtrace:
    █
 1. ├─`%>%`(...)
 2. ├─tidyr::pivot_wider(., names_from = Indiv, values_from = gt_GT)
 3. └─tidyr:::pivot_wider.data.frame(., names_from = Indiv, values_from = gt_GT)
 4.   └─tidyr::build_wider_spec(...)
 5.     └─vctrs::vec_unique(data[names_from])
 6.       ├─vctrs::vec_slice(x, vec_unique_loc(x))
 7.       └─vctrs::vec_unique_loc(x)
 8.         └─(function () ...
Execution halted
Too big
```

### Filter VCF
```
sbatch scripts/fltrVCF.sbatch \
	fltrVCF_MiSeq_CHYA \
	mkVCF_MiSeq_CHYA/TotalRawSNPs.2.1.vcf \
	config_files/fltrVCF_chya_A.config \
	chyaA
50581

sbatch scripts/fltrVCF.sbatch \
	fltrVCF_MiSeq_CHYA2 \
	mkVCF_MiSeq_CHYA/TotalRawSNPs.2.1.vcf \
	config_files/fltrVCF_chya_B.config \
	chyaB
50582

sbatch scripts/fltrVCF.sbatch \
	fltrVCF_MiSeq_CHYA3 \
	mkVCF_MiSeq_CHYA/TotalRawSNPs.2.1.vcf \
	config_files/fltrVCF_chya_C.config \
	chyaC
50583

mv fltrVCF_MiSeq_CHYA2/* fltrVCF_MiSeq_CHYA/; rmdir fltrVCF_MiSeq_CHYA2
mv fltrVCF_MiSeq_CHYA3/* fltrVCF_MiSeq_CHYA/; rmdir fltrVCF_MiSeq_CHYA3

ls fltrVCF_MiSeq_CHYA/MiSeq*MostInformativeSNP.vcf

module load R/gcc/64/3.5.1
Rscript scripts/summarizeVCF.R fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaA.2.1.Fltr21.26.MostInformativeSNP.vcf

Rscript scripts/summarizeVCF.R fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaB.2.1.Fltr21.28.MostInformativeSNP.vcf

Rscript scripts/summarizeVCF.R fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaC.2.1.Fltr21.30.MostInformativeSNP.vcf
```

Genotyping Stats
| Metric | [chyaA](config_files/fltrVCF_chya_A.config) | [chyaB](config_files/fltrVCF_chya_B.config) | [chyaC](config_files/fltrVCF_chya_C.config) | [chyaD](config_files/fltrVCF_chya_D.config) | [chyaE](config_files/fltrVCF_chya_E.config) | [chyaF](config_files/fltrVCF_chya_F.config) |  |
| --- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
| JobID | [`50250`](SLURM_out/fltrVCF-50250.out) | [`50251`](SLURM_out/fltrVCF-50251.out) | [`50252`](SLURM_out/fltrVCF-50252.out) | [`50665`](SLURM_out/fltrVCF-50665.out) | [`50677`](SLURM_out/fltrVCF-50677.out) | [`50683`](SLURM_out/fltrVCF-50683.out) |  |
| Summary Graph | [chyaA](fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaA.fltrStats2.plots.pdf) | [chyaB](fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaB.fltrStats2.plots.pdf) | [chyaC](fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaC.fltrStats2.plots.pdf) | [chyaD](fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaD.fltrStats2.plots.pdf) | [chyaE](fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaE.fltrStats2.plots.pdf) | [chyaF](fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaF.fltrStats2.plots.pdf) |  |
| Number Individuals | 473 | 473 | 490 | 495 | 489 | 458 |  |
| Number SNPs | 2,405 | 2,365 | 2,412 | 502 | 513 | 2,371 |  |
| Number Contigs | 2,405 | 2,365 | 2,412 | 502 | 513 | 2,371 |  |
| Mean Coverage | 54,131 ± 17,251 SD | 53,856 ± 16,843 SD | 52,881 ± 16,676 SD | 49,486 ± 17,499 SD | 48,460 ± 17,991 SD | 45,519 ± 17,242 SD | ±  SD |
| Range Coverage | 14,423 - 93,710 | 22,233 - 93,587 | 18,516 - 90,213 | 16,809 - 86,397 | 15,491 - 88,578 | 15,341 - 91,224 |  -  |
| Mean PHRED | 311,290 ± 315,927 SD | 311,922 ± 315,777 SD | 308,786 ± 314,114 SD | 347,059 ± 275,062 SD | 340,234 ± 267,528 SD | 339,880 ± 269,139 SD | ±  SD |
| Range PHRED | 3,708 - 2,196,180 | 5,882 - 2,196,180 | 7,287 - 2,196,180 | 30,008 - 1,563,150 | 24,755 - 1,563,150 | 20,243 - 1,563,150 |  -  |
| Mean Missing (Ind) | 6.5% ± 8.5% | 6.5% ± 8.5% | 8% ± 11% | 5.6% ± 9% | 5.8% ± 9% | 6.0% ± 6.5% | % ± % |
| Range Missing (Ind) | 0% - 39% | 0% - 40% | 0% - 49% | 0% - 40% | 0% - 41% | 0% - 30% | % - % |
| Mean Missing (Loci) | 6.5% ± 3.6% | 6.5% ± 3.6% | 8% ± 4% | 5.6% ± 3.7% | 5.8% ± 4% | 6.0% ± 3.9% | % ± % |
| Range Missing (Loci) | 0% - 15% | 0% - 15% | 0% - 15% | 0% - 14% | 0% - 15% | 0% - 15% | % - % |

A and B are for all intents and purposes identical. C is better than either. Try going harder earlier on the SNPs?? - but still not keeping as many as I'd like....

Change for D.
- Increase required PHRED to be >200 (from 40 - fltr03)
- Require average depth across all individuals to keep a locus to be >25
- Reduce filter 14 to call an individual with 3 reads (from 5)
- Remove entire contigs if loci fail (using fltr86) since at the end I'll only be using one contig per SNP
  - this happens after a variety of steps
- Keep only loci with MAF > 0.1 (previously was 0.005)
- Remove contigs where any loci out of HWE (p < 0.01)

Notes
- Can't run filter 86 after filter 14 as it removes all contigs

Changes for E.
- D was far too restrictive eliminating too many SNPs
- swap the missing data individuals to SNPs after making NA loci/individuals with less than 3 reads
- Swap tail of individual/loci missing data

Changes for F.
- Remove HWE filter - far too restrictive and based solely on p-value (unadjusted)
  - can incorporate into a later filter with p-values adjusted so as to not whily nilly remove all the loci
- add final filter for missing data in individuals @ 30%

Changes for G.
- Increase stringency of first missing data filter for individuals (85 to 75) to get rid of them right at the beginning so maybe save some on loci on the back
- Increase stingency of second missing data filter for loci (50 to 60) to get rid of more which became NA when individuals were set to NA with fewer than 3 loci

```
sbatch scripts/fltrVCF.sbatch \
	fltrVCF_MiSeq_CHYA \
	mkVCF_MiSeq_CHYA/TotalRawSNPs.2.1.vcf \
	config_files/fltrVCF_chya_G.config \
	chyaG
50683

ls fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaF*MostInformativeSNP.vcf

module load R/gcc/64/3.5.1

Rscript scripts/summarizeVCF.R fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaF.2.1.Fltr21.37.MostInformativeSNP.vcf
```
