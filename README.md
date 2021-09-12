# Coryphopterus_RADseq
 Process Coryphopterus NovaSeq samples with dDocent

## Preprocessing Tracking
- MiSeq - Finished
- NovaSeq - Finished

## dDocent Tracking
- MiSeq only - Finished
- MiSeq Assembly, NovaSeq Mapping -
- NovaSeq only -

## ToDo
-

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
| Number Samples | 8 |  |
| Mean Number Reads | 835,009 ± 613,265 SD |  ±  SD |
| Range Number Reads | 6018 - 1,636,205 |  -  |

## Step 9. Get dDocent
I copied [dDocentHPC](https://github.com/cbirdlab/dDocentHPC) to `/work/hobi/jselwyn/Coryphopterus_RAD/scripts`, and added it to `.gitignore`.

## Step 10. Assemble *de novo* reference genomes
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
6. Cutoff1 = 15, Cutoff2 = 1

Must edit config file each time and then run below line
```
sbatch scripts/mkRef.sbatch mkREF_MiSeq config_files/mkREF_MiSeq.config
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
```
MiSeq Reference Stats
| Metric | 1.1 | 2.2 | 10.1 | 5.1 | 5.2 | 10.2 | 15.1 |
| --- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
| Number Contigs | 88,540 | 13,330 | 33,685 | 47,443 | 8,097 | 5,459 | 26,910 |
| Mean Length | 470 ± 65 SD | 476 ± 56 SD | 497 ± 27 SD | 493 ± 33 SD | 492 ± 29 SD | 496 ± 25 SD | 499 ± 25 SD |
| Range Length | 282 - 586 | 282 - 585 | 285 - 575 | 282 - 581 | 284 - 575 | 285 - 575 | 285 - 562 |
| Total Length | 41,644,725 | 6,339,890 | 16,733,979 | 23,379,362 | 3,986,179 | 2,706,419 | 13,416,533 |
| Contigs with Central Ns | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### NovaSeq
```
mkdir mkREF_NovaSeq
mv NovaSeq/fq_fp1_fp2_fqscrn_repaired/*gz mkREF_NovaSeq


```
Move poorly sequenced samples back to the preprocessing area so they aren't used moving forward

Keep all others since they have >10k unique reads. Also exclude the two blanks. Command to move all samples generated by `utils/investigate_preprocessing.R`

| Metric | Remaining Samples |
| --- | ----- |
| Number Samples |  |
| Mean Number Unique Reads |  ±  |
| Range Number Unique Reads |  -  |

Run multiple times with different cutoffs
Cutoff1 is the minimum coverage required to keep a contig
Cutoff2 is the minimum number of individuals a contig must be present in to keep

1. Cutoff1 = , Cutoff2 =
2. Cutoff1 = , Cutoff2 =
3. Cutoff1 = , Cutoff2 =
4. Cutoff1 = , Cutoff2 =
5. Cutoff1 = , Cutoff2 =


Must edit config file each time and then run below line
```
sbatch scripts/mkRef.sbatch mkREF_NovaSeq config_files/mkREF_NovaSeq.config #PE
sbatch scripts/mkRef.sbatch mkREF_NovaSeq_OL config_files/mkREF_NovaSeq.config #OL
```

Check Reference Genome Stats
```
module load R/gcc/64/3.5.1
Rscript scripts/checkContigs.R  mkREF_NovaSeq/reference.2.2.fasta
```
NovaSeq Reference Stats
| Metric |  |  |  |  |  |
| --- | ----- | ----- | ----- | ----- | ----- |
| Number Contigs |  |  |  |  |  |
| Mean Length |  ±  SD |  ±  SD |  ±  SD |  ±  SD |  ±  SD |
| Range Length |  -  |  -  |  -  |  -  |  -  |
| Total Length |  |  |  |  |  |
| Contigs with Central Ns |  |  |  |  |  |

## Step 11. Map reads to *de novo* reference genomes
### Choose Reference Genomes to use
1. MiSeq - mkREF/reference.10.1.fasta
  - Balance confidence in existence of locus with not removing too many.
  - Discuss with others to see about choice
  - 33,685 contigs
2. NovaSeq -

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
  config_files/mkBAM_MiSeq.config \
  mkREF_MiSeq/reference.5.1.fasta

sbatch -o SLURM_out/bam_summary-%j.out \
  scripts/runRscript.sbatch \
  scripts/checkBAM.R \
  mkBAM_MiSeq RAW
```
Mapping Stats
| Metric | # Reads Mapped |
| --- | ----- |
| Mean |  ±  SD |
| Range |  -  |


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

sbatch scripts/fltrBAM.sbatch \
  fltrBAM_test \
  config_files/fltrBAM_MiSeq.config \
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
| Mean Reads Mapped |  ±  SD |
| Range Reads Mapped |  -  |

### Filter reads mapped to NovaSeq reference

## Step 13. Genotyping
### Test scripts by genotyping MiSeq reads mapped to MiSeq reference Genome
```
mkdir mkVCF_test
mv fltrBAM_test/*RG* mkVCF_test

sbatch scripts/mkVCF.sbatch \
  mkVCF_test \
  config_files/mkVCF_MiSeq.config \
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


### Genotype reads mapped to NovaSeq reference
