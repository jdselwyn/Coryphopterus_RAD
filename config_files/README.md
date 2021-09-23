# Config Files
## dDocent Configs
- [MiSeq.config](MiSeq.config)
  - Setting used for making reference through genotyping using the MiSeq sequencing to build reference genome with NovaSeq reads mapped and genotyped
- [NovaSeq.config](NovaSeq.config)
  - Setting used for making reference through genotyping using the NovaSeq sequencing to build reference genome and also mapping/genotyping

## Filter VCF Configs
- [A](fltrVCF_A.config)
  - Default settings
  - Change to keep only most Informative SNP per contig
- [initialTwo](fltrVCF_initialTwo.config)
  - Just the first two filter steps
	- Keep only Biallelic SNPs (infinite sites)
	- Remove Indels
  - Can use before anything else since these will always be first two steps
- [lightSpecies](fltrVCF_lightSpecies.config)
  - Light filtering used to get a SNP dataset with loci in most individuals which can be used for DAPC to split the species
