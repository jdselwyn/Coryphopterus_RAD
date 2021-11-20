## Takes a VCF file and calculates Relatedness for all pairs using 1 locus per contig with the least amount of missing data present that is shared by both individuals in the pair.
## For each pair also simulate UNREL number of unrelated pairs using the loci shared by that pair
## For each pair bootstrap the loci (after filtering out those with missing data) and calculate relatedness for BS interval
## Write rds file with relatedness, unrel simulation, and bootstrap.

##TODO - many issues with accumulating memory - issue with how parallelized??m


if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  vcf_file <- args[1]
  gds_path <- args[2]
  UNREL <- as.integer(args[3])
  BOOT <- as.integer(args[4])
  SUBGROUPS <- as.integer(args[5])
  
} else {
  
  if(Sys.info()['sysname'] == 'Linux'){
    vcf_file <- '~/Documents/Coryphopterus/Bioinformatics/Coryphopterus_RAD/fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaK.2.1.Fltr041.22.vcf'
    gds_path <- '~/Documents/Coryphopterus/Bioinformatics/Coryphopterus_RAD/tmp_dir/relatedness_results'
  } else {
    vcf_file <- '~/Coryphopterus/Bioinformatics/Coryphopterus_RAD/fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaK.2.1.Fltr041.22.vcf'
    gds_path <- '~/../Desktop/tmp_relatedness'
  }
  
  UNREL <- 100
  BOOT <- 100
  SUBGROUPS <- 10
}

library(magrittr)


# library(tidyverse)
# library(SNPRelate)
# library(GWASTools)
# library(furrr)
# 


library(future)
library(doFuture)
library(doRNG)
library(progressr)
library(batchtools)


dir.create(gds_path, showWarnings = FALSE, recursive = TRUE)
storage_dir <- purrr::map_lgl(paste(gds_path, c('unrel_sim', 'boot_rel'), sep = '/'), 
                       ~dir.create(.x, showWarnings = FALSE, recursive = TRUE))


# filter_steps <- tibble(filter = c(41, 16, 5, 16, 5, 16, 5, 16, 5, 16, 21),
#                        step = c(22, 23, 24, 26, 27, 29, 30, 32, 33, 35, 37),
#                        perc_ind = c(NA_real_, 0.7, NA_real_, 0.6, NA_real_, 0.5, NA_real_, 0.4, NA_real_, 0.3, NA_real_),
#                        perc_snp = c(NA_real_, NA_real_, 0.7, NA_real_, 0.75, NA_real_, 0.8, NA_real_, 0.85, NA_real_, NA_real_)) %>%
#   mutate(perc_ind = 1 - perc_ind)

#### Functions ####
convert_vcf_gds <- function(vcf, out_path){
  
  out_file <- stringr::str_extract(vcf, '/[A-Za-z0-9\\._]+vcf') %>%
    stringr::str_replace('vcf', 'gds') %>%
    stringr::str_c(out_path, ., sep = '')
  
  if(!file.exists(out_file)){
    SNPRelate::snpgdsVCF2GDS(vcf.fn = vcf,
                  out.fn = out_file,
                  method = 'biallelic.only',
                  verbose = TRUE)
  }
  
  out_file
}

get_gds <- function(gds, column){
  gdsfmt::index.gdsn(gds, index = column) %>%
    gdsfmt::read.gdsn()
}

# summarize_gds <- function(gds_file){
#   gds <- snpgdsOpen(gds_file)
#   
#   genotypes <- get_gds(gds, 'genotype')
#   
#   ind_df <- tibble(sample = get_gds(gds, 'sample.id'),
#                    percent_missing = rowMeans(genotypes == 3))
#   
#   snp_df <- tibble(snp = get_gds(gds, 'snp.id'),
#                    contig = get_gds(gds, 'snp.chromosome'),
#                    # percent_missing = colMeans(genotypes == 3),
#                    # p = (2 * colSums(genotypes == 0) + colSums(genotypes == 1)) / (2 * colSums(genotypes != 3)),
#                    # q = (2 * colSums(genotypes == 2) + colSums(genotypes == 1)) / (2 * colSums(genotypes != 3)),
#                    hwe_p = snpgdsHWE(gds, with.id = FALSE)) %>%
#     bind_cols(snpgdsSNPRateFreq(gds) %>%
#                 bind_cols()) %>%
#     rename(MajorFreq = AlleleFreq)
#   
#   snpgdsClose(gds)
#   
#   rownames(genotypes) <- ind_df$sample
#   
#   tibble(n_ind = nrow(ind_df), n_snp = nrow(snp_df), 
#          n_contig = n_distinct(snp_df$contig),
#          individuals = list(ind_df), loci = list(snp_df),
#          genotypes = list(genotypes))
# }
# 
# shared_missing <- function(geno, inds, loci){
#   geno_sub <- geno[rownames(geno) %in% inds,loci]
#   
#   1 - ((geno_sub != 3) %*% t(geno_sub != 3)) / ncol(geno_sub)
# }

relatedness_genotypes <- function(geno_pair, snps, rel_method, allele_freq){
  SNPRelate::snpgdsPairIBD(geno1 = geno_pair[,1],
                geno2 = geno_pair[,2],
                allele.freq = allele_freq,
                coeff.correct = FALSE,
                method = rel_method,
                verbose = FALSE) %>%
    # dplyr::mutate(kinship = (1 - k0 - k1) * 0.5 + k1 * 0.25) %>%
    dplyr::mutate(number_loci = sum(rowSums(is.na(geno_pair)) == 0))
}

sim_ind <- function(genos, parents, snps){
  p0 <- SNPRelate::snpgdsGetGeno(genos, sample.id = parents, snp.id = snps, snpfirstdim = TRUE, verbose = FALSE) 
  
  p1 <- matrix(c(as.integer(!p0[, 1] %in% c(0, 1)), 
                 as.integer(p0[, 1] %in% c(1, 2))),
               ncol = 2)
  p1[which(is.na(p0[,1])), ] <- NA
  
  
  p2 <- matrix(c(as.integer(!p0[, 2] %in% c(0, 1)), 
                 as.integer(p0[, 2] %in% c(1, 2))),
               ncol = 2)
  p2[which(is.na(p0[,2])), ] <- NA
  
  p1_allele <- sample(c(1,2), nrow(p1), replace = TRUE)
  p2_allele <- sample(c(1,2), nrow(p1), replace = TRUE)
  
  cbind(p1[cbind(seq_along(p1_allele), p1_allele)], 
        p2[cbind(seq_along(p2_allele), p2_allele)]) %>%
    rowSums()
  
}

sim_dyad <- function(genos, ids, snps, relationship = 'UR', af){
  # af <- snpgdsSNPRateFreq(genos, snp.id = snps, with.id = TRUE)$AlleleFreq
  
  parents1 <- sample(ids, 2)
  
  if(relationship == 'UR'){
    parents2 <- sample(ids[!ids %in% parents1], 2)
  } else if(relationship == 'PO'){
    1+1 #don't need 2 parents for the PO relationship
    
  } else if(relationship == 'FS'){
    parents2 <- parents1
    
  } else if(relationship == 'HS'){
    parents2 <- c(sample(parents1, 1), sample(ids[!ids %in% parents1], 1))
    
  } else {
    print('Only PO/FS/HS/UR relationships are implemented')
    break
  }
  
  o1 <- sim_ind(genos, parents1, snps)
  
  if(relationship != 'PO'){
    o2 <- sim_ind(genos, parents2, snps)
  } else {
    o2 <- snpgdsGetGeno(genos, sample.id = sample(parents1, 1), snp.id = snps, snpfirstdim = TRUE, verbose = FALSE) %>%
      as.numeric()
  }
  
  cbind(o1, o2) %>%
    relatedness_genotypes(., snps = snps, allele_freq = af, rel_method)
}

sim_ur <- function(genos, ids, snps, af, rel_method){
  
  parents1 <- sample(ids, 2)
  parents2 <- sample(ids[!ids %in% parents1], 2)
  
  o1 <- sim_ind(genos, parents1, snps)
  o2 <- sim_ind(genos, parents2, snps)
  
  cbind(o1, o2) %>%
    relatedness_genotypes(., snps = snps, allele_freq = af, rel_method = rel_method)
}

bootstrap_rel <- function(pair_geno, snps, af, ...){
  boot_snps <- sort(sample(length(snps), size = length(snps), replace = TRUE))
  
  relatedness_genotypes(pair_geno[boot_snps,], snps = snps[boot_snps], allele_freq = af[boot_snps], rel_method = 'EM')
  
}

calculate_relatedness_pair <- function(ind1, ind2, gds_file, rel_method, N_unrel = 0, N_boot = 0){
  gds <- SNPRelate::snpgdsOpen(gds_file, readonly = TRUE, allow.duplicate = TRUE, allow.fork = TRUE)
  
  #Get sample pair and find loci which are shared
  all_samples <- get_gds(gds, 'sample.id')
  sample_choice <- which(all_samples %in% c(ind1, ind2))
  genotypes <- get_gds(gds, 'genotype')[sample_choice,]
  loci_shared <- which(colSums(genotypes == 3) == 0)
  
  #Step 1 - Filter to only loci both have
  snp_use <- tibble::tibble(snp = get_gds(gds, 'snp.id'),
                   contig = get_gds(gds, 'snp.chromosome')) %>%
    dplyr::bind_cols(SNPRelate::snpgdsSNPRateFreq(gds) %>%
                       dplyr::bind_cols()) %>%
    dplyr::rename(MajorFreq = AlleleFreq) %>%
    dplyr::filter(snp %in% loci_shared) %>%
    
    #Step 2 - filter 1 snp per contig. Based on best estimate of population allele frequency
    # - may want to reevaluate as random snp each round of a bootstrap if possible - also need to figure out how to do bootstrapping here
    dplyr::group_by(contig) %>%
    dplyr::filter(MissingRate == min(MissingRate)) %>%
    
    #Step 3 - if multiple snps per contig then pick most informative
    dplyr::filter(MinorFreq == max(MinorFreq)) %>%
    
    #Step 4 - if still multiple snps per contig sample 1 randomly
    dplyr::sample_n(1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(snp)
  
  if(nrow(snp_use) > 0){
    out <- relatedness_genotypes(t(genotypes[,snp_use$snp]), snps = snp_use$snp, allele_freq = snp_use$MajorFreq, rel_method = 'EM')
    
    if(N_unrel > 0){
      unrel_file <- stringr::str_remove(gds_file, '/[A-Za-z0-9_\\-\\.]+$') %>%
        stringr::str_c('/unrel_sim/', ind1, '-', ind2, '-unrel_sim.csv')
      
      unrel_sim <- replicate(N_unrel, 
                             sim_ur(gds, ids = all_samples, snps = snp_use$snp, af = snp_use$MajorFreq, rel_method = 'EM'), 
                             simplify = FALSE) %>%
        dplyr::bind_rows() %>%
        # as_tibble %>%
        readr::write_csv(unrel_file)

      out$unrel <- unrel_file
    }
    
    if(N_boot > 0){
      boot_file <- stringr::str_remove(gds_file, '/[A-Za-z0-9_\\-\\.]+$') %>%
        stringr::str_c('/boot_rel/', ind1, '-', ind2, '-bootstrap_relatedness.csv')
      
      boot_rel <- replicate(N_boot, 
                            bootstrap_rel(t(genotypes[,snp_use$snp]), snps = snp_use$snp, af = snp_use$MajorFreq, rel_method = 'EM'),
                            simplify = FALSE) %>%
        dplyr::bind_rows() %>%
        # as_tibble %>%
        readr::write_csv(boot_file)
 
      out$boot_rel <- boot_file
    }
    
  } else {
    out <- tibble::tibble(k0 = NA_real_, k1 = NA_real_, loglik = NA_real_, niter = NA_integer_, kinship = NA_real_, number_loci = NA_integer_, unrel = NA_character_, boot_rel = NA_character_)
  }
  
  SNPRelate::snpgdsClose(gds)
  
  out$ind1 <- ind1
  out$ind2 <- ind2
  
  #cleanup
  rm(list = stringr::str_subset(ls(), 'out', negate = TRUE))
  # rm(list = ls())
  gc(verbose = FALSE)
  
  out
  # list(k0 = out$k0, k1 = out$k1, loglik = out$loglik, niter = out$niter, kinship = out$kinship, number_loci = out$number_loci, unrel = out$unrel, boot_rel = out$boot_rel)
}


#### Get summary stats about each vcf ####
gds_file <- convert_vcf_gds(vcf_file, gds_path)

gds <- SNPRelate::snpgdsOpen(gds_file)
all_pairs <- tidyr::expand_grid(sample1 = get_gds(gds, 'sample.id'), 
            sample2 = get_gds(gds, 'sample.id')) %>%
  dplyr::filter(sample1 < sample2) %>%
  # dplyr::sample_n(50) %>%
  identity() %>%
  dplyr::group_by(groupings = dplyr::row_number() %% SUBGROUPS) %>%
  dplyr::group_split()
SNPRelate::snpgdsClose(gds)

#### Set up HPC Jobs ####
reg <- makeRegistry(file.dir = paste0(gds_path, '/batch_files'), 
                    packages = c('magrittr', 'future',
                                 'doFuture', 'doRNG',
                                 'progressr'))

reg$cluster.functions <- makeClusterFunctionsSlurm(template = "slurm_template.tmpl",
                                                   array.jobs = TRUE,
                                                   nodename = "localhost",
                                                   scheduler.latency = 1,
                                                   fs.latency = 65)

send_to_node <- function(in_pairs){
  
  registerDoFuture()
  registerDoRNG()
  
  plan('multicore', gc = TRUE)
  
  handlers(global = TRUE)
  handlers("progress")
  
  run_relatedness <- function(the_pairs){
    p <- progressor(steps = nrow(the_pairs))
    y <- foreach(x = the_pairs$sample1, y = the_pairs$sample2, 
                 .export = c("gds_file", 'UNREL', 'BOOT'), 
                 .packages = c('magrittr'),
                 .combine = 'rbind',
                 .inorder = FALSE) %dopar% {
                   
                   out <- calculate_relatedness_pair(ind1 = x, 
                                                     ind2 = y, 
                                                     gds_file = gds_file,
                                                     N_unrel = UNREL,
                                                     N_boot = BOOT,
                                                     rel_method = 'EM')
                   p()
                   out
                 }
    return(y)
  }
  
  run_relatedness(in_pairs)
}

batchMap(fun = send_to_node, in_pairs = all_pairs)
batchExport(list(gds_file = gds_file, 
                 UNREL = UNREL, 
                 BOOT = BOOT, 
                 SUBGROUPS = SUBGROUPS,
                 get_gds = get_gds, 
                 relatedness_genotypes = relatedness_genotypes,
                 sim_ind = sim_ind, 
                 sim_ur = sim_ur, 
                 bootstrap_rel = bootstrap_rel,
                 calculate_relatedness_pair = calculate_relatedness_pair))


#### Calculate Relatedness ####
submitJobs(resources = list(max.concurrent.jobs = 20))
waitForJobs()

relatedness <- purrr::map_dfr(1:SUBGROUPS, loadResult)

readr::write_csv(relatedness,
                 stringr::str_replace(gds_file, '\\.gds$', '_relatedness.csv'))

#### Read in sim/boot files and output as single rds
# 
full_relatedness <- relatedness %>%
  dplyr::rowwise() %>%
  dplyr::mutate(dplyr::across(c(unrel, boot_rel), 
                              ~list(readr::read_csv(., 
                                                    col_types = readr::cols(k0 = readr::col_double(),
                                                                            k1 = readr::col_double(),
                                                                            loglik = readr::col_double(),
                                                                            niter = readr::col_double(),
                                                                            number_loci = readr::col_double()))
                              ))) %>%
  dplyr::ungroup() %>%
  readr::write_rds(stringr::str_replace(gds_file, '\\.gds$', '_relatedness.rds'), compress = 'xz')
