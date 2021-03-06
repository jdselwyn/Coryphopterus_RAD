## Takes a VCF file and calculates Relatedness for all pairs using 1 locus per contig with the least amount of missing data present that is shared by both individuals in the pair.
## For each pair also simulate UNREL number of unrelated pairs using the loci shared by that pair
## For each pair bootstrap the loci (after filtering out those with missing data) and calculate relatedness for BS interval
## Write rds file with relatedness, unrel simulation, and bootstrap.

if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  vcf_file <- args[1]
  gds_path <- args[2]
  UNREL <- as.integer(args[3])
  BOOT <- as.integer(args[4])
  SUBGROUPS <- as.integer(args[5])
  
} else {
  
  if(Sys.info()['sysname'] == 'Linux'){
    vcf_file <- 'fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaK.2.1.Fltr041.22.vcf'
    gds_path <- 'tmp_dir/relatedness_results'
  } else {
    vcf_file <- '~/Coryphopterus/Bioinformatics/Coryphopterus_RAD/fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaK.2.1.Fltr041.22.vcf'
    gds_path <- '~/../Desktop/tmp_relatedness'
  }
  
  UNREL <- 10
  BOOT <- 50
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

#### Functions ####
convert_vcf_gds <- function(vcf, out_path){
  
  out_file <- stringr::str_extract(vcf, '/[A-Za-z0-9\\._]+vcf$') %>%
    stringr::str_replace('vcf$', 'gds') %>%
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

relatedness_genotypes <- function(geno_pair, allele_freq, rel_method){
  SNPRelate::snpgdsPairIBD(geno1 = geno_pair[,1],
                geno2 = geno_pair[,2],
                allele.freq = allele_freq,
                coeff.correct = FALSE,
                method = rel_method,
                verbose = FALSE) %>%
    dplyr::mutate(kinship = (1 - k0 - k1) * 0.5 + k1 * 0.25) %>%
    dplyr::mutate(number_loci = sum(rowSums(is.na(geno_pair)) == 0))
}

simulate_unrelated_pair <- function(af, rel_method){
  
  ur_pair <- purrr::map(af, ~rbinom(n = 2, size = 2, prob = .x)) %>% 
    purrr::transpose() %>%
    purrr::simplify_all() %>%
    do.call(cbind, .)
  
  relatedness_genotypes(ur_pair, allele_freq = af, rel_method = rel_method)
}

bootstrap_rel <- function(pair_geno, snps, af, ...){
  boot_snps <- sort(sample(length(snps), size = length(snps), replace = TRUE))
  
  boot_val <- relatedness_genotypes(pair_geno[boot_snps,], allele_freq = af[boot_snps], rel_method = 'EM')
  
  boot_val <- dplyr::bind_cols(boot_val, calc_likelihoods(pair_geno[boot_snps,], af = af[boot_snps], logLik = boot_val$loglik))
  boot_val
}

calc_logLik <- function(geno1, geno2, allele.freq, k0, k1){
  make_pr_table <- function(g1, g2, p){
    #recreation of PrIBDTable C++ function
    
    q <- 1 - p
    
    if(g1 == 0){
      if(g2 == 0){
        t2 = q*q; t1 = t2*q; t0 = t1*q
      }
      
      if(g2 == 1){
        t1 = p*q*q; t0 = 2*t1*q; t2 = 0
      }
      
      if(g2 == 2){
        t0 = p*p*q*q; t1 = t2 = 0
      }
    }
    
    if(g1 == 1){
      if(g2 == 0){
        t1 = p*q*q; t0 = 2*t1*q; t2 = 0
      }
      
      if(g2 == 1){
        t1 = p*q; t0 = 4*t1*t1; t2 = 2*t1
      }
      
      if(g2 == 2){
        t1 = p*p*q; t0 = 2*p*t1; t2 = 0
      }
    }
    
    if(g1 == 2){
      if(g2 == 0){
        t0 = p*p*q*q; t1 = t2 = 0
      }
      
      if(g2 == 1){
        t1 = p*p*q; t0 = 2*p*t1; t2 = 0
      }
      
      if(g2 == 2){
        t2 = p*p; t1 = t2*p; t0 = t1*p
      }
    }
    
    c(t0, t1, t2)
  }
  
  n <- length(geno1)
  
  pr_tab <- sapply(1:n, function(x) make_pr_table(geno1[x], geno2[x], allele.freq[x]))
  k <- c(k0, k1, 1 - k0 - k1)
  
  sum(log(t(k) %*% pr_tab))
}

calc_likelihoods <- function(genos, af, logLik){
  
  #Coefficients from SNPrelate, http://faculty.washington.edu/tathornt/BIOST551/lectures_2012/Lecture7_Coefficients_of_Identity_and_Coancestry.pdf
  #avuncular is same as halfsib
  
  k0 <- c(0, 0.25, 0, 0.5, 0.75, 1, 9/16, 15/16)
  k1 <- c(0, 0.5, 1, 0.5, 0.25, 0, 6/16, 1/16)
  relationships <- c("self", "fullsib", "offspring", "halfsib", "cousin", "unrelated", 'double.first.cousin', 'second.cousin')
  
  likelihoods <- purrr::map2_dbl(k0, k1,
                                 ~calc_logLik(geno1 = genos[, 1], 
                                              geno2 = genos[, 2],
                                              allele.freq = af, 
                                              k0 = .x, k1 = .y)) %>%
    magrittr::set_names(stringr::str_c('logLik', relationships, sep = '_'))
  
  p_values <- pchisq(logLik - likelihoods, 1, lower.tail = FALSE) %>%
    magrittr::set_names(stringr::str_c('logLik', relationships, 'pValue', sep = '_'))
  
  lik_out <- dplyr::bind_rows(c(likelihoods, p_values))
  
  lik_out$most_likely <- relationships[likelihoods == max(likelihoods)]
  lik_out
}

permutation_EMD_pvalue<- function(x, y, N_permutation = 1000){
  #https://divingintogeneticsandgenomics.rbind.io/post/how-to-test-if-two-distributions-are-different/
  
  permutation_EMD<- function(d){
    d$group<- sample(d$group)
    calculate_EMD(d)
  }
  
  calculate_EMD <- function(df){
    
    #My input to theirs
    
    num<- 1:nrow(df)
    exp_data<- df$value
    names(exp_data)<- glue::glue("sample_{num}")
    labels<- df$group
    names(labels)<- names(exp_data)
    
    EMDomics:::calculate_emd_gene(exp_data, labels, names(exp_data))
  }
  
  
  dat <- tibble::tibble(x = x, y = y) %>%
    tidyr::pivot_longer(cols = 1:2, names_to = "group", values_to = "value")
  
  obs <- calculate_EMD(dat)
  
  permutation_EMDs<- replicate(N_permutation - 1, permutation_EMD(dat))
  
  ### p-value
  p <- mean(obs <= c(permutation_EMDs, obs))
  tibble::tibble(emd_stat = obs, emd_p = p)
}

permute_ks <- function(x, y, N_permutation){
  ind_permutation <- function(d){
    d$group<- sample(d$group)
    ks.test(d$value[d$group == 'x'], d$value[d$group == 'y'])$statistic
  }
  
  dat <- tibble::tibble(x = x, y = y) %>%
    tidyr::pivot_longer(cols = 1:2, names_to = "group", values_to = "value")
  
  obs <- ks.test(dat$value[dat$group == 'x'], dat$value[dat$group == 'y'])$statistic
  
  permutation_ks<- replicate(N_permutation - 1, ind_permutation(dat))
  
  ### p-value
  p <- mean(obs <= c(permutation_ks, obs))
  tibble::tibble(ks_stat = obs, ks_p = p)
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
    out <- relatedness_genotypes(t(genotypes[,snp_use$snp]), allele_freq = snp_use$MajorFreq, rel_method = 'EM')
    
    out <- dplyr::bind_cols(out, calc_likelihoods(t(genotypes[,snp_use$snp]), af = snp_use$MajorFreq, logLik = out$loglik))
    
    if(N_unrel > 0){
      unrel_file <- stringr::str_remove(gds_file, '/[A-Za-z0-9_\\-\\.]+$') %>%
        stringr::str_c('/unrel_sim/', ind1, '-', ind2, '-unrel_sim.csv')
      
      unrel_sim <- replicate(N_unrel, 
                             simulate_unrelated_pair(af = snp_use$MajorFreq, 
                                                     rel_method = 'EM'), 
                             simplify = FALSE) %>%
        dplyr::bind_rows() %>%
        tibble::as_tibble()
      # readr::write_csv(unrel_file)
      
      # out$unrel <- unrel_file
      out$unrel <- list(unrel_sim)
      out$unrel_cutoff_999 <- quantile(unrel_sim$kinship, 0.999, na.rm = TRUE)
      out$unrel_cutoff_99 <- quantile(unrel_sim$kinship, 0.99, na.rm = TRUE)
      out$unrel_cutoff_95 <- quantile(unrel_sim$kinship, 0.95, na.rm = TRUE)
    }
    
    if(N_boot > 0){
      boot_file <- stringr::str_remove(gds_file, '/[A-Za-z0-9_\\-\\.]+$') %>%
        stringr::str_c('/boot_rel/', ind1, '-', ind2, '-bootstrap_relatedness.csv')
      
      boot_rel <- replicate(N_boot, 
                            bootstrap_rel(t(genotypes[,snp_use$snp]), snps = snp_use$snp, 
                                          af = snp_use$MajorFreq, rel_method = 'EM'),
                            simplify = FALSE) %>%
        dplyr::bind_rows() %>%
        tibble::as_tibble()
      # readr::write_csv(boot_file)
      
      # out$boot_rel <- boot_file
      out$boot_rel <- list(boot_rel)
      out$lwr_kinship_95 <- quantile(boot_rel$kinship, 0.025, na.rm = TRUE)
      out$upr_kinship_95 <- quantile(boot_rel$kinship, 0.975, na.rm = TRUE)
      
      pct_boot_types <- (table(c(boot_rel$most_likely, "self", "fullsib", 
                                 "offspring", "halfsib", "cousin", "unrelated", 
                                 'double.first.cousin', 'second.cousin')) - 1) / N_boot
      
      pct_boot_types <- set_names(pct_boot_types, stringr::str_c('pctBoot_', names(pct_boot_types))) %>%
        rbind()
      
      out <- cbind(out, pct_boot_types)
    } 
    
    if(N_boot > 0 & N_unrel > 0 & N_boot == N_unrel){
      out <- cbind(out, 
                   permutation_EMD_pvalue(out$unrel[[1]]$kinship, out$boot_rel[[1]]$kinship, 10000),
                   permute_ks(out$unrel[[1]]$kinship, out$boot_rel[[1]]$kinship, 10000))
    }
    
  } else {
    out <- tibble::tibble(k0 = NA_real_, k1 = NA_real_, loglik = NA_real_, niter = NA_integer_, 
                          kinship = NA_real_, number_loci = NA_integer_, 
                          logLik_self = NA_real_, logLik_fullsib = NA_real_, logLik_offspring = NA_real_, 
                          logLike_halfsib = NA_real_, logLik_cousin = NA_real_, log_Lik_unrelated = NA_real_, 
                          logLike_double.first.cousin = NA_real_, logLik_second.cousin = NA_real_,
                          most_likely = NA_character_,
                          unrel = NA_character_, boot_rel = NA_character_,
                          unrel_cutoff_999 = NA_real_, unrel_cutoff_99 = NA_real_, unrel_cutoff_95 = NA_real_,
                          lwr_kinship_95 = NA_real_, upr_kinship_95 = NA_real_,
                          ind1 = ind1, ind2 = ind2,
                          pctBoot_cousin = NA_real_, pctBoot_double.first.cousin = NA_real_,
                          pctBoot_fullsib = NA_real_, pctBoot_halfsib = NA_real_, pctBoot_offspring = NA_real_, 
                          pctBoot_second.cousin = NA_real_, pctBoot_self = NA_real_, pctBoot_unrelated = NA_real_)
  }
  
  SNPRelate::snpgdsClose(gds)
  
  out$ind1 <- ind1
  out$ind2 <- ind2
  
  #cleanup
  rm(list = stringr::str_subset(ls(), 'out', negate = TRUE))
  # rm(list = ls())
  gc(verbose = FALSE)
  
  dplyr::select(out, ind1, ind2, dplyr::everything())
  # list(k0 = out$k0, k1 = out$k1, loglik = out$loglik, niter = out$niter, kinship = out$kinship, number_loci = out$number_loci, unrel = out$unrel, boot_rel = out$boot_rel)
}

create_individual <- function(gds){
  
  tibble::tibble(snp = get_gds(gds, 'snp.id'),
                 contig = get_gds(gds, 'snp.chromosome')) %>%
    dplyr::bind_cols(SNPRelate::snpgdsSNPRateFreq(gds) %>%
                       dplyr::bind_cols()) %>%
    dplyr::rename(MajorFreq = AlleleFreq) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(genotype = rbinom(1, 2, MajorFreq)) %>%
    dplyr::ungroup() %>%
    dplyr::select(snp, contig, MajorFreq, genotype)
  
}

create_offspring <- function(parent_genotypes){
  p0 <- parent_genotypes
  
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

simulate_relationship <- function(gds_file, n_loci, relationship, rel_method){
  
  gds <- SNPRelate::snpgdsOpen(gds_file, readonly = TRUE, allow.duplicate = TRUE, allow.fork = TRUE)
  
  parents <- dplyr::full_join(create_individual(gds),
                              create_individual(gds),
                              by = c('snp', 'contig', 'MajorFreq')) %>%
    dplyr::group_by(contig) %>%
    dplyr::sample_n(1) %>%
    dplyr::ungroup() %>%
    dplyr::sample_n(n_loci) %>%
    dplyr::mutate(genos = cbind(genotype.x, genotype.y)) %>%
    dplyr::select(-genotype.x, -genotype.y)
  
  if(relationship == 'UR'){ #unrelated
    sim_genotypes <- parents$genos
    
  } else if(relationship == 'PO'){ #parent-offspring
    
    parent <- parents$genos[,sample(2, 1)]
    offspring <- create_offspring(parents$genos)
    
    sim_genotypes <- cbind(parent, offspring)
    
  } else if(relationship == 'FS'){ #full sib
    
    sim_genotypes <- cbind(create_offspring(parents$genos), create_offspring(parents$genos))
    
  } else if(relationship == 'HS'){ #half sib
    
    trio <- dplyr::inner_join(parents, create_individual(gds), 
                              by = c('snp', 'contig', 'MajorFreq')) %>%
      dplyr::mutate(genos = cbind(genos, genotype)) %>%
      dplyr::select(-genotype)
    
    set1 <- sample(3, 2)
    set2 <- c(which(!1:3 %in% set1), sample(set1, 1))
    
    sim_genotypes <- cbind(create_offspring(trio$genos[,set1]), create_offspring(trio$genos[,set2]))
    
  } else if(relationship == 'GG'){ #Grandparent-grandchild
    
    #two pairs of grandparents
    grandparents <- dplyr::inner_join(dplyr::rename(parents, gp1 = genos), 
                                      dplyr::full_join(create_individual(gds),
                                                       create_individual(gds),
                                                       by = c('snp', 'contig', 'MajorFreq')), 
                                      by = c('snp', 'contig', 'MajorFreq')) %>%
      dplyr::mutate(gp2 = cbind(genotype.x, genotype.y)) %>%
      dplyr::select(-genotype.x, -genotype.y)
    
    #one offspring each = parents
    parent1 <- create_offspring(grandparents$gp1)
    parent2 <- create_offspring(grandparents$gp2)
    
    #one kid
    kid <- create_offspring(cbind(parent1, parent2))
    
    #relatedness of random grandparent and the kid
    sim_genotypes <- cbind(cbind(grandparents$gp1, grandparents$gp2)[,sample(4, 1)], kid)
    
  } else if(relationship == 'AV'){ #Avuncular
    
    #two pairs of grandparents
    grandparents <- dplyr::inner_join(dplyr::rename(parents, gp1 = genos), 
                                      dplyr::full_join(create_individual(gds),
                                                       create_individual(gds),
                                                       by = c('snp', 'contig', 'MajorFreq')), 
                                      by = c('snp', 'contig', 'MajorFreq')) %>%
      dplyr::mutate(gp2 = cbind(genotype.x, genotype.y)) %>%
      dplyr::select(-genotype.x, -genotype.y)
    
    
    #one offspring one pair, two offspring one pair
    parent1 <- create_offspring(grandparents$gp1)
    parent2 <- create_offspring(grandparents$gp2)
    aunt_uncle <- cbind(create_offspring(grandparents$gp1), 
                        create_offspring(grandparents$gp2))[,sample(2, 1)]
    
    #one kid
    kid <- create_offspring(cbind(parent1, parent2))
    
    #relatedness btw non-parent & kid
    sim_genotypes <- cbind(aunt_uncle, kid)
    
  } else if(relationship == 'C1'){ #First Cousins
    
    #3 pairs of grandparents
    grandparents <- dplyr::inner_join(dplyr::rename(parents, gp1 = genos), 
                                      dplyr::full_join(create_individual(gds),
                                                       create_individual(gds),
                                                       by = c('snp', 'contig', 'MajorFreq')), 
                                      by = c('snp', 'contig', 'MajorFreq')) %>%
      dplyr::mutate(gp2 = cbind(genotype.x, genotype.y)) %>%
      dplyr::select(-genotype.x, -genotype.y) %>%
      dplyr::inner_join(dplyr::full_join(create_individual(gds),
                                         create_individual(gds),
                                         by = c('snp', 'contig', 'MajorFreq')), 
                        by = c('snp', 'contig', 'MajorFreq')) %>%
      dplyr::mutate(gp3 = cbind(genotype.x, genotype.y)) %>%
      dplyr::select(-genotype.x, -genotype.y) 
    
    #one with 2 kids, two with one kid each
    parent1 <- create_offspring(grandparents$gp1)
    
    parent2 <- create_offspring(grandparents$gp2)
    parent3 <- create_offspring(grandparents$gp2)
    
    parent4 <- create_offspring(grandparents$gp3)
    
    #one kid from each non-inbred pair
    kid1 <- create_offspring(cbind(parent1, parent2))
    kid2 <- create_offspring(cbind(parent3, parent4))
    
    #relatedness btw two kids
    sim_genotypes <- cbind(kid1, kid2)
    
  } else if(relationship == 'C2'){ #Double First Cousins
    
    #two pairs of grandparents
    grandparents <- dplyr::inner_join(dplyr::rename(parents, gp1 = genos), 
                                      dplyr::full_join(create_individual(gds),
                                                       create_individual(gds),
                                                       by = c('snp', 'contig', 'MajorFreq')), 
                                      by = c('snp', 'contig', 'MajorFreq')) %>%
      dplyr::mutate(gp2 = cbind(genotype.x, genotype.y)) %>%
      dplyr::select(-genotype.x, -genotype.y)
    
    #two pairs of kids each
    parent1 <- create_offspring(grandparents$gp1)
    parent2 <- create_offspring(grandparents$gp1)
    
    parent3 <- create_offspring(grandparents$gp2)
    parent4 <- create_offspring(grandparents$gp2)
    
    #cross the families and produce two kids
    kid1 <- create_offspring(cbind(parent1, parent3))
    kid2 <- create_offspring(cbind(parent2, parent4))
    
    #relatedness btw two kids
    sim_genotypes <- cbind(kid1, kid2)
    
  } else if(relationship == 'SC'){
    
    great_grand_parents <- parents
    
    
    grandparent0 <- dplyr::inner_join(dplyr::select(great_grand_parents, -genos),
                                      create_individual(gds),
                                      by = c('snp', 'contig', 'MajorFreq')) %>%
      dplyr::pull(genotype)
    grandparent1 <- create_offspring(great_grand_parents$genos)
    grandparent2 <- create_offspring(great_grand_parents$genos)
    grandparent3 <- dplyr::inner_join(dplyr::select(great_grand_parents, -genos),
                                      create_individual(gds),
                                      by = c('snp', 'contig', 'MajorFreq')) %>%
      dplyr::pull(genotype)
    
    
    
    parent0 <- dplyr::inner_join(dplyr::select(great_grand_parents, -genos),
                                 create_individual(gds),
                                 by = c('snp', 'contig', 'MajorFreq')) %>%
      dplyr::pull(genotype)
    parent1 <- create_offspring(cbind(grandparent0, grandparent1))
    parent2 <- create_offspring(cbind(grandparent2, grandparent3))
    parent3 <- dplyr::inner_join(dplyr::select(great_grand_parents, -genos),
                                 create_individual(gds),
                                 by = c('snp', 'contig', 'MajorFreq')) %>%
      dplyr::pull(genotype)
    
    kid1 <- create_offspring(cbind(parent0, parent1))
    kid2 <- create_offspring(cbind(parent2, parent3))
    
    
    #relatedness btw two kids
    sim_genotypes <- cbind(kid1, kid2)
    
    
    } else {
    message('Only PO/FS/HS/UR relationships are implemented')
    break
  }
  
  sim_rel <- relatedness_genotypes(sim_genotypes, parents$MajorFreq, rel_method = rel_method)
  sim_rel <- dplyr::bind_cols(sim_rel, calc_likelihoods(sim_genotypes, af = parents$MajorFreq, logLik = sim_rel$loglik))
  sim_rel
}

#### Initialize Data ####
gds_file <- convert_vcf_gds(vcf_file, gds_path)

gds <- SNPRelate::snpgdsOpen(gds_file)

max_loci <- dplyr::n_distinct(get_gds(gds, 'snp.chromosome'))

all_pairs <- tidyr::expand_grid(sample1 = get_gds(gds, 'sample.id'), 
            sample2 = get_gds(gds, 'sample.id')) %>%
  dplyr::filter(sample1 < sample2) %>%
  # dplyr::sample_n(50) %>%
  identity() %>%
  dplyr::group_by(groupings = dplyr::row_number() %% SUBGROUPS) %>%
  dplyr::group_split()
SNPRelate::snpgdsClose(gds)

if(Sys.info()['sysname'] == 'Windows'){
  all_pairs <- all_pairs[[sample(length(all_pairs), 1)]] %>%
    dplyr::sample_n(500)
}

#### Set up HPC Jobs ####
relatedness_to_node <- function(in_pairs){
  
  registerDoFuture()
  registerDoRNG()
  
  if(Sys.info()['sysname'] != 'Windows' & !interactive()){
    plan('multicore', gc = TRUE)
  } else {
    plan('multisession', gc = TRUE)
  }
  
  
  # handlers(global = TRUE)
  handlers("progress")
  
  run_relatedness <- function(the_pairs){
    p <- progressor(steps = nrow(the_pairs))
    y <- foreach(x = the_pairs$sample1, y = the_pairs$sample2, 
                 .export = c("gds_file", 'UNREL', 'BOOT'), 
                 .packages = c('magrittr'),
                 .combine = 'rbind',
                 .inorder = FALSE) %dorng% {
                   
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
  
  with_progress({
    out <- run_relatedness(in_pairs)
  })
  out
}

if(Sys.info()['sysname'] != 'Windows'){
  reg <- makeRegistry(file.dir = paste0(gds_path, '/batch_files'), 
                      packages = c('magrittr', 'future',
                                   'doFuture', 'doRNG',
                                   'progressr'))
  
  reg$cluster.functions <- makeClusterFunctionsSlurm(template = "slurm_template.tmpl",
                                                     array.jobs = TRUE,
                                                     nodename = "localhost",
                                                     scheduler.latency = 1,
                                                     fs.latency = 65)
}


#### Calculate Relatedness ####
if(Sys.info()['sysname'] != 'Windows'){
  batchMap(fun = relatedness_to_node, in_pairs = all_pairs)
  batchExport(list(gds_file = gds_file, 
                   UNREL = UNREL, 
                   BOOT = BOOT, 
                   SUBGROUPS = SUBGROUPS,
                   get_gds = get_gds, 
                   relatedness_genotypes = relatedness_genotypes,
                   simulate_unrelated_pair = simulate_unrelated_pair, 
                   bootstrap_rel = bootstrap_rel,
                   calc_logLik = calc_logLik,
                   calc_likelihoods = calc_likelihoods,
                   permutation_EMD_pvalue = permutation_EMD_pvalue,
                   permute_ks = permute_ks,
                   calculate_relatedness_pair = calculate_relatedness_pair))
  
  submitJobs(resources = list(max.concurrent.jobs = 20))
  waitForJobs()
  
  relatedness <- purrr::map_dfr(1:SUBGROUPS, loadResult)
  
  clearRegistry()
} else {
  relatedness <- relatedness_to_node(all_pairs)
}


readr::write_rds(relatedness, stringr::str_replace(gds_file, '\\.gds$', '_relatedness.rds'), compress = 'xz')

readr::write_csv(dplyr::select(relatedness, -unrel, -boot_rel), 
                 stringr::str_replace(gds_file, '\\.gds$', '_relatedness.csv'))

unlink(paste0(gds_path, '/batch_files'), recursive = TRUE)



#### Make a few plots before closing ####
suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))

make_plot <- function(unrel, boot_rel, kinship, lwr_kinship_95, upr_kinship_95, ind1, ind2){
  ind1 = str_remove(ind1, '.fp2.repr.2.1')
  ind2 = str_remove(ind2, '.fp2.repr.2.1')
  
  ggplot() +
    geom_histogram(data = unrel,
                   aes(x = kinship,
                       after_stat(count/sum(count))),
                   bins = 50) +
    geom_histogram(data = boot_rel,
                   aes(x = kinship, after_stat(count/sum(count))),
                   bins = 50, fill = 'red') +
    geom_vline(xintercept = c(lwr_kinship_95,
                              upr_kinship_95),
               linetype = 'dashed') +
    geom_vline(xintercept = kinship) +
    labs(title = str_c(ind1, ind2, sep = '-'),
         x = 'Kinship',
         y = 'Percent Simulated Pairs') +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_classic()
}

N_samples <- min(table(relatedness$most_likely))

various_plots <- relatedness %>%
  as_tibble %>%
  group_by(most_likely) %>%
  sample_n(N_samples) %>%
  ungroup %>%
  
  # rowwise %>%
  # mutate(unrel = list(tibble(kinship = rbeta(1000, 1, 10))),
  #        boot_rel = list(tibble(kinship = rbeta(1000, 1, 3)))) %>%
  # ungroup %>%
  
  select(ind1, ind2, kinship, lwr_kinship_95, upr_kinship_95, most_likely, contains('unrel'), contains('boot_rel')) %>%
  select(-contains('logLik'), -contains('cutoff')) %>%
  rowwise %>%
  mutate(plot = list(make_plot(unrel, boot_rel, kinship, lwr_kinship_95, upr_kinship_95, ind1, ind2))) %>%
  group_by(most_likely) %>%
  
  summarise(group_plots = list(wrap_plots(plot) + plot_annotation(title = most_likely))) %>%
  mutate(out_name = str_c(str_remove(gds_file, '\\.gds$'), most_likely, 'plots.png', sep = '_'))

walk2(various_plots$group_plots, various_plots$out_name, ~ggsave(.y, plot = .x, height = 10, width = 10))
