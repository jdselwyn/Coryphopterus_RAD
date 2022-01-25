
if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  vcf_file <- args[1]
  gds_path <- args[2]
  NSIM <- as.integer(args[3])
  
} else {
  
  if(Sys.info()['sysname'] == 'Linux'){
    vcf_file <- 'fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaK.2.1.Fltr041.22.vcf'
    gds_path <- 'tmp_dir/relatedness_results'
  } else {
    vcf_file <- '~/Coryphopterus/Bioinformatics/Coryphopterus_RAD/fltrVCF_MiSeq_CHYA/MiSeq_CHYA_chyaK.2.1.Fltr041.22.vcf'
    gds_path <- '~/../Desktop/tmp_relatedness'
  }
  
  NSIM <- 500
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
sims_to_node <- function(param_set){
  
  registerDoFuture()
  registerDoRNG()
  
  if(Sys.info()['sysname'] != 'Windows' & !interactive()){
    plan('multicore', gc = TRUE)
  } else {
    plan('multisession', gc = TRUE)
  }
  
  
  # handlers(global = TRUE)
  handlers("progress")
  
  n_loci <- param_set$n_loci
  relationship <- param_set$relationship
  rel_method <- param_set$rel_method
  
  run_simulation <- function(NSIM){
    p <- progressor(steps = NSIM)
    
    
    
    y <- foreach(x = 1:NSIM, 
                 .export = c("gds_file", 'n_loci', 'relationship', 'rel_method'), 
                 .packages = c('magrittr'),
                 .combine = 'rbind',
                 .inorder = FALSE) %dorng% {
                   
                   out <- simulate_relationship(gds_file = gds_file, 
                                                n_loci, 
                                                relationship, 
                                                rel_method)
                   p()
                   out
                 }
    return(y)
  }
  
  with_progress({
    out <- run_simulation(NSIM)
  })
  
  cbind(param_set, out)
  
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


#### Simulate Relatedness ####
simulation_settings <- tidyr::expand_grid(relationship = c('PO', 'FS', 'HS', 'UR',
                                                           'GG', 'AV', 'C1', 'C2',
                                                           'SC'),
                                          rel_method = c('EM', 'MoM'),
                                          n_loci = unique(c(floor(seq(1, 100, length.out = 5)),
                                                            floor(seq(100, 1000, length.out = 5)),
                                                            round(seq(1000, max_loci, 
                                                                      length.out = 5))))) %>%
  dplyr::mutate(n_loci = dplyr::if_else(n_loci == 1, 2, n_loci)) %>%
  dplyr::sample_frac(ifelse(Sys.info()['sysname'] == 'Windows', 0.05, 1)) %>%
  dplyr::group_by(groupings = dplyr::row_number()) %>%
  dplyr::group_split()

if(Sys.info()['sysname'] != 'Windows'){
  batchMap(fun = sims_to_node, param_set = simulation_settings)
  batchExport(list(gds_file = gds_file, 
                   get_gds = get_gds, 
                   relatedness_genotypes = relatedness_genotypes,
                   calc_logLik = calc_logLik,
                   calc_likelihoods = calc_likelihoods,
                   create_individual = create_individual,
                   create_offspring = create_offspring,
                   simulate_relationship = simulate_relationship,
                   NSIM = NSIM))
  
  submitJobs(resources = list(max.concurrent.jobs = 20))
  waitForJobs()
  
  
  successful_jobs <- dplyr::as_tibble(getJobTable()) %>%
    dplyr::filter(is.na(error)) %>%
    dplyr::pull(job.id)
  
  message(length(successful_jobs), ' successful jobs\n',
          nrow(getJobTable()) - length(successful_jobs), ' failed jobs')
  
  if(nrow(getJobTable()) - length(successful_jobs) > 0){
    message('Failure Reasons:')
    dplyr::as_tibble(getErrorMessages()) %>%
      dplyr::select(job.id, message) %>%
      print()
    
    message('Failure Parameter Sets:')
    failed_jobs <- dplyr::as_tibble(getJobTable()) %>%
      dplyr::filter(!is.na(error)) %>%
      dplyr::select(job.id, job.pars) %>%
      tidyr::unnest(job.pars) %>%
      tidyr::unnest(job.pars) 
    
    print(failed_jobs)
  }
  
  simulated_relatedness <- purrr::map_dfr(successful_jobs, loadResult)
  
  clearRegistry()
  
  z <- 0
  while(length(simulation_settings) != length(successful_jobs)){
    z <- z + 1
    message('Rerunning failed jobs round ', z)
    
    simulation_settings2 <- dplyr::bind_rows(simulation_settings) %>%
      dplyr::filter(groupings %in% failed_jobs$groupings) %>%
      dplyr::group_by(groupings) %>%
      dplyr::group_split()
    
    batchMap(fun = sims_to_node, param_set = simulation_settings2)
    batchExport(list(gds_file = gds_file, 
                     get_gds = get_gds, 
                     relatedness_genotypes = relatedness_genotypes,
                     calc_logLik = calc_logLik,
                     calc_likelihoods = calc_likelihoods,
                     create_individual = create_individual,
                     create_offspring = create_offspring,
                     simulate_relationship = simulate_relationship,
                     NSIM = NSIM))
    
    submitJobs(resources = list(max.concurrent.jobs = 20))
    waitForJobs()
    
    successful_jobs2 <- dplyr::as_tibble(getJobTable()) %>%
      dplyr::filter(is.na(error)) %>%
      dplyr::pull(job.id)
    
    if(length(successful_jobs2) != nrow(dplyr::as_tibble(getJobTable()))){
      failed_jobs <- dplyr::as_tibble(getJobTable()) %>%
        dplyr::filter(!is.na(error)) %>%
        dplyr::select(job.id, job.pars) %>%
        tidyr::unnest(job.pars) %>%
        tidyr::unnest(job.pars) 
    }
    
    message('Jobs failed after refitting round ', z, ' = ', length(simulation_settings) - length(successful_jobs) - length(successful_jobs2))
    
    simulated_relatedness2 <- purrr::map_dfr(successful_jobs2, loadResult)
    simulated_relatedness <- dplyr::bind_rows(simulated_relatedness, simulated_relatedness2)
    
    clearRegistry()
    
    successful_jobs <- c(successful_jobs, successful_jobs2)
  }
  
} else {
  simulated_relatedness <- purrr::map_dfr(simulation_settings, sims_to_node)
}

readr::write_csv(simulated_relatedness, stringr::str_replace(gds_file, '\\.gds$', '_simulated_relationships.csv'))

# Correlation between simulation and expectation
method_comparison_initial <- simulated_relatedness %>%
  dplyr::left_join(tibble::tibble(relationship = c('PO', 'FS', 'HS', 'UR',
                                                   'GG', 'AV', 'C1', 'C2',
                                                   'SC'),
                                  mean_rel = c(0.5, 0.5, 0.25, 0,
                                               0.25, 0.25, 0.125, 0.25,
                                               0.03125) / 2),
                   by = 'relationship') %>%
  dplyr::group_by(rel_method, n_loci) %>%
  dplyr::summarise(broom::tidy(cor.test(mean_rel, kinship)),
                   n = dplyr::n(),
                   .groups = 'drop') %>%
  dplyr::rename(correlation = estimate)

readr::write_csv(method_comparison_initial, stringr::str_replace(gds_file, '\\.gds$', '_relatedness_method_comparison_intervals.csv'))


method_comparison <- method_comparison_initial %>%
  dplyr::select(-statistic:-alternative) %>%
  
  tidyr::pivot_wider(names_from = 'rel_method',
                     values_from = 'correlation') %>%
  
  #https://stats.stackexchange.com/questions/278751/how-do-i-determine-whether-two-correlations-are-significantly-different
  dplyr::mutate(dplyr::across(c(EM, MoM), ~0.5 * log((1 + .) / (1 - .)), .names = '{.col}_prime'),
                S = sqrt(2 * (1 / (n - 3))),
                z = (EM_prime - MoM_prime) / S,
                p = 2 * pnorm(abs(z), lower.tail = FALSE),
                p_adj = p.adjust(p, 'holm'))

readr::write_csv(method_comparison, stringr::str_replace(gds_file, '\\.gds$', '_relatedness_method_comparison.csv'))


correlation_plot <- ggplot2::ggplot(method_comparison_initial, 
                                    ggplot2::aes(x = n_loci, y = correlation, ymin = conf.low, ymax = conf.high, colour = rel_method)) +
  # ggplot2::geom_line(show.legend = FALSE) +
  ggplot2::geom_point() +
  ggplot2::geom_text(data = dplyr::filter(method_comparison, p_adj < 0.05), colour = 'black', y = 1,
                     ggplot2::aes(x = n_loci, label = '*'), inherit.aes = FALSE) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::scale_colour_discrete(labels = c('EM' = 'ML', 'MoM')) +
  # ggplot2::scale_x_log10() +
  ggplot2::labs(x = 'Number of Loci Shared',
                y = 'Pearson Correlation Coefficient',
                colour = 'Relatedness\nMethod') +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = c(1, 0),
                 legend.justification = c(1,0))

ggplot2::ggsave(stringr::str_replace(gds_file, '\\.gds$', '_correlation_plot.svg'), 
                plot = correlation_plot, height = 7, width = 7)

ggplot2::ggsave(stringr::str_replace(gds_file, '\\.gds$', '_correlation_plot.png'), 
                plot = correlation_plot, height = 7, width = 7)


# Distingish Relationship Types
pointwise_equivilence <- simulated_relatedness %>%
  dplyr::left_join(tibble::tibble(relationship = c('PO', 'FS', 'HS', 'UR',
                                                   'GG', 'AV', 'C1', 'C2',
                                                   'SC'),
                                  mean_rel = c(0.5, 0.5, 0.25, 0,
                                               0.25, 0.25, 0.125, 0.25,
                                               0.03125) / 2),
                   by = 'relationship') %>%
  dplyr::group_by(relationship, rel_method, n_loci, mean_rel) %>%
  dplyr::summarise(broom::tidy(t.test(kinship, mu = unique(mean_rel))),
                   dplyr::bind_rows(fitdistrplus::fitdist(kinship + 1e-15, 'beta')$estimate),
                   .groups = 'drop')

readr::write_csv(pointwise_equivilence, stringr::str_replace(gds_file, '\\.gds$', '_relatedness_pointwise_equivilence.csv'))


simulation_plot <- simulated_relatedness %>%
  tibble::as_tibble() %>%
  # sample_frac(0.1) %>%
  ggplot2::ggplot(ggplot2::aes(x = n_loci, y = kinship, colour = rel_method, 
                               group = interaction(n_loci, rel_method))) +
  ggplot2::geom_hline(data = tibble::tibble(relationship = c('PO', 'FS', 'HS', 'UR',
                                                             'GG', 'AV', 'C1', 'C2',
                                                             'SC'),
                                            mean_rel = c(0.5, 0.5, 0.25, 0,
                                                         0.25, 0.25, 0.125, 0.25,
                                                         0.03125) / 2),
                      ggplot2::aes(yintercept = mean_rel)) +
  # ggbeeswarm::geom_beeswarm() +
  # ggplot2::geom_boxplot(position = ggplot2::position_dodge(5)) +
  ggplot2::stat_summary(position = ggplot2::position_dodge(5), fun.data = ggplot2::median_hilow) + #
  ggplot2::geom_text(data = pointwise_equivilence %>%
                       dplyr::mutate(p.adj = p.adjust(p.value, method = 'holm')) %>%
                       dplyr::filter(p.adj < 0.05) %>%
                       dplyr::mutate(y = dplyr::if_else(rel_method == 'EM', 0.5, -Inf)),
                     ggplot2::aes(y = y, label = '*'), show.legend = FALSE) +
  ggplot2::scale_y_continuous(limits = c(0, 0.5)) +
  ggplot2::labs(x = 'Number of Loci Shared',
                y = 'Kinship Coefficient',
                colour = 'Relatedness\nMethod') +
  ggplot2::facet_wrap(~relationship) +
  ggplot2::theme_classic()


ggplot2::ggsave(stringr::str_replace(gds_file, '\\.gds$', '_simulated_relationships.png'), 
                plot = simulation_plot, height = 15, width = 15)


simulation_plot_beta <- pointwise_equivilence %>%
  dplyr::mutate(relationship = forcats::fct_reorder(relationship, -mean_rel),
                p_adj = p.adjust(p.value, 'holm')) %>%
  dplyr::mutate(lwr = qbeta(0.025, shape1, shape2),
                upr = qbeta(0.975, shape1, shape2)) %>%
  ggplot2::ggplot(ggplot2::aes(x = n_loci, y = estimate, ymin = lwr, ymax = upr, colour = rel_method)) +
  ggplot2::geom_pointrange() +
  ggplot2::geom_hline(ggplot2::aes(yintercept = mean_rel)) +
  ggplot2::geom_text(ggplot2::aes(y = 1, label = dplyr::if_else(p_adj < 0.05, '*', "")), show.legend = FALSE) +
  ggplot2::facet_grid(relationship ~ rel_method) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(x = 'Number of Loci Shared',
                y = 'Kinship Coefficient',
                colour = 'Relatedness\nMethod') +
  ggplot2::theme_classic()

ggplot2::ggsave(stringr::str_replace(gds_file, '\\.gds$', '_simulated_relationships_beta.png'), 
                plot = simulation_plot_beta, height = 15, width = 15)

#Unrelated Miscalssification
ur_miscategorization <- simulated_relatedness %>%
  dplyr::group_by(relationship, rel_method, n_loci) %>%
  dplyr::summarise(kinship = list(kinship),
                   .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = 'relationship',
                     values_from = 'kinship') %>%
  dplyr::rowwise(rel_method, n_loci) %>%
  dplyr::summarise(dplyr:::across(tidyselect:::where(is.list), ~sum(UR >= quantile(., 0.025)) / NSIM), 
                   .groups = 'drop') %>%
  dplyr::select(-UR) %>%
  tidyr::pivot_longer(cols = -rel_method:-n_loci,
                      names_to = 'relationship',
                      values_to = 'percent_ur_misidentified') %>%
  dplyr::mutate(relationship = stringr::str_replace_all(relationship, c('AV' = 'Avuncular',
                                                                        'C1' = 'First Cousin',
                                                                        'C2' = 'Double First Cousin',
                                                                        'FS' = 'Full-sib',
                                                                        'GG' = 'Grandparent - Grandchild',
                                                                        'HS' = 'Half-sib',
                                                                        'PO' = 'Parent - Offspring',
                                                                        'SC' = 'Second Cousin'))) %>%
  dplyr::mutate(relationship = forcats::fct_reorder(relationship, percent_ur_misidentified)) 

readr::write_csv(ur_miscategorization, stringr::str_replace(gds_file, '\\.gds$', '_unrelated_miscategorization.csv'))


## Find cutoff loci
gam_model <- mgcv::gam(percent_ur_misidentified ~ relationship + rel_method + 
                         s(n_loci, interaction(relationship, rel_method), bs = 'fs', 
                           k = dplyr::n_distinct(ur_miscategorization$n_loci) - 1),
                       family = 'betar',
                       data = ur_miscategorization)

gam_model_log <- mgcv::gam(percent_ur_misidentified ~ relationship + rel_method + 
                             s(log(n_loci, base = 10), interaction(relationship, rel_method), bs = 'fs', 
                               k = dplyr::n_distinct(ur_miscategorization$n_loci) - 1),
                           family = 'betar',
                           data = ur_miscategorization)

the_aic <- AIC(gam_model, gam_model_log)

gam_model_use <- list(gam_model, gam_model_log)[[which.min(the_aic$AIC)]]

predictions <- tidyr::expand_grid(relationship = unique(ur_miscategorization$relationship),
                                  rel_method = unique(ur_miscategorization$rel_method),
                                  n_loci = modelr::seq_range(ur_miscategorization$n_loci, 1000)) %>%
  predict(gam_model_use, newdata = ., se.fit = TRUE) %>%
  purrr::map(as.numeric) %>%
  dplyr::bind_cols() %>%
  dplyr::bind_cols(tidyr::expand_grid(relationship = unique(ur_miscategorization$relationship),
                                      rel_method = unique(ur_miscategorization$rel_method),
                                      n_loci = modelr::seq_range(ur_miscategorization$n_loci, 1000)), 
                   .) %>%
  dplyr::mutate(lwr = fit - se.fit,
                upr = fit + se.fit) %>%
  dplyr::mutate(dplyr:::across(c(fit, lwr, upr), ~binomial()$linkinv(.)))

ur_miscategorization_graph <- ur_miscategorization%>%
  ggplot2::ggplot(ggplot2::aes(x = n_loci, y = percent_ur_misidentified, colour = rel_method)) +
  ggplot2::geom_point() +
  ggplot2::geom_ribbon(data = predictions, ggplot2::aes(y = fit, ymin = lwr, ymax = upr, fill = rel_method), 
                       alpha = 0.5, colour = NA, show.legend = FALSE) +
  ggplot2::geom_line(data = predictions, ggplot2::aes(y = fit)) +
  ggplot2::facet_wrap(~relationship) +
  ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(1)) +
  ggplot2::scale_colour_discrete(labels = c('EM' = 'ML', 'MoM')) +
  ggplot2::labs(x = 'Number of Loci Shared',
                y = 'Percent Unrelated Simulations in 95% CI',
                colour = 'Relatedness\nMethod') +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = c(1, 0),
                 legend.justification = c(1, 0))

ggplot2::ggsave(stringr::str_replace(gds_file, '\\.gds$', '_unrelated_miscategorization.svg'), 
                plot = ur_miscategorization_graph, height = 7, width = 7)

ggplot2::ggsave(stringr::str_replace(gds_file, '\\.gds$', '_unrelated_miscategorization.png'), 
                plot = ur_miscategorization_graph, height = 7, width = 7)

cutoff_n_loci <- predictions %>%
  dplyr::bind_rows(dplyr::distinct(predictions, relationship, rel_method) %>%
                     dplyr::mutate(fit = 0.01,
                                   n_loci = Inf)) %>%
  dplyr::group_by(relationship, rel_method) %>%
  dplyr::filter(fit < 0.05) %>%
  dplyr::filter(n_loci == min(n_loci)) %>%
  dplyr::ungroup()

readr::write_csv(cutoff_n_loci, stringr::str_replace(gds_file, '\\.gds$', '_min_number_loci.csv'))
