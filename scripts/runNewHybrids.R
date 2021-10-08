##TODO - better way to pick loci to use??
if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  DIR <- args[1]
  VCF <- args[2]
  ADMIXTURE <- args[3]
  burnin <- as.integer(args[4])
  iterations <- as.integer(args[5])
  thin <- as.integer(args[6])
  chains <- as.integer(args[7])
  
} else {
  setwd('..')
  DIR <- 'splitSpecies/newHybrids'
  VCF <- 'fltrVCF_MiSeq/MiSeq_lightSpecies2.10.1.Fltr20.8.randSNPperLoc.vcf'
  ADMIXTURE <- 'splitSpecies/ADMIXTURE/MiSeq_lightSpecies2.10.1.2.results.csv'
  burnin <- 100
  iterations <- 1000
  thin <- 100
  chains <- 5
}




suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(vcfR))
suppressMessages(library(furrr))
suppressMessages(library(posterior))
suppressMessages(library(patchwork))

#### Functions ####
write_newhybrid_data <- function(data, path, loci, missing_cutoff){
  if(!is.na(loci)){
    out <- data[,c(1,2,sample(ncol(data) - 2, loci) + 2)] %>%
      rowwise %>%
      mutate(missing = rowSums(across(starts_with('L'), ~. == '0')) / loci) %>%
      ungroup %>%
      filter(missing <= missing_cutoff) %>%
      select(-missing) %>%
      mutate(ID = str_c(row_number(), ID, sep = '\t'))
  } else {
    out <- data %>%
      mutate(ID = str_c(row_number(), ID, sep = '\t'))
  }
  
  
  write_csv(out, str_c(path, 'newhybrid_data.csv', sep = '/'))
  
  out_file <- str_c(path, 'cope.dat', sep = '/')
  
  tibble(setting = c('NumIndivs', 'NumLoci', 'Digits', 'Format'), 
         value = c(nrow(out), ncol(out) - 2, '1', 'Lumped')) %>%
    write_delim(out_file, delim = ' ', col_names = FALSE)
  write_lines("", out_file, append = TRUE)
  write_lines(str_c(c('LocusNames', colnames(select(out, starts_with('L')))), collapse = ' '),
              out_file, append = TRUE)
  write_lines("", out_file, append = TRUE)
  write_delim(select(out, -Indiv), out_file, delim = ' ', col_names = FALSE, append = TRUE)
}

run_newHybrids <- function(prior, chain, data, newHybrids_path, path, burn, iter, thin, loci, missing_cutoff){
  start_dir <- getwd()
  working_dir <- str_c(path, prior, chain, sep = '/')
  dir.create(working_dir, recursive = TRUE, showWarnings = FALSE)
  
  if(Sys.info()['sysname'] == 'Windows'){
    file.copy(from = newHybrids_path,
              to = paste0(working_dir, '/', "NewHybrids_PC_1_1_WOG.exe"))
  }
  
  
  write_lines("6", str_c(working_dir, "freqTwoGeneration.dat", sep = '/'), append = FALSE)
  tibble(type = c('pure_a', 'pure_b', 'f1', 'f2', 'a_bx', 'b_bx'), 
         a = c(1, 0, 0, 0.25, 0.5, 0), 
         b = c(0, 0, 0.5, 0.25, 0.25, 0.25), 
         c = c(0, 0, 0.5, 0.25, 0.25, 0.25), 
         d = c(0, 1, 0, 0.25, 0, 0.5)) %>%
    write_delim(str_c(working_dir, "freqTwoGeneration.dat", sep = '/'), 
                delim = ' ', col_names = FALSE, append = TRUE)
  
  write_newhybrid_data(data, working_dir, loci, missing_cutoff)

  setwd(working_dir)

  if(Sys.info()['sysname'] == 'Windows'){
    out <- system2("NewHybrids_PC_1_1_WOG.exe",
                   input = c("cope.dat",
                             "freqTwoGeneration.dat",
                             "0", #no priors for allele freq
                             sample(1:20, 2), #rng
                             as.integer(prior != 'Jefferies'), #jefferies priors
                             as.integer(prior != 'Jefferies'), #jefferies priors
                             format(burn, scientific = FALSE),
                             format(iter, scientific = FALSE),
                             'X'),
                   stdout = 'terminal.out',
                   wait = TRUE)
    
    #Clean up unnecessary files
    str_c('./', 
          c('NewHybrids_PC_1_1_WOG.exe',
            'cope.dat', 'freqTwoGeneration.dat', 
            'EchoedGtypData.txt', 'EchoedGtypFreqCats.txt',
            'aa-ProcessedPriors.txt', 'aa-LociAndAlleles.txt'),
          sep = '/') %>%
      file.remove()
  } else {
    
    the_command <- str_c(
      'module load newhybrids/gcc7/2.0; newhybs',
      '-d cope.dat',
      '-c freqTwoGeneration.dat',
      str_c('--theta-prior', prior, sep = ' '),
      str_c('--pi-prior', prior, sep = ' '),
      str_c('--burn-in', burn, sep = ' '),
      str_c('--num-sweeps', iter, sep = ' '),
      str_c('--seeds', str_c(sample(1:20, 2), collapse = ' '), sep = ' '),
      '--no-gui',
      str_c('--print-traces Pi', thin, sep = ' '),
      '> terminal.out',
      sep = ' '
    )
    
    out <- system(the_command,
                  wait = TRUE)
  }
  
  setwd(start_dir)
  working_dir
}

read_newhybrids_trace <- function(file){
  out_file <- read_lines(file)
  
  if(Sys.info()['sysname'] == 'Windows'){
    out <- tibble(.iteration = out_file[str_which(out_file, 'Currently on Sweep')],
                  pi_values = out_file[str_which(out_file, 'Currently on Sweep') + 5]) %>%
      mutate(sweep = str_extract(sweep, '[0-9]+') %>% parse_integer,
             pi_values = str_trim(pi_values)) %>%
      separate(pi_values, sep = '[ ]+', into = c('pure_a', 'pure_b', 'F1', 'F2', 'a_bx', 'b_bx')) %>%
      mutate(across(where(is.character), as.numeric))
  } else {
    out <- out_file[str_which(out_file, 'PI_TRACE:')] %>%
      tibble(lines = .) %>%
      slice(-1:-2) %>%
      separate(lines, sep = '\t', 
               into = c('sweep', 'pure_a', 'pure_b', 
                        'F1', 'F2', 'a_bx', 'b_bx')) %>%
      mutate(sweep = str_remove(sweep, 'PI_TRACE:')) %>%
      mutate(across(everything(), parse_number),
             sweep = as.integer(sweep))
  }
  
}

read_newhybrids_results <- function(path, results_file, data_file){
  
  if(Sys.info()['sysname'] == 'Windows'){
    results <- str_c(path, results_file, sep = '/') %>%
      read_delim(delim = '\t', 
                 col_types = cols(.default = col_double()),
                 skip = 1,
                 col_names = c('id', 'pure_a', 'pure_b', 'F1', 'F2', 'a_bx', 'b_bx')) 
  } else {
    results <- str_c(path, results_file, sep = '/') %>%
      read_delim(delim = '\t', 
                 col_types = cols(.default = col_double(),
                                  tmp = col_character()),
                 skip = 1,
                 col_names = c('id', 'tmp', 'pure_a', 'pure_b', 'F1', 'F2', 'a_bx', 'b_bx')) %>%
      select(-tmp)
  }
  
  
  data <- str_c(path, data_file, sep = '/') %>%
    read_csv(col_types = cols(.default = col_double(),
                              Indiv = col_character(),
                              ID = col_character())) %>%
    select(-starts_with('L')) %>%
    mutate(id = row_number() - if_else(Sys.info()['sysname'] == 'Windows', 1L, 0L))
  
  full_join(data, results, by = 'id') %>% 
    mutate(Z = str_extract(ID, 'z[0-1]')) %>%
    select(Indiv, Z, everything()) %>%
    select(-id, -ID)
}

#### Read in Data ####
raw_vcf <- read.vcfR(VCF)

full_vcf <- vcfR2tidy(raw_vcf)

admix_data <- read_csv(ADMIXTURE,
                       col_types = cols(ID = col_character(),
                                        Cluster1 = col_double(),
                                        Cluster2 = col_double(),
                                        admix_assignment = col_character())) %>%
  select(-admix_assignment) %>%
  rowwise %>%
  mutate(max_prob = max(c_across(starts_with('Cluster'))),
         which_clust = which.max(c_across(starts_with('Cluster')))) %>%
  ungroup %>%
  mutate(admix_id = if_else(round(max_prob, 4) >= 1, 
                            as.character(which_clust - 1), 'U')) %>%
  select(ID, admix_id) %>%
  rename(Indiv = ID)

#### Setup NewHybrids ####
new_hybrids_data <- full_vcf$gt %>%
  select(ChromKey, POS, Indiv, gt_GT) %>%
  unite(locus, ChromKey, POS, sep = '_') %>%
  mutate(locus = str_c('L', locus)) %>%
  separate(gt_GT, into = c('A1', 'A2'), sep = '/') %>%
  mutate(across(c(A1, A2), ~as.integer(.) + 1L)) %>%
  unite(gt_GT, A1, A2, sep = '') %>%
  mutate(gt_GT = if_else(gt_GT == 'NANA', NA_character_, gt_GT)) %>%
  pivot_wider(names_from = 'locus', 
              values_from = 'gt_GT') %>%
  mutate(across(starts_with('L'), ~replace_na(., '0'))) %>%
  full_join(admix_data, tst, by = 'Indiv') %>%
  arrange(Indiv) %>%
  mutate(ID = if_else(admix_id != 'U', str_c('z', admix_id), '')) %>%
  select(Indiv, ID, starts_with('L')) 

#### Run NewHybrids ####
strat <- if_else(Sys.info()['sysname'] == 'Windows', 'multisession', 'multicore')
plan(strategy = strat)


new_hybrids_results <- expand_grid(.chain = 1:chains,
                                   prior = c('Jeffreys', 'uniform')) %>% 
  mutate(new_hybrids = future_map2_chr(prior, .chain, run_newHybrids,
                                       data = new_hybrids_data,
                                       newHybrids_path = if_else(Sys.info()['sysname'] == 'Windows', 'scripts/NewHybrids_PC_1_1_WOG.exe', NULL),
                                       path = DIR,
                                       burn = burnin,
                                       iter = iterations,
                                       thin = thin,
                                       loci = 200,
                                       missing_cutoff = 0.99,
                                       .progress = Sys.info()['sysname'] == 'Windows',
                                       .options = furrr_options(seed = TRUE))) %>%
  mutate(trace = str_c(new_hybrids, 'terminal.out', sep = '/'),
         trace = map(trace, read_newhybrids_trace)) %>%
  mutate(results = map(new_hybrids, read_newhybrids_results,
                       results_file = 'aa-PofZ.txt', 
                       data_file = 'newhybrid_data.csv'))

#### Diagnostics ####
traces <- new_hybrids_results %>%
  select(.chain, prior, trace) %>%
  unnest(trace) %>%
  group_by(.chain, prior) %>%
  mutate(.iteration = 1:n(),
         max_iter = max(.iteration)) %>%
  ungroup %>%
  filter(.iteration <= min(max_iter)) %>%
  select(-max_iter) 

trace_plot <- traces %>%
  pivot_longer(cols = c(pure_a, pure_b, F1, F2, a_bx, b_bx)) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = sweep, y = value, colour = as.character(.chain))) +
  geom_line(show.legend = FALSE) +
  facet_grid(name ~ prior, scales = 'free_y') +
  labs(x = 'iteration',
       y = 'Proportion of Population',
       colour = 'Chain') +
  theme_classic()
ggsave(str_c(DIR, 'newHybrids_trace.png', sep = '/'), plot = trace_plot,
       height = 14, width = 27)

rhats <- traces %>%
  select(-sweep) %>%
  nest(data = -c(prior)) %>%
  rowwise(prior) %>%
  mutate(data = list(select(data, where(~!all(is.na(.x))))),
         data = list(as_draws(data))) %>%
  summarise(summarise_draws(data),
            .groups = 'drop') %T>% 
  print %T>%
  write_csv(str_c(DIR, 'newHybrids_rhat.csv', sep = '/'))

#### Results ####
processed_results <- new_hybrids_results %>%
  select(.chain, prior, results) %>%
  unnest(results) %>%
  pivot_longer(cols = -c(.chain, prior, Indiv, Z),
               names_to = 'hybrid_class',
               values_to = 'posterior_probability') %>%
  pivot_wider(names_from = '.chain',
              values_from = 'posterior_probability',
              names_prefix = 'C') %>%
  
  rowwise %>%
  mutate(mean_posterior_probability = mean(c_across(starts_with('C')), na.rm = TRUE),
         sd_posterior_probability = sd(c_across(starts_with('C')), na.rm = TRUE),
         cv_posterior_probability = sd_posterior_probability / mean_posterior_probability) %>%
  ungroup %>%
  nest(data = -c(prior)) %>%
  rowwise(prior) %>%
  summarise(full_join(data, 
                      select(new_hybrids_data, Indiv) %>%
                        expand_grid(hybrid_class = c('pure_a', 'pure_b',
                                                     'F1', 'F2', 
                                                     'a_bx', 'b_bx')),
                      by = c('Indiv', 'hybrid_class')),
            .groups = 'drop') %>%
  arrange(Indiv) %T>%
  write_csv(str_c(DIR, 'newHybrids_fullResults.csv', sep = '/'))



# processed_results %>%
#   filter(!is.nan(cv_posterior_probability)) %>%
#   ggplot(aes(x = hybrid_class, y = cv_posterior_probability)) +
#   geom_boxplot()

simplified_results <- processed_results %>%
  select(-starts_with('C'), -sd_posterior_probability, -cv_posterior_probability) %>%
  group_by(prior, Indiv) %>%
  mutate(mean_posterior_probability = mean_posterior_probability / sum(mean_posterior_probability)) %>%
  ungroup %>%
  pivot_wider(names_from = hybrid_class,
              values_from = mean_posterior_probability) %>%
  
  nest(data = -prior) %T>%
  mutate(file_out = map2(prior, data, ~write_csv(.y, str_c(DIR, '/', 'newHybrids_', .x, '.csv')))) %>%
  unnest(data)


#### Individual Results #### 
all_fish_hybridization <- simplified_results %>%
  drop_na(where(is.numeric)) %>%
  nest(data = -c(prior)) %>%
  mutate(plot = map2(data, prior, ~.x %>%
                      rowwise() %>%
                      mutate(max = max(c_across(where(is.numeric))),
                             match = which.max(c_across(where(is.numeric)))) %>%
                      ungroup %>%
                      arrange(match, 
                              -max) %>%
                      mutate(id = row_number()) %>%
                      select(-match, -max) %>%
                      pivot_longer(cols = -c(Indiv, Z, id)) %>%
                      ggplot(aes(x = id, y = value, fill = name)) +
                      geom_bar(position="stack", stat="identity") +
                      scale_x_continuous(expand = c(0, 0)) +
                      scale_y_continuous(expand = c(0, 0)) +
                      scale_fill_manual(values = c('a_bx' = 'red', 'b_bx' = 'green', 
                                                   'F1' = 'blue', 'F2' = 'purple',
                                                   'pure_a' = 'gray10', 'pure_b' = 'gray70')) +
                      labs(x = NULL,
                           y = 'Assignment Probability',
                           fill = NULL,
                           title = str_to_sentence(.y)) +
                      theme_classic() +
                      theme(legend.direction = "horizontal", 
                            legend.position = "bottom",
                            legend.box = "vertical",
                            axis.title = element_text(colour = 'black', size = 20),
                            strip.text = element_text(colour = 'black', size = 20),
                            axis.text = element_text(colour = 'black', size = 16),
                            legend.text = element_text(colour = 'black', size = 16),
                            axis.text.x = element_blank(),
                            panel.border = element_rect(colour = 'black', fill = 'transparent')))) %>%
  pull(plot) %>%
  wrap_plots(ncol = 1) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
  
ggsave(str_c(DIR, 'newHybrids_all_prob.png', sep = '/'), all_fish_hybridization,
       height = 7, width = 21)
