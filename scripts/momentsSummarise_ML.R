
if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  in_dir <- args[1]
  out_dir <- args[2]
  
} else {
  setwd("~/Coryphopterus/Bioinformatics/Coryphopterus_RAD")
  in_dir <- 'Moments/p15.p15/momentsOptimize'
  out_dir <- 'Moments/p15.p15'
}



suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(ggbeeswarm))


#### Functions ####
get_param_names <- function(inFILE){
  # param_names <-
  read_tsv(inFILE, col_types = cols(.default = col_character())) %>%
    clean_names() %>%
    colnames() %>%
    tibble() %>%
    rename(col_name = 1) %>%
    filter(str_detect(col_name,
                      "optimized")) %>%
    mutate(col_names = str_remove(col_name,
                                  "optimized_params_"),
           col_names = str_replace_all(col_names,
                                       "_",
                                       ",")) %>%
    pull(col_names) %>%
    str_split(pattern = ",") %>%
    as_vector()
}

get_model_results <- function(inFILE){
  # data_summaries <-
  param_names <-
    get_param_names(inFILE)
  
  read_tsv(inFILE,
           na = c("",
                  "NA",
                  "nan",
                  "--"),
           col_types = cols(
             .default = col_character(),
             `log-likelihood` = col_double(),
             AIC = col_double(),
             `chi-squared` = col_double(),
             theta = col_double()
           )) %>%
    clean_names() %>%
    #if there are multiple headers, remove additional
    filter(model != "Model") %>%
    #if there are identical rows, remove them, this is from multiple runs in same dir without deleting
    distinct() %>%
    mutate(summary_file_name = inFILE) %>%
    separate(7,
             into = param_names,
             sep = ",",
             remove = FALSE) %>%
    separate(replicate,
             into = c("round",
                      "replicate"),
             sep = "_Replicate_") %>%
    mutate(round = as_factor(str_remove(round,
                                        "Round_")),
           replicate = as_factor(replicate)) %>%
    arrange(model,
            desc(log_likelihood))
}

get_all_model_results <- function(inFILES){
  inFILES %>%
    purrr::map(get_model_results) %>%
    bind_rows() %>%
    select(model:theta,
           starts_with("nu"),
           starts_with("m"),
           starts_with("t"),
           contains("optimized_params"))
}

#### Read in Data ####
file_names_optimized <- list.files(in_dir, full.names = TRUE) %>%
  str_subset("optimized.txt")

data <- get_all_model_results(file_names_optimized) %>%
  mutate(across(c(round, replicate, nu1:t3), as.numeric))

setwd(out_dir)

#### Visualize model logLike & AIC change over rounds ####
reorder_function <- function(x, y){
  median(x[y == max(y)], na.rm = TRUE)
}

models_through_rounds_ll <- data %>%
  filter(! is.na(chi_squared) & chi_squared >= 0) %>%
  #Reorder models based on which has the highest median logLik in final round
  mutate(model = fct_reorder2(model, log_likelihood, round,
                              .fun = reorder_function)) %>%
  
  ggplot(aes(x = round, y = log_likelihood)) +
  geom_beeswarm(groupOnX = TRUE, show.legend = FALSE, size = 0.5) +
  stat_summary(fun.data = mean_cl_boot) +
  facet_wrap(~model) +
  labs(x = 'Round',
       y = 'log(Likelihood)') +
  theme_classic() +
  theme(axis.title = element_text(colour = 'black', size = 20),
        axis.text = element_text(colour = 'black', size = 18),
        strip.text = element_text(colour = 'black', size = 20),
        strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'))
ggsave('model_logLikelihood_overRounds.png', plot = models_through_rounds_ll,
       height = 15, width = 15)

models_through_rounds_aic <- data %>%
  filter(! is.na(chi_squared) & chi_squared >= 0) %>%
  #Reorder models based on which has the lowest median aic in final round
  mutate(model = fct_reorder2(model, -aic, round,
                              .fun = reorder_function)) %>%
  
  ggplot(aes(x = round, y = aic)) +
  geom_beeswarm(groupOnX = TRUE, show.legend = FALSE, size = 0.5) +
  stat_summary(fun.data = mean_cl_boot) +
  facet_wrap(~model) +
  labs(x = 'Round',
       y = 'AIC') +
  theme_classic() +
  theme(axis.title = element_text(colour = 'black', size = 20),
        axis.text = element_text(colour = 'black', size = 18),
        strip.text = element_text(colour = 'black', size = 20),
        strip.background = element_blank(),
        panel.background = element_rect(colour = 'black'))
ggsave('model_AIC_overRounds.png', plot = models_through_rounds_aic,
       height = 15, width = 15)


#### Just loglikelihood of final round ####
model_logLike <- data %>%
  filter(! is.na(chi_squared) & chi_squared >= 0) %>%
  group_by(model) %>%
  filter(round == max(round)) %>%
  mutate(max_loglikelihood = max(log_likelihood)) %>%
  ungroup %>%
  mutate(model = fct_reorder(model, log_likelihood)) %>%
  
  ggplot(aes(y = model,
             x = log_likelihood,
             colour = max_loglikelihood)) +
  geom_beeswarm(groupOnX = FALSE, size = 0.5, show.legend = FALSE) +
  stat_summary(fun.data = mean_cl_boot, colour = 'black') +
  labs(y = NULL,
       x = 'log(Likelihood)') +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black', size = 16),
        axis.title = element_text(size = 20))
ggsave('model_logLikelihood_finalRound.png', plot = model_logLike,
       height = 12, width = 8)

model_aic <- data %>%
  filter(! is.na(chi_squared) & chi_squared >= 0) %>%
  group_by(model) %>%
  filter(round == max(round)) %>%
  mutate(min_aic = min(aic)) %>%
  ungroup %>%
  mutate(model = fct_reorder(model, -aic)) %>%
  
  ggplot(aes(y = model,
             x = aic,
             colour = -min_aic)) +
  geom_beeswarm(groupOnX = FALSE, size = 0.5, show.legend = FALSE) +
  stat_summary(fun.data = mean_cl_boot, colour = 'black') +
  labs(y = NULL,
       x = 'AIC') +
  theme_classic() +
  theme(axis.text = element_text(colour = 'black', size = 16),
        axis.title = element_text(size = 20))
ggsave('model_AIC_finalRound.png', plot = model_aic,
       height = 12, width = 8)

#### Check Model Failure Rate ####
failure_plot <- data %>%
  filter(! is.na(chi_squared) & chi_squared >= 0) %>%
  group_by(model,
           round) %>%
  summarize(n = n(), .groups = 'drop') %>%
  ggplot(aes(x = model,
             y = n)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Number of Successful Replicates") +
  facet_grid(round ~ .,
             scales ="free_y")
ggsave('model_successRate.png', plot = failure_plot,
       height = 10, width = 10)
