args <- commandArgs(trailingOnly = TRUE)
outNAME <- args[1]
# outNAME <- 'MiSeq_lightSpecies2.10.1'

suppressWarnings(suppressMessages(library(tidyverse)))

cv_data <- read_lines(str_c(outNAME, '.CV_error')) %>%
  tibble(line = .) %>%
  mutate(K = str_extract(line, '(K=)[0-9]+') %>%
           str_remove('K=') %>% as.integer,
         Error = str_extract(line, '[0-9\\.]+$') %>%
           as.numeric) %>%
  arrange(K)

cv_plot <- cv_data %>%
  ggplot(aes(x = K, y = Error)) +
  geom_line() +
  geom_point(aes(colour = Error == min(Error))) +
  labs(x = "K Clusters",
       y = 'CV - Error Rate') +
  scale_colour_manual(values = c('TRUE' = 'red', 'FALSE' = 'black')) +
  theme_classic() +
  theme(legend.position = 'none')
ggsave(str_c(outNAME, '.cv_plot.png'), plot = cv_plot, width = 7, height = 7)

write(cv_data$K[cv_data$Error == min(cv_data$Error)][1], stdout())