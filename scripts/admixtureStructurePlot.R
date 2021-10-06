args <- commandArgs(trailingOnly = TRUE)
outNAME <- args[1]
K <- as.integer(args[2])

# outNAME <- 'MiSeq_lightSpecies2.10.1'
# K <- 2

suppressWarnings(suppressMessages(library(tidyverse)))


sample_info <- str_c(outNAME, '.nosex') %>%
  read_delim(delim = '\t', 
             col_names = c('ID1', 'ID2'),
             col_types = c(.default = col_character())) %>%
  unite(ID, ID1, ID2, sep = '_')

cluster_assignment <- str_c(outNAME, K, 'Q', sep = '.') %>%
  read_delim(delim = ' ',
             col_names = str_c('Cluster', 1:K, sep = ''),
             col_types = c(.default = col_number()))

full_results <- bind_cols(sample_info, cluster_assignment) %>%
  rowwise %>%
  mutate(admix_assignment = which.max(c_across(starts_with('Cluster'))) %>% 
           str_c('Cluster', .)) %>%
  ungroup
write_csv(full_results, str_c(outNAME, K, 'results.csv', sep = '.'))

struct_plot <- full_results %>%
  rowwise() %>%
  mutate(max = max(c_across(starts_with('Cluster'))),
         match = which.max(c_across(starts_with('Cluster')))) %>%
  ungroup %>%
  arrange(match, 
          -max) %>%
  mutate(number = row_number()) %>%
  pivot_longer(cols = starts_with('Cluster')) %>%
  ggplot(aes(x = number, y = value, fill = name)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(~ admix_assignment, scales = 'free_x', space = 'free_x') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL,
       y = 'Assignment Probability',
       fill = NULL) +
  theme_classic() +
  theme(legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "vertical",
        axis.text = element_text(colour = 'black'),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour = 'black', fill = 'transparent'))
ggsave(str_c(outNAME, K, 'structure_plot.png', sep = '.'), 
       plot = struct_plot, height = 7, width = 21)
