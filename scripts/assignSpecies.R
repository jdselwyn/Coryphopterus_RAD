#TODO - set up args
#TODO - better method of sorting out cutting off plot
#TODO - no a priori guess

args <- commandArgs(trailingOnly = TRUE)
print(args)

mitoID_file <- args[1]
vcf_file <- args[2]
NBOOT <- as.numeric(args[3])
prefix <- args[4]

# mitoID_file <- 'Mitochondrial_Mapping/blast_speciesID.csv'
# vcf_file <- 'fltrVCF_MiSeq/MiSeq_lightSpecies.10.1.Fltr20.7.randSNPperLoc.vcf'
# NBOOT <- 10
# prefix <- 'MiSeq_lightSpecies'

outDir <- 'splitSpecies'

suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(adegenet))
suppressMessages(library(vcfR))
suppressMessages(library(ggbeeswarm))

parallel_cluster <- parallel::makeCluster(parallel::detectCores())

dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

#### Read in Data ####
vcf <- read.vcfR(vcf_file, verbose = TRUE)

vcf_genind <- vcfR2genind(vcf)

mito_id <- read_csv(mitoID_file, col_types = cols(.default = col_character())) %>%
  select(ID, species) %>%
  mutate(ID = str_c(ID, '.fp2.repr'),
         id_index = which(rownames(vcf_genind@tab) %in% ID))

#### DAPC Just Known ####
just_mito_id <- vcf_genind[mito_id$id_index, ]
strata(just_mito_id) <- as.data.frame(mito_id)
setPop(just_mito_id) <- ~species

cv_mito_dapc <- xvalDapc(tab(just_mito_id, NA.method = "mean"), grp = pop(just_mito_id),
                         n.pca = 1:(nLoc(just_mito_id) - 1), n.rep = NBOOT,
                         center = TRUE, scale = TRUE,
                         parallel = "snow", 
                         cl = parallel_cluster)

message('Number of PCs retained: ', cv_mito_dapc$`Number of PCs Achieving Highest Mean Success`, 
        '\n, RMSE = ', 
        round(cv_mito_dapc$`Root Mean Squared Error by Number of PCs of PCA`[as.integer(cv_mito_dapc$`Number of PCs Achieving Highest Mean Success`)], 3), '\n')

cv_dapc_plot <- as_tibble(cv_mito_dapc$`Cross-Validation Results`) %T>%
  write_csv(str_c(outDir, '/', prefix, '_dapc_mitoCV_results.csv')) %>%
  ggplot(aes(x = n.pca, y = success)) +
  geom_rect(xmin = -Inf, xmax = Inf,
            ymin = cv_mito_dapc$`Median and Confidence Interval for Random Chance`[1],
            ymax = cv_mito_dapc$`Median and Confidence Interval for Random Chance`[3]) +
  geom_hline(yintercept = cv_mito_dapc$`Median and Confidence Interval for Random Chance`[2]) +
  geom_beeswarm(size = 0.5, colour = 'gray50', groupOnX = TRUE) +
  stat_summary(aes(colour = if_else(as.character(n.pca) %in% cv_mito_dapc$`Number of PCs Achieving Highest Mean Success`, 'Best', 'Other')),
               fun.data = mean_cl_boot) +
  geom_smooth(method = 'gam',  formula = y ~ s(x, bs = "cs")) +
  scale_color_manual(values = c('Best' = 'red', 'Other' = 'black')) +
  scale_y_continuous(labels = scales::percent_format(1), limits = c(0, 1)) +
  labs(x = 'Number of Principle Components',
       y = 'Percent Successful Assignments') +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text(colour = 'black'),
        panel.border = element_rect(colour = 'black', fill = 'transparent'))
ggsave(str_c(outDir, '/', prefix, '_dapc_mitoCV_results.png'), plot = cv_dapc_plot, height = 7, width = 15)

#### Plot Species PCA - just mito ID'd - DAPC labels ####
mito_dapc <- cv_mito_dapc$DAPC

mito_pca <- dudi.pca(tab(just_mito_id, NA.method = "mean"), scannf=FALSE, 
                     center = TRUE, scale=TRUE, nf = nLoc(just_mito_id) - 1)

get_range <- function(vect, scal = 3){
  the_range <- range(vect)
  outlier_range <- 0 * median(vect) + c(-1 * scal, scal) * (IQR(vect)/2)
  
  x_option <- c(the_range[1], outlier_range[1])
  y_option <- c(the_range[2], outlier_range[2])

  use <- expand.grid(x_option, y_option) %>%
    mutate(range = abs(Var1 - Var2)) %>%
    filter(range == min(range)) %$%
    c(Var1, Var2)
  
  use[use %in% the_range] <- NA
  use
}

x_range <- get_range(mito_pca$li$Axis1)
y_range <- get_range(mito_pca$li$Axis2)

cut_off_left <- sum(mito_pca$li$Axis1 <= x_range[1], na.rm = TRUE)
cut_off_right <- sum(mito_pca$li$Axis1 >= x_range[2], na.rm = TRUE)
cut_off_bot <- sum(mito_pca$li$Axis2 <= y_range[1], na.rm = TRUE)
cut_off_top <- sum(mito_pca$li$Axis2 >= y_range[2], na.rm = TRUE)
the_caption <- str_c('# Cut-off outlier points:\nup = ',  cut_off_top, '; down = ', cut_off_bot,
                     '\nleft = ', cut_off_left, '; right = ', cut_off_right)

mito_pca_plot <- bind_cols(arrange(mito_id, id_index), 
          dapc_assign = as.character(mito_dapc$assign),
          mito_dapc$posterior) %>% 
  mutate(across(c(species, dapc_assign), ~str_remove(., 'Coryphopterus ')),
         id_match = if_else(species == dapc_assign, species, 
                            str_c(species, dapc_assign, sep = '_'))) %>%
  rowwise %>%
  mutate(posterior_max = max(c_across(starts_with('Coryphopterus')))) %>%
  ungroup %>%
  bind_cols(mito_pca$li) %T>%
  
  # write a csv to same dir as plot
  write_csv(str_c(outDir, '/', prefix, '_dapc_mito_pca.csv')) %>%
  
  ggplot(aes(x = Axis1, y = Axis2, colour = id_match, alpha = posterior_max)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_point() +
  scale_y_continuous(limits = y_range) +
  scale_x_continuous(limits = x_range) +
  labs(x = str_c('PC1 (', round(100 * mito_pca$eig/sum(mito_pca$eig), 1)[1], '%)'),
       y = str_c('PC2 (', round(100 * mito_pca$eig/sum(mito_pca$eig), 1)[2], '%)'),
       colour = 'ID',
       alpha = 'Posterior Probability',
       caption = the_caption) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text(colour = 'black'),
        panel.border = element_rect(colour = 'black', fill = 'transparent'))
  
ggsave(str_c(outDir, '/', prefix, '_dapc_mito_pca.png'), mito_pca_plot, height = 7, width = 7)

#### Structure Plot ####


#### Predict Assignments ####
full_pca <- dudi.pca(tab(vcf_genind, NA.method = "mean"), scannf=FALSE, 
                     center = TRUE, scale=TRUE, nf = nLoc(vcf_genind) - 1)

x_range <- get_range(full_pca$li$Axis1, 5)
y_range <- get_range(full_pca$li$Axis2, 5)

cut_off_left <- sum(full_pca$li$Axis1 <= x_range[1], na.rm = TRUE)
cut_off_right <- sum(full_pca$li$Axis1 >= x_range[2], na.rm = TRUE)
cut_off_bot <- sum(full_pca$li$Axis2 <= y_range[1], na.rm = TRUE)
cut_off_top <- sum(full_pca$li$Axis2 >= y_range[2], na.rm = TRUE)
# the_caption <- str_c('# Cut-off outlier points:\nup = ',  cut_off_top, '; down = ', cut_off_bot,
#                      '\nleft = ', cut_off_left, '; right = ', cut_off_right)
the_caption <- str_c('# Cut-off outlier points:\nup = ',  cut_off_top, '; down = ', cut_off_bot)

full_pca_plot <- predict(mito_dapc, vcf_genind) %$%
  bind_cols(dapc_prediction = assign,
            as_tibble(posterior, rownames = 'ID'),
            as_tibble(ind.scores)) %>%
  full_join(select(mito_id, -id_index), by = 'ID') %>%
  mutate(ID = str_remove(ID, '.fp2.repr')) %>%
  select(ID, species, everything()) %>%
  
  mutate(across(c(species, dapc_prediction), ~str_remove(., 'Coryphopterus ')),
         id_match = case_when(is.na(species) ~ str_c('pred', dapc_prediction, sep = '_'),
                              species == dapc_prediction ~ species, 
                              TRUE ~ str_c(species, dapc_prediction, sep = '_'))) %>%
  rowwise %>%
  mutate(posterior_max = max(c_across(starts_with('Coryphopterus')))) %>%
  ungroup %>%
  
  bind_cols(full_pca$li) %T>%
  #output csv
  write_csv(str_c(outDir, '/', prefix, '_dapc_all_pca.csv')) %>%
  
  ggplot(aes(x = Axis1, y = Axis2, colour = id_match, alpha = posterior_max)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_point() +
  scale_y_continuous(limits = y_range) +
  # scale_x_continuous(limits = x_range) +
  labs(x = str_c('PC1 (', round(100 * mito_pca$eig/sum(mito_pca$eig), 1)[1], '%)'),
       y = str_c('PC2 (', round(100 * mito_pca$eig/sum(mito_pca$eig), 1)[2], '%)'),
       colour = 'ID',
       alpha = 'Posterior Probability',
       caption = the_caption) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text(colour = 'black'),
        panel.border = element_rect(colour = 'black', fill = 'transparent'))
ggsave(str_c(outDir, '/', prefix, '_dapc_all_pca.png'), full_pca_plot, height = 7, width = 7)

#### Structure Plot Everyone ####


#### DAPC 2 Groups ####
# grp <- find.clusters(vcf_genind, max.n.clust = 15, 
#                      n.pca = nLoc(vcf_genind) - 1, n.clust = 2)
# 
# cv_mito_dapc <- xvalDapc(tab(vcf_genind, NA.method = "mean"), grp = grp$grp,
#                          n.pca = 1:(nLoc(just_mito_id) - 1), n.rep = NBOOT,
#                          center = TRUE, scale = TRUE,
#                          parallel = "snow", 
#                          cl = parallel_cluster)

