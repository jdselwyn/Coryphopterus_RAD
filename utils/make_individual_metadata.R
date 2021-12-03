#### Make CSV with translation of specimen ID to site/shoal/species etc. ####
library(sf)
library(tidyverse)
library(magrittr)
library(janitor)
library(readxl)

distribute_fish <- function(number_fish, centroid, spread_m = 0.25){
  #Distribute fish around centroid within a shoal. Assumes crs is in "m" and will distribute with stdev of spread_m (in m)
  if(number_fish == 1){
    out <- as.numeric(centroid) %>%
      as_tibble() %>%
      mutate(name = c('x', 'y')) %>%
      pivot_wider(names_from = 'name', values_from = 'value')
  } else {
    out <- MASS::mvrnorm(n = number_fish, mu = as.numeric(centroid), 
                         Sigma = diag(spread_m, nrow = 2, ncol = 2)) %>%
      set_colnames(c('x', 'y')) %>%
      as_tibble()
  }
  out
}

shoal_data <- list.files("~/Coryphopterus/Maps/COPE_Sites/Shoals", pattern = '.shp$', recursive = TRUE, full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(Site = str_extract(file, 'BZ17-[0-9ABNSK]+'),
         Site = factor(Site, levels = str_c('BZ17', c('10KN', '5KN', '1KN', '500N', '100N', 
                                                      '0A', '0B', '60S', '100S', '500S', '1KS', '5KS'),
                                            sep = '-'))) %>%
  mutate(shoals = map(file, st_read, quiet = TRUE)) %>%
  select(-file) %>%
  unnest(shoals) %>%
  filter(Shoal.Size > 0) %>%
  mutate(Nut = map_chr(Nut, ~str_split(.x, '-', simplify = TRUE) %>% sort %>% str_to_lower %>% str_c(collapse = '-')),
         Nut = str_replace(Nut, 'orage', 'orange')) %>%
  select(Site, Nut, Shoal.Size, geometry) %>%
  st_as_sf %>%
  st_transform(crs = "+proj=utm +zone=16Q +ellps=WGS84 +datum=WGS84 +units=m +no_defs") %>%
  as_tibble

individual_data <- unlist(str_split('cdcccccccnnnciiiccclcccccicccc','')) %>%
  tibble(col_type = .) %>%
  mutate(col_type = case_when(col_type == 'c' ~ 'text',
                              col_type == 'd' ~ 'numeric',
                              col_type == 'n' ~ 'numeric',
                              col_type == 'i' ~ 'numeric',
                              col_type == 'l' ~ 'text')) %$%
  read_excel('~/Coryphopterus/Master COPE collection database.xlsx',sheet = 1, 
             col_types = col_type, na = c('N/A','NA', "?")) %>% 
  clean_names %>%
  dplyr::select(tube_label, year, site, cloud, sl_mm, tl_mm) %>%
  dplyr::rename(ID = tube_label) %>%
  left_join(read_csv('~/Coryphopterus/Bioinformatics/Coryphopterus_RAD/splitSpecies/hybrid_types.csv', 
                     col_types = cols(.default = col_character())) %>%
              mutate(ID = str_remove(ID, '.fp2.repr') %>% str_replace('_', '-')),
            by = 'ID') %>%
  mutate(cloud = map_chr(cloud, ~str_split(.x, '-', simplify = TRUE) %>% sort %>% str_to_lower %>% str_c(collapse = '-')),
         cloud = str_remove(cloud, '^[ab]-')) %>%
  select(ID, year, site, cloud, sl_mm, tl_mm, hybrid) %>%
  inner_join(shoal_data,
             by = c('site' = 'Site', 'cloud' = 'Nut')) %>%
  rename(shoal = cloud) %>%
  nest(fish_data = -c(site, shoal, geometry)) %>%
  mutate(n_fish = map_int(fish_data, nrow)) %>%
  mutate(points = map2(n_fish, geometry, distribute_fish, spread_m = 0.06)) %>% #Shuffle around fish so none on top of each other
  select(-geometry, -n_fish) %>%
  unnest(c(fish_data, points)) %>%
  st_as_sf(coords = c('x', 'y'),
           crs = "+proj=utm +zone=16Q +ellps=WGS84 +datum=WGS84 +units=m +no_defs") %>%
  st_transform("+proj=utm +zone=16Q +ellps=WGS84 +datum=WGS84 +units=km +no_defs") %>%
  as_tibble %>%
  select(ID, year, site, shoal, sl_mm, tl_mm, hybrid, Shoal.Size, geometry) %>%
  # filter(hybrid == 'CHYA') %>%
  mutate(ID = str_c(ID, '.fp2.repr.2.1') %>% str_replace('-', '_')) %>%
  arrange(ID) 

#### Impute Missing Length ####
just_non_hybrids <- filter(individual_data, hybrid %in% c('CHYA', 'CPER'))

imputed_lengths <- Hmisc::aregImpute(~sl_mm + tl_mm + site + hybrid, data = just_non_hybrids, n.impute = 20)

impLengthMean <- imputed_lengths$imputed[1:2] %>%
  map(rowMeans)

just_non_hybrids$sl_mm[as.numeric(names(impLengthMean$sl_mm))] <- unname(impLengthMean$sl_mm)
just_non_hybrids$tl_mm[as.numeric(names(impLengthMean$tl_mm))] <- unname(impLengthMean$tl_mm)

bind_rows(just_non_hybrids,
          filter(individual_data, !hybrid %in% c('CHYA', 'CPER'))) %>%
  arrange(ID) %>%
  st_as_sf %>%
  write_sf('../individual_metadata.shp')

