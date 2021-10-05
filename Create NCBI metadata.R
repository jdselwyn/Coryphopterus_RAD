library(tidyverse)
library(sf)
library(readxl)
library(janitor)
library(lubridate)

#### Gather all the data together ####
reads <- read_delim('sample_list.txt', delim = '\t', col_names = 'file_name', col_types = cols(file_name = col_character())) %>%
  mutate(ID = str_extract(file_name, '[0-9]+'),
         species = str_extract(file_name, '^C[A-Z]+'),
         sequencing_type = str_extract(file_name, 'miseq|novaseq'),
         read = str_extract(file_name, 'r[12]')) 

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
  mutate(ID = str_remove(ID, 'COPE-')) %>%
  mutate(cloud = map_chr(cloud, ~str_split(.x, '-', simplify = TRUE) %>% sort %>% str_to_lower %>% str_c(collapse = '-')),
         cloud = str_remove(cloud, '^[ab]-'))

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
  select(Site, Nut, geometry) %>%
  filter(Nut != '') %>%
  st_as_sf

dive_log <- read_excel('~/../Google Drive/Diving/Dive Log.xlsx') %>%
  select(Date, `Dive Site`, `Country/Region`) %>%
  mutate(Date = as_date(Date),
         year = year(Date)) %>%
  filter(`Country/Region` == 'Belize') %>%
  filter(year == 2017) %>%
  group_by(`Dive Site`) %>%
  filter(Date == max(Date)) %>%
  sample_n(1) %>%
  ungroup %>%
  select(-`Country/Region`)

#### Join Data ####
get_tl <- lm(tl_mm ~ sl_mm, data = individual_data)

ncbi_data <- right_join(individual_data, reads, by = 'ID') %>%
  nest(files = c(file_name, read)) %>%
  nest(sequencing = c(sequencing_type, files)) %>%
  nest(fish = -c(site, cloud, year)) %>%
  inner_join(shoal_data, by = c('site' = 'Site', 'cloud' = 'Nut')) %>%
  inner_join(dive_log, by = c('site' = 'Dive Site', 'year')) %>%
  
  unnest(fish) %>%
  mutate(sample_name = str_c(species, ID, sep = '_')) %>%
  select(-cloud, -ID) %>%
  mutate(st_coordinates(geometry) %>% as_tibble) %>%
  rename(lon = X, lat = Y) %>%
  mutate(lat_lon = str_c(lat, ' N ', abs(lon), ' W'),
         species = case_when(species == 'CHYA' ~ 'Coryphopterus hyalinus',
                             species == 'CPER' ~ 'Coryphopterus personatus',
                             species == 'CSP' ~ 'Coryphopterus sp.'),
         site = str_c('Belize: Turneffe Atoll, ', site),
         tl_mm = case_when(is.na(tl_mm) & !is.na(sl_mm) ~ predict(get_tl, newdata = data.frame(sl_mm)),
                           TRUE ~ tl_mm),
         age = str_c('~', round(-0.85 + 2.99 * tl_mm, 1), ' days post hatch'),
         description = str_c(str_replace_na(sl_mm, 'unmeasured'), ' mm Standard Length specimen')) %>%
  select(sample_name, species, age, Date, site, lat_lon, description, sequencing)

write_csv(select(ncbi_data, -sequencing), 'for_ncbi.csv')

range(ncbi_data$age, na.rm = TRUE)


ncbi_data %>%
  mutate(age = str_extract(age, '[0-9\\.]+') %>% as.numeric()) %>%
  arrange(-age) %>%
  ggplot(aes(x = age)) +
  geom_histogram()

#Make SRA Metadata 
sra_info <- ncbi_data %>%
  select(sample_name, species, sequencing) %>%
  unnest(sequencing) %>%
  unnest(files) %>%
  pivot_wider(names_from = 'read',
              values_from = 'file_name') %>%
  mutate(title = str_c('ddRAD of ', species, ': Caudal fin clip'),
         library_id = str_c(sample_name, sequencing_type, sep = '.')) %>%
  arrange(sequencing_type, sample_name) %>%
  select(sample_name, library_id, title, r1, r2)
write_csv(sra_info, 'for_ncbi.csv')
