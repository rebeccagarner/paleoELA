# ELA phytoplankton species cell count, density, and biomass data

# Load libraries
library(tidyverse)
library(lubridate)
library(scales)
library(ggforce)
library(ggnewscale)
library(patchwork)

# Load palettes
source("scripts/00-palettes.R")


#### Import and format data ####
# Import phytoplankton species count and biomass data
phytoplankton_species <- read_tsv("output/ela/ela_monitoring_phytoplankton_species.tsv", col_names = TRUE)

# Import and format experimental manipulation data
ela_manipulations <- read_tsv("data/ela/ela_manipulations.txt", col_names = TRUE, comment = "#")


#### Format seasons ####
# Define function to assign season from date
assignSeason <- function(date) {
  month_date_numeric <- 100 * month(date) + day(date)
  season_cuts <- cut(month_date_numeric, breaks = c(0, 319, 620, 921, 1220, 1231))
  levels(season_cuts) <- c("Winter","Spring","Summer","Fall","Winter")
  return(season_cuts)
}

# Assign season to sampling dates
phytoplankton_species <- phytoplankton_species %>%
  mutate(season = assignSeason(date))

# Summarize season start and end dates of phytoplankton monitoring
phytoplankton_species %>%
  mutate(year = year(date)) %>%
  group_by(lake_id, stratum, year, season) %>%
  summarize(season_start_date = min(date),
            season_end_date = max(date)) %>%
  ungroup()

# Plot sampling frequency
season_days <- tibble(season = c("Winter", "Spring", "Summer", "Fall", "Winter"),
                      season_start_day = c(1, 79, 172, 265, 355),
                      season_end_day = c(78, 171, 264, 354, 365))

(phytoplankton_species_freq_plot <- phytoplankton_species %>%
  group_by(lake_id, date, start_depth, end_depth, stratum, season) %>%
  summarize(total_cell_count = sum(cell_count),
            total_biomass = sum(biomass)) %>%
  ungroup() %>%
  mutate(year = year(date),
         year_day = yday(date)) %>%
  ggplot() +
  facet_grid(factor(stratum, levels = c("Epilimnion", "Metalimnion", "Hypolimnion"))~
               factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = season_start_day, ymax = season_end_day,
                fill = season),
            data = season_days,
            alpha = 0.1) +
  scale_fill_manual(values = palette_season, name = "Season") +
  new_scale_fill() +
  geom_tile(aes(x = year, y = year_day, fill = total_biomass)) +
  scale_y_reverse(expand = c(0,0)) +
  scale_fill_gradientn(labels = comma,
                       colours = topo.colors(5), na.value = "#EBEBEB",
                       trans = "log10",
                       name = "Biomass (\u00b5g/L)") +
  labs(y = "Day") +
  theme_bw() %+replace%
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        panel.grid = element_blank()))
#ggsave("figures/ela_phytoplankton_species_freq_plot.pdf", phytoplankton_species_freq_plot, width = 18, height = 8, units = "in", device = cairo_pdf)


#### Assess water column sampling coverage by stratum and depth ####
hypolimnion_depths <- phytoplankton_species %>%
  filter(stratum == "Hypolimnion") %>%
  distinct(lake_id, date, stratum, start_depth, end_depth) %>%
  pivot_longer(c(start_depth, end_depth), names_to = "location", values_to = "depth") %>%
  pivot_wider(names_from = c(stratum, location), values_from = depth)

epiliminon_metalimnion_depths <- phytoplankton_species %>%
  filter(stratum == "Epilimnion" | stratum == "Metalimnion") %>%
  distinct(lake_id, date, stratum, start_depth, end_depth) %>%
  pivot_longer(c(start_depth, end_depth), names_to = "location", values_to = "depth") %>%
  pivot_wider(names_from = c(stratum, location), values_from = depth)

strata_depths <- epiliminon_metalimnion_depths %>%
  full_join(hypolimnion_depths, by = c("lake_id", "date"))

strata_depths %>%
  filter(!is.na(Hypolimnion_start_depth) | !is.na(Hypolimnion_end_depth))

# Remove hypolimnion samples
phytoplankton_species <- phytoplankton_species %>%
  filter(stratum == "Epilimnion" | stratum == "Metalimnion")


#### Summarize data by taxonomic groups ####
# Assess non-phytoplankton taxa (e.g., Choanoflagellates)
phytoplankton_species %>%
  distinct(supergroup_pr2, division_pr2, subdivision_pr2, class_pr2, trophic_mode) %>%
  arrange(supergroup_pr2, division_pr2, subdivision_pr2, class_pr2) %>%
  filter(trophic_mode != "phototrophs")


#### Format phytoplankton species monitoring data ####
# Format site by taxon table (cell counts)
sample_by_taxon_cellcount <- phytoplankton_species %>%
  mutate(sample_id = str_c(lake_id, date, stratum, sep = "_")) %>%
  dplyr::select(sample_id, taxon_code, cell_count) %>%
  pivot_wider(names_from = taxon_code, values_from = cell_count, values_fill = 0)

# Format site by taxon table (cell density)
sample_by_taxon_celldensity <- phytoplankton_species %>%
  mutate(sample_id = str_c(lake_id, date, stratum, sep = "_")) %>%
  dplyr::select(sample_id, taxon_code, cell_density) %>%
  pivot_wider(names_from = taxon_code, values_from = cell_density, values_fill = 0)

# Format site by taxon table (biomass)
sample_by_taxon_biomass <- phytoplankton_species %>%
  mutate(sample_id = str_c(lake_id, date, stratum, sep = "_")) %>%
  dplyr::select(sample_id, taxon_code, biomass) %>%
  pivot_wider(names_from = taxon_code, values_from = biomass, values_fill = 0)

# Format phytoplankton taxonomy (and function)
phytoplankton_taxonomy <- phytoplankton_species %>%
  distinct(taxon_code, supergroup_pr2, division_pr2, subdivision_pr2, class_pr2, order_pr2, family_pr2, genus_pr2, trophic_mode) %>%
  arrange(supergroup_pr2, division_pr2, subdivision_pr2, class_pr2, order_pr2, family_pr2, genus_pr2, trophic_mode)

# Parse sample information
samples <- phytoplankton_species %>%
  mutate(sample_id = str_c(lake_id, date, stratum, sep = "_")) %>%
  distinct(sample_id, lake_id, date, stratum)

# Convert site by taxon tables to long format (to fill in missing zeros)
site_by_taxon_long <- sample_by_taxon_cellcount %>%
  pivot_longer(!sample_id, names_to = "taxon_code", values_to = "cell_count") %>%
  left_join(sample_by_taxon_celldensity %>%
              pivot_longer(!sample_id, names_to = "taxon_code", values_to = "cell_density"),
            by = c("sample_id", "taxon_code")) %>%
  left_join(sample_by_taxon_biomass %>%
              pivot_longer(!sample_id, names_to = "taxon_code", values_to = "biomass"),
            by = c("sample_id", "taxon_code")) %>%
  left_join(samples, by = "sample_id") %>%
  dplyr::select(lake_id, date, stratum, taxon_code, cell_count, cell_density, biomass)


#### Calculate time-weighted means for phytoplankton biomass ####
# Summarize data by month for water column strata (based on summarizeByMonthWC)
phytoplankton_month_means <- site_by_taxon_long %>%
  mutate(year = year(date),
         month = month(date)) %>%
  group_by(lake_id, year, month, stratum, taxon_code) %>%
  summarize(mean_cellcount = mean(cell_count),
            sd_cellcount = sd(cell_count),
            
            mean_celldensity = mean(cell_density),
            sd_celldensity = sd(cell_density),
            
            mean_biomass = mean(biomass),
            sd_biomass = sd(biomass),
            
            ndays = n()) %>%
  ungroup() %>%
  arrange(lake_id, year, month,
          factor(stratum, levels = c("Epilimnion", "Metalimnion")),
          taxon_code)

# Summarize data by month-season for water column strata (based on summarizeByMonthSeasonWC)
phytoplankton_monthseason_means <- phytoplankton_month_means %>%
  mutate(month_season = case_when(month == 1 |
                                    month == 2 |
                                    month == 3 ~ "Winter",
                                  month == 4 |
                                    month == 5 |
                                    month == 6 ~ "Spring",
                                  month == 7 |
                                    month == 8 |
                                    month == 9 ~ "Summer",
                                  month == 10 |
                                    month == 11 |
                                    month == 12 ~ "Fall")) %>%
  group_by(lake_id, year, month_season, stratum, taxon_code) %>%
  summarize(mean_cellcount = mean(mean_cellcount),
            mean_celldensity = mean(mean_celldensity),
            mean_biomass = mean(mean_biomass),
            nmonths = n(),
            ndays = sum(ndays)) %>%
  ungroup() %>%
  arrange(lake_id, year, factor(month_season, levels = names(palette_season)),
          factor(stratum, levels = c("Epilimnion", "Metalimnion")),
          taxon_code)

# Define functions to interpolate variable mean for missing month-seasons
interpolateMonthSeasonsTaxonWC <- function(data, names, method = "approx") {
  data_tmp <- data %>%
    arrange(year)
  
  lakeid_tmp <- unique(data_tmp$lake_id)
  monthseason_tmp <- unique(data_tmp$month_season)
  stratum_tmp <- unique(data_tmp$stratum)
  taxon_tmp <- unique(data_tmp$taxon_code)
  
  print(paste0(lakeid_tmp, "_", monthseason_tmp, "_", stratum_tmp, "_", taxon_tmp))
  
  year_range <- seq(data_tmp$year[1], tail(data_tmp$year, n = 1), 1)
  
  if (length(year_range) > 1) {
    if (method == "approx") {
      interpolated <- approx(data_tmp$year, data_tmp$mean, year_range) %>%
        as_tibble() %>%
        rename(year = x,
               mean = y) %>%
        filter(!year %in% unique(data_tmp$year)) %>%
        mutate(lake_id = lakeid_tmp,
               month_season = monthseason_tmp,
               stratum = stratum_tmp,
               taxon_code = taxon_tmp)
    } else if (method == "spline") {
      interpolated <- spline(data_tmp$year, data_tmp$mean,
                             xout = year_range,
                             method = "fmm") %>%
        as_tibble() %>%
        rename(year = x,
               mean = y) %>%
        filter(!year %in% unique(data_tmp$year)) %>%
        mutate(lake_id = lakeid_tmp,
               month_season = monthseason_tmp,
               stratum = stratum_tmp,
               taxon_code = taxon_tmp)
    }
    
    data_tmp %>%
      bind_rows(interpolated) %>%
      arrange(year) %>%
      mutate(interpolated = case_when(is.na(nmonths) ~ TRUE,
                                      TRUE ~ FALSE)) %>%
      mutate(nmonths = case_when(is.na(nmonths) ~ 0L,
                                 TRUE ~ nmonths),
             ndays = case_when(is.na(ndays) ~ 0L,
                               TRUE ~ ndays)) %>%
      return()
  }
}


#### Interpolate missing cell count data ####
# Interpolate missing cell count data by lake, year, season, stratum (based on interpolateMonthSeasonsAllWC)
cellcount_split <- phytoplankton_monthseason_means %>%
  dplyr::select(-mean_biomass) %>%
  rename(mean = mean_cellcount) %>%
  filter(month_season %in% c("Spring", "Summer", "Fall")) %>%  # Remove winter sampling dates
  group_split(lake_id, month_season, stratum, taxon_code)

cellcount_split %>%
  map(~mutate(., lakeid_monthseason_stratum_taxon = str_c(lake_id, "_", month_season, "_", stratum, "_", taxon_code))) %>%
  map(~pull(.,lakeid_monthseason_stratum_taxon)) %>%  # Pull out variable
  map(~as.character(.)) %>%  # Convert factor to character
  map(~unique(.)) -> names(cellcount_split)

cellcount_monthseason_interpolated <- list(data = cellcount_split, names = names(cellcount_split), method = "approx") %>%
  pmap(interpolateMonthSeasonsTaxonWC) %>%
  map_dfr(`[`, c("lake_id", "year", "month_season", "stratum", "taxon_code", "mean", "nmonths", "ndays", "interpolated"))

(cellcount_monthseason_interpolated_approx_plot <- cellcount_monthseason_interpolated %>%
    group_by(lake_id, year, month_season, stratum, nmonths, ndays, interpolated) %>%
    summarize(mean = sum(mean)) %>%
    ungroup() %>%
    mutate(monthseason_floor = case_when(month_season == "Spring" ~ ymd(paste0(year, "-04-01")),
                                         month_season == "Summer" ~ ymd(paste0(year, "-07-01")),
                                         month_season == "Fall" ~ ymd(paste0(year, "-10-01")))) %>%
    ggplot() +
    facet_col(lake_id~factor(stratum, levels = c("Epilimnion", "Metalimnion")),
              scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = monthseason_floor, y = mean,
                   colour = as.character(interpolated),
                   shape = factor(month_season, levels = c("Winter", "Spring", "Summer", "Fall")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(labels = comma) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "dodgerblue"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_shape_manual(values = c("Spring" = 50, "Summer" = 51, "Fall" = 52)) +
    labs(y = "No. cells",
         colour = "Interpolated",
         shape = "Season (month)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))

# Assess month-seasons with high incidences of data interpolation
cellcount_monthseason_interpolated %>%
  distinct(lake_id, year, month_season, stratum, nmonths, ndays, interpolated) %>%
  filter(interpolated) %>%
  group_by(lake_id, stratum, month_season, interpolated) %>%
  count() %>%
  ungroup() %>%
  arrange(-n)

# Summarize interpolated data by year (based on summarizeByYearInterpolatedWC)
# Remove all fall samples for all lakes
cellcount_byyear_interpolated <- cellcount_monthseason_interpolated %>%
  filter(month_season != "Fall") %>%
  group_by(lake_id, year, stratum, taxon_code) %>%
  summarize(mean = mean(mean),
            nseasons = n(),
            nmonths = sum(nmonths),
            ndays = sum(ndays)) %>%
  ungroup() %>%
  arrange(lake_id, year,
          factor(stratum, levels = c("Epilimnion", "Metalimnion")),
          taxon_code)

# Plot interpolated mean annual phytoplankton cell count
(cellcount_species_byyear_interpolated_plot <- cellcount_byyear_interpolated %>%
    group_by(lake_id, year, stratum) %>%
    summarize(mean = sum(mean)) %>%
    ungroup() %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    geom_point(aes(x = year_floor, y = mean,
                   colour = factor(stratum, levels = c("Epilimnion", "Metalimnion")))) +
    geom_line(aes(x = year_floor, y = mean,
                  colour = factor(stratum, levels = c("Epilimnion", "Metalimnion")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(labels = comma,
                       breaks = seq(0, sum(cellcount_byyear_interpolated$mean), 1000)) +
    scale_colour_manual(values = palette_stratum[which(names(palette_stratum) %in% unique(cellcount_byyear_interpolated$stratum))]) +
    scale_fill_manual(values = palette_manipulation) +
    labs(y = "No. cells",
         colour = "Stratum",
         fill = "Experiment") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_phytoplankton_species_cellcount_byyear_interpolated_plot.pdf", cellcount_species_byyear_interpolated_plot, width = 12, height = 10, units = "in", device = cairo_pdf)


#### Interpolate missing cell density data ####
# Interpolate missing cell density data by lake, year, season, stratum (based on interpolateMonthSeasonsAllWC)
celldensity_split <- phytoplankton_monthseason_means %>%
  dplyr::select(-mean_biomass) %>%
  rename(mean = mean_celldensity) %>%
  filter(month_season %in% c("Spring", "Summer", "Fall")) %>%  # Remove winter sampling dates
  group_split(lake_id, month_season, stratum, taxon_code)

celldensity_split %>%
  map(~mutate(., lakeid_monthseason_stratum_taxon = str_c(lake_id, "_", month_season, "_", stratum, "_", taxon_code))) %>%
  map(~pull(.,lakeid_monthseason_stratum_taxon)) %>%  # Pull out variable
  map(~as.character(.)) %>%  # Convert factor to character
  map(~unique(.)) -> names(celldensity_split)

celldensity_monthseason_interpolated <- list(data = celldensity_split, names = names(celldensity_split), method = "approx") %>%
  pmap(interpolateMonthSeasonsTaxonWC) %>%
  map_dfr(`[`, c("lake_id", "year", "month_season", "stratum", "taxon_code", "mean", "nmonths", "ndays", "interpolated"))

(celldensity_monthseason_interpolated_approx_plot <- celldensity_monthseason_interpolated %>%
    group_by(lake_id, year, month_season, stratum, nmonths, ndays, interpolated) %>%
    summarize(mean = sum(mean)) %>%
    ungroup() %>%
    mutate(monthseason_floor = case_when(month_season == "Spring" ~ ymd(paste0(year, "-04-01")),
                                         month_season == "Summer" ~ ymd(paste0(year, "-07-01")),
                                         month_season == "Fall" ~ ymd(paste0(year, "-10-01")))) %>%
    ggplot() +
    facet_col(lake_id~factor(stratum, levels = c("Epilimnion", "Metalimnion")),
              scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = monthseason_floor, y = mean,
                   colour = as.character(interpolated),
                   shape = factor(month_season, levels = c("Winter", "Spring", "Summer", "Fall")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(labels = comma) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "dodgerblue"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_shape_manual(values = c("Spring" = 50, "Summer" = 51, "Fall" = 52)) +
    labs(y = "Cell density (L-1)",
         colour = "Interpolated",
         shape = "Season (month)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))

# Assess month-seasons with high incidences of data interpolation
celldensity_monthseason_interpolated %>%
  distinct(lake_id, year, month_season, stratum, nmonths, ndays, interpolated) %>%
  filter(interpolated) %>%
  group_by(lake_id, stratum, month_season, interpolated) %>%
  count() %>%
  ungroup() %>%
  arrange(-n)

# Summarize interpolated data by year (based on summarizeByYearInterpolatedWC)
# Remove all fall samples for all lakes
celldensity_byyear_interpolated <- celldensity_monthseason_interpolated %>%
  filter(month_season != "Fall") %>%
  group_by(lake_id, year, stratum, taxon_code) %>%
  summarize(mean = mean(mean),
            nseasons = n(),
            nmonths = sum(nmonths),
            ndays = sum(ndays)) %>%
  ungroup() %>%
  arrange(lake_id, year,
          factor(stratum, levels = c("Epilimnion", "Metalimnion")),
          taxon_code)

# Plot interpolated mean annual phytoplankton cell density
(celldensity_species_byyear_interpolated_plot <- celldensity_byyear_interpolated %>%
    group_by(lake_id, year, stratum) %>%
    summarize(mean = sum(mean)) %>%
    ungroup() %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    geom_point(aes(x = year_floor, y = mean,
                   colour = factor(stratum, levels = c("Epilimnion", "Metalimnion")))) +
    geom_line(aes(x = year_floor, y = mean,
                  colour = factor(stratum, levels = c("Epilimnion", "Metalimnion")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(labels = comma,
                       breaks = seq(0, sum(celldensity_byyear_interpolated$mean), 20000000)) +
    scale_colour_manual(values = palette_stratum[which(names(palette_stratum) %in% unique(celldensity_byyear_interpolated$stratum))]) +
    scale_fill_manual(values = palette_manipulation) +
    labs(y = "Cell density (L-1)",
         colour = "Stratum",
         fill = "Experiment") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_phytoplankton_species_celldensity_byyear_interpolated_plot.pdf", celldensity_species_byyear_interpolated_plot, width = 12, height = 10, units = "in", device = cairo_pdf)


#### Interpolate missing biomass data ####
# Interpolate missing biomass data by lake, year, season, stratum (based on interpolateMonthSeasonsAllWC)
biomass_split <- phytoplankton_monthseason_means %>%
  dplyr::select(-mean_cellcount) %>%
  rename(mean = mean_biomass) %>%
  filter(month_season %in% c("Spring", "Summer", "Fall")) %>%  # Remove winter sampling dates
  group_split(lake_id, month_season, stratum, taxon_code)

biomass_split %>%
  map(~mutate(., lakeid_monthseason_stratum_taxon = str_c(lake_id, "_", month_season, "_", stratum, "_", taxon_code))) %>%
  map(~pull(.,lakeid_monthseason_stratum_taxon)) %>%  # Pull out variable
  map(~as.character(.)) %>%  # Convert factor to character
  map(~unique(.)) -> names(biomass_split)

biomass_monthseason_interpolated <- list(data = biomass_split, names = names(biomass_split), method = "approx") %>%
  pmap(interpolateMonthSeasonsTaxonWC) %>%
  map_dfr(`[`, c("lake_id", "year", "month_season", "stratum", "taxon_code", "mean", "nmonths", "ndays", "interpolated"))

(biomass_monthseason_interpolated_approx_plot <- biomass_monthseason_interpolated %>%
    group_by(lake_id, year, month_season, stratum, nmonths, ndays, interpolated) %>%
    summarize(mean = sum(mean)) %>%
    ungroup() %>%
    mutate(monthseason_floor = case_when(month_season == "Spring" ~ ymd(paste0(year, "-04-01")),
                                         month_season == "Summer" ~ ymd(paste0(year, "-07-01")),
                                         month_season == "Fall" ~ ymd(paste0(year, "-10-01")))) %>%
    ggplot() +
    facet_col(lake_id~factor(stratum, levels = c("Epilimnion", "Metalimnion")),
              scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = monthseason_floor, y = mean,
                   colour = as.character(interpolated),
                   shape = factor(month_season, levels = c("Winter", "Spring", "Summer", "Fall")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(labels = comma) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "dodgerblue"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_shape_manual(values = c("Spring" = 50, "Summer" = 51, "Fall" = 52)) +
    labs(y = "Biomass (\u00b5g/L)",
         colour = "Interpolated",
         shape = "Season (month)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))

# Assess month-seasons with high incidences of data interpolation
biomass_monthseason_interpolated %>%
  distinct(lake_id, year, month_season, stratum, nmonths, ndays, interpolated) %>%
  filter(interpolated) %>%
  group_by(lake_id, stratum, month_season, interpolated) %>%
  count() %>%
  ungroup() %>%
  arrange(-n)

# Summarize interpolated data by year (based on summarizeByYearInterpolatedWC)
# Remove all fall samples for all lakes
biomass_byyear_interpolated <- biomass_monthseason_interpolated %>%
  filter(month_season != "Fall") %>%
  group_by(lake_id, year, stratum, taxon_code) %>%
  summarize(mean = mean(mean),
            nseasons = n(),
            nmonths = sum(nmonths),
            ndays = sum(ndays)) %>%
  ungroup() %>%
  arrange(lake_id, year,
          factor(stratum, levels = c("Epilimnion", "Metalimnion")),
          taxon_code)

# Plot interpolated mean annual phytoplankton biomass
(biomass_species_byyear_interpolated_plot <- biomass_byyear_interpolated %>%
    group_by(lake_id, year, stratum) %>%
    summarize(mean = sum(mean)) %>%
    ungroup() %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    geom_point(aes(x = year_floor, y = mean,
                   colour = factor(stratum, levels = c("Epilimnion", "Metalimnion")))) +
    geom_line(aes(x = year_floor, y = mean,
                  colour = factor(stratum, levels = c("Epilimnion", "Metalimnion")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(labels = comma,
                       breaks = seq(0, max(biomass_byyear_interpolated$mean), 3000)) +
    scale_colour_manual(values = palette_stratum[which(names(palette_stratum) %in% unique(biomass_byyear_interpolated$stratum))]) +
    scale_fill_manual(values = palette_manipulation) +
    labs(y = "Biomass (\u00b5g/L)",
         colour = "Stratum",
         fill = "Experiment") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_phytoplankton_species_biomass_byyear_interpolated_plot.pdf", biomass_species_byyear_interpolated_plot, width = 12, height = 10, units = "in", device = cairo_pdf)


#### Area plots ####
# Order taxa by decreasing total relative sequence abundance
phytoplankton_tax_order <- phytoplankton_species %>%
  dplyr::select(domain_pr2, supergroup_pr2, division_pr2, subdivision_pr2, class_pr2, order_pr2, family_pr2, genus_pr2,
                cell_count, cell_density, biomass) %>%
  group_by(supergroup_pr2, division_pr2, subdivision_pr2) %>%
  summarize(biomass = sum(biomass)) %>%
  ungroup() %>%
  arrange(-biomass)

phytoplankton_supergroup_order <- phytoplankton_tax_order %>%
  distinct(supergroup_pr2) %>%
  mutate(supergroup_order = row_number())

phytoplankton_division_order <- phytoplankton_tax_order %>%
  distinct(division_pr2) %>%
  mutate(division_order = row_number())

phytoplankton_subdivision_order <- phytoplankton_tax_order %>%
  distinct(subdivision_pr2) %>%
  mutate(subdivision_order = row_number())

phytoplankton_tax_order <- phytoplankton_tax_order %>%
  left_join(phytoplankton_supergroup_order, by = "supergroup_pr2") %>%
  left_join(phytoplankton_division_order, by = "division_pr2") %>%
  left_join(phytoplankton_subdivision_order, by = "subdivision_pr2") %>%
  arrange(supergroup_order, division_order, subdivision_order)

# Plot interpolated mean annual phytoplankton cell count
(cellcount_species_byyear_interpolated_areaplot <- cellcount_byyear_interpolated %>%
    left_join(phytoplankton_taxonomy, by = "taxon_code") %>%
    group_by(lake_id, year, stratum, subdivision_pr2) %>%
    summarize(mean = sum(mean)) %>%
    ungroup() %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id+stratum~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    scale_fill_manual(values = palette_manipulation, name = "Experiment") +
    new_scale("fill") +
    geom_area(aes(x = year_floor, y = mean,
                  fill = factor(subdivision_pr2, levels = phytoplankton_tax_order$subdivision_pr2)),
              stat = "identity") +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(labels = comma,
                       breaks = seq(0, max(cellcount_byyear_interpolated$mean), 3000)) +
    scale_fill_manual(values = palette_subdivision, name = "Taxonomy") +
    labs(y = "No. cells",
         colour = "Stratum") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_phytoplankton_species_cellcount_byyear_interpolated_areaplot.pdf", cellcount_species_byyear_interpolated_areaplot, width = 12, height = 8, units = "in", device = cairo_pdf)

# Plot interpolated mean annual phytoplankton cell density
(celldensity_species_byyear_interpolated_areaplot <- celldensity_byyear_interpolated %>%
    left_join(phytoplankton_taxonomy, by = "taxon_code") %>%
    group_by(lake_id, year, stratum, subdivision_pr2) %>%
    summarize(mean = sum(mean)) %>%
    ungroup() %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id+stratum~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    scale_fill_manual(values = palette_manipulation, name = "Experiment") +
    new_scale("fill") +
    geom_area(aes(x = year_floor, y = mean,
                  fill = factor(subdivision_pr2, levels = phytoplankton_tax_order$subdivision_pr2)),
              stat = "identity") +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(labels = comma,
                       breaks = seq(0, max(celldensity_byyear_interpolated$mean), 20000000)) +
    scale_fill_manual(values = palette_subdivision, name = "Taxonomy") +
    labs(y = "Cell density (L-1)",
         colour = "Stratum") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_phytoplankton_species_celldensity_byyear_interpolated_areaplot.pdf", celldensity_species_byyear_interpolated_areaplot, width = 12, height = 8, units = "in", device = cairo_pdf)

# Plot interpolated mean annual phytoplankton biomass
(biomass_species_byyear_interpolated_areaplot <- biomass_byyear_interpolated %>%
    left_join(phytoplankton_taxonomy, by = "taxon_code") %>%
    group_by(lake_id, year, stratum, subdivision_pr2) %>%
    summarize(mean = sum(mean)) %>%
    ungroup() %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id+stratum~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    scale_fill_manual(values = palette_manipulation, name = "Experiment") +
    new_scale("fill") +
    geom_area(aes(x = year_floor, y = mean,
                  fill = factor(subdivision_pr2, levels = phytoplankton_tax_order$subdivision_pr2)),
              stat = "identity") +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(labels = comma,
                       breaks = seq(0, max(biomass_byyear_interpolated$mean), 3000)) +
    scale_fill_manual(values = palette_subdivision, name = "Taxonomy") +
    labs(y = "Biomass (\u00b5g/L)",
         colour = "Stratum") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_phytoplankton_species_biomass_byyear_interpolated_areaplot.pdf", biomass_species_byyear_interpolated_areaplot, width = 12, height = 8, units = "in", device = cairo_pdf)


#### Write data summary to file ####
# cellcount_byyear_interpolated %>%
#   dplyr::select(lake_id, year, stratum, taxon_code, mean) %>%
#   rename(cell_count = mean) %>%
#   left_join(celldensity_byyear_interpolated %>%
#               dplyr::select(lake_id, year, stratum, taxon_code, mean) %>%
#               rename(cell_density = mean), by = c("lake_id", "year", "stratum", "taxon_code")) %>%
#   left_join(biomass_byyear_interpolated %>%
#               dplyr::select(lake_id, year, stratum, taxon_code, mean) %>%
#               rename(biomass = mean), by = c("lake_id", "year", "stratum", "taxon_code")) %>%
#   filter(!(cell_count == 0 & cell_density == 0 & biomass == 0)) %>%
#   left_join(phytoplankton_taxonomy, by = "taxon_code") %>%
#   mutate(stratum = str_to_lower(stratum)) %>%
#   write_tsv("output/ela/ela_phytoplankton_species_byyear_interpolated.tsv", col_names = TRUE)
