# Summarized monitoring data

# setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
# setwd("~/Dropbox/projects/ela18s/")

# Load libraries
library(tidyverse)


#### Import and format data ####
# Import and format experimental manipulation data
ela_manipulations <- read_tsv("data/ela/ela_manipulations.txt", col_names = TRUE, comment = "#")

# Import and format climate monitoring data annual summaries
climate_eccc_byyear <- read_tsv("output/environmental/eccc_climate_byyear.tsv", col_names = TRUE) %>%
  crossing(tibble(lake_id = c("L224", "L226N", "L226S", "L227", "L373")))
climate_ela_byyear <- read_tsv("output/ela/ela_climate_byyear.tsv", col_names = TRUE) %>%
  crossing(tibble(lake_id = c("L224", "L226N", "L226S", "L227", "L373")))

# Import monitoring data annual summaries
watertemp_byyear_interpolated <- read_tsv("output/ela/ela_watertemp_strata_byyear_interpolated.tsv", col_names = TRUE)
sdd_byyear_interpolated <- read_tsv("output/ela/ela_sdd_byyear_interpolated.tsv", col_names = TRUE)
chemistry_byyear_interpolated <- read_tsv("output/ela/ela_chemistry_strata_byyear_interpolated.tsv", col_names = TRUE)

# Import phytoplankton species data
phytoplankton_species <- read_tsv("output/ela/ela_phytoplankton_species_byyear_interpolated.tsv", col_names = TRUE)


#### Summarize climate and limnological monitoring data ####
# Define function to extract environmental monitoring data for specified lake
extractEnvYear <- function(lakeid, env_data) {
  env_data %>%
    filter(lake_id == lakeid) %>%
    return()
}

# Define function to summarize environmental variables across sediment time windows
syncEnv <- function(lakeid, env_data, env_var) {
  metadata_lake <- metadata %>%
    filter(lake_id == lakeid)
  
  data_round <- metadata_lake %>%
    mutate(lowerpt_year_round = round(lowerpt_year),
           upperpt_year_round = round(upperpt_year)) %>%
    select(layer_order, lowerpt_year_round, upperpt_year_round) %>%
    mutate(year_diff_round = upperpt_year_round - lowerpt_year_round)
  
  year_range <- seq(round(tail(metadata_lake$lowerpt_year, n = 1)), 2018, 1)
  
  internal_years <- tibble(year = year_range) %>%
    arrange(-year) %>%
    filter(!year %in% data_round$lowerpt_year_round)
  
  data_sync <- data_round %>%
    mutate(year = lowerpt_year_round) %>%
    bind_rows(internal_years) %>%
    arrange(-year) %>%
    mutate(layer_sync = layer_order) %>%
    fill(layer_sync, .direction = "up")
  
  # Monitoring years correspond to which sediment layers?
  monitoring <- env_data %>%
    filter(lake_id == lakeid) %>%
    arrange(-year)
  
  data_sync_monitoring <- data_sync %>%
    left_join(monitoring, by = "year") 
  
  env_bylayer_mean <- data_sync_monitoring %>%
    select(layer_sync, as.name(paste(env_var))) %>%
    filter(!is.na(eval(as.name(paste(env_var))))) %>%
    group_by(layer_sync) %>%
    summarize({{env_var}} := mean(eval(as.name(paste(env_var)))),
              nyear = n()) %>%
    ungroup() %>%
    arrange(-layer_sync)
  
  return(env_bylayer_mean)
}
# syncEnv("L224", sdd_byyear_interpolated, "sdd")
# syncEnv("L224", climate_eccc_byyear, "meanairtemp")


#### Summarize phytoplankton species monitoring data ####
# Remove heterotrophs
phytoplankton_phototrophs <- phytoplankton_species %>%
  filter(trophic_mode %in% c("phototrophs", "phototrophs/mixotrophs", "mixotrophs"))

# Define function to sum annual phytoplankton biomass
sumPhytoplanktonAnnual <- function(phytoplankton, taxon) {
  monitoring <- phytoplankton %>%
    group_by(lake_id, year, stratum, eval(as.name(paste(taxon)))) %>%
    summarize(biomass = sum(biomass)) %>%
    ungroup() %>%
    rename(tax = `eval(as.name(paste(taxon)))`) %>%
    filter(!is.na(tax)) %>%
    arrange(lake_id, stratum, -year)

  return(monitoring)
}

# Define function to summarize phytoplankton monitoring across sediment time windows
syncPhytoplankton <- function(phytoplankton, taxon, lakeid, strat) {
  metadata_lake <- metadata %>%
    filter(lake_id == lakeid)
  
  data_round <- metadata_lake %>%
    mutate(lowerpt_year_round = round(lowerpt_year),
           upperpt_year_round = round(upperpt_year)) %>%
    select(layer_order, lowerpt_year_round, upperpt_year_round) %>%
    mutate(year_diff_round = upperpt_year_round - lowerpt_year_round)
  
  year_range <- seq(round(tail(metadata_lake$lowerpt_year, n = 1)), 2018, 1)
  
  internal_years <- tibble(year = year_range) %>%
    arrange(-year) %>%
    filter(!year %in% data_round$lowerpt_year_round)
  
  data_sync <- data_round %>%
    mutate(year = lowerpt_year_round) %>%
    bind_rows(internal_years) %>%
    arrange(-year) %>%
    mutate(layer_sync = layer_order) %>%
    fill(layer_sync, .direction = "up")
  
  # Monitoring years correspond to which sediment layers?
  monitoring <- phytoplankton %>%
    group_by(lake_id, year, stratum, eval(as.name(paste(taxon)))) %>%
    summarize(biomass = sum(biomass)) %>%
    ungroup() %>%
    rename(tax = `eval(as.name(paste(taxon)))`) %>%
    filter(lake_id == lakeid & stratum == strat) %>%
    filter(!is.na(tax))
  
  monitoring_wide <- monitoring %>%
    pivot_wider(names_from = tax, values_from = biomass, values_fill = 0) %>%
    arrange(-year)
  
  data_sync_monitoring <- data_sync %>%
    left_join(monitoring_wide, by = "year") %>%
    filter(!is.na(lake_id))  # Delete empty rows (layers/years without monitoring data)
  
  phytoplankton_bylayer_mean <- data_sync_monitoring %>%
    select(layer_sync, (unique(monitoring$tax))) %>%
    pivot_longer(!layer_sync, names_to = "tax", values_to = "biomass") %>%
    group_by(layer_sync, tax) %>%
    summarize(mean = mean(biomass),
              nyear = n()) %>%
    ungroup() %>%
    arrange(-layer_sync)
  
  layer_sync_nyear <- phytoplankton_bylayer_mean %>%
    distinct(layer_sync, nyear)
  
  phytoplankton_bylayer_mean_wide <- phytoplankton_bylayer_mean %>%
    select(-nyear) %>%
    pivot_wider(names_from = tax, values_from = mean, values_fill = 0) %>%
    arrange(-layer_sync) %>%
    replace(is.na(.), 0) %>%
    left_join(layer_sync_nyear, by = "layer_sync") %>%
    relocate(nyear, .after = "layer_sync")
  
  return(phytoplankton_bylayer_mean_wide)
}
# syncPhytoplankton(phytoplankton_phototrophs, "taxon_code", "L226N", "epilimnion")
# syncPhytoplankton(phytoplankton_phototrophs, "taxon_code", "L226S", "epilimnion")
# syncPhytoplankton(phytoplankton_phototrophs, "taxon_code", "L227", "epilimnion")
# syncPhytoplankton(phytoplankton_phototrophs, "taxon_code", "L224", "epilimnion")
# syncPhytoplankton(phytoplankton_phototrophs, "taxon_code", "L373", "epilimnion")
