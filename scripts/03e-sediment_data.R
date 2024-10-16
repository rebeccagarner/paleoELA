# Sediment data

# setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
# setwd("~/Dropbox/projects/ela18s/")

# Load libraries
library(tidyverse)

# Load palettes
source("scripts/00-palettes.R")


#### Import and format data ####
# Import melted ASV table
asv_melt <- read_tsv("output/melted/ela18s_microeuks_melt_function.tsv", col_names = TRUE,
                     col_types = cols(date = col_date()))

# Import sediment dating data
dating <- read_tsv("output/environmental/ela_gravitycore_dating.tsv", col_names = TRUE) %>%
  dplyr::select(-lakepulse_id)


#### Format sediment and dating data ####
# Filter sediment samples and join dating data
sediment_melt <- asv_melt %>%
  filter(sample_type == "sediment")

# Join lowerpt and upperpt dating data
dating_split <- dating %>%
  group_split(lake_id)

dating_split %>%
  map(~pull(.,lake_id)) %>%  # Pull out variable
  map(~as.character(.)) %>%  # Convert factor to character
  map(~unique(.)) -> names(dating_split)

datePts <- function(data) {
  data_all <- data %>%
    arrange(midpt)
  
  lakeid_tmp <- unique(data_all$lake_id)
  
  midpt_range <- seq(0.25, tail(data_all$midpt, n = 1), 0.25)
  
  interpolated <- approx(data_all$midpt, data_all$year, midpt_range) %>%
    as_tibble() %>%
    rename(midpt = x,
           year = y) %>%
    filter(!midpt %in% unique(data_all$midpt)) %>%
    mutate(lake_id = lakeid_tmp)
  
  data_all <- data_all %>%
    bind_rows(interpolated) %>%
    arrange(midpt)
  
  data %>%
    filter(midpt %in% unique(sediment_melt$midpt[which(sediment_melt$lake_id == unique(data$lake_id))])) %>%
    arrange(midpt) %>%
    rename(midpt_year = year) %>%
    mutate(upperpt = midpt - 0.25,
           lowerpt = midpt + 0.25) %>%
    relocate(midpt_year, .after = lowerpt) %>%
    left_join(data_all, by = c("lake_id", "upperpt" = "midpt")) %>%
    rename(upperpt_year = year) %>%
    mutate(upperpt_year = case_when(upperpt == 0 ~ 2018,
                                    TRUE ~ upperpt_year)) %>%
    left_join(data_all, by = c("lake_id", "lowerpt" = "midpt")) %>%
    rename(lowerpt_year = year) %>%
    arrange(-midpt) %>%
    mutate(layer_order = row_number()) %>%
    arrange(midpt) %>%
    return()
}

dating <- list(data = dating_split) %>%
  pmap(datePts) %>%
  map_dfr(`[`, c("lake_id", "midpt", "upperpt", "lowerpt",
                 "midpt_year", "upperpt_year", "lowerpt_year",
                 "layer_order"))

# Filter sediment samples and join dating data
sediment_melt <- sediment_melt %>%
  left_join(dating, by = c("lake_id", "midpt", "upperpt", "lowerpt")) %>%
  mutate(interval = paste0(format(upperpt, nsmall = 1), " \U2012 ", format(lowerpt, nsmall = 1))) %>%
  dplyr::select(-date, -depth_unit)

# Format taxonomy
taxonomy <- sediment_melt %>%
  distinct(asv_code, sequence, domain, supergroup, division, subdivision, class, order, family, genus, species)

# Format metadata
# Calculate time intervals (in years) captured by sediment layers
metadata <- sediment_melt %>%
  distinct(sample_id, lake_id, interval, midpt, upperpt, lowerpt,
           midpt_year, upperpt_year, lowerpt_year,
           layer_order) %>%
  mutate(year_diff = upperpt_year - lowerpt_year) %>%
  arrange(lake_id, midpt)

# Isolate unique samples
samples_unique <- sediment_melt %>%
  distinct(sample_id, lake_id, interval, upperpt)

# Define sample order
sample_order <- metadata %>%
  arrange(lake_id, midpt) %>%
  pull(sample_id)


#### Format taxonomic and trophic groups ####
# Extract phototroph ASVs
phototrophs_melt <- sediment_melt %>%
  filter(trophic_mode %in% c("phototrophs", "phototrophs/mixotrophs", "mixotrophs"))

# Extract diatom ASVs
diatoms_melt <- sediment_melt %>%
  filter(class == "Bacillariophyceae" | class == "Coscinodiscophyceae" | class == "Mediophyceae")

# Extract Chrysophyceae ASVs
chrysophyceae_melt <- sediment_melt %>%
  filter(class == "Chrysophyceae")

# Order taxa by decreasing total relative sequence abundance
tax_order <- sediment_melt %>%
  group_by(supergroup, division, subdivision) %>%
  summarize(relseqs = sum(relseqs)) %>%
  ungroup() %>%
  arrange(-relseqs) %>%
  filter(!is.na(subdivision))

supergroup_order <- tax_order %>%
  distinct(supergroup) %>%
  mutate(supergroup_order = row_number())

division_order <- tax_order %>%
  distinct(division) %>%
  mutate(division_order = row_number())

subdivision_order <- tax_order %>%
  distinct(subdivision) %>%
  mutate(subdivision_order = row_number())

tax_order <- tax_order %>%
  left_join(supergroup_order, by = "supergroup") %>%
  left_join(division_order, by = "division") %>%
  left_join(subdivision_order, by = "subdivision") %>%
  arrange(supergroup_order, division_order, subdivision_order)


#### Remove top sediment layers ####
# Define function to remove sediment layers from the last 10 years
removeTopSediment <- function(melt, lakeid) {
  melt %>%
    filter(upperpt_year < 2008) %>%
    filter(lake_id %in% lakeid) %>%
    return()
}


#### Format sediment ASV data ####
# Define function to format sediment layer by taxon matrix
formatLayerByTaxon <- function(melt, taxon) {
  # Recalculate relative abundance
  layer_nseqs <- melt %>%
    group_by(layer_order) %>%
    dplyr::summarize(nseqs_total = sum(nseqs)) %>%
    ungroup()
  
  # Format site (layer) by taxon table
  melt %>%
    dplyr::select(-relseqs) %>%
    left_join(layer_nseqs, by = "layer_order") %>%
    mutate(relseqs = nseqs/nseqs_total) %>%
    select(-nseqs_total) %>%
    select(layer_order, taxon, relseqs) %>%
    filter(!is.na(eval(as.name(paste(taxon))))) %>%
    group_by(layer_order, eval(as.name(paste(taxon)))) %>%
    dplyr::summarize(relseqs = sum(relseqs)) %>%
    pivot_wider(names_from = `eval(as.name(paste(taxon)))`, values_from = relseqs, values_fill = 0) %>%
    arrange(-layer_order) %>%
    return()
}
