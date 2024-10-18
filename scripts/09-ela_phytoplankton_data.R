# Format ELA phytoplankton monitoring data

# Load libraries
library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)

# Load palettes
source("scripts/00-palettes.R")


#### Import and format data ####
# Define function for filtering study lakes up until date of core sampling
filterObs <- function(data) {
  data %>%
    filter(lake_id %in% c(224, 226, 227, 373)) %>%
    filter(date < "2018-09-01") %>%
    return()
}

# Define function to format lake IDs
formatLakeID <- function(data) {
  data_tmp <- data
  if ("lake_id" %in% names(data_tmp) & "sublocation" %in% names(data_tmp)) {
    data_tmp <- data_tmp %>%
      unite("lake_id", c(lake_id, sublocation), sep = "_") %>%
      mutate(lake_id = case_when(lake_id == "224_LA" ~ "L224",
                                 lake_id == "226_NEB" ~ "L226N",
                                 lake_id == "226_SWB" ~ "L226S",
                                 lake_id == "227_LA" ~ "L227",
                                 lake_id == "373_LA" ~ "L373"))
  } else if ("lake_id" %in% names(data_tmp) & !("sublocation" %in% names(data_tmp))) {
    data_tmp <- data_tmp %>%
      mutate(lake_id = str_c("L", lake_id))
  }
  return(data_tmp)
}

# Summarize lake maximum depths
lake_depths <- c("L226N" = 14.7,
                 "L226S" = 11.6,
                 "L227" = 10.0,
                 "L224" = 27.3,
                 "L373" = 21.2) %>%
  enframe("lake_id", "max_depth")

# Define function for sampling depth within lake maximum depth
filterDepths <- function(data) {
  data %>%
    left_join(lake_depths, by = "lake_id") %>%
    filter(!(start_depth > max_depth | end_depth > max_depth)) %>%
    dplyr::select(-max_depth) %>%
    return()
}

# Import and format phytoplankton biomass data until 2012
phytoplankton_biomass_til2012 <- read_xls("data/ela/phyto_daily_group_data_to_2012.xls") %>%
  dplyr::select(-Year, -`Station ID`) %>%
  rename(lake_id = `Location ID`,
         sublocation = `Sublocation ID`,
         start_depth = `Start Depth (m)`,
         date = `Sample Date`,
         stratum = Stratum,
         layer_volume = `Layer Volume (m³)`,
         end_depth = `End Depth (m)`,
         Cyanobacteria = `Cyanophyte Biomass (mg·m¯³)`,
         Chlorophytes = `Chlorophyte Biomass (mg·m¯³)`,
         Euglenophytes = `Euglenophyte Biomass (mg·m¯³)`,
         Chrysophytes = `Chrysophyte Biomass (mg·m¯³)`,
         Diatoms = `Diatom Biomass (mg·m¯³)`,
         Cryptophytes = `Cryptophyte Biomass (mg·m¯³)`,
         Dinoflagellates = `Dinophyte Biomass (mg·m¯³)`,
         `Total phytoplankton` = `Total Biomass Phyto (mg·m¯³)`,
         dataset_id = `Data Set ID`,
         update_date = `Update date`) %>%
  mutate(date = date(date),
         update_date = date(update_date)) %>%
  filterObs() %>%
  formatLakeID() %>%
  filterDepths()

# Import and format phytoplankton species count data until 2012
phytoplankton_species_til2012 <- read_xlsx("data/ela/L224_L226_L227_L373_phyto_species_counts_to_2012.xlsx") %>%
  rename(lake_id = `Location ID`, sublocation = `Sublocation ID`,
         date = `Start Date`) %>%
  dplyr::select(-Year, -`Start Time`, -`Station ID`) %>%
  mutate(date = date(date)) %>%
  filterObs() %>%
  formatLakeID() %>%
  rename(start_depth = `Start Depth (m)`,
         end_depth = `End Depth (m)`,
         stratum = Stratum,
         stratum_volume_m3 = `Stratum Volume (m³)`,
         sequence_n = `Sequence Number`,
         taxon_code = `Taxon Code`,
         old_cell_volume_um3 = `Old Cell Volume (µm³)`,
         cell_volume_um3 = `Cell Volume (µm³)`,
         correction_factor = `Correction Factor`,
         cell_count1 = `Number Cells-count 1`,
         cell_count2 = `Number Cells-count 2`,
         chloroplast_length_um = `Chloroplast Length (µm)`,
         chloroplast_width_um = `Chloroplast Width (µm)`,
         cell_density_per_l = `Cell·L¯¹`,
         biomass = `Biomass (µg·L¯¹)`,
         sample_type = `Sample Type`,
         collection_method = `Collection Method`,
         collection_authority = `Collection Authority`,
         field_reference = `Field Reference`,
         dataset_id = `Data Set ID`,
         update_date = `Update date`) %>%
  filterDepths() %>%
  filter(!(biomass == 0 | (cell_count1 == 0 & cell_count2 == 0)))

# Import and format phytoplankton biomass data 2013-2020
phytoplankton_biomass_20132020 <- read_xlsx("data/ela/IISD-ELA_Phyto_2013-20_Garner.xlsx",
                                            sheet = "Biomass summary",
                                            skip = 8,
                                            col_names = c("lake_id", "date", "stratum", "start_depth", "end_depth",
                                                          "Cyanobacteria", "Chlorophytes", "Euglenophytes",
                                                          "Chrysophytes", "Diatoms", "Cryptophytes",
                                                          "Dinoflagellates", "Total biomass")) %>%
  rename(`Total phytoplankton` = `Total biomass`) %>%
  mutate(date = date(date)) %>%
  filterObs() %>%
  formatLakeID() %>%
  filterDepths()

# Import and format phytoplankton species count data 2013-2020
# Remove entries for Cyanobacteria
phytoplankton_species_20132020 <- read_xlsx("data/ela/IISD-ELA_Phyto_2013-20_Garner.xlsx",
                                            sheet = "Species data",
                                            skip = 20,
                                            col_names = c("lake_id", "date",
                                                          "Strata", "Start depth (m)", "End depth (m)", 
                                                          "Species code", "Species name",
                                                          "Subsample volume (mls)", "Raw cell count",
                                                          "microscope factor", "length µ", "width µ",
                                                          "cell volume µ3", "density cells L-1", "biomass mg m-3", "X16")) %>%
  dplyr::select(-X16) %>%  # Removes column producing the warnings 
  mutate(date = date(date)) %>%
  filterObs() %>%
  mutate(tax_group = case_when(substr(as.character(`Species code`), 1, 1) == "1" ~ "Cyanophytes",
                               substr(as.character(`Species code`), 1, 1) == "2" ~ "Chlorophytes",
                               substr(as.character(`Species code`), 1, 1) == "3" ~ "Euglenophytes",
                               substr(as.character(`Species code`), 1, 1) == "4" ~ "Chrysophytes",
                               substr(as.character(`Species code`), 1, 1) == "5" ~ "Diatoms",
                               substr(as.character(`Species code`), 1, 1) == "6" ~ "Cryptophytes",
                               substr(as.character(`Species code`), 1, 1) == "7" ~ "Dinoflagellates")) %>%
  formatLakeID() %>%
  rename(stratum = Strata,
         start_depth = `Start depth (m)`,
         end_depth = `End depth (m)`,
         taxon_code = `Species code`,
         taxon = `Species name`,
         subsample_volume_ml = `Subsample volume (mls)`,
         cell_count_raw = `Raw cell count`,
         microscope_factor = `microscope factor`,
         length_um = `length µ`,
         width_um = `width µ`,
         cell_volume_um3 = `cell volume µ3`,
         cell_density_per_l = `density cells L-1`,
         biomass_mg_per_m3 = `biomass mg m-3`) %>%
  rename(biomass = biomass_mg_per_m3) %>%
  mutate(taxon_code = str_c("P", taxon_code)) %>%
  filterDepths() %>%
  filter(!(biomass == 0 | cell_count_raw == 0)) %>%
  filter(tax_group != "Cyanophytes")

# Import and format species codes
phytoplankton_species_codes <- read_xls("data/ela/phyto species codes.xls")
phytoplankton_species_codes[phytoplankton_species_codes == "."] <- NA
phytoplankton_species_codes <- phytoplankton_species_codes %>%
  rename(taxon_code = `Access Code`) %>%
  clean_names() %>%
  dplyr::select(where(~any(!is.na(.))), -common_name, -preferred) %>%
  mutate(species = case_when(species == "sp" |
                               species == "sp." |
                               species == "spp." |
                               is.na(species) ~ "species",
                             TRUE ~ species))
phytoplankton_species_codes$species[phytoplankton_species_codes$species == "species"] <- NA

# Import PR2 v. 5.0.0 taxonomy
pr2 <- read_tsv("data/taxonomy/pr2_version_5.0.0_SSU_dada2_taxonomy_nr.tsv", col_names = TRUE)


#### Join phytoplankton data across monitoring years ####
# Join and format phytoplankton biomass data ----
phytoplankton_biomass <- phytoplankton_biomass_til2012 %>%
  bind_rows(phytoplankton_biomass_20132020) %>%
  # mutate(total_phytoplankton = rowSums(across(c(Cyanobacteria,
  #                                               Chlorophytes,
  #                                               Euglenophytes,
  #                                               Chrysophytes,
  #                                               Diatoms,
  #                                               Cryptophytes,
  #                                               Dinoflagellates)),
  #                                      na.rm = TRUE)) %>%
  dplyr::select(-`Total phytoplankton`) %>%  # Max difference of 0.2 ug/L between original and recalculated total phytoplankton biomass
  mutate(stratum = case_when(grepl("^EPI", stratum, ignore.case = TRUE) ~ "Epilimnion",
                             grepl("^MET", stratum, ignore.case = TRUE) ~ "Metalimnion",
                             grepl("^HYP", stratum, ignore.case = TRUE) ~ "Hypolimnion")) %>%
  mutate(across(c(Cyanobacteria,
                  Chlorophytes,
                  Euglenophytes,
                  Chrysophytes,
                  Diatoms,
                  Cryptophytes,
                  Dinoflagellates), ~case_when(is.na(.) ~ 0,
                                                   TRUE ~ .))) %>%
  select(lake_id, date, stratum, start_depth, end_depth, Cyanobacteria:Dinoflagellates)

# Write phytoplankton biomass data to file
# phytoplankton_biomass %>%
#   arrange(lake_id, date, start_depth) %>%
#   write_tsv("output/ela/ela_monitoring_phytoplankton_biomass.tsv", col_names = TRUE)

# Join and format phytoplankton species data ----
phytoplankton_species <- phytoplankton_species_til2012 %>%
  bind_rows(phytoplankton_species_20132020) %>%
  mutate(stratum = case_when(grepl("^EPI", stratum, ignore.case = TRUE) ~ "Epilimnion",
                             grepl("^MET", stratum, ignore.case = TRUE) ~ "Metalimnion",
                             grepl("^HYP", stratum, ignore.case = TRUE) ~ "Hypolimnion"))

# Correct putative errors in data
phytoplankton_species <- phytoplankton_species %>%
  mutate(start_depth = case_when(lake_id == "L227" & date == "2017-06-05" & stratum == "Metalimnion" ~ 1,
                                 TRUE ~ start_depth),
         end_depth = case_when(lake_id == "L227" & date == "2017-06-05" & stratum == "Metalimnion" ~ 4,
                               TRUE ~ end_depth)) %>%
  mutate(stratum = case_when(lake_id == "L373" & date == "2007-01-09" & start_depth == 5 & end_depth == 5 ~ "Metalimnion",
                             TRUE ~ stratum))

# Summarize cell count data
phytoplankton_species <- phytoplankton_species %>%
  mutate(cell_count = case_when(!is.na(cell_count1) & !is.na(cell_count2) ~ ceiling((cell_count1 + cell_count2)/2),
                              !is.na(cell_count1) & is.na(cell_count2) ~ cell_count1,
                              !is.na(cell_count_raw) ~ cell_count_raw)) %>%
  dplyr::select(-c(cell_count1, cell_count2, cell_count_raw))


#### Resolve duplicate entries ####
# Remove duplicated entries
phytoplankton_species <- phytoplankton_species %>%
  dplyr::select(c(lake_id, date, start_depth, end_depth, stratum,
                  taxon_code,
                  old_cell_volume_um3, cell_volume_um3, correction_factor,
                  cell_count, cell_density_per_l, biomass)) %>%
  distinct()

# Collapse data for taxon codes with multiple entries per lake/stratum/date
phytoplankton_species <- phytoplankton_species %>%
  group_by(lake_id, date, stratum, start_depth, end_depth, taxon_code) %>%
  summarize(cell_count = sum(cell_count),
            cell_density = sum(cell_density_per_l),
            biomass = sum(biomass)) %>%
  ungroup()


#### Curate ELA taxonomy ####
# Join taxonomy
phytoplankton_species <- phytoplankton_species %>%
  left_join(phytoplankton_species_codes, by = "taxon_code") %>%
  dplyr::select(-description_code)

# Summarize phytoplankton species taxonomy
phytoplankton_taxonomy <- phytoplankton_species %>%
  distinct(taxon_code, phylum, class, order, family, genus, species, variety)

# Assess taxon codes with missing taxonomy
phytoplankton_taxonomy_nogenus <- phytoplankton_taxonomy %>%
  filter(is.na(genus))

# Correct missing taxonomy
phytoplankton_taxonomy_genus <- phytoplankton_taxonomy_nogenus %>%
  left_join(phytoplankton_species_20132020 %>%
              dplyr::select(taxon_code, taxon), by = "taxon_code") %>%
  mutate(taxon = str_remove(taxon, ".*\\(")) %>%
  mutate(taxon = str_remove(taxon, "\\).*")) %>%
  mutate(taxon = sub("^(\\S*\\s+\\S+).*", "\\1", taxon)) %>%
  mutate(taxon = case_when(taxon == "Chrysophyte cyst" ~ "Chrysophyte",
                           TRUE ~ taxon)) %>%
  mutate(class = case_when(taxon == "Chrysophyte" ~ "Chrysophyceae"),
         genus = case_when(str_count(taxon, " ") == 1 ~ str_remove(taxon, " .*")),
         species = case_when(str_count(taxon, " ") == 1 ~ str_remove(taxon, ".* "))) %>%
  mutate(phylum = case_when(class == "Chrysophyceae" ~ "Chrysophyta",
                            genus == "Monomastix" |
                              genus == "Selenastrum" |
                              genus == "Phacotus" |
                              genus == "Scenedesmus" ~ "Chlorophyta"),
         class = case_when(genus == "Monomastix" |
                             genus == "Selenastrum" |
                             genus == "Phacotus" |
                             genus == "Scenedesmus" ~ "Chlorophyceae",
                           TRUE ~ class)) %>%
  dplyr::select(-taxon)

phytoplankton_taxonomy_nr <- phytoplankton_taxonomy %>%
  filter(!taxon_code %in% unique(phytoplankton_taxonomy_nogenus$taxon_code)) %>%
  bind_rows(phytoplankton_taxonomy_genus) %>%
  distinct(taxon_code, phylum, class, order, family, genus, species) %>%
  arrange(phylum, class, order, family, genus, species) %>%
  mutate(tax_group = case_when(class == "Bacillariophyceae" ~ "Diatoms",
                               class == "Chloromonadophyceae" ~ "Raphidophytes",
                               class == "Chlorophyceae" ~ "Chlorophytes",
                               class == "Chrysophyceae" ~ "Chrysophytes",
                               class == "Cryptophyceae" ~ "Cryptophytes",
                               class == "Dinophyceae" ~ "Dinoflagellates",
                               class == "Euglenophyceae" ~ "Euglenophytes",
                               phylum == "Cyanophyta" ~ "Cyanobacteria")) %>%
  relocate(tax_group, .after = taxon_code)
rm(phytoplankton_taxonomy)

phytoplankton_taxonomy_nr <- phytoplankton_taxonomy_nr %>%
  mutate(species = case_when(grepl("sp\\.", species) ~ NA_character_,
                             species == "cysts" ~ NA_character_,
                             TRUE ~ species)) %>%
  mutate(genus = case_when(grepl("small|large|blue green|penate", genus, ignore.case = TRUE) ~ NA_character_,
                           genus == "Dinoflagellate" ~ NA_character_,
                           TRUE ~ genus))

# Update taxonomy
phytoplankton_species <- phytoplankton_species %>%
  dplyr::select(-c(phylum, class, order, family, genus, species, variety)) %>%
  left_join(phytoplankton_taxonomy_nr, by = "taxon_code")
rm(phytoplankton_taxonomy_nr)

# Remove remaining Cyanobacteria
phytoplankton_species <- phytoplankton_species %>%
  filter(phylum != "Cyanophyta")


#### Harmonize ELA taxonomy with PR2 taxonomy ####
# Summarize ELA taxonomy
phytoplankton_taxonomy_nr <- phytoplankton_species %>%
  distinct(taxon_code, tax_group, phylum, class, order, family, genus, species) %>%
  arrange(tax_group, phylum, class, order, family, genus, species)

# Assess ELA phytoplankton genera not found in PR2
phytoplankton_taxonomy_nr %>%
  filter(!genus %in% unique(pr2$genus)) %>%
  distinct(phylum, class, order, family, genus)

# Correct ELA phytoplankton taxonomy spelling and add taxonomic synonyms in PR2
phytoplankton_taxonomy_pr2 <- phytoplankton_taxonomy_nr %>%
  mutate(genus = case_when(genus == "Cormarium" ~ "Cosmarium",
                           genus == "Ankrya" ~ "Ankyra",
                           genus == "Kichnerella" ~ "Kirchneriella",
                           genus == "Nephrocytiun" ~ "Nephrocytium",
                           genus == "Pyramidomonas" ~ "Pyramimonas",
                           genus == "Chrysolkos" ~ "Chrysolykos",
                           genus == "Bicoeca" ~ "Bicosoeca",
                           genus == "Chrysostephanospaera" ~ "Chrysostephanosphaera",
                           genus == "Katablepharis" ~ "Kathablepharis",
                           genus == "Gonyostomun" ~ "Gonyostomum",
                           TRUE ~ genus),
         species = case_when(species == "granutatum" ~ "granulatum",
                             species == "judai" ~ "judayi",
                             TRUE ~ species)) %>%
  mutate(genus_synonym = case_when(genus == "Teilingia" ~ "Sphaerozosma",
                                   genus == "Chodatella" ~ "Lagerheimia",
                                   genus == "Chrysococcus" ~ "Chrysastrella",
                                   genus == "Crucigeniella" ~ "Crucigenia",
                                   genus == "Gloeobotrys" ~ "Chlorobotrys",
                                   genus == "Isthmochloron" ~ "Arthrodesmus",
                                   genus == "Pseudokephyrion" ~ "Kephyrion",
                                   TRUE ~ genus))

# Identify distinct PR2 genera
pr2_genera <- pr2 %>%
  dplyr::select(-species) %>%
  filter(!is.na(genus)) %>%
  distinct(domain, supergroup, division, subdivision, class, order, family, genus)

# Identify distinct PR2 families
pr2_families <- pr2 %>%
  dplyr::select(-c(genus, species)) %>%
  filter(!is.na(family)) %>%
  distinct(domain, supergroup, division, subdivision, class, order, family)

# Identify distinct PR2 orders
pr2_orders <- pr2 %>%
  dplyr::select(-c(family, genus, species)) %>%
  filter(!is.na(order)) %>%
  distinct(domain, supergroup, division, subdivision, class, order)

# Identify distinct PR2 classes
pr2_classes <- pr2 %>%
  dplyr::select(domain,supergroup, division, subdivision, class) %>%
  filter(!is.na(class)) %>%
  distinct(domain, supergroup, division, subdivision, class)

# Harmonize with PR2 taxonomy
phytoplankton_taxonomy_pr2$domain_pr2 <- NA
phytoplankton_taxonomy_pr2$supergroup_pr2 <- NA
phytoplankton_taxonomy_pr2$division_pr2 <- NA
phytoplankton_taxonomy_pr2$subdivision_pr2 <- NA
phytoplankton_taxonomy_pr2$class_pr2 <- NA
phytoplankton_taxonomy_pr2$order_pr2 <- NA
phytoplankton_taxonomy_pr2$family_pr2 <- NA
phytoplankton_taxonomy_pr2$genus_pr2 <- NA
for (i in 1:nrow(phytoplankton_taxonomy_pr2)) {
  genus_tmp <- phytoplankton_taxonomy_pr2$genus[i]
  genus_synonym_tmp <- phytoplankton_taxonomy_pr2$genus_synonym[i]
  family_tmp <- phytoplankton_taxonomy_pr2$family[i]
  order_tmp <- phytoplankton_taxonomy_pr2$order[i]
  class_tmp <- phytoplankton_taxonomy_pr2$class[i]
  
  if (genus_tmp %in% pr2_genera$genus) {
    phytoplankton_taxonomy_pr2$genus_pr2[i] <- genus_tmp
    phytoplankton_taxonomy_pr2$family_pr2[i] <- pr2_genera$family[which(pr2_genera$genus == genus_tmp)]
    phytoplankton_taxonomy_pr2$order_pr2[i] <- pr2_genera$order[which(pr2_genera$genus == genus_tmp)]
    phytoplankton_taxonomy_pr2$class_pr2[i] <- pr2_genera$class[which(pr2_genera$genus == genus_tmp)]
    phytoplankton_taxonomy_pr2$subdivision_pr2[i] <- pr2_genera$subdivision[which(pr2_genera$genus == genus_tmp)]
    phytoplankton_taxonomy_pr2$division_pr2[i] <- pr2_genera$division[which(pr2_genera$genus == genus_tmp)]
    phytoplankton_taxonomy_pr2$supergroup_pr2[i] <- pr2_genera$supergroup[which(pr2_genera$genus == genus_tmp)]
    phytoplankton_taxonomy_pr2$domain_pr2[i] <- pr2_genera$domain[which(pr2_genera$genus == genus_tmp)]
  } else if (genus_synonym_tmp %in% pr2_genera$genus) {
    phytoplankton_taxonomy_pr2$genus_pr2[i] <- genus_synonym_tmp
    phytoplankton_taxonomy_pr2$family_pr2[i] <- pr2_genera$family[which(pr2_genera$genus == genus_synonym_tmp)]
    phytoplankton_taxonomy_pr2$order_pr2[i] <- pr2_genera$order[which(pr2_genera$genus == genus_synonym_tmp)]
    phytoplankton_taxonomy_pr2$class_pr2[i] <- pr2_genera$class[which(pr2_genera$genus == genus_synonym_tmp)]
    phytoplankton_taxonomy_pr2$subdivision_pr2[i] <- pr2_genera$subdivision[which(pr2_genera$genus == genus_synonym_tmp)]
    phytoplankton_taxonomy_pr2$division_pr2[i] <- pr2_genera$division[which(pr2_genera$genus == genus_synonym_tmp)]
    phytoplankton_taxonomy_pr2$supergroup_pr2[i] <- pr2_genera$supergroup[which(pr2_genera$genus == genus_synonym_tmp)]
    phytoplankton_taxonomy_pr2$domain_pr2[i] <- pr2_genera$domain[which(pr2_genera$genus == genus_synonym_tmp)]
  } else if (is.na(family_tmp) & is.na(order_tmp) & class_tmp %in% pr2_classes$class) {
    phytoplankton_taxonomy_pr2$class_pr2[i] <- class_tmp
    phytoplankton_taxonomy_pr2$subdivision_pr2[i] <- pr2_classes$subdivision[which(pr2_classes$class == class_tmp)]
    phytoplankton_taxonomy_pr2$division_pr2[i] <- pr2_classes$division[which(pr2_classes$class == class_tmp)]
    phytoplankton_taxonomy_pr2$supergroup_pr2[i] <- pr2_classes$supergroup[which(pr2_classes$class == class_tmp)]
    phytoplankton_taxonomy_pr2$domain_pr2[i] <- pr2_classes$domain[which(pr2_classes$class == class_tmp)]
  }
}

# Harmonize remaining taxonomy
phytoplankton_taxonomy_pr2 %>%
  filter(is.na(genus_pr2) & is.na(class_pr2)) %>%
  distinct(genus, .keep_all = TRUE)

phytoplankton_taxonomy_pr2 <- phytoplankton_taxonomy_pr2 %>%
  mutate(family_pr2 = case_when(genus == "Nephrocytium" ~ "Chlorellales_X",
                                genus == "Sphaerocystis" ~ "Chlamydomonadales_X",
                                genus == "Scourfieldia" ~ "Pedinophyceae_XX",
                                genus == "Bitrichia" ~ "Ochromonadaceae",
                                genus == "Chrysolykos" ~ "Ochromonadaceae",
                                genus == "Spiniferomonas" ~ "Synuraceae",
                                genus == "Desmarella" ~ "Salpingoecidae",
                                genus == "Chrysostephanosphaera" ~ "Chrysamoebaceae",
                                genus == "Rhizochrysis" ~ "Chrysamoebaceae",
                                genus == "Stelexomonas" ~ "Salpingoecidae",
                                genus == "Cryptaulax" ~ "Cryptomonadales_X",
                                genus == "Glenodinium" ~ "Peridiniales_X",
                                is.na(genus) & family == "Peridiniaceae" ~ "Peridiniaceae",
                                TRUE ~ family_pr2),
         order_pr2 = case_when(genus == "Nephrocytium" ~ "Chlorellales",
                               genus == "Sphaerocystis" ~ "Chlamydomonadales",
                               genus == "Scourfieldia" ~ "Pedinophyceae_X",
                               genus == "Desmarella" ~ "Craspedida",
                               genus == "Stelexomonas" ~ "Craspedida",
                               genus == "Cryptaulax" ~ "Cryptomonadales",
                               genus == "Glenodinium" ~ "Peridiniales",
                               TRUE ~ order_pr2),
         class_pr2 = case_when(genus == "Scourfieldia" ~ "Pedinophyceae",
                               TRUE ~ class_pr2))

for (i in 1:nrow(phytoplankton_taxonomy_pr2)) {
  if (is.na(phytoplankton_taxonomy_pr2$genus_pr2[i])) {
    if (phytoplankton_taxonomy_pr2$family_pr2[i] %in% pr2_families$family) {
      phytoplankton_taxonomy_pr2$genus_pr2[i] <- phytoplankton_taxonomy_pr2$genus[i]
      phytoplankton_taxonomy_pr2$order_pr2[i] <- pr2_families$order[which(pr2_families$family == phytoplankton_taxonomy_pr2$family_pr2[i])]
      phytoplankton_taxonomy_pr2$class_pr2[i] <- pr2_families$class[which(pr2_families$family == phytoplankton_taxonomy_pr2$family_pr2[i])]
      phytoplankton_taxonomy_pr2$subdivision_pr2[i] <- pr2_families$subdivision[which(pr2_families$family == phytoplankton_taxonomy_pr2$family_pr2[i])]
      phytoplankton_taxonomy_pr2$division_pr2[i] <- pr2_families$division[which(pr2_families$family == phytoplankton_taxonomy_pr2$family_pr2[i])]
      phytoplankton_taxonomy_pr2$supergroup_pr2[i] <- pr2_families$supergroup[which(pr2_families$family == phytoplankton_taxonomy_pr2$family_pr2[i])]
      phytoplankton_taxonomy_pr2$domain_pr2[i] <- pr2_families$domain[which(pr2_families$family == phytoplankton_taxonomy_pr2$family_pr2[i])]
    } else if (!(phytoplankton_taxonomy_pr2$family_pr2[i] %in% pr2_families$family) & phytoplankton_taxonomy_pr2$order_pr2[i] %in% pr2_orders$order) {
      phytoplankton_taxonomy_pr2$genus_pr2[i] <- phytoplankton_taxonomy_pr2$genus[i]
      phytoplankton_taxonomy_pr2$class_pr2[i] <- pr2_orders$class[which(pr2_orders$order == phytoplankton_taxonomy_pr2$order_pr2[i])]
      phytoplankton_taxonomy_pr2$subdivision_pr2[i] <- pr2_orders$subdivision[which(pr2_orders$order == phytoplankton_taxonomy_pr2$order_pr2[i])]
      phytoplankton_taxonomy_pr2$division_pr2[i] <- pr2_orders$division[which(pr2_orders$order == phytoplankton_taxonomy_pr2$order_pr2[i])]
      phytoplankton_taxonomy_pr2$supergroup_pr2[i] <- pr2_orders$supergroup[which(pr2_orders$order == phytoplankton_taxonomy_pr2$order_pr2[i])]
      phytoplankton_taxonomy_pr2$domain_pr2[i] <- pr2_orders$domain[which(pr2_orders$order == phytoplankton_taxonomy_pr2$order_pr2[i])]
    }
  }
}

# Join harmonized PR2 taxonomy to ELA phytoplankton data
phytoplankton_species <- phytoplankton_species %>%
  dplyr::select(-c(tax_group, phylum, class, order, family, genus, species)) %>%
  left_join(phytoplankton_taxonomy_pr2, by = "taxon_code")


#### Link taxonomy and function ####
# Create non-redundant taxonomy table
taxonomy_nr <- phytoplankton_species %>%
  distinct(across(ends_with("_pr2"))) %>%
  rename_with(~str_remove(.x, "_pr2")) %>%
  arrange(supergroup, division, subdivision, class, order, family, genus)

# Link taxonomy and trophic function summarized from Adl et al. 2019 "Revisions to
# the Classification, Nomenclature, and Diversity of Eukaryotes"
taxonomy_function <- taxonomy_nr %>%
  mutate(trophic_mode = case_when(
    # Amoebozoa
    #species == "Difflugia_bacillariarum" ~ "mixotrophs",  # To verify
    supergroup == "Amoebozoa" ~ "heterotrophs",
    
    # Apusomonada
    division == "Apusomonada" ~ "heterotrophs",
    
    # Breviatea
    division == "Breviatea" ~ "heterotrophs",
    
    # Opisthokonta
    genus == "Archaeorhizomyces" ~ "heterotrophs",
    genus == "Kionochaeta" ~ "saprotrophs",
    genus == "Teratosphaeria" ~ "parasites",
    order == "Gromochytriales" ~ "parasites",
    class == "Aphelidiomycota" ~ "parasites",
    class == "Blastocladiomycota" ~ "heterotrophs",
    class == "Chytridiomycota" ~ "saprotrophs/parasites",
    class == "Opisthosporidia" ~ "parasites",
    class == "Rozellomycota" ~ "parasites",
    subdivision == "Choanoflagellata" ~ "heterotrophs",
    division == "Opisthokonta" ~ "heterotrophs",
    
    supergroup == "Obazoa" ~ "heterotrophs",
    
    # Archaeplastida
    genus == "Helicosporidium" ~ "parasites",
    subdivision == "Chlorophyta" ~ "phototrophs",
    subdivision == "Streptophyta" ~ "phototrophs",
    supergroup == "Archaeplastida" ~ "phototrophs",
    
    # Cryptista
    genus == "Chilomonas" ~ "heterotrophs",  # ELA monitoring
    genus == "Chroomonas" ~ "phototrophs",  # ELA monitoring
    genus == "Cryptaulax" ~ "heterotrophs",  # ELA monitoring
    genus == "Cryptomonas" ~ "mixotrophs",
    genus == "Geminigera" ~ "mixotrophs",
    genus == "Goniomonas" ~ "heterotrophs",
    genus == "Komma" ~ "mixotrophs",
    genus == "Plagioselmis" ~ "mixotrophs",
    genus == "Rhodomonas" ~ "phototrophs/mixotrophs",  # ELA monitoring
    division == "Kathablepharidacea" ~ "heterotrophs",
    
    # Haptista
    genus == "Diacronema"~ "mixotrophs",
    genus == "Chrysochromulina" ~ "mixotrophs",
    genus == "Prymnesium" ~ "mixotrophs",
    division == "Centroplasthelida" ~ "heterotrophs",
    
    # Alveolata
    genus == "Alphamonas" ~ "heterotrophs",
    genus == "Amphidinium" ~ "phototrophs/mixotrophs/heterotrophs",  # ELA monitoring
    genus == "Anurofeca" ~ "parasites",
    genus == "Apocalathium" ~ "phototrophs",
    genus == "Asulcocephalium" ~ "phototrophs",
    genus == "Biecheleria" ~ "heterotrophs/mixotrophs",
    genus == "Blastodinium" ~ "mixotrophs",
    genus == "Borghiella" ~ "phototrophs",
    genus == "Ceratium" ~ "phototrophs/mixotrophs/heterotrophs",
    genus == "Chromera" ~ "phototrophs",
    genus == "Chytriodinium" ~ "parasites",
    genus == "Cladocopium" ~ "phototrophs",
    genus == "Gymnodinium" ~ "phototrophs/mixotrophs",
    genus == "Gyrodinium" ~ "heterotrophs",
    genus == "Ichthyophthirius" ~ "parasites",
    genus == "Karenia" ~ "heterotrophs/mixotrophs",
    genus == "Nusuttodinium" ~ "heterotrophs/mixotrophs",
    genus == "Paramoebidium" ~ "parasites",
    genus == "Peridinium" ~ "phototrophs",
    genus == "Tetrahymena" ~ "parasites",
    family == "Bursariomorphida" ~ "heterotrophs",
    family == "Symbiodiniaceae" ~ "phototrophs",
    order == "Peridiniales" ~ "phototrophs/mixotrophs/heterotrophs",
    class == "Phyllopharyngea" ~ "parasites",
    subdivision == "Apicomplexa" ~ "parasites",
    subdivision == "Ciliophora" ~ "heterotrophs",
    subdivision == "Perkinsea" ~ "parasites",
    
    # Stramenopiles
    genus == "Acrispumella" ~ "heterotrophs",
    genus == "Dinobryon" ~ "mixotrophs",
    genus == "Epipyxis" ~ "phototrophs",
    genus == "Hydrurus" ~ "phototrophs",
    order == "Chromulinales" ~ "phototrophs",
    order == "Hyphochytriales" ~ "saprotrophs",
    order == "Paraphysomonadales" ~ "heterotrophs",
    order == "Synurales" ~ "phototrophs",
    class == "Bacillariophyceae" ~ "phototrophs",
    class == "Bolidophyceae" ~ "phototrophs",
    class == "Chrysophyceae" ~ "phototrophs",
    class == "Coscinodiscophyceae" ~ "phototrophs",
    class == "Dictyochophyceae" ~ "phototrophs",
    class == "Eustigmatophyceae" ~ "phototrophs",
    class == "Mediophyceae" ~ "phototrophs",
    class == "Peronosporomycetes" ~ "parasites",
    class == "Pirsoniales" ~ "parasites",
    class == "Phaeothamniophyceae" ~ "phototrophs",
    class == "Raphidophyceae" ~ "phototrophs",
    class == "Xanthophyceae" ~ "phototrophs",
    subdivision == "Bigyra" ~ "heterotrophs",
    
    # Rhizaria
    genus == "Aquavolon" ~ "heterotrophs",
    order == "Cercomonadida" ~ "heterotrophs",
    order == "Euglyphida" ~ "heterotrophs",
    order == "Glissomonadida" ~ "heterotrophs",
    order == "Novel-clade-10" ~ "heterotrophs",
    order == "Pansomonadida" ~ "heterotrophs",
    order == "Paracercomonadida" ~ "heterotrophs",
    order == "Plasmodiophorida" ~ "parasites",
    order == "Thaumatomonadida" ~ "heterotrophs",
    order == "Vampyrellida" ~ "heterotrophs",
    
    # Telonemia
    division == "Telonemia" ~ "heterotrophs",
    
    # Excavata
    genus == "Hexamita" ~ "parasites",
    order == "Aphagea" ~ "heterotrophs",
    order == "Euglenophyceae" ~ "phototrophs",
    order == "Hibberdiales" ~ "phototrophs",
    order == "Trypanosomatida" ~ "parasites",
    class == "Heterolobosea" ~ "heterotrophs",
    class == "Kinetoplastea" ~ "heterotrophs",
    division == "Metamonada" ~ "heterotrophs",
    
    # CRuMs
    division == "Collodictyonidae" ~ "heterotrophs",  # ELA monitoring only
    division == "Rigifilida" ~ "heterotrophs",
    
    # Other
    division == "Ancyromonadida" ~ "heterotrophs"
  ))

# Join trophic function information with phytoplankton data
phytoplankton_species <- phytoplankton_species %>%
  left_join(taxonomy_function, by = c("supergroup_pr2" = "supergroup",
                                      "division_pr2" = "division",
                                      "subdivision_pr2" = "subdivision",
                                      "class_pr2" = "class",
                                      "order_pr2" = "order",
                                      "family_pr2" = "family",
                                      "genus_pr2" = "genus"))

# Write phytoplankton species biomass data to file
# phytoplankton_species %>%
#   write_tsv("output/ela/ela_monitoring_phytoplankton_species.tsv", col_names = TRUE)


#### Compare total species-level and group-level biomass ####
# Calculate total biomass across group- and species-level datasets
biomass_total_phytoplankton <- phytoplankton_biomass %>%
  mutate(eukaryotes = Chlorophytes + Euglenophytes + Chrysophytes + Diatoms + Cryptophytes + Dinoflagellates) %>%
  dplyr::select(lake_id, date, stratum, eukaryotes) %>%
  arrange(lake_id, date, stratum)

biomass_species <- phytoplankton_species %>%
  group_by(lake_id, date, stratum) %>%
  summarize(biomass_species = sum(biomass)) %>%
  ungroup() %>%
  arrange(lake_id, date, stratum)

# Assess entries in species data that do not appear in group data
# There are no lake/date/stratum entries in the group data that do not also appear
# in the species-level data
biomass_species %>%
  left_join(biomass_total_phytoplankton, by = c("lake_id", "date", "stratum")) %>%
  filter(is.na(eukaryotes)) %>%
  filter(stratum %in% c("Epilimnion", "Metalimnion")) %>%
  left_join(phytoplankton_species, by = c("lake_id", "date", "stratum"))

# Assess entries that appear in both group- and species-level data
biomass_total_phytoplankton %>%
  left_join(biomass_species, by = c("lake_id", "date", "stratum")) %>%
  filter(!is.na(eukaryotes) | !is.na(biomass_species))

# Verify agreement between phylum- and species-level biomass totals
biomass_total_phytoplankton %>%
  inner_join(biomass_species, by = c("lake_id", "date", "stratum")) %>%
  mutate(biomass_diff = abs(eukaryotes - biomass_species)) %>%
  #filter(biomass_diff > 1) %>%
  arrange(-biomass_diff) %>%
  mutate(`group > spp` = case_when(eukaryotes > biomass_species ~ TRUE,
                                   TRUE ~ FALSE)) %>%
  ggplot() +
  geom_histogram(aes(x = biomass_diff, y = ..count.., fill = `group > spp`), binwidth = 100) +
  theme_bw()
