# ENA sample submission template

setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
setwd("~/Dropbox/projects/ela18s/")

# Load libraries
library(tidyverse)

# Load sediment data
source("scripts/03e-sediment_data.R")


#### Format ENA sample template ####
metadata_format <- asv_melt %>%
  filter(sample_type == "filter") %>%
  distinct(sample_id, lake_id, date, depth_unit) %>%
  mutate(sample_id = str_c("ELA-", str_replace_all(sample_id, "_", "-"))) %>%
  mutate(depth = as.numeric(str_remove(depth_unit, "m$")))

# Format ENA group sample template
sample_template <- metadata_format %>%
  mutate(tax_id = 449393,
         scientific_name = "freshwater metagenome",
         sample_alias = sample_id,
         sample_title = str_c("IISD-ELA ", lake_id, " ", depth_unit, " sampled on ", date),
         sample_description = "18S rRNA gene V7 region amplicons",
         isolation_source = "lake water",
         `collection date` = date,
         `geographic location (country and/or sea)` = "Canada",
         `geographic location (region and locality)` = "International Institute for Sustainable Development Experimental Lakes Area (IISD-ELA)",
         depth = depth,
         `environmental_sample` = "Yes") %>%
  select(tax_id,
         scientific_name,
         sample_alias,
         sample_title,
         sample_description,
         isolation_source,
         `collection date`,
         `geographic location (country and/or sea)`,
         `geographic location (region and locality)`,
         depth,
         `environmental_sample`)

# Write sample template to file
# ** Note: first three rows of TSV file to be added manually:
# Checklist	ERC000011	ENA default sample checklist
# tax_id	scientific_name	sample_alias	sample_title	sample_description	isolation_source	collection date	geographic location (country and/or sea)	geographic location (region and locality)	depth	environmental_sample
# #units	 	 	 	 	 	 	 	 
# sample_template %>%
#   write_tsv("~/Desktop/ela18s_sample_checklist-all_filters.tsv", col_names = FALSE)
