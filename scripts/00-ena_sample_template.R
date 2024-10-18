# ENA sample submission template

# Load libraries
library(tidyverse)

# Load sediment data
source("scripts/03e-sediment_data.R")


#### Format ENA sample template ####
metadata_format <- metadata %>%
  mutate(upperpt_str = case_when(grepl(".5", upperpt) ~ as.character(upperpt),
                                 TRUE ~ str_c(upperpt, ".0")),
         lowerpt_str = case_when(grepl(".5", lowerpt) ~ as.character(lowerpt),
                                 TRUE ~ str_c(lowerpt, ".0"))) %>%
  mutate(sample_id = str_c(lake_id, "_", upperpt_str, "-", lowerpt_str, "cm"))

# Format ENA group sample template
sample_template <- metadata_format %>%
  mutate(tax_id = 556182,
         scientific_name = "freshwater sediment metagenome",
         sample_alias = sample_id,
         sample_title = str_c("IISD-ELA ", lake_id, " ", upperpt_str, "-", lowerpt_str, " cm"),
         sample_description = "18S rRNA gene V7 region amplicons",
         isolation_source = "lake sediment",
         `collection date` = "2018-08",
         `geographic location (country and/or sea)` = "Canada",
         `geographic location (region and locality)` = "International Institute for Sustainable Development Experimental Lakes Area (IISD-ELA)",
         depth = format(upperpt/100, nsmall = 3),
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
#   write_tsv("~/Desktop/ela18s_sample_checklist-all.tsv", col_names = FALSE)
