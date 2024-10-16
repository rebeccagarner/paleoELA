# Sediment dating information

setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
setwd("~/Dropbox/projects/ela18s/")
setwd("~/Library/CloudStorage/Dropbox-Personal/projects/ela18s/")

# Load libraries
library(tidyverse)
library(janitor)


#### Import and format data ####
# Import and format dating data
dating <- read_csv("data/dating/ELA_InterDates_0.25.csv", col_names = TRUE) %>%
  clean_names() %>%
  rename(lakepulse_id = lake_id, midpt = midpoint_depth) %>%
  dplyr::select(lakepulse_id, midpt, year)


### Join and curate dating data ####
# Filter study lakes
dating <- dating %>%
  mutate(lake_id = case_when(lakepulse_id == "06-311" ~ "L224",
                             lakepulse_id == "06-308" ~ "L373",
                             lakepulse_id == "06-304S" ~ "L226S",
                             lakepulse_id == "06-304N" ~ "L226N",
                             lakepulse_id == "06-307" ~ "L227")) %>%
  filter(!is.na(lake_id)) %>%
  relocate(lake_id) %>%
  arrange(lake_id, midpt)

# Write dating information to file
# dating %>%
#   write_tsv("output/environmental/ela_gravitycore_dating.tsv", col_names = TRUE, na = "")
