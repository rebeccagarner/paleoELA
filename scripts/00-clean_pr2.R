# Clean PR2 v. 5.0.0 taxonomy

setwd("C:/Users/Gandalf/Downloads/")
setwd("~/Desktop/trophic_function/")

# Load libraries
library(tidyverse)


#### Import and format data ####
# Import and format PR2 v. 5.0.0 fasta sequence headers
# pr2 <- read_delim("pr2_version_5.0.0_SSU_dada2_taxonomy.txt",
#                  col_names = c("domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species")) %>%
#   mutate(domain = str_remove(domain, "^>")) %>%
#   select(domain, supergroup, division, subdivision, class, order, family, genus, species) %>%
#   arrange(domain, supergroup, division, subdivision, class, order, family, genus, species) %>%
#   filter(grepl("^Eukaryota", domain))

# Write cleaned headers to file
# pr2 %>%
#   write_csv("pr2_version_5.0.0_SSU_dada2_taxonomy.csv", col_names = TRUE)

# Import and format cleaned headers
pr2_clean <- read_csv("pr2_version_5.0.0_SSU_dada2_taxonomy.csv", col_names = TRUE)
unique(pr2_clean$domain)


#### Summarize PR2 taxonomy ####
# Create non-redundant PR2 v. 5.0.0 taxonomy
pr2_nr <- pr2_clean %>%
  mutate(across(c(domain, supergroup, division, subdivision, class, order, family, genus, species),
                ~str_remove_all(.x, ":.*"))) %>%
  mutate(across(c(supergroup, division, subdivision), ~str_remove_all(.x, "_.*"))) %>%
  distinct(domain, supergroup, division, subdivision, class, order, family, genus, species) %>%
  arrange(domain, supergroup, division, subdivision, class, order, family, genus, species)
  
# Write non-redundant taxonomy to file
# pr2_nr %>%
#   write_csv("pr2_version_5.0.0_SSU_dada2_taxonomy_nr.csv", col_names = TRUE)
