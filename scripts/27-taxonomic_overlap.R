# Taxonomic overlap between datasets

setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
setwd("~/Dropbox/projects/ela18s/")
setwd("~/Library/CloudStorage/Dropbox-Personal/projects/ela18s/")

# Load libraries
library(tidyverse)
library(patchwork)

# Load palettes
source("scripts/00-palettes.R")

# Load sediment data
source("scripts/03e-sediment_data.R")

# Load monitoring data
source("scripts/11-ela_monitoring.R")


#### Import data ####
# Import LakePulse ELA filter ASV sequence table
load("data/lakepulse/lp2017-2019watercolumn18s_ela_seqtab_nochim.rda")

# Import LakePulse ELA filter ASV taxonomy table
load("data/lakepulse/lp2017-2019watercolumn18s_ela_taxonomy_pr2v5_tryrc.rda")


#### Format phytoplankton species monitoring data ####
# Format phytoplankton taxonomy
phytoplankton_taxonomy <- phytoplankton_species %>%
  distinct(taxon_code, supergroup_pr2, division_pr2, subdivision_pr2, class_pr2, order_pr2, family_pr2, genus_pr2) %>%
  arrange(supergroup_pr2, division_pr2, subdivision_pr2, class_pr2, order_pr2, family_pr2, genus_pr2)

# Sum total phytoplankton biomass by year
phytoplankton_phototrophs_annual_taxoncode <- sumPhytoplanktonAnnual(phytoplankton_phototrophs, "taxon_code")


#### Format LakePulse ELA filter ASV data ####
# Format LakePulse ELA filter ASV taxonomy
lakepulse_taxonomy <- taxonomy %>%
  as_tibble(rownames = "sequence") %>%
  mutate(across(c(domain, supergroup, division, subdivision, class, order, family, genus, species),
                ~str_remove_all(.x, ":.*"))) %>%
  mutate(across(c(supergroup, division, subdivision), ~str_remove_all(.x, "_.*")))
rm(taxonomy)

# Format LakePulse ELA filter ASV data
lakepulse_melt <- seqtab.nochim %>%
  as_tibble(rownames = "lakepulse_id") %>%
  pivot_longer(!lakepulse_id, names_to = "sequence", values_to = "nseqs") %>%
  filter(nseqs > 0) %>%
  mutate(lake_id = case_when(lakepulse_id == "06-311" ~ "L224",
                             lakepulse_id == "06-308" ~ "L373",
                             lakepulse_id == "06-304" ~ "L226S",
                             lakepulse_id == "06-307" ~ "L227")) %>%
  dplyr::select(-lakepulse_id) %>%
  relocate(lake_id, 1) %>%
  left_join(lakepulse_taxonomy, by = "sequence")

# Identify ASVs assigned to Metazoa
lakepulse_metazoa_asvs <- lakepulse_taxonomy %>%
  filter(subdivision == "Metazoa") %>%
  pull(sequence)

# Identify ASVs assigned to land plants
lakepulse_plant_asvs <- lakepulse_taxonomy %>%
  filter(class == "Embryophyceae") %>%
  pull(sequence)

# Identify ASVs not assigned at the domain level
lakepulse_unassigned_domain_asvs <- lakepulse_taxonomy %>%
  filter(is.na(domain)) %>%
  pull(sequence)

# Remove ASVs assigned to non-target taxa or not assigned at the domain level
lakepulse_melt <- lakepulse_melt %>%
  filter(!sequence %in% lakepulse_metazoa_asvs & !sequence %in% lakepulse_plant_asvs & !sequence %in% lakepulse_unassigned_domain_asvs)
length(unique(lakepulse_melt$sequence)) == length(lakepulse_taxonomy$sequence) - length(c(lakepulse_metazoa_asvs, lakepulse_plant_asvs, lakepulse_unassigned_domain_asvs))
# Should evaluate to TRUE

# Calculate relative abundance of sequences
lakepulse_samples_nseqs <- lakepulse_melt %>%
  group_by(lake_id) %>%
  summarize(nseqs_total = sum(nseqs))

lakepulse_melt <- lakepulse_melt %>%
  left_join(lakepulse_samples_nseqs, by = "lake_id") %>%
  mutate(relseqs = nseqs/nseqs_total) %>%
  select(-nseqs_total)
sum(lakepulse_melt$relseqs) == length(unique(lakepulse_melt$lake_id))  # Should evaluate to TRUE


#### Format GF/C filter data ####
gfc_melt <- asv_melt %>%
  filter(sample_type == "filter")
sum(gfc_melt$relseqs)  # Should evaluate to n of samples


#### Format assemblage data ####
# Calculate total biomass across epilimnion dataset
biomass_epilimnion <- phytoplankton_phototrophs_annual_taxoncode %>%
  filter(stratum == "epilimnion") %>%
  summarize(biomass = sum(biomass)) %>%
  pull(biomass)

# Assess phytoplankton monitoring diversity by group
phytoplankton_biomass <- phytoplankton_phototrophs_annual_taxoncode %>%
  filter(stratum == "epilimnion") %>%
  left_join(phytoplankton_taxonomy, by = c("tax" = "taxon_code")) %>%
  group_by(subdivision_pr2) %>%
  summarize(biomass = sum(biomass)) %>%
  ungroup() %>%
  mutate(pct_biomass = biomass/biomass_epilimnion * 100) %>%
  arrange(-pct_biomass) %>%
  mutate(dataset = "microscopy") %>%
  rename(subdivision = subdivision_pr2)

# Assess sediment ASV diversity by group
sum(sediment_melt$relseqs)  # Should evaluate to n of samples
sediment_seqs <- sediment_melt %>%
  group_by(subdivision) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  mutate(pct_seqs = nseqs/sum(sediment_melt$nseqs) * 100) %>%
  arrange(-pct_seqs) %>%
  mutate(dataset = "sediment")
sum(sediment_seqs$pct_seqs) == 100  # Should evaluate to TRUE

# Assess LakePulse filter ASV diversity by group
lakepulse_seqs <- lakepulse_melt %>%
  group_by(subdivision) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  mutate(pct_seqs = nseqs/sum(lakepulse_melt$nseqs) * 100) %>%
  arrange(-pct_seqs) %>%
  mutate(dataset = "lakepulse")
sum(lakepulse_seqs$pct_seqs) == 100  # Should evaluate to TRUE

# Assess ELA GF/C filter ASV diversity by group
gfc_seqs <- gfc_melt %>%
  group_by(subdivision) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  mutate(pct_seqs = nseqs/sum(gfc_melt$nseqs) * 100) %>%
  arrange(-pct_seqs) %>%
  mutate(dataset = "gfc")
sum(gfc_seqs$pct_seqs) == 100  # Should evaluate to TRUE

# Join taxonomic summaries across datasets
subdivision_dataset <- bind_rows(phytoplankton_biomass,
                                 sediment_seqs,
                                 lakepulse_seqs,
                                 gfc_seqs) %>%
  mutate(pct_dataset = case_when(!is.na(pct_biomass) ~ pct_biomass,
                                 !is.na(pct_seqs) ~ pct_seqs)) %>%
  select(subdivision, dataset, pct_dataset) %>%
  pivot_wider(names_from = dataset, values_from = pct_dataset, values_fill = 0) %>%
  pivot_longer(!subdivision, names_to = "dataset", values_to = "pct_dataset")


#### Visualize datasets by group ####
# Identify groups based on abundance across datasets
threshold_wc <- 1
threshold_sediment <- 1
subdivision_dataset_abundance <- subdivision_dataset %>%
  pivot_wider(names_from = dataset, values_from = pct_dataset, values_fill = 0) %>%
  select(subdivision, microscopy, gfc, lakepulse, sediment) %>%
  mutate(dataset_abundance = case_when((microscopy >= threshold_wc |
                                          gfc >= threshold_wc |
                                          lakepulse >= threshold_wc) &
                                         sediment < threshold_wc ~ "watercolumn",
                                       microscopy < threshold_sediment &
                                         gfc < threshold_sediment &
                                         lakepulse < threshold_sediment &
                                         sediment >= threshold_sediment ~ "sediment",
                                       (microscopy >= threshold_wc |
                                          gfc >= threshold_wc |
                                          lakepulse >= threshold_wc) &
                                         sediment >= threshold_wc ~ "watercolumn+sediment")) %>%
  filter(!is.na(dataset_abundance)) %>%
  mutate(wc_max = case_when(microscopy > gfc &
                              microscopy > lakepulse ~ microscopy,
                            gfc > microscopy &
                              gfc > lakepulse ~ gfc,
                            lakepulse > microscopy &
                              lakepulse > gfc ~ lakepulse)) %>%
  mutate(wc_sed_diff = abs(wc_max - sediment)) %>%
  arrange(wc_sed_diff)

# Plot percentage of dataset assigned to each group
(subdivision_dataset_plot <- subdivision_dataset %>%
    right_join(subdivision_dataset_abundance, by = "subdivision") %>%
    ggplot() +
    geom_bar(aes(x = fct_reorder(subdivision, desc(wc_sed_diff)),
                 y = pct_dataset,
                 fill = factor(dataset, levels = c("microscopy", "gfc", "lakepulse", "sediment"))),
             stat = "identity", position = "dodge") +
    # geom_text(aes(x = fct_reorder(subdivision, desc(wc_sed_diff)),
    #               y = pct_dataset,
    #               group = factor(dataset, levels = c("microscopy", "gfc", "lakepulse", "sediment")),
    #               label = round(pct_dataset, 4)), position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = c("microscopy" = "#07A0C3",
                                 "gfc" = "#73FBD3",
                                 "lakepulse" = "#F0C808",
                                 "sediment" = "#503D3F"), name = "Dataset") +
    scale_y_continuous(breaks = seq(0, 100, 1)) +
    labs(y = "% of dataset") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          axis.ticks = element_line(colour = "black"),
          axis.title.x = element_blank(),
          panel.grid = element_blank()))
#ggsave("figures/ela_taxonomic_overlap.pdf", subdivision_dataset_plot, width = 8, height = 4, units = "in", device = cairo_pdf)
