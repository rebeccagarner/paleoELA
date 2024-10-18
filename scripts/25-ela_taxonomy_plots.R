# Phytoplankton taxonomic diversity and biomass

# Load libraries
library(tidyverse)
library(ggforce)
library(scales)
library(ggnewscale)
library(patchwork)

# Load palettes
source("scripts/00-palettes.R")

# Load sediment data
source("scripts/03e-sediment_data.R")

# Load environmental data
source("scripts/08-ela_env_data.R")

# Load monitoring data
source("scripts/11-ela_monitoring.R")


#### Format phytoplankton species monitoring data ####
# Format phytoplankton taxonomy
phytoplankton_taxonomy <- phytoplankton_species %>%
  distinct(taxon_code, supergroup_pr2, division_pr2, subdivision_pr2, class_pr2, order_pr2, family_pr2, genus_pr2) %>%
  arrange(supergroup_pr2, division_pr2, subdivision_pr2, class_pr2, order_pr2, family_pr2, genus_pr2)

# Sum total phytoplankton biomass by year
phytoplankton_phototrophs_annual_taxoncode <- sumPhytoplanktonAnnual(phytoplankton_phototrophs, "taxon_code")


#### Phytoplankton phototroph diversity by year ####
# Calculate total biomass across epilimnion dataset
biomass_epilimnion <- phytoplankton_phototrophs_annual_taxoncode %>%
  filter(stratum == "epilimnion") %>%
  summarize(biomass = sum(biomass)) %>%
  pull(biomass)

# Assess phytoplankton monitoring diversity by group
phytoplankton_biomass <- phytoplankton_phototrophs_annual_taxoncode %>%
  filter(stratum == "epilimnion") %>%
  left_join(phytoplankton_taxonomy, by = c("tax" = "taxon_code")) %>%
  group_by(class_pr2) %>%
  summarize(biomass = sum(biomass)) %>%
  ungroup() %>%
  mutate(pct_biomass = biomass/biomass_epilimnion * 100) %>%
  arrange(-pct_biomass) %>%
  rename_with(~str_remove(.x, "_pr2$"), everything()) %>%
  select(-biomass)

# Sediment phototrophs
phototrophs_seqs <- phototrophs_melt %>%
  group_by(class) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  mutate(pct_seqs = nseqs/sum(phototrophs_melt$nseqs) * 100) %>%
  arrange(-pct_seqs) %>%
  select(-nseqs)

class_datasets <- phytoplankton_biomass %>%
  full_join(phototrophs_seqs,
            by = "class") %>%
  mutate(pct_biomass = case_when(is.na(pct_biomass) ~ 0,
                                 TRUE ~ pct_biomass),
         pct_seqs = case_when(is.na(pct_seqs) ~ 0,
                              TRUE ~ pct_seqs)) %>%
  mutate(pct_datasets = pct_biomass + pct_seqs) %>%
  mutate(pct_datasets_all = pct_datasets/(sum(phytoplankton_biomass$pct_biomass) + sum(phototrophs_seqs$pct_seqs)) * 100)

phytoplankton_taxonomy_class <- phytoplankton_taxonomy %>%
  rename_with(~str_remove(.x, "_pr2$"), everything()) %>%
  distinct(supergroup, division, subdivision, class)
phototrophs_taxonomy_class <- phototrophs_melt %>%
  distinct(supergroup, division, subdivision, class)
taxonomy_all <- bind_rows(phytoplankton_taxonomy_class,
                          phototrophs_taxonomy_class) %>%
  distinct(supergroup, division, subdivision, class)

# Order taxa by decreasing total dataset(s) percentage
phototrophs_tax_order <- class_datasets %>%
  left_join(taxonomy_all, by = "class")

phototrophs_supergroup_order <- phototrophs_tax_order %>%
  distinct(supergroup) %>%
  mutate(supergroup_order = row_number())

phototrophs_division_order <- phototrophs_tax_order %>%
  distinct(division) %>%
  mutate(division_order = row_number())

phototrophs_subdivision_order <- phototrophs_tax_order %>%
  distinct(subdivision) %>%
  mutate(subdivision_order = row_number())

phototrophs_class_order <- phototrophs_tax_order %>%
  distinct(class) %>%
  mutate(class_order = row_number())

phototrophs_tax_order <- phototrophs_tax_order %>%
  left_join(phototrophs_supergroup_order, by = "supergroup") %>%
  left_join(phototrophs_division_order, by = "division") %>%
  left_join(phototrophs_subdivision_order, by = "subdivision") %>%
  left_join(phototrophs_class_order, by = "class") %>%
  arrange(supergroup_order, division_order, subdivision_order, class_order)

# Plot interpolated mean annual epilimnion phytoplankton biomass
biomass_epilimnion <- phytoplankton_phototrophs_annual_taxoncode %>%
  filter(stratum == "epilimnion") %>%
  group_by(lake_id, year) %>%
  summarize(biomass = sum(biomass)) %>%
  ungroup() %>%
  mutate(year_floor = ymd(paste0(year, "-01-01")))

(biomass_phototrophs_byyear_epilimnion_interpolated_areaplot <- phytoplankton_phototrophs_annual_taxoncode %>%
    filter(stratum == "epilimnion") %>%
    pivot_wider(names_from = tax, values_from = biomass, values_fill = 0) %>%
    pivot_longer(!c(lake_id, year, stratum), names_to = "tax", values_to = "biomass") %>%
    left_join(phytoplankton_taxonomy, by = c("tax" = "taxon_code")) %>%
    group_by(lake_id, year, stratum, class_pr2) %>%
    summarize(biomass = sum(biomass)) %>%
    ungroup() %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_wrap(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               ncol = 1) +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations) +
    scale_fill_manual(values = palette_manipulation, name = "Experiment") +
    new_scale("fill") +
    geom_area(aes(x = year_floor, y = biomass/1000,
                  fill = factor(class_pr2, levels = phototrophs_tax_order$class)),
              stat = "identity",
              colour = "black", linewidth = 0.05) +
    geom_line(aes(x = year_floor, y = biomass/1000),
              data = biomass_epilimnion,
              linewidth = 0.01) +
    scale_x_date(date_breaks = "10 years", date_labels = "%Y",
                 limits = c(ymd("1870-01-01"), ymd("2020-01-01")), expand = c(0,0)) +
    scale_y_continuous(labels = comma,
                       breaks = seq(0, sum(phytoplankton_phototrophs_annual_taxoncode$biomass)/1000, 1)) +
    scale_fill_manual(values = palette_phototroph, name = "Class") +
    labs(y = "Biomass (mg/L)",
         colour = "Stratum") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "#00B0F0")))

# Sediment
phototrophs_top <- class_datasets %>%
  slice_max(order_by = pct_seqs, n = 15)

phototrophs_nseqs <- phototrophs_melt %>%
  group_by(sample_id) %>%
  summarize(nseqs_total = sum(nseqs)) %>%
  ungroup()

# Sediment phototroph taxonomic composition
(sediment_barplot <- phototrophs_melt %>%
    left_join(phototrophs_nseqs, by = "sample_id") %>%
    mutate(relseqs = nseqs/nseqs_total) %>%
    select(sample_id, asv_code, relseqs) %>%
    pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
    pivot_longer(!sample_id, names_to = "asv_code", values_to = "relseqs") %>%
    left_join(taxonomy, by = "asv_code") %>%
    left_join(metadata, by = "sample_id") %>%
    mutate(class = case_when(class %in% phototrophs_top$class ~ class,
                             TRUE ~ "Other")) %>%
    group_by(lake_id, sample_id, lowerpt_year, upperpt_year,
             class) %>%
    summarize(relseqs = sum(relseqs)) %>%
    ungroup() %>%
    arrange(sample_id, factor(class, levels = phototrophs_tax_order$class)) %>%
    group_by(sample_id) %>%
    mutate(cumsum_relseqs = cumsum(relseqs)) %>%
    mutate(ymax = cumsum_relseqs) %>%
    mutate(ymin = lag(cumsum_relseqs, default = 0)) %>%
    ungroup() %>%
    mutate(lowerpt_year = as.Date(date_decimal(lowerpt_year)),
           upperpt_year = as.Date(date_decimal(upperpt_year))) %>%
    ggplot() +
    facet_wrap(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               ncol = 1) +
    geom_rect(aes(xmin = lowerpt_year, xmax = upperpt_year,
                  ymin = ymin, ymax = ymax,
                  fill = factor(class, levels = phototrophs_tax_order$class)),
              colour = "black", linewidth = 0.05) +
    scale_fill_manual(values = c(palette_phototroph, "Other" = "black"), na.value = "black") +
    scale_x_date(date_breaks = "10 years", date_labels = "%Y",
                 limits = c(ymd("1870-01-01"), ymd("2020-01-01")), expand = c(0,0)) +
    labs(x = "Year",
         y = "Relative sequence abundance",
         fill = "Class") +
    theme_bw() %+replace%
    theme(#axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      #axis.title.x = element_blank(),
      axis.text = element_text(colour = "black"),
      #axis.line.y = element_line(colour = "black"),
      strip.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      panel.spacing = unit(2.3, "lines"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))

# Combine plots
(ela_taxonomy_plots_phototrophs <- biomass_phototrophs_byyear_epilimnion_interpolated_areaplot /
    sediment_barplot)
#ggsave("~/Desktop/ela_taxonomy_plots_phototrophs.pdf", ela_taxonomy_plots_phototrophs, width = 12, height = 18, units = "in", device = cairo_pdf)
