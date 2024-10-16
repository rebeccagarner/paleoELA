# PCA sediment DNA and annually resolved phytoplankton biomass monitoring

setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
setwd("~/Dropbox/projects/ela18s/")

# Load libraries
library(tidyverse)
library(vegan)
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
phytoplankton_phototrophs_annual_genus <- sumPhytoplanktonAnnual(phytoplankton_phototrophs, "genus_pr2")


#### Perform PCA ####
plotPCATimeseries <- function(dna_melt, phytoplankton_melt, lakeid) {
  # Format site by taxon tables for DNA and phytoplankton assemblages
  asvtable <- dna_melt %>%
    filter(lake_id == lakeid) %>%
    select(sample_id, genus, nseqs) %>%
    rename(diversity_unit = nseqs)
  
  phytoplankton_table <- phytoplankton_melt %>%
    filter(lake_id == lakeid) %>%
    mutate(sample_id = str_c(lake_id, year, stratum, sep = "_")) %>%
    dplyr::select(sample_id, tax, biomass) %>%
    rename(genus = tax,
           diversity_unit = biomass)
  
  sample_by_taxon <- bind_rows(asvtable,
                               phytoplankton_table) %>%
    group_by(sample_id, genus) %>%
    summarize(diversity_unit = sum(diversity_unit)) %>%
    ungroup() %>%
    pivot_wider(names_from = genus, values_from = diversity_unit, values_fill = 0) %>%
    column_to_rownames("sample_id")
  
  # Transform assemblage data
  sample_by_taxon_hellinger <- decostand(sample_by_taxon, "hellinger")
  
  # Compute PCA
  pca <- rda(sample_by_taxon_hellinger)
  
  pca_sites <- scores(pca)$sites %>%
    as_tibble(rownames = "sample_id") %>%
    mutate(timeseries = case_when(grepl("cm$", sample_id) ~ "sediment",
                                  grepl("\\d\\d\\d\\d", sample_id) ~ "monitoring")) %>%
    mutate(lake_id = str_remove_all(sample_id, "_.*")) %>%
    mutate(sample_id_abbr = case_when(timeseries == "sediment" ~ str_remove(sample_id, str_c(lakeid, "_")),
                                      timeseries == "monitoring" ~ str_extract(sample_id, "\\d\\d\\d\\d")))
  pca_species <- scores(pca)$species %>%
    as_tibble(rownames = "genus")
  
  pc1_percent_var <- round(x = (eigenvals(x = pca)[1])/(sum(eigenvals(x = pca)))*100, digits = 1)
  pc2_percent_var <- round(x = (eigenvals(x = pca)[2])/(sum(eigenvals(x = pca)))*100, digits = 1)
  
  # # PC contributions
  # contributors <- pca_species %>%
  #   mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
  #   slice_max(magnitude, n = 10) %>%
  #   mutate(label = case_when(!is.na(genus_pr2) ~ str_c(genus_pr2, " (", class_pr2, ")"),
  #                            !is.na(family_pr2) ~ str_c(family_pr2, " (", class_pr2, ")"),
  #                            !is.na(order_pr2) ~ str_c(order_pr2, " (", class_pr2, ")"),
  #                            !is.na(class_pr2) ~ str_c(class_pr2)))
  
  # scores(pca, display = "species", choices = c(1), scaling = 1) %>%
  #   as_tibble(rownames = "taxon_code") %>%
  #   left_join(phytoplankton_taxonomy, by = "taxon_code") %>%
  #   slice_max(abs(PC1), n = 5)
  # scores(pca, display = "species", choices = c(2), scaling = 1) %>%
  #   as_tibble(rownames = "taxon_code") %>%
  #   left_join(phytoplankton_taxonomy, by = "taxon_code") %>%
  #   slice_max(abs(PC2), n = 5)
  
  # Plot PCA
  (pca_plot <- ggplot() +
      # geom_text(aes(x = x*PC1, y = x*PC2, label = label),
      #           data = contributors,
      #           colour = "#A0A0A0", size = 4) +
      # geom_segment(aes(x = 0, y = 0, xend = (x - 0.1*x)*PC1, yend = (x - 0.1*x)*PC2, colour = subdivision_pr2),
      #              data = contributors,
      #              arrow = arrow(length = unit(x = 0.2, units = "cm")),
      #              alpha = 1) +
      # scale_colour_manual(values = palette_subdivision, guide = "none") +
      # new_scale_colour() +
      # geom_path(aes(x = PC1, y = PC2, group = lake_id,
      #               colour = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))),
    #           data = pca_sites,
    #           alpha = 0.7) +
    geom_point(aes(x = PC1, y = PC2,
                   colour = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
                   shape = timeseries),
               data = pca_sites, alpha = 0.7, stroke = 0) +
      scale_colour_manual(values = palette_lake) +
      # scale_size(range = c(1, 8), guide = guide_legend(reverse = TRUE)) +
      geom_text(aes(x = PC1, y = PC2,
                    label = sample_id_abbr),
                data = pca_sites,
                colour = "black",
                size = 2, fontface = "bold", check_overlap = TRUE) +
      # annotate("text", x = Inf, y = Inf,
      #          label = str_c(min(pca_sites$year), "-", max(pca_sites$year)),
      #          hjust = 1, vjust = 1) +
      labs(x = paste("PC1 (", pc1_percent_var, "%)", sep = ""),
           y = paste("PC2 (", pc2_percent_var, "%)", sep = ""),
           colour = "Lake",
           size = "Year") +
      theme(axis.text = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
            legend.key = element_blank(),
            legend.position = "none"))
}

(pca_plot_l226n_timeseries <- plotPCATimeseries(phototrophs_melt, phytoplankton_phototrophs_annual_genus %>%
                                                  filter(stratum == "epilimnion"), "L226N"))
(pca_plot_l226s_timeseries <- plotPCATimeseries(phototrophs_melt, phytoplankton_phototrophs_annual_genus %>%
                                                  filter(stratum == "epilimnion"), "L226S"))
(pca_plot_l227_timeseries <- plotPCATimeseries(phototrophs_melt, phytoplankton_phototrophs_annual_genus %>%
                                                 filter(stratum == "epilimnion"), "L227"))
(pca_plot_l224_timeseries <- plotPCATimeseries(phototrophs_melt, phytoplankton_phototrophs_annual_genus %>%
                                                 filter(stratum == "epilimnion"), "L224"))
(pca_plot_l373_timeseries <- plotPCATimeseries(phototrophs_melt, phytoplankton_phototrophs_annual_genus %>%
                                                 filter(stratum == "epilimnion"), "L373"))
(pca_plot_timeseries <- plotPCATimeseries(phototrophs_melt, phytoplankton_phototrophs_annual_genus %>%
                                            filter(stratum == "epilimnion"), c("L226N", "L226S", "L227", "L224", "L373")))
(pca_plots_timeseries_all <- pca_plot_l226n_timeseries +
    pca_plot_l226s_timeseries +
    pca_plot_l227_timeseries +
    pca_plot_l224_timeseries +
    pca_plot_l373_timeseries +
    pca_plot_timeseries +
    plot_layout(guides = "collect"))
