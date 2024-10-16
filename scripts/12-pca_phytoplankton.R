# PCA annually resolved phytoplankton biomass monitoring

setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
setwd("~/Dropbox/projects/ela18s/")
setwd("~/Library/CloudStorage/Dropbox-Personal/projects/ela18s/")

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
# Define function to plot PCA
plotPCA <- function(melt, legend = "none", x = 1) {
  # Format site by taxon table
  sample_by_taxon <- melt %>%
    mutate(sample_id = str_c(lake_id, year, stratum, sep = "_")) %>%
    dplyr::select(sample_id, tax, biomass) %>%
    pivot_wider(names_from = tax, values_from = biomass, values_fill = 0) %>%
    column_to_rownames("sample_id")
  
  # Transform assemblage data
  sample_by_taxon_hellinger <- decostand(sample_by_taxon, "hellinger")
  
  # Compute PCA
  pca <- rda(sample_by_taxon_hellinger)
  
  pca_sites <- scores(pca)$sites %>%
    as_tibble(rownames = "sample_id") %>%
    separate(sample_id, c("lake_id", "year", "stratum"), sep = "_") %>%
    mutate(year = as.numeric(year))
  
  # Flip PC1 axis (left to right: old assemblages to new assemblages)
  pca_sites <- pca_sites %>%
    mutate(PC1 = -PC1)
  
  pca_species <- scores(pca)$species %>%
    as_tibble(rownames = "tax")
  
  if (all(unique(pca_species$tax) %in% unique(phytoplankton_taxonomy$taxon_code))) {
    pca_species <- pca_species %>%
      rename(taxon_code = tax) %>%
      left_join(phytoplankton_taxonomy, by = "taxon_code")
  } else if (all(unique(pca_species$tax) %in% unique(phytoplankton_taxonomy$genus_pr2))) {
    pca_species <- pca_species %>%
      rename(genus_pr2 = tax) %>%
      left_join(phytoplankton_taxonomy %>%
                  distinct(supergroup_pr2, division_pr2, subdivision_pr2, class_pr2, order_pr2, family_pr2, genus_pr2),
                by = "genus_pr2")
  }
  
  pc1_percent_var <- round(x = (eigenvals(x = pca)[1])/(sum(eigenvals(x = pca)))*100, digits = 1)
  pc2_percent_var <- round(x = (eigenvals(x = pca)[2])/(sum(eigenvals(x = pca)))*100, digits = 1)
  
  # PC contributions
  contributors <- pca_species %>%
    mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
    slice_max(magnitude, n = 10) %>%
    mutate(label = case_when(!is.na(genus_pr2) ~ str_c(genus_pr2, " (", class_pr2, ")"),
                             !is.na(family_pr2) ~ str_c(family_pr2, " (", class_pr2, ")"),
                             !is.na(order_pr2) ~ str_c(order_pr2, " (", class_pr2, ")"),
                             !is.na(class_pr2) ~ str_c(class_pr2)))
  
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
      geom_path(aes(x = PC1, y = PC2, group = lake_id,
                    colour = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))),
                data = pca_sites %>%
                  arrange(year),
                alpha = 0.7) +
      geom_point(aes(x = PC1, y = PC2,
                     colour = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
                     size = as.numeric(year)),
                 data = pca_sites, alpha = 0.7, stroke = 0) +
      scale_colour_manual(values = palette_lake) +
      scale_size(limits = c(1880, 2020),
                 range = c(1, 8),
                 guide = guide_legend(reverse = TRUE)) +
      geom_text(aes(x = PC1, y = PC2,
                    label = year),
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
            legend.position = legend))
}
# Plot PCAs for phytoplankton phototroph taxon_code epilimnion monitoring assemblages
(pca_plot_l226n_phototrophs_taxoncode_epilimnion <- plotPCA(phytoplankton_phototrophs_annual_taxoncode %>%
                                                              filter(lake_id == "L226N" & stratum == "epilimnion"), x = 1))
(pca_plot_l226s_phototrophs_taxoncode_epilimnion <- plotPCA(phytoplankton_phototrophs_annual_taxoncode %>%
                                                              filter(lake_id == "L226S" & stratum == "epilimnion"), x = 1))
(pca_plot_l227_phototrophs_taxoncode_epilimnion <- plotPCA(phytoplankton_phototrophs_annual_taxoncode %>%
                                                             filter(lake_id == "L227" & stratum == "epilimnion"), x = 0.9))
(pca_plot_l224_phototrophs_taxoncode_epilimnion <- plotPCA(phytoplankton_phototrophs_annual_taxoncode %>%
                                                             filter(lake_id == "L224" & stratum == "epilimnion"), x = 1))
(pca_plot_l373_phototrophs_taxoncode_epilimnion <- plotPCA(phytoplankton_phototrophs_annual_taxoncode %>%
                                                             filter(lake_id == "L373" & stratum == "epilimnion"), x = 1.2))
(pca_plot_phototrophs_taxoncode_epilimnion <- plotPCA(phytoplankton_phototrophs_annual_taxoncode %>%
                                                        filter(stratum == "epilimnion"), "right", x = 0.6))
(pca_plots_phototrophs_taxoncode_epilimnion_all <- pca_plot_l226n_phototrophs_taxoncode_epilimnion +
    pca_plot_l226s_phototrophs_taxoncode_epilimnion +
    pca_plot_l227_phototrophs_taxoncode_epilimnion +
    pca_plot_l224_phototrophs_taxoncode_epilimnion +
    pca_plot_l373_phototrophs_taxoncode_epilimnion +
    pca_plot_phototrophs_taxoncode_epilimnion +
    plot_layout(guides = "collect"))
#ggsave("figures/ela_phytoplankton_phototrophs_annual_taxoncode_epilimnion_pca.pdf", pca_plots_phototrophs_taxoncode_epilimnion_all, width = 18, height = 11, units = "in", device = cairo_pdf)

# Plot PCAs for phytoplankton phototroph genus epilimnion monitoring assemblages
(pca_plot_l226n_phototrophs_genus_epilimnion <- plotPCA(phytoplankton_phototrophs_annual_genus %>%
                                                          filter(lake_id == "L226N" & stratum == "epilimnion"), x = 1))
(pca_plot_l226s_phototrophs_genus_epilimnion <- plotPCA(phytoplankton_phototrophs_annual_genus %>%
                                                          filter(lake_id == "L226S" & stratum == "epilimnion"), x = 1))
(pca_plot_l227_phototrophs_genus_epilimnion <- plotPCA(phytoplankton_phototrophs_annual_genus %>%
                                                         filter(lake_id == "L227" & stratum == "epilimnion"), x = 0.9))
(pca_plot_l224_phototrophs_genus_epilimnion <- plotPCA(phytoplankton_phototrophs_annual_genus %>%
                                                         filter(lake_id == "L224" & stratum == "epilimnion"), x = 1))
(pca_plot_l373_phototrophs_genus_epilimnion <- plotPCA(phytoplankton_phototrophs_annual_genus %>%
                                                         filter(lake_id == "L373" & stratum == "epilimnion"), x = 1.2))
(pca_plot_phototrophs_genus_epilimnion <- plotPCA(phytoplankton_phototrophs_annual_genus %>%
                                                    filter(stratum == "epilimnion"), "right", x = 0.6))
(pca_plots_phototrophs_genus_epilimnion_all <- pca_plot_l226n_phototrophs_genus_epilimnion +
    pca_plot_l226s_phototrophs_genus_epilimnion +
    pca_plot_l227_phototrophs_genus_epilimnion +
    pca_plot_l224_phototrophs_genus_epilimnion +
    pca_plot_l373_phototrophs_genus_epilimnion +
    pca_plot_phototrophs_genus_epilimnion +
    plot_layout(guides = "collect"))

# Define function to plot PCA with both epilimnion and metalimnion assemblages
plotPCAStrat <- function(melt, legend = "none", x = 1) {
  # Format site by taxon table
  sample_by_taxon <- melt %>%
    mutate(sample_id = str_c(lake_id, year, stratum, sep = "_")) %>%
    dplyr::select(sample_id, tax, biomass) %>%
    pivot_wider(names_from = tax, values_from = biomass, values_fill = 0) %>%
    column_to_rownames("sample_id")
  
  # Transform assemblage data
  sample_by_taxon_hellinger <- decostand(sample_by_taxon, "hellinger")
  
  # Compute PCA
  pca <- rda(sample_by_taxon_hellinger)
  
  pca_sites <- scores(pca)$sites %>%
    as_tibble(rownames = "sample_id") %>%
    separate(sample_id, c("lake_id", "year", "stratum"), sep = "_") %>%
    mutate(year = as.numeric(year)) %>%
    mutate(stratum = str_to_sentence(stratum))
  
  # Flip PC1 axis (left to right: old assemblages to new assemblages)
  pca_sites <- pca_sites %>%
    mutate(PC1 = -PC1)
  
  pca_species <- scores(pca)$species %>%
    as_tibble(rownames = "tax")
  
  if (all(unique(pca_species$tax) %in% unique(phytoplankton_taxonomy$taxon_code))) {
    pca_species <- pca_species %>%
      rename(taxon_code = tax) %>%
      left_join(phytoplankton_taxonomy, by = "taxon_code")
  } else if (all(unique(pca_species$tax) %in% unique(phytoplankton_taxonomy$genus_pr2))) {
    pca_species <- pca_species %>%
      rename(genus_pr2 = tax) %>%
      left_join(phytoplankton_taxonomy %>%
                  distinct(supergroup_pr2, division_pr2, subdivision_pr2, class_pr2, order_pr2, family_pr2, genus_pr2),
                by = "genus_pr2")
  }
  
  pc1_percent_var <- round(x = (eigenvals(x = pca)[1])/(sum(eigenvals(x = pca)))*100, digits = 1)
  pc2_percent_var <- round(x = (eigenvals(x = pca)[2])/(sum(eigenvals(x = pca)))*100, digits = 1)
  
  # PC contributions
  contributors <- pca_species %>%
    mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
    slice_max(magnitude, n = 10) %>%
    mutate(label = case_when(!is.na(genus_pr2) ~ str_c(genus_pr2, " (", class_pr2, ")"),
                             !is.na(family_pr2) ~ str_c(family_pr2, " (", class_pr2, ")"),
                             !is.na(order_pr2) ~ str_c(order_pr2, " (", class_pr2, ")"),
                             !is.na(class_pr2) ~ str_c(class_pr2)))
  
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
      new_scale_colour() +
      geom_path(aes(x = PC1, y = PC2, group = interaction(lake_id, year),
                    colour = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))),
                data = pca_sites %>%
                  arrange(year),
                alpha = 0.7) +
      geom_point(aes(x = PC1, y = PC2,
                     colour = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
                     shape = stratum,
                     size = as.numeric(year)),
                 data = pca_sites, alpha = 0.7, stroke = 0) +
      scale_colour_manual(values = palette_lake) +
      scale_shape_manual(values = palette_stratum_shape) +
      scale_size(range = c(1, 8), guide = guide_legend(reverse = TRUE)) +
      geom_text(aes(x = PC1, y = PC2,
                    label = year),
                data = pca_sites,
                colour = "black",
                size = 2, fontface = "bold", check_overlap = TRUE) +
      annotate("text", x = Inf, y = Inf,
               label = str_c(min(pca_sites$year), "-", max(pca_sites$year)),
               hjust = 1, vjust = 1) +
      labs(x = paste("PC1 (", pc1_percent_var, "%)", sep = ""),
           y = paste("PC2 (", pc2_percent_var, "%)", sep = ""),
           colour = "Lake",
           shape = "Stratum",
           size = "Year") +
      theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black"),
            panel.grid = element_blank(), panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
            legend.key = element_blank(), legend.position = legend))
}

# Plot PCAs for phytoplankton phototroph taxon_code epilimnion and metalimnion monitoring assemblages
(pca_plot_l226n_phototrophs_taxoncode <- plotPCAStrat(phytoplankton_phototrophs_annual_taxoncode %>%
                                                        filter(lake_id == "L226N"), x = 1))
(pca_plot_l226s_phototrophs_taxoncode <- plotPCAStrat(phytoplankton_phototrophs_annual_taxoncode %>%
                                                        filter(lake_id == "L226S"), x = 1))
(pca_plot_l227_phototrophs_taxoncode <- plotPCAStrat(phytoplankton_phototrophs_annual_taxoncode %>%
                                                       filter(lake_id == "L227"), x = 0.9))
(pca_plot_l224_phototrophs_taxoncode <- plotPCAStrat(phytoplankton_phototrophs_annual_taxoncode %>%
                                                       filter(lake_id == "L224"), x = 1))
(pca_plot_l373_phototrophs_taxoncode <- plotPCAStrat(phytoplankton_phototrophs_annual_taxoncode %>%
                                                       filter(lake_id == "L373"), x = 1.2))
(pca_plot_phototrophs_taxoncode <- plotPCAStrat(phytoplankton_phototrophs_annual_taxoncode, "right", x = 0.6))
(pca_plots_phototrophs_taxoncode_all <- pca_plot_l226n_phototrophs_taxoncode +
    pca_plot_l226s_phototrophs_taxoncode +
    pca_plot_l227_phototrophs_taxoncode +
    pca_plot_l224_phototrophs_taxoncode +
    pca_plot_l373_phototrophs_taxoncode +
    pca_plot_phototrophs_taxoncode +
    plot_layout(guides = "collect"))


#### Save PCA axes to file ####
# Define function to compute PCA
computePCA <- function(melt) {
  # Format site by taxon table
  sample_by_taxon <- melt %>%
    mutate(sample_id = str_c(lake_id, year, stratum, sep = "_")) %>%
    dplyr::select(sample_id, tax, biomass) %>%
    pivot_wider(names_from = tax, values_from = biomass, values_fill = 0) %>%
    column_to_rownames("sample_id")
  
  # Transform assemblage data
  sample_by_taxon_hellinger <- decostand(sample_by_taxon, "hellinger")
  
  # Compute PCA
  pca <- rda(sample_by_taxon_hellinger)
  
  return(pca)
}

# Define function to extract PCA axes 1 and 2
extractPCAAxes <- function(pca, dataset) {
  # Extract site coordinates
  pca_sites <- scores(pca)$sites %>%
    as_tibble(rownames = "sample_id") %>%
    separate(sample_id, c("lake_id", "year", "stratum"), sep = "_") %>%
    mutate(year = as.numeric(year)) %>%
    mutate(dataset = dataset) %>%
    arrange(-year)
  
  # Flip PC1 axis (left to right: old assemblages to new assemblages)
  pca_sites <- pca_sites %>%
    mutate(PC1 = -PC1)
  
  return(pca_sites)
}

# Define function to compute PCAs for epilimnion and metalimnion assemblages in all lakes
computePCAAll <- function(melt, dataset) {
  pca_all <- bind_rows(computePCA(melt %>%
                                    filter(lake_id == "L226N" & stratum == "epilimnion")) %>%
                         extractPCAAxes(dataset),
                       computePCA(melt %>%
                                    filter(lake_id == "L226S" & stratum == "epilimnion")) %>%
                         extractPCAAxes(dataset),
                       computePCA(melt %>%
                                    filter(lake_id == "L227" & stratum == "epilimnion")) %>%
                         extractPCAAxes(dataset),
                       computePCA(melt %>%
                                    filter(lake_id == "L224" & stratum == "epilimnion")) %>%
                         extractPCAAxes(dataset),
                       computePCA(melt %>%
                                    filter(lake_id == "L373" & stratum == "epilimnion")) %>%
                         extractPCAAxes(dataset),
                       
                       computePCA(melt %>%
                                    filter(lake_id == "L226N" & stratum == "metalimnion")) %>%
                         extractPCAAxes(dataset),
                       computePCA(melt %>%
                                    filter(lake_id == "L226S" & stratum == "metalimnion")) %>%
                         extractPCAAxes(dataset),
                       computePCA(melt %>%
                                    filter(lake_id == "L227" & stratum == "metalimnion")) %>%
                         extractPCAAxes(dataset),
                       computePCA(melt %>%
                                    filter(lake_id == "L224" & stratum == "metalimnion")) %>%
                         extractPCAAxes(dataset),
                       computePCA(melt %>%
                                    filter(lake_id == "L373" & stratum == "metalimnion")) %>%
                         extractPCAAxes(dataset))
  
  return(pca_all)
}

pca_all <- bind_rows(computePCAAll(phytoplankton_phototrophs_annual_taxoncode, "Phototroph taxon_codes"),
                     computePCAAll(phytoplankton_phototrophs_annual_genus, "Phototroph genera")) %>%
  mutate(dataset = str_c(dataset, " ", stratum))

# Write PCA axis data to file
# pca_all %>%
#   write_tsv("output/ordinations/ela_phytoplankton_pca_axes.tsv", col_names = TRUE)
