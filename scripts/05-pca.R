# Principal component analyses (PCAs)

# Load libraries
library(tidyverse)
library(vegan)
library(ggnewscale)
library(patchwork)

# Load palettes
source("scripts/00-palettes.R")

# Load sediment data
source("scripts/03e-sediment_data.R")


#### Plot individual lake stratigraphy ordinations ####
# Define function to plot PCA
plotPCA <- function(melt, legend = "none") {
  asvtable <- melt %>%
    select(sample_id, asv_code, nseqs) %>%
    pivot_wider(names_from = asv_code, values_from = nseqs, values_fill = 0) %>%
    column_to_rownames("sample_id")
  
  asvtable_hellinger <- decostand(asvtable, "hellinger")
  
  pca <- rda(asvtable_hellinger)
  
  pca_sites <- scores(pca)$sites %>%
    as_tibble(rownames = "sample_id") %>%
    left_join(metadata, by = "sample_id")
  
  # Flip PC1 axis (left to right: old assemblages to new assemblages)
  pca_sites <- pca_sites %>%
    mutate(PC1 = -PC1)
  
  pca_species <- scores(pca)$species %>%
    as_tibble(rownames = "asv_code")
  (pc1_percent_var <- round(x = (eigenvals(x = pca)[1])/(sum(eigenvals(x = pca)))*100, digits = 1))
  (pc2_percent_var <- round(x = (eigenvals(x = pca)[2])/(sum(eigenvals(x = pca)))*100, digits = 1))
  
  # PC contributions
  contributors <- pca_species %>%
    mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
    slice_max(magnitude, n = 10) %>%
    pull(asv_code)
  
  # Curate taxon contributor labels in PCA plot
  contribution <- pca_species %>%
    filter(asv_code %in% contributors) %>%
    left_join(taxonomy, by = "asv_code") %>%
    mutate(label = case_when(!is.na(species) ~ str_c(asv_code, " ", species, " (", subdivision, ")", sep = ""),
                             is.na(species) & !is.na(genus) ~ str_c(asv_code, " ", genus, " (", subdivision, ")", sep = ""),
                             is.na(genus) & !is.na(family) ~ str_c(asv_code, " ", family, " (", subdivision, ")", sep = ""),
                             is.na(family) & !is.na(order) ~ str_c(asv_code, " ", order, " (", subdivision, ")", sep = ""),
                             is.na(order) & !is.na(class) ~ str_c(asv_code, " ", class, " (", subdivision, ")", sep = ""),
                             is.na(class) & !is.na(subdivision) ~ str_c(asv_code, " ", subdivision),
                             is.na(subdivision) & !is.na(division) ~ str_c(asv_code, " ", division),
                             is.na(subdivision) & is.na(division) & !is.na(supergroup) ~ str_c(asv_code, " ", supergroup))) %>%
    mutate(label = str_replace_all(label, "_", " "))
  
  # Plot PCA
  (pca_plot <- ggplot() +
      geom_path(aes(x = PC1, y = PC2, group = lake_id,
                    colour = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))),
                data = pca_sites %>%
                  arrange(midpt),
                alpha = 0.7) +
      geom_point(aes(x = PC1, y = PC2,
                     colour = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
                     size = midpt_year),
                 data = pca_sites, alpha = 0.7, stroke = 0) +
      scale_colour_manual(values = palette_lake) +
      scale_size(limits = c(1880, 2020),
                 range = c(1, 8),
                 guide = guide_legend(reverse = TRUE)) +
      geom_text(aes(x = PC1, y = PC2,
                    label = round(midpt_year)),
                data = pca_sites,
                colour = "black",
                size = 2, fontface = "bold", check_overlap = TRUE) +
      labs(x = paste("PC1 (", pc1_percent_var, "%)", sep = ""),
           y = paste("PC2 (", pc2_percent_var, "%)", sep = ""),
           colour = "Lake",
           size = "Year") +
      theme(axis.text = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
            legend.key = element_blank(),
            legend.position = legend))
}
# Plot PCAs for entire microeukaryotic assemblages in individual lakes
(pca_plot_l226n <- plotPCA(sediment_melt %>%
                             filter(lake_id == "L226N")))
(pca_plot_l226s <- plotPCA(sediment_melt %>%
                             filter(lake_id == "L226S")))
(pca_plot_l227 <- plotPCA(sediment_melt %>%
                            filter(lake_id == "L227")))
(pca_plot_l224 <- plotPCA(sediment_melt %>%
                            filter(lake_id == "L224")))
(pca_plot_l373 <- plotPCA(sediment_melt %>%
                            filter(lake_id == "L373")))
(pca_plot2 <- plotPCA(sediment_melt, "right"))
(pca_plots_all <- pca_plot_l226n +
    pca_plot_l226s +
    pca_plot_l227 +
    pca_plot_l224 +
    pca_plot_l373 +
    pca_plot2 +
    plot_layout(guides = "collect"))
#ggsave("figures/ela18s_pca_all.pdf", pca_plots_all, width = 20, height = 12, units = "in")

# Plot PCAs for phototroph assemblages in individual lakes
(pca_plot_l226n_phototrophs <- plotPCA(phototrophs_melt %>%
                                         filter(lake_id == "L226N")))
(pca_plot_l226s_phototrophs <- plotPCA(phototrophs_melt %>%
                                         filter(lake_id == "L226S")))
(pca_plot_l227_phototrophs <- plotPCA(phototrophs_melt %>%
                                        filter(lake_id == "L227")))
(pca_plot_l224_phototrophs <- plotPCA(phototrophs_melt %>%
                                        filter(lake_id == "L224")))
(pca_plot_l373_phototrophs <- plotPCA(phototrophs_melt %>%
                                        filter(lake_id == "L373")))
(pca_plot_phototrophs <- plotPCA(phototrophs_melt, "right"))
(pca_plots_phototrophs_all <- pca_plot_l226n_phototrophs +
    pca_plot_l226s_phototrophs +
    pca_plot_l227_phototrophs +
    pca_plot_l224_phototrophs +
    pca_plot_l373_phototrophs +
    pca_plot_phototrophs +
    plot_layout(guides = "collect"))
#ggsave("figures/ela18s_pca_phototrophs_all.pdf", pca_plots_phototrophs_all, width = 18, height = 11, units = "in")


#### Save PCA axes to file ####
# Define function to compute PCA
computePCA <- function(melt) {
  asvtable <- melt %>%
    select(sample_id, asv_code, nseqs) %>%
    pivot_wider(names_from = asv_code, values_from = nseqs, values_fill = 0) %>%
    column_to_rownames("sample_id")
  
  # Transform community data
  asvtable_hellinger <- decostand(asvtable, "hellinger")
  
  # Compute PCA
  pca <- rda(asvtable_hellinger)
  
  return(pca)
}

# Define function to extract PCA axes 1 and 2
extractPCAAxes <- function(pca, dataset) {
  # Extract site coordinates
  pca_sites <- scores(pca)$sites %>%
    as_tibble(rownames = "sample_id") %>%
    left_join(metadata, by = "sample_id") %>%
    mutate(dataset = dataset) %>%
    arrange(midpt)
  
  # Flip PC1 axis (left to right: old assemblages to new assemblages)
  pca_sites <- pca_sites %>%
    mutate(PC1 = -PC1)
  
  return(pca_sites)
}

# Define function to compute PCAs for all lakes
computePCAAll <- function(melt, dataset) {
  pca_all <- bind_rows(computePCA(melt %>%
                                    filter(lake_id == "L226N")) %>%
                         extractPCAAxes(dataset),
                       computePCA(melt %>%
                                    filter(lake_id == "L226S")) %>%
                         extractPCAAxes(dataset),
                       computePCA(melt %>%
                                    filter(lake_id == "L227")) %>%
                         extractPCAAxes(dataset),
                       computePCA(melt %>%
                                    filter(lake_id == "L224")) %>%
                         extractPCAAxes(dataset),
                       computePCA(melt %>%
                                    filter(lake_id == "L373")) %>%
                         extractPCAAxes(dataset))
  
  return(pca_all)
}

pca_all <- bind_rows(computePCAAll(sediment_melt, "Sediment ASVs"),
                     computePCAAll(phototrophs_melt, "Phototroph ASVs"),
                     computePCAAll(chrysophyceae_melt, "Chrysophyceae ASVs"))

# Write PCA axis data to file
# pca_all %>%
#   write_tsv("output/ordinations/ela18s_pca_axes.tsv", col_names = TRUE)
