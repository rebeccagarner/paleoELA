# PCA of ELA phytoplankton monitoring assemblages

# Load libraries
library(tidyverse)
library(vegan)
library(lubridate)
library(ggnewscale)
library(patchwork)

# Load palettes
source("scripts/00-palettes.R")

# Load environmental data
source("scripts/08-ela_env_data.R")


#### Import and format data ####
# Import phytoplankton species count and biomass data
phytoplankton_species <- read_tsv("output/ela/ela_monitoring_phytoplankton_species.tsv", col_names = TRUE)


#### Format phytoplankton species monitoring data ####
# Format site by taxon table
sample_by_taxon <- phytoplankton_species %>%
  mutate(sample_id = str_c(lake_id, date, stratum, sep = "_")) %>%
  dplyr::select(sample_id, taxon_code, biomass) %>%
  pivot_wider(names_from = taxon_code, values_from = biomass, values_fill = 0) %>%
  column_to_rownames("sample_id")

# Format phytoplankton taxonomy
phytoplankton_taxonomy <- phytoplankton_species %>%
  distinct(taxon_code, supergroup_pr2, division_pr2, subdivision_pr2, class_pr2, order_pr2, family_pr2, genus_pr2) %>%
  arrange(supergroup_pr2, division_pr2, subdivision_pr2, class_pr2, order_pr2, family_pr2, genus_pr2)

# Parse sample information
samples <- phytoplankton_species %>%
  mutate(sample_id = str_c(lake_id, date, stratum, sep = "_")) %>%
  distinct(sample_id, lake_id, date, stratum) %>%
  mutate(year = year(date)) %>%
  assignSeason()


#### Perform PCA ####
# Define function to plot PCA
plotPCA <- function(melt, legend = "none", x = 1) {
  # Format site by taxon table
  sample_by_taxon <- melt %>%
    mutate(sample_id = str_c(lake_id, date, stratum, sep = "_")) %>%
    dplyr::select(sample_id, taxon_code, biomass) %>%
    pivot_wider(names_from = taxon_code, values_from = biomass, values_fill = 0) %>%
    column_to_rownames("sample_id")
  
  # Transform assemblage data
  sample_by_taxon_hellinger <- decostand(sample_by_taxon, "hellinger")
  
  # Compute PCA
  pca <- rda(sample_by_taxon_hellinger)
  
  pca_sites <- scores(pca)$sites %>%
    as_tibble(rownames = "sample_id") %>%
    left_join(samples, by = "sample_id")
  pca_species <- scores(pca)$species %>%
    as_tibble(rownames = "taxon_code") %>%
    left_join(phytoplankton_taxonomy, by = "taxon_code")
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
  
  scores(pca, display = "species", choices = c(1), scaling = 1) %>%
    as_tibble(rownames = "taxon_code") %>%
    left_join(phytoplankton_taxonomy, by = "taxon_code") %>%
    slice_max(abs(PC1), n = 5)
  scores(pca, display = "species", choices = c(2), scaling = 1) %>%
    as_tibble(rownames = "taxon_code") %>%
    left_join(phytoplankton_taxonomy, by = "taxon_code") %>%
    slice_max(abs(PC2), n = 5)
  
  # Plot PCA
  (pca_plot <- ggplot() +
      geom_text(aes(x = x*PC1, y = x*PC2, label = label),
                data = contributors,
                colour = "#A0A0A0", size = 4) +
      geom_segment(aes(x = 0, y = 0, xend = (x - 0.1*x)*PC1, yend = (x - 0.1*x)*PC2, colour = subdivision_pr2),
                   data = contributors,
                   arrow = arrow(length = unit(x = 0.2, units = "cm")),
                   alpha = 1) +
      scale_colour_manual(values = palette_subdivision, guide = "none") +
      new_scale_colour() +
      geom_point(aes(x = PC1, y = PC2,
                     colour = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
                     shape = season,
                     size = as.numeric(year)),
                 data = pca_sites, alpha = 0.4, stroke = 0) +
      scale_colour_manual(values = palette_lake) +
      scale_shape_manual(values = c("Winter" = 49, "Spring" = 50, "Summer" = 51, "Fall" = 52)) +
      scale_size(range = c(1, 8), guide = guide_legend(reverse = TRUE)) +
      annotate("text", x = Inf, y = Inf,
               label = str_c(min(pca_sites$year), "-", max(pca_sites$year)),
               hjust = 1, vjust = 1) +
      labs(x = paste("PC1 (", pc1_percent_var, "%)", sep = ""),
           y = paste("PC2 (", pc2_percent_var, "%)", sep = ""),
           colour = "Lake",
           shape = "Season",
           size = "Year") +
      theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black"),
            panel.grid = element_blank(), panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
            legend.key = element_blank(), legend.position = legend))
}

# Plot PCAs for all phytoplankton epilimnion monitoring assemblages
(pca_plot_l226n_epilimnion <- plotPCA(phytoplankton_species %>%
                             filter(lake_id == "L226N" & stratum == "Epilimnion"), x = 0.8))
(pca_plot_l226s_epilimnion <- plotPCA(phytoplankton_species %>%
                             filter(lake_id == "L226S" & stratum == "Epilimnion"), x = 0.8))
(pca_plot_l227_epilimnion <- plotPCA(phytoplankton_species %>%
                            filter(lake_id == "L227" & stratum == "Epilimnion"), x = 0.7))
(pca_plot_l224_epilimnion <- plotPCA(phytoplankton_species %>%
                            filter(lake_id == "L224" & stratum == "Epilimnion"), x = 0.6))
(pca_plot_l373_epilimnion <- plotPCA(phytoplankton_species %>%
                            filter(lake_id == "L373" & stratum == "Epilimnion"), x = 0.9))
(pca_epilimnion_plot <- plotPCA(phytoplankton_species %>%
                                  filter(stratum == "Epilimnion"), "right", x = 0.3))
(pca_epilimnion_plots_all <- pca_plot_l226n_epilimnion +
    pca_plot_l226s_epilimnion +
    pca_plot_l227_epilimnion +
    pca_plot_l224_epilimnion +
    pca_plot_l373_epilimnion +
    pca_epilimnion_plot +
    plot_layout(guides = "collect"))
#ggsave("figures/ela_phytoplankton_species_epilimnion_pca.pdf", pca_epilimnion_plots_all, width = 20, height = 12, units = "in", device = cairo_pdf)

# Plot PCAs for all phytoplankton metalimnion monitoring assemblages
(pca_plot_l226n_metalimnion <- plotPCA(phytoplankton_species %>%
                                         filter(lake_id == "L226N" & stratum == "Metalimnion"), x = 0.8))
(pca_plot_l226s_metalimnion <- plotPCA(phytoplankton_species %>%
                                         filter(lake_id == "L226S" & stratum == "Metalimnion"), x = 0.6))
(pca_plot_l227_metalimnion <- plotPCA(phytoplankton_species %>%
                                        filter(lake_id == "L227" & stratum == "Metalimnion"), x = 0.5))
(pca_plot_l224_metalimnion <- plotPCA(phytoplankton_species %>%
                                        filter(lake_id == "L224" & stratum == "Metalimnion"), x = 0.7))
(pca_plot_l373_metalimnion <- plotPCA(phytoplankton_species %>%
                                        filter(lake_id == "L373" & stratum == "Metalimnion"), x = 0.9))
(pca_metalimnion_plot <- plotPCA(phytoplankton_species %>%
                                   filter(stratum == "Metalimnion"), "right", x = 0.3))
(pca_metalimnion_plots_all <- pca_plot_l226n_metalimnion +
    pca_plot_l226s_metalimnion +
    pca_plot_l227_metalimnion +
    pca_plot_l224_metalimnion +
    pca_plot_l373_metalimnion +
    pca_metalimnion_plot +
    plot_layout(guides = "collect"))
#ggsave("figures/ela_phytoplankton_species_metalimnion_pca.pdf", pca_metalimnion_plots_all, width = 20, height = 12, units = "in", device = cairo_pdf)
