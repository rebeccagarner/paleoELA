# PC taxon contributions

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


#### Extract PC1 taxon loadings ####
extractPCl1oadings <- function(melt, timeseries) {
  if (timeseries == "sediment") {
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
    
    # Flip PC1 axis (left to right: old assemblages to new assemblages)
    pca_species <- pca_species %>%
      mutate(PC1 = -PC1)
    
    (pc1_percent_var <- round(x = (eigenvals(x = pca)[1])/(sum(eigenvals(x = pca)))*100, digits = 1))
    (pc2_percent_var <- round(x = (eigenvals(x = pca)[2])/(sum(eigenvals(x = pca)))*100, digits = 1))
    
    # PC contributions
    sum(abs(pca_species$PC1))
    
    pc1_loadings <- pca_species %>%
      left_join(taxonomy, by = "asv_code") %>%
      slice_max(abs(PC1), n = 20) %>%
      mutate(direction = case_when(PC1 > 0 ~ "positive",
                                   PC1 < 0 ~ "negative")) %>%
      mutate(label = case_when(!is.na(species) ~ str_c(asv_code, " s__", species, " (", subdivision, ")", sep = ""),
                               is.na(species) & !is.na(genus) ~ str_c(asv_code, " g__", genus, " (", subdivision, ")", sep = ""),
                               is.na(genus) & !is.na(family) ~ str_c(asv_code, " f__", family, " (", subdivision, ")", sep = ""),
                               is.na(family) & !is.na(order) ~ str_c(asv_code, " o__", order, " (", subdivision, ")", sep = ""),
                               is.na(order) & !is.na(class) ~ str_c(asv_code, " c__", class, " (", subdivision, ")", sep = ""),
                               is.na(class) & !is.na(subdivision) ~ str_c(asv_code, " sd__", subdivision),
                               is.na(subdivision) & !is.na(division) ~ str_c(asv_code, " d__", division),
                               is.na(subdivision) & is.na(division) & !is.na(supergroup) ~ str_c(asv_code, " sg__", supergroup))) %>%
      mutate(timeseries = timeseries,
             lake_id = unique(melt$lake_id))
    
    return(pc1_loadings)
    
  } else if (timeseries == "monitoring") {
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
    
    # Flip PC1 axis (left to right: old assemblages to new assemblages)
    pca_species <- pca_species %>%
      mutate(PC1 = -PC1)
    
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
    
    (pc1_percent_var <- round(x = (eigenvals(x = pca)[1])/(sum(eigenvals(x = pca)))*100, digits = 1))
    (pc2_percent_var <- round(x = (eigenvals(x = pca)[2])/(sum(eigenvals(x = pca)))*100, digits = 1))
    
    # PC contributions
    sum(abs(pca_species$PC1))
    
    pc1_loadings <- pca_species %>%
      rename_with(~str_remove(.x, "_pr2"), ends_with("_pr2")) %>%
      slice_max(abs(PC1), n = 20) %>%
      mutate(direction = case_when(PC1 > 0 ~ "positive",
                                   PC1 < 0 ~ "negative")) %>%
      mutate(label = case_when(!is.na(genus) ~ str_c(taxon_code, " g__", genus, " (", subdivision, ")", sep = ""),
                               is.na(genus) & !is.na(family) ~ str_c(taxon_code, " f__", family, " (", subdivision, ")", sep = ""),
                               is.na(family) & !is.na(order) ~ str_c(taxon_code, " o__", order, " (", subdivision, ")", sep = ""),
                               is.na(order) & !is.na(class) ~ str_c(taxon_code, " c__", class, " (", subdivision, ")", sep = ""),
                               is.na(class) & !is.na(subdivision) ~ str_c(taxon_code, " sd__", subdivision),
                               is.na(subdivision) & !is.na(division) ~ str_c(taxon_code, " d__", division),
                               is.na(subdivision) & is.na(division) & !is.na(supergroup) ~ str_c(taxon_code, " sg__", supergroup))) %>%
      mutate(timeseries = timeseries,
             lake_id = unique(melt$lake_id))
    
    return(pc1_loadings)
  }
}

# Extract PC1 loadings
(pc1_loadings_l226n_monitoring <- extractPCl1oadings(phytoplankton_phototrophs_annual_taxoncode %>%
                                                       filter(lake_id == "L226N" & stratum == "epilimnion"), timeseries = "monitoring"))
(pc1_loadings_l226s_monitoring <- extractPCl1oadings(phytoplankton_phototrophs_annual_taxoncode %>%
                                                       filter(lake_id == "L226S" & stratum == "epilimnion"), timeseries = "monitoring"))
(pc1_loadings_l227_monitoring <- extractPCl1oadings(phytoplankton_phototrophs_annual_taxoncode %>%
                                                      filter(lake_id == "L227" & stratum == "epilimnion"), timeseries = "monitoring"))
(pc1_loadings_l224_monitoring <- extractPCl1oadings(phytoplankton_phototrophs_annual_taxoncode %>%
                                                      filter(lake_id == "L224" & stratum == "epilimnion"), timeseries = "monitoring"))
(pc1_loadings_l373_monitoring <- extractPCl1oadings(phytoplankton_phototrophs_annual_taxoncode %>%
                                                      filter(lake_id == "L373" & stratum == "epilimnion"), timeseries = "monitoring"))

(pc1_loadings_l226n_sediment <- extractPCl1oadings(phototrophs_melt %>%
                                                     filter(lake_id == "L226N"), timeseries = "sediment"))
(pc1_loadings_l226s_sediment <- extractPCl1oadings(phototrophs_melt %>%
                                                     filter(lake_id == "L226S"), timeseries = "sediment"))
(pc1_loadings_l227_sediment <- extractPCl1oadings(phototrophs_melt %>%
                                                    filter(lake_id == "L227"), timeseries = "sediment"))
(pc1_loadings_l224_sediment <- extractPCl1oadings(phototrophs_melt %>%
                                                    filter(lake_id == "L224"), timeseries = "sediment"))
(pc1_loadings_l373_sediment <- extractPCl1oadings(phototrophs_melt %>%
                                                    filter(lake_id == "L373"), timeseries = "sediment"))

(pc1_loadings_all <- bind_rows(pc1_loadings_l226n_monitoring,
                               pc1_loadings_l226s_monitoring,
                               pc1_loadings_l227_monitoring,
                               pc1_loadings_l224_monitoring,
                               pc1_loadings_l373_monitoring,
                               
                               pc1_loadings_l226n_sediment,
                               pc1_loadings_l226s_sediment,
                               pc1_loadings_l227_sediment,
                               pc1_loadings_l224_sediment,
                               pc1_loadings_l373_sediment) %>%
    mutate(tax = case_when(!is.na(taxon_code) ~ taxon_code,
                           !is.na(asv_code) ~ asv_code)) %>%
    select(timeseries, lake_id, tax,
           supergroup, division, subdivision, class, order, family, genus, species,
           label, direction) %>%
    pivot_wider(names_from = lake_id, values_from = direction) %>%
    arrange(timeseries, supergroup, division, subdivision, class, order, family, genus, species, tax))


#### Visualize PC1 taxon loadings ####
# Define function to plot top PC1 loadings
plotPC1loadings <- function(melt, timeseries) {
  if (timeseries == "sediment") {
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
    
    # Flip PC1 axis (left to right: old assemblages to new assemblages)
    pca_species <- pca_species %>%
      mutate(PC1 = -PC1)
    
    (pc1_percent_var <- round(x = (eigenvals(x = pca)[1])/(sum(eigenvals(x = pca)))*100, digits = 1))
    (pc2_percent_var <- round(x = (eigenvals(x = pca)[2])/(sum(eigenvals(x = pca)))*100, digits = 1))
    
    # PC contributions
    sum(abs(pca_species$PC1))
    
    pc1_loadings_accumulation_plot <- pca_species %>%
      mutate(direction = case_when(PC1 > 0 ~ "positive",
                                   PC1 < 0 ~ "negative")) %>%
      mutate(proportion_contrib = abs(PC1)/sum(abs(PC1))) %>%
      arrange(-proportion_contrib) %>%
      mutate(n = row_number()) %>%
      mutate(cum_contrib = cumsum(proportion_contrib)) %>%
      ggplot() +
      geom_bar(aes(x = fct_reorder(asv_code, desc(proportion_contrib)), y = cum_contrib,
                   fill = direction),
               stat = "identity") +
      geom_text(aes(x = fct_reorder(asv_code, desc(proportion_contrib)), y = cum_contrib,
                    label = n),
                alpha = 0.2) +
      scale_fill_manual(values = c("positive" = "dodgerblue", "negative" = "red")) +
      labs(fill = "Loading") +
      theme_bw() %+replace%
      theme(axis.text = element_text(colour = "black"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title.x = element_blank(),
            panel.grid = element_blank())
    
    (pc1_loadings_plot <- pca_species %>%
        left_join(taxonomy, by = "asv_code") %>%
        slice_max(abs(PC1), n = 20) %>%
        mutate(direction = case_when(PC1 > 0 ~ "positive",
                                     PC1 < 0 ~ "negative")) %>%
        mutate(label = case_when(!is.na(species) ~ str_c(asv_code, " s__", species, " (", subdivision, ")", sep = ""),
                                 is.na(species) & !is.na(genus) ~ str_c(asv_code, " g__", genus, " (", subdivision, ")", sep = ""),
                                 is.na(genus) & !is.na(family) ~ str_c(asv_code, " f__", family, " (", subdivision, ")", sep = ""),
                                 is.na(family) & !is.na(order) ~ str_c(asv_code, " o__", order, " (", subdivision, ")", sep = ""),
                                 is.na(order) & !is.na(class) ~ str_c(asv_code, " c__", class, " (", subdivision, ")", sep = ""),
                                 is.na(class) & !is.na(subdivision) ~ str_c(asv_code, " sd__", subdivision),
                                 is.na(subdivision) & !is.na(division) ~ str_c(asv_code, " d__", division),
                                 is.na(subdivision) & is.na(division) & !is.na(supergroup) ~ str_c(asv_code, " sg__", supergroup))) %>%
        ggplot() +
        geom_bar(aes(y = fct_reorder(label, PC1), x = PC1, fill = direction),
                 stat = "identity") +
        scale_fill_manual(values = c("positive" = "dodgerblue", "negative" = "red")) +
        labs(fill = "Loading") +
        ggtitle(str_c(unique(melt$lake_id), " sediment")) +
        theme_bw() %+replace%
        theme(axis.text = element_text(colour = "black"),
              #axis.text.x = element_text(angle = 90, hjust = 1),
              axis.ticks = element_line(colour = "black"),
              axis.title.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              legend.position = "bottom"))
  } else if (timeseries == "monitoring") {
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
    
    # Flip PC1 axis (left to right: old assemblages to new assemblages)
    pca_species <- pca_species %>%
      mutate(PC1 = -PC1)
    
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
    
    (pc1_percent_var <- round(x = (eigenvals(x = pca)[1])/(sum(eigenvals(x = pca)))*100, digits = 1))
    (pc2_percent_var <- round(x = (eigenvals(x = pca)[2])/(sum(eigenvals(x = pca)))*100, digits = 1))
    
    # PC contributions
    sum(abs(pca_species$PC1))
    
    pc1_loadings_accumulation_plot <- pca_species %>%
      mutate(direction = case_when(PC1 > 0 ~ "positive",
                                   PC1 < 0 ~ "negative")) %>%
      mutate(proportion_contrib = abs(PC1)/sum(abs(PC1))) %>%
      arrange(-proportion_contrib) %>%
      mutate(n = row_number()) %>%
      mutate(cum_contrib = cumsum(proportion_contrib)) %>%
      ggplot() +
      geom_bar(aes(x = fct_reorder(taxon_code, desc(proportion_contrib)), y = cum_contrib,
                   fill = direction),
               stat = "identity") +
      geom_text(aes(x = fct_reorder(taxon_code, desc(proportion_contrib)), y = cum_contrib,
                    label = n),
                alpha = 0.2) +
      scale_fill_manual(values = c("positive" = "dodgerblue", "negative" = "red")) +
      labs(fill = "Loading") +
      theme_bw() %+replace%
      theme(axis.text = element_text(colour = "black"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title.x = element_blank(),
            panel.grid = element_blank())
    
    (pc1_loadings_plot <- pca_species %>%
        rename_with(~str_remove(.x, "_pr2"), ends_with("_pr2")) %>%
        slice_max(abs(PC1), n = 20) %>%
        mutate(direction = case_when(PC1 > 0 ~ "positive",
                                     PC1 < 0 ~ "negative")) %>%
        mutate(label = case_when(!is.na(genus) ~ str_c(taxon_code, " g__", genus, " (", subdivision, ")", sep = ""),
                                 is.na(genus) & !is.na(family) ~ str_c(taxon_code, " f__", family, " (", subdivision, ")", sep = ""),
                                 is.na(family) & !is.na(order) ~ str_c(taxon_code, " o__", order, " (", subdivision, ")", sep = ""),
                                 is.na(order) & !is.na(class) ~ str_c(taxon_code, " c__", class, " (", subdivision, ")", sep = ""),
                                 is.na(class) & !is.na(subdivision) ~ str_c(taxon_code, " sd__", subdivision),
                                 is.na(subdivision) & !is.na(division) ~ str_c(taxon_code, " d__", division),
                                 is.na(subdivision) & is.na(division) & !is.na(supergroup) ~ str_c(taxon_code, " sg__", supergroup))) %>%
        ggplot() +
        geom_bar(aes(y = fct_reorder(label, PC1), x = PC1, fill = direction),
                 stat = "identity") +
        scale_fill_manual(values = c("positive" = "dodgerblue", "negative" = "red")) +
        labs(fill = "Loading") +
        ggtitle(str_c(unique(melt$lake_id), " monitoring")) +
        theme_bw() %+replace%
        theme(axis.text = element_text(colour = "black"),
              #axis.text.x = element_text(angle = 90, hjust = 1),
              axis.ticks = element_line(colour = "black"),
              axis.title.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              legend.position = "bottom"))
  }
}

# Plot PC1 loadings for monitoring assemblages in individual lakes
(pc1_loadings_plot_l226n_monitoring <- plotPC1loadings(phytoplankton_phototrophs_annual_taxoncode %>%
                                                         filter(lake_id == "L226N" & stratum == "epilimnion"), timeseries = "monitoring"))
(pc1_loadings_plot_l226s_monitoring <- plotPC1loadings(phytoplankton_phototrophs_annual_taxoncode %>%
                                                         filter(lake_id == "L226S" & stratum == "epilimnion"), timeseries = "monitoring"))
(pc1_loadings_plot_l227_monitoring <- plotPC1loadings(phytoplankton_phototrophs_annual_taxoncode %>%
                                                        filter(lake_id == "L227" & stratum == "epilimnion"), timeseries = "monitoring"))
(pc1_loadings_plot_l224_monitoring <- plotPC1loadings(phytoplankton_phototrophs_annual_taxoncode %>%
                                                        filter(lake_id == "L224" & stratum == "epilimnion"), timeseries = "monitoring"))
(pc1_loadings_plot_l373_monitoring <- plotPC1loadings(phytoplankton_phototrophs_annual_taxoncode %>%
                                                        filter(lake_id == "L373" & stratum == "epilimnion"), timeseries = "monitoring"))

# Plot PC1 loadings for sediment assemblages in individual lakes
(pc1_loadings_plot_l226n_sediment <- plotPC1loadings(phototrophs_melt %>%
                                                       filter(lake_id == "L226N"), timeseries = "sediment"))
(pc1_loadings_plot_l226s_sediment <- plotPC1loadings(phototrophs_melt %>%
                                                       filter(lake_id == "L226S"), timeseries = "sediment"))
(pc1_loadings_plot_l227_sediment <- plotPC1loadings(phototrophs_melt %>%
                                                      filter(lake_id == "L227"), timeseries = "sediment"))
(pc1_loadings_plot_l224_sediment <- plotPC1loadings(phototrophs_melt %>%
                                                      filter(lake_id == "L224"), timeseries = "sediment"))
(pc1_loadings_plot_l373_sediment <- plotPC1loadings(phototrophs_melt %>%
                                                      filter(lake_id == "L373"), timeseries = "sediment"))

# Combine plots
(pc1_loadings_plots_all <- pc1_loadings_plot_l226n_monitoring +
    pc1_loadings_plot_l226n_sediment +
    pc1_loadings_plot_l226s_monitoring +
    pc1_loadings_plot_l226s_sediment +
    pc1_loadings_plot_l227_monitoring +
    pc1_loadings_plot_l227_sediment +
    pc1_loadings_plot_l224_monitoring +
    pc1_loadings_plot_l224_sediment +
    pc1_loadings_plot_l373_monitoring +
    pc1_loadings_plot_l373_sediment +
    plot_layout(guides = "collect", ncol = 2) &
    theme(legend.position = "bottom"))
#ggsave("figures/ela_pc1_loadings_all.pdf", pc1_loadings_plots_all, width = 14, height = 20, units = "in")
