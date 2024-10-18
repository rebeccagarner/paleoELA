# Correlation/s between sediment and water column monitoring assemblages

# Load libraries
library(tidyverse)
library(vegan)
library(FactoMineR)
library(factoextra)

# Load palettes
source("scripts/00-palettes.R")

# Load sediment data
source("scripts/03e-sediment_data.R")

# Load monitoring data
source("scripts/11-ela_monitoring.R")


#### Calculate RV coefficients ####
# Define function to calculate RV coefficient between sediment and monitoring assemblages
calculateRV <- function(melt, taxon1, x_data = "Matrix1",
                        monitoring_data, taxon2, y_data = "Matrix2",
                        lakeid, strat) {
  # Format X (sediment) and Y (monitoring) matrices
  sediment <- formatLayerByTaxon(melt %>%
                                   filter(lake_id == lakeid), taxon1)
  
  monitoring <- syncPhytoplankton(monitoring_data, taxon2, lakeid, strat) %>%
    select(-nyear)
  
  sediment_df <- sediment %>%
    filter(layer_order %in% unique(monitoring$layer_sync)) %>%
    column_to_rownames("layer_order")
  
  monitoring_df <- monitoring %>%
    column_to_rownames("layer_sync")
  
  # Transform assemblage data
  sediment_hellinger <- decostand(sediment_df, "hellinger")
  monitoring_hellinger <- decostand(monitoring_df, "hellinger")
  
  # Calculate RV coefficient
  comp_rv <- coeffRV(sediment_hellinger, monitoring_hellinger)
  rv <- comp_rv$rv
  nobs <- nrow(sediment_hellinger)
  pvalue <- comp_rv$p.value
  
  rv_table <- tibble(sediment_data = x_data,
                     monitoring_data = y_data,
                     lake_id = lakeid,
                     stratum = strat,
                     rv_coefficient = rv,
                     n = nobs,
                     p_value = pvalue) %>%
    mutate(signif_level = case_when(p_value >= 0.05 ~ "",
                                    p_value < 0.05 & p_value > 0.01 ~ "*",
                                    p_value <= 0.01 & p_value > 0.001 ~ "**",
                                    p_value <= 0.001 ~ "***"))
  return(rv_table)
}

# Define function to calculate RV coefficients for all lakes and strata
calculateRVAll <- function(melt, taxon1, x_data = "Matrix1",
                           monitoring_data, taxon2, y_data = "Matrix2") {
  rv_all <- bind_rows(calculateRV(melt, taxon1, x_data,
                                  monitoring_data, taxon2, y_data,
                                  "L226N", "epilimnion"),
                      
                      calculateRV(melt, taxon1, x_data,
                                  monitoring_data, taxon2, y_data,
                                  "L226S", "epilimnion"),
                      
                      calculateRV(melt, taxon1, x_data,
                                  monitoring_data, taxon2, y_data,
                                  "L227", "epilimnion"),
                      
                      calculateRV(melt, taxon1, x_data,
                                  monitoring_data, taxon2, y_data,
                                  "L224", "epilimnion"),
                      
                      calculateRV(melt, taxon1, x_data,
                                  monitoring_data, taxon2, y_data,
                                  "L373", "epilimnion"))
  
  return(rv_all)
}

(rv_all <- bind_rows(
  calculateRVAll(phototrophs_melt, "asv_code", "Phototroph ASVs",
                 phytoplankton_phototrophs, "taxon_code", "Phytoplankton (phototrophs) taxon codes"),
  
  calculateRVAll(phototrophs_melt, "genus","Phototroph genera",
                 phytoplankton_phototrophs, "genus_pr2", "Phytoplankton (phototrophs) genera")))

# Write RV coefficients to file
# rv_all %>%
#   write_tsv("output/rv/ela18s_rv_coefficients.tsv", col_names = TRUE)

# Plot RV coefficients for different compared matrices
rv_all %>%
  mutate(comparison = str_c(sediment_data, " vs. ", monitoring_data)) %>%
  ggplot(aes(x = comparison, y = rv_coefficient)) +
  facet_wrap(~stratum) +
  geom_bar(aes(fill = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))),
           position = "dodge", stat = "identity") +
  geom_text(aes(label = signif_level,
                group = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))),
            position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = palette_lake, name = "Lake") +
  labs(y = "RV coefficient") +
  ylim(c(0, 1)) +
  theme_bw() %+replace%
  theme(axis.title.x = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank())

rv_all %>%
    filter(grepl("ASV", sediment_data) & grepl("taxon", monitoring_data)) %>%
    mutate(comparison = str_c(sediment_data, " vs. ", monitoring_data)) %>%
    ggplot(aes(x = rv_coefficient, y = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")))) +
    facet_wrap(~stratum, ncol = 1, strip.position = "left") +
    geom_bar(aes(fill = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))),
             stat = "identity") +
    geom_text(aes(label = signif_level,
                  group = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))),
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = palette_lake, name = "Lake") +
    scale_x_continuous(position = "bottom", limits = c(0, 1)) +
    scale_y_discrete(limits = rev(c("L226N", "L226S", "L227", "L224", "L373"))) +
    labs(x = "RV coefficient") +
    theme_bw() %+replace%
    theme(axis.title.y = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.background = element_blank()))

(rv_plot <- rv_all %>%
    filter(grepl("ASV", sediment_data) & grepl("taxon", monitoring_data)) %>%
    filter(stratum == "epilimnion") %>%
    ggplot(aes(x = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               y = rv_coefficient)) +
    geom_bar(aes(fill = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))),
             stat = "identity") +
    geom_text(aes(label = signif_level,
                  group = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))),
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = palette_lake, name = "Lake") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(y = "RV coefficient") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position = "none"))
#ggsave("figures/ela18s_rv_coefficients.pdf", rv_plot, width = 5, height = 5, units = "in", device = cairo_pdf)
