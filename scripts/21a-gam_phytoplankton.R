# Generalized Additive Models of phytoplankton monitoring
# https://r.qcbs.ca/workshops/r-workshop-08/

setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
setwd("~/Dropbox/projects/ela18s/")
setwd("~/Library/CloudStorage/Dropbox-Personal/projects/ela18s/")

# Load GAM functions
source("scripts/21-gam_functions_phytoplankton.R")

# Load sediment data
source("scripts/03e-sediment_data.R")

# Load environmental data
source("scripts/08-ela_env_data.R")

# Load monitoring data
source("scripts/11-ela_monitoring.R")


#### Import and format data ####
# Import PCA axes 1 and 2 data
pca <- read_tsv("output/ordinations/ela_phytoplankton_pca_axes.tsv", col_names = TRUE)


#### Compute Generalized Additive Models (GAMs) ####
# Plot relationships between response and explanatory variables
# (var_scatterplot_l226n_pc1_phototrophs_epilimnion <- plotVars("L226N", pca %>%
#                                                                  filter(dataset == "Phototroph taxon_codes epilimnion"), "PC1"))
# (var_scatterplot_l226s_pc1_phototrophs_epilimnion <- plotVars("L226S", pca %>%
#                                                                  filter(dataset == "Phototroph taxon_codes epilimnion"), "PC1"))
# (var_scatterplot_l227_pc1_phototrophs_epilimnion <- plotVars("L227", pca %>%
#                                                                 filter(dataset == "Phototroph taxon_codes epilimnion"), "PC1"))
# (var_scatterplot_l224_pc1_phototrophs_epilimnion <- plotVars("L224", pca %>%
#                                                                 filter(dataset == "Phototroph taxon_codes epilimnion"), "PC1"))
# (var_scatterplot_l373_pc1_phototrophs_epilimnion <- plotVars("L373", pca %>%
#                                                                 filter(dataset == "Phototroph taxon_codes epilimnion"), "PC1"))
# 
# (var_scatterplot_pc1_phototrophs_epilimnion_all <- var_scatterplot_l226n_pc1_phototrophs_epilimnion /
#     var_scatterplot_l226s_pc1_phototrophs_epilimnion /
#     var_scatterplot_l227_pc1_phototrophs_epilimnion /
#     var_scatterplot_l224_pc1_phototrophs_epilimnion /
#     var_scatterplot_l373_pc1_phototrophs_epilimnion)
# #ggsave("figures/ela_phytoplankton_pc1_phototrophs_epilimnion_vars.pdf", var_scatterplot_pc1_phototrophs_epilimnion_all, width = 8, height = 36, units = "in", device = cairo_pdf)

# Clear up environment
rm(list = ls(pattern = "^var_scatterplot_"))


#### Fit GAMs based on nutrients, climate, and nutrients + climate ####
# Compute GAMs based on nutrients, climate, and nutrients + climate
gams_nutriclim_phototrophs_taxoncode_epilimnion_pc1 <- computeGAMNutriClimAll(pca %>%
                                                                                filter(dataset == "Phototroph taxon_codes epilimnion"), "PC1")

# Visualize deviance explained by GAMs based on explanatory variables across environmental categories
gams_nutriclim_phototrophs_taxoncode_epilimnion_pc_all <- gams_nutriclim_phototrophs_taxoncode_epilimnion_pc1 %>%
  mutate(dataset_response_var = str_c(dataset, response_var, sep = " "))

# gams_nutriclim_phototrophs_taxoncode_epilimnion_pc_all %>%
#   write_csv("output/gams/ela_phototrophs_taxoncode_epilimnion_pc_gams_nutriclim.csv", col_names = TRUE)

gams_nutriclim_phototrophs_taxoncode_epilimnion_pc_barplot <- gams_nutriclim_phototrophs_taxoncode_epilimnion_pc_all %>%
  ggplot() +
  facet_grid(dataset_response_var~factor(explanatory_category, levels = c("nutrients", "climate", "nutrients.climate"))) +
  geom_bar(aes(x = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               y = gam_dev_expl * 100,
               fill = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))),
           position = position_dodge(preserve = "single"), stat = "identity") +
  scale_fill_manual(values = palette_lake, name = "Lake") +
  labs(y = "Deviance explained (%)") +
  ylim(c(0, 100)) +
  theme_bw() %+replace%
  theme(axis.title.x = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank())

gams_nutriclim_phototrophs_taxoncode_epilimnion_pc_heatmap <- gams_nutriclim_phototrophs_taxoncode_epilimnion_pc_all %>%
  mutate(gam_var_signif = strsplit(gam_vars_p_value, ", ")) %>%
  unnest(gam_var_signif) %>%
  separate(gam_var_signif, c("gam_var", "var_p_value"), "\\|") %>%
  mutate(var_p_value = as.numeric(var_p_value)) %>%
  mutate(var_signif = case_when(var_p_value >= 0.05 ~ "",
                                var_p_value < 0.05 & var_p_value > 0.01 ~ "*",
                                var_p_value <= 0.01 & var_p_value > 0.001 ~ "**",
                                var_p_value <= 0.001 ~ "***")) %>%
  ggplot(aes(x = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
             y = factor(gam_var, levels = rev(c(vars_nutrients, vars_climate))))) +
  facet_grid(response_var~factor(explanatory_category, levels = c("nutrients", "climate", "nutrients.climate"))) +
  geom_tile(fill = "#E0E0E0", colour = "white") +
  geom_text(aes(label = var_signif)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() %+replace%
  theme(axis.title = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank())

(gams_nutriclim_phototrophs_taxoncode_epilimnion_pc_plots <- gams_nutriclim_phototrophs_taxoncode_epilimnion_pc_barplot /
    gams_nutriclim_phototrophs_taxoncode_epilimnion_pc_heatmap +
    plot_layout(heights = c(3, 2)))
#ggsave("figures/ela_phototrophs_taxoncode_epilimnion_pc_gams_nutriclim.pdf", gams_nutriclim_phototrophs_taxoncode_epilimnion_pc_plots, width = 8, height = 7, units = "in", device = cairo_pdf)

# gams_nutriclim_phototrophs_taxoncode_epilimnion_pc_all %>%
#   mutate(gam_var_signif = strsplit(gam_vars_p_value, ", ")) %>%
#   unnest(gam_var_signif) %>%
#   separate(gam_var_signif, c("gam_var", "var_p_value"), "\\|") %>%
#   mutate(var_p_value = as.numeric(var_p_value)) %>%
#   mutate(var_signif = case_when(var_p_value >= 0.05 ~ "",
#                                 var_p_value < 0.05 & var_p_value > 0.01 ~ "*",
#                                 var_p_value <= 0.01 & var_p_value > 0.001 ~ "**",
#                                 var_p_value <= 0.001 ~ "***")) %>%
#   write_csv("output/gams/ela_phototrophs_taxoncode_epilimnion_pc_gams_nutriclim_vars.csv", col_names = TRUE)

# Combine contribution plots
(gams_nutriclim_phototrophs_taxoncode_epilimnion_pc1_contribplot <- plotGAMNutriClimContrib(pca %>%
                                                                                              filter(dataset == "Phototroph taxon_codes epilimnion"), "PC1"))
#ggsave("figures/ela_phototrophs_taxoncode_epilimnion_pc1_gams_nutriclim_contribplot.pdf", gams_nutriclim_phototrophs_taxoncode_epilimnion_pc1_contribplot, width = 16, height = 12, units = "in", device = cairo_pdf)
