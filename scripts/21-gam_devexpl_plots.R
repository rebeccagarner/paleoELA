# Generalized Additive Model plots

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
# Import GAM summaries
gams_nutriclim_phototrophs_taxoncode_epilimnion_pc_all <- read_csv("output/gams/ela_phototrophs_taxoncode_epilimnion_pc_gams_nutriclim.csv", col_names = TRUE)
gams_nutriclim_phototrophs_asvs_pc_all <- read_csv("output/gams/ela18s_phototrophs_asvs_pc_gams_nutriclim.csv", col_names = TRUE)
gams_nutriclim_phototrophs_sedmonitor_asvs_pc_all <- read_csv("output/gams/ela18s_phototrophs_sedmonitor_asvs_pc_gams_nutriclim.csv", col_names = TRUE) %>%
  mutate(dataset = str_c(dataset, " monitoring"),
         dataset_response_var = str_c(dataset_response_var, " monitoring")) %>%
  filter(explanatory_category == "climate")

# Combine GAM summaries from phytoplankton and sediment DNA response data
gams_nutriclim_phototrophs_pc_all <- bind_rows(gams_nutriclim_phototrophs_taxoncode_epilimnion_pc_all,
                                               gams_nutriclim_phototrophs_sedmonitor_asvs_pc_all,
                                               gams_nutriclim_phototrophs_asvs_pc_all) %>%
  filter(response_var == "PC1") %>%
  mutate(label = case_when(grepl("taxon_codes", dataset_response_var) ~ "1monitoring",
                           grepl("asv", dataset_response_var, ignore.case = TRUE) & grepl("monitoring", dataset_response_var) ~ "2paleo in monitoring",
                           grepl("asv", dataset_response_var, ignore.case = TRUE) ~ "3paleo full"))  %>%
  mutate(label = case_when(dataset == "Phototroph ASVs" &
                             label == "3paleo full" &
                             grepl("nutrients", explanatory_category) ~ "2paleo in monitoring",
                           TRUE ~ label))


#### Visualize GAMs ####
(gams_nutriclim_phototrophs_pc_barplot <- gams_nutriclim_phototrophs_pc_all %>%
   ggplot() +
   facet_grid(factor(explanatory_category, levels = c("nutrients", "climate", "nutrients.climate"))~label,
              switch = "y") +
   geom_bar(aes(x = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
                y = gam_dev_expl * 100,
                fill = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))),
            position = position_dodge(preserve = "single"), stat = "identity") +
   scale_fill_manual(values = palette_lake, name = "Lake") +
   labs(y = "Deviance explained (%)") +
   ylim(c(0, 100)) +
   scale_y_continuous(position = "right") +
   theme_bw() %+replace%
   theme(axis.title.x = element_blank(),
         axis.text = element_text(colour = "black"),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         strip.background = element_blank(),
         strip.placement = "outside",
         panel.grid = element_blank(),
         panel.background = element_blank(),
         legend.position = "right"))

(gams_nutriclim_phototrophs_pc_heatmap <- gams_nutriclim_phototrophs_pc_all %>%
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
    facet_grid(factor(explanatory_category, levels = c("nutrients", "climate", "nutrients.climate"))~label) +
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
          strip.text = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank()))

(gams_nutriclim_phototrophs_pc_plots <- gams_nutriclim_phototrophs_pc_barplot /
    gams_nutriclim_phototrophs_pc_heatmap +
    plot_layout(heights = c(2.5, 1)))
#ggsave("figures/ela_pc1_gams_nutriclim.pdf", gams_nutriclim_phototrophs_pc_plots, width = 12, height = 12, units = "in", device = cairo_pdf)
