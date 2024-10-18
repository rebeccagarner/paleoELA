# PC1 plot

# Load sediment data
source("scripts/03e-sediment_data.R")

# Load environmental data
source("scripts/08-ela_env_data.R")

# Load monitoring data
source("scripts/11-ela_monitoring.R")


#### Import and format data ####
# Import PCA axes 1 and 2 data
pca_phytoplankton <- read_tsv("output/ordinations/ela_phytoplankton_pca_axes.tsv", col_names = TRUE) %>%
  filter(dataset == "Phototroph taxon_codes epilimnion")
pca_sediment <- read_tsv("output/ordinations/ela18s_pca_axes.tsv", col_names = TRUE) %>%
  filter(dataset == "Phototroph ASVs")


#### Plot PC1 ####
# Combine PCA data from
pca_all <- pca_phytoplankton %>%
  bind_rows(pca_sediment) %>%
  mutate(timeseries = case_when(dataset == "Phototroph taxon_codes epilimnion" ~ "monitoring",
                                dataset == "Phototroph ASVs" ~ "sediment")) %>%
  mutate(year = case_when(timeseries == "monitoring" ~ year,
                          timeseries == "sediment" ~ midpt_year))

# Plot PC1 (point colour is year)
(pc1_plot_rainbow_timeseries <- pca_all %>%
    mutate(year_label = as.character(round(year))) %>%
    mutate(year_label = substr(year_label, 3, 4)) %>%
    ggplot(aes(x = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               y = PC1)) +
    facet_wrap(~timeseries, scales = "free_x",
               ncol = 1, strip.position = "right") +
    # geom_text(aes(label = year_label),
    #           nudge_x = 0.2, size = 1) +
    geom_point(aes(fill = year),
               colour = "black", pch = 21,
               alpha = 0.7, stroke = 0.1) +
    scale_fill_gradientn(colours = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"),
                         na.value = "#EBEBEB",
                         trans = "sqrt",
                         breaks = seq(1880, 2020, 20),
                         name = "Year") +
    scale_x_discrete(limits = rev) +
    coord_flip() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_line(linetype = "dashed", color = "#E6E6E6"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(colour = "black"),
          legend.key = element_blank(),
          legend.position = "bottom"))
#ggsave("figures/ela_pc1.pdf", pc1_plot_rainbow_timeseries, width = 7, height = 4, units = "in", device = cairo_pdf)

# Plot time intervals (point colour is PC1 loading)
(pc1_plot_timeseries <- pca_all %>%
    mutate(year_label = as.character(round(year))) %>%
    mutate(year_label = substr(year_label, 3, 4)) %>%
    ggplot(aes(x = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               y = year)) +
    facet_wrap(~timeseries, scales = "free_x",
               ncol = 1, strip.position = "right") +
    geom_point(aes(fill = PC1),
               colour = "black", pch = 21,
               size = 3, alpha = 1, stroke = 0.1) +
    scale_fill_gradient2(low = "red", mid = "white", high = "dodgerblue",
                         midpoint = 0) +
    scale_x_discrete(limits = rev) +
    scale_y_continuous(breaks = seq(1880, 2020, 10)) +
    coord_flip() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_line(linetype = "dashed", color = "#E6E6E6"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(colour = "black"),
          legend.key = element_blank(),
          legend.position = "bottom"))
#ggsave("figures/ela_pc1_years.pdf", pc1_plot_timeseries, width = 7, height = 4, units = "in", device = cairo_pdf)
