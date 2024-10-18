# Visualize nutrient and climate monitoring

# Load libraries
library(tidyverse)
library(lubridate)
library(ggforce)
library(ggnewscale)
library(patchwork)
library(scales)

# Load palettes
source("scripts/00-palettes.R")

# Load monitoring data
source("scripts/11-ela_monitoring.R")


#### Format environmental monitoring data ####
# Join all annual environmental variable summaries
env_byyear_all <- sdd_byyear_interpolated %>%
  full_join(watertemp_byyear_interpolated, by = c("lake_id", "year")) %>%
  full_join(chemistry_byyear_interpolated, by = c("lake_id", "year")) %>%
  full_join(climate_eccc_byyear %>%
              rename(precip_eccc = precip), by = c("lake_id", "year")) %>%
  full_join(climate_ela_byyear %>%
              rename(precip_ela = precip), by = c("lake_id", "year"))

# Format climate (air temperature) data
climate <- env_byyear_all %>%
  select(lake_id, year, meanairtemp) %>%
  distinct(year, meanairtemp) %>%
  arrange(year)

# Format nutrient data
nutrients <- env_byyear_all %>%
  select(lake_id, year, epilimnion_tp, epilimnion_tn) %>%
  filter(!is.na(epilimnion_tp) | !is.na(epilimnion_tn))


#### Plot climate and nutrient time series ####
# Plot climate (air temperature) time series
(climate_plot <- climate %>%
   mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
   ggplot(aes(x = year_floor, y = meanairtemp)) +
   geom_line(colour = "steelblue", alpha = 0.7) +
   geom_smooth(colour = NA, alpha = 0.2) +
   geom_point(aes(fill = meanairtemp),
              colour = "black", size = 2.6, pch = 21, stroke = 0.1) +
   geom_smooth(colour = "steelblue", fill = NA) +
   labs(y = "Air temperature (\u00B0C)") +
   scale_fill_distiller(palette = "RdBu",
                        direction = -1,
                        limits = c(min(climate$meanairtemp), max(climate$meanairtemp)),
                        name = "Temp (\u00B0C)") +
   scale_x_date(limits = as.Date(c("1890-01-01", "2029-12-31")),
                date_breaks = "10 years", date_labels = "%Y",
                expand = c(0,0)) +
   scale_y_continuous(limits = c(-1, 6)) +
   theme_bw() %+replace%
   theme(axis.text = element_text(colour = "black"),
         axis.title.x = element_blank(),
         panel.grid = element_blank()))

# Plot TP and TN with two axes
max(nutrients$epilimnion_tp)/max(nutrients$epilimnion_tn)
axis2_coeff <- 50

(nutrients_plot <- nutrients %>%
    bind_rows(climate) %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_wrap(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               ncol = 1) +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations) +
    scale_fill_manual(values = palette_manipulation) +
    new_scale_fill() +
    geom_line(aes(x = year_floor, y = epilimnion_tp), colour = "green") +
    geom_line(aes(x = year_floor, y = epilimnion_tn*axis2_coeff), colour = "orange") +
    geom_point(aes(x = year_floor, y = epilimnion_tp, fill = epilimnion_tp),
               colour = "black", size = 2.6, pch = 21, stroke = 0.1) +
    scale_fill_distiller(palette = "Greens",
                         direction = 1,
                         trans = "log10",
                         limits = c(min(nutrients$epilimnion_tp), max(nutrients$epilimnion_tp)),
                         labels = comma,
                         name = "TP (\u00b5g/L)") +
    new_scale_fill() +
    geom_point(aes(x = year_floor, y = epilimnion_tn*axis2_coeff, fill = epilimnion_tn),
               colour = "black", size = 2.6, pch = 21, stroke = 0.1) +
    scale_fill_distiller(palette = "Oranges",
                         direction = 1,
                         trans = "sqrt",
                         limits = c(min(nutrients$epilimnion_tn), max(nutrients$epilimnion_tn)),
                         #breaks = c(0, 1, 2),
                         name = "TN (mg/L)") +
    scale_y_continuous(limits = c(0,80),
                       name = "TP (ug/L)",
                       sec.axis = sec_axis(~./axis2_coeff, name = "TN (mg/L)")) +
    scale_x_date(limits = as.Date(c("1890-01-01", "2029-12-31")),
                 date_breaks = "10 years", date_labels = "%Y",
                 expand = c(0,0)) +
    expand_limits(y = 0) +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "#00B0F0")))

# Combine climate and nutrient plots
(nutriclim_plots <- climate_plot +
    nutrients_plot +
    plot_layout(ncol = 1, heights = c(1, 6)))
#ggsave("figures/ela_nutriclim.pdf", nutriclim_plots, width = 16, height = 10, units = "in", device = cairo_pdf)
