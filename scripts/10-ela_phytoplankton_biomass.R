# ELA phytoplankton biomass data

setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
setwd("~/Dropbox/projects/ela18s/")

# Load libraries
library(tidyverse)
library(lubridate)
library(ggnewscale)
library(patchwork)
library(ggforce)
library(scales)

# Load palettes
source("scripts/00-palettes.R")


#### Import and format data ####
# Import phytoplankton biomass data
# There is only ever a single observation for each lake, date, stratum, and taxon
phytoplankton_biomass <- read_tsv("output/ela/ela_monitoring_phytoplankton_biomass.tsv", col_names = TRUE)

# Format long-form matrix for phytoplankton biomass
phytoplankton_melt <- phytoplankton_biomass %>%
  pivot_longer(!c(lake_id, date, stratum, start_depth, end_depth),
               names_to = "taxon", values_to = "biomass") %>%
  filter(!is.na(biomass))

# Import and format experimental manipulation data
ela_manipulations <- read_tsv("data/ela/ela_manipulations.txt", col_names = TRUE, comment = "#")


#### Format seasons ####
# Define function to assign season from date
assignSeason <- function(date) {
  month_date_numeric <- 100 * month(date) + day(date)
  season_cuts <- cut(month_date_numeric, breaks = c(0, 319, 620, 921, 1220, 1231))
  levels(season_cuts) <- c("Winter","Spring","Summer","Fall","Winter")
  return(season_cuts)
}

# Assign season to sampling dates
phytoplankton_melt <- phytoplankton_melt %>%
  mutate(season = assignSeason(date))

# Summarize season start and end dates of phytoplankton monitoring
phytoplankton_melt %>%
  filter(taxon != "total_phytoplankton") %>%
  mutate(year = year(date)) %>%
  group_by(lake_id, stratum, year, season) %>%
  summarize(season_start_date = min(date),
            season_end_date = max(date)) %>%
  ungroup()


#### Visualize phytoplankton biomass time series ####
# Define function to plot biomass time series
plotTimeSeries <- function(melt) {
  # Extract sampling start and end dates
  start_date <- min(melt$date)
  end_date <- max(melt$date)
  
  # Create vector of sampling years
  sampling_years <- seq(from = year(start_date), to = year(end_date), by = 1)
  
  # Format seasons
  season_dates <- tibble(year = rep(sampling_years, 5)) %>%
    arrange(year) %>%
    mutate(season = rep(c("Winter", "Spring", "Summer", "Fall", "Winter"), length(sampling_years)),
           season_start_month_day = rep(c("01-01", "03-20", "06-21", "09-22", "12-21"), length(sampling_years)),
           season_end_month_day = rep(c("03-19", "06-20", "09-21", "12-20", "12-31"), length(sampling_years))) %>%
    mutate(season_start_date = date(str_c(year, "-", season_start_month_day)),
           season_end_date = date(str_c(year, "-", season_end_month_day))) %>%
    filter(season_start_date >= start_date,
           season_end_date <= end_date)
  
  season_dates <- bind_rows(season_dates %>%
                              mutate(stratum = "Epilimnion"),
                            season_dates %>%
                              mutate(stratum = "Metalimnion"))
  
  # Remove total phytoplankton biomass data from melted data set
  melt_tmp <- melt %>%
    filter(taxon != "total_phytoplankton")
  
  # Create a separate melted data set for total phytoplankton biomass
  total_phytoplankton_melt <- melt %>%
    filter(taxon == "total_phytoplankton")
  
  # Plot time series
  (timeseries_plot <- melt_tmp %>%
      ggplot() +
      facet_grid(stratum~lake_id, space = "free_y", scales = "free_y") +
      geom_rect(aes(xmin = season_start_date, xmax = season_end_date,
                    ymin = 0, ymax = Inf,
                    fill = season),
                data = season_dates,
                alpha = 0.1) +
      scale_fill_manual(values = palette_season, name = "Season") +
      new_scale_fill() +
      geom_area(aes(x = date, y = biomass,
                    fill = factor(taxon, levels = c("Cyanobacteria", "Chlorophytes", "Chrysophytes", "Dinoflagellates", "Cryptophytes",
                                                    "Diatoms", "Euglenophytes"))),
                stat = "identity") +
      geom_point(aes(x = date, y = biomass),
                 data = total_phytoplankton_melt,
                 size = 0.1) +
      geom_line(aes(x = date, y = biomass),
                data = total_phytoplankton_melt,
                size = 0.1) +
      scale_x_date(date_breaks = "1 year", date_labels = "%y",
                   #limits = c(start_date, end_date),
                   expand = c(0,0)) +
      scale_y_continuous(labels = scales::comma) +
      scale_fill_manual(values = palette_phytoplankton[which(names(palette_phytoplankton) != "total_phytoplankton")]) +
      labs(x = "Date",
           y = "Biomass (\u00b5g/L)",
           fill = "Taxon") +
      theme_bw() %+replace%
      theme(axis.text = element_text(colour = "black"),
            #axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            axis.title.x = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            #panel.spacing = unit(2.3, "lines"),
            axis.line = element_line(colour = "black"),
            plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))
}

# Plot phytoplankton biomass time series
(phytoplankton_timeseries_l224 <- plotTimeSeries(phytoplankton_melt %>%
                                                   filter(lake_id == "L224")))
(phytoplankton_timeseries_l226n <- plotTimeSeries(phytoplankton_melt %>%
                                                    filter(lake_id == "L226N")))
(phytoplankton_timeseries_l226s <- plotTimeSeries(phytoplankton_melt %>%
                                                    filter(lake_id == "L226S")))
(phytoplankton_timeseries_l227 <- plotTimeSeries(phytoplankton_melt %>%
                                                   filter(lake_id == "L227")))
(phytoplankton_timeseries_l373 <- plotTimeSeries(phytoplankton_melt %>%
                                                   filter(lake_id == "L373")))
# ggsave("figures/ela_phytoplankton_biomass_timeseries_l224.pdf", phytoplankton_timeseries_l224, width = 16, height = 8, units = "in", device = cairo_pdf)
# ggsave("figures/ela_phytoplankton_biomass_timeseries_l226n.pdf", phytoplankton_timeseries_l226n, width = 16, height = 8, units = "in", device = cairo_pdf)
# ggsave("figures/ela_phytoplankton_biomass_timeseries_l226s.pdf", phytoplankton_timeseries_l226s, width = 16, height = 8, units = "in", device = cairo_pdf)
# ggsave("figures/ela_phytoplankton_biomass_timeseries_l227.pdf", phytoplankton_timeseries_l227, width = 16, height = 8, units = "in", device = cairo_pdf)
# ggsave("figures/ela_phytoplankton_biomass_timeseries_l373.pdf", phytoplankton_timeseries_l373, width = 16, height = 8, units = "in", device = cairo_pdf)
# 

#### Calculate time-weighted means for phytoplankton biomass ####
# Summarize data by month for water column strata (based on summarizeByMonthWC)
phytoplankton_month_means <- phytoplankton_melt %>%
  mutate(year = year(date),
         month = month(date)) %>%
  group_by(lake_id, year, month, stratum, taxon) %>%
  summarize(mean = mean(biomass),
            sd = sd(biomass),
            ndays = n()) %>%
  ungroup() %>%
  arrange(lake_id, year, month,
          factor(stratum, levels = c("Epilimnion", "Metalimnion")),
          taxon)

# Summarize data by month-season for water column strata (based on summarizeByMonthSeasonWC)
phytoplankton_monthseason_means <- phytoplankton_month_means %>%
  mutate(month_season = case_when(month == 1 |
                                    month == 2 |
                                    month == 3 ~ "Winter",
                                  month == 4 |
                                    month == 5 |
                                    month == 6 ~ "Spring",
                                  month == 7 |
                                    month == 8 |
                                    month == 9 ~ "Summer",
                                  month == 10 |
                                    month == 11 |
                                    month == 12 ~ "Fall")) %>%
  group_by(lake_id, year, month_season, stratum, taxon) %>%
  summarize(mean = mean(mean),
            nmonths = n(),
            ndays = sum(ndays)) %>%
  ungroup() %>%
  arrange(lake_id, year, factor(month_season, levels = names(palette_season)),
          factor(stratum, levels = c("Epilimnion", "Metalimnion")),
          taxon)

# # Summarize data by year for water column strata (based on summarizeByYearWC)
# phytoplankton_year_means <- phytoplankton_monthseason_means %>%
#   filter(month_season %in% c("Spring", "Summer", "Fall")) %>%  # Remove winter sampling dates
#   group_by(lake_id, year, stratum, taxon) %>%
#   summarize(mean = mean(mean),
#             nseasons = n(),
#             nmonths = sum(nmonths),
#             ndays = sum(ndays)) %>%
#   ungroup() %>%
#   arrange(lake_id, year,
#           factor(stratum, levels = c("Epilimnion", "Metalimnion")),
#           taxon)

# Define functions to interpolate variable mean for missing month-seasons
interpolateMonthSeasonsTaxonWC <- function(data, names, method = "approx") {
  data_tmp <- data %>%
    arrange(year)
  
  lakeid_tmp <- unique(data_tmp$lake_id)
  monthseason_tmp <- unique(data_tmp$month_season)
  stratum_tmp <- unique(data_tmp$stratum)
  taxon_tmp <- unique(data_tmp$taxon)
  
  print(paste0(lakeid_tmp, "_", monthseason_tmp, "_", stratum_tmp, "_", taxon_tmp))
  
  year_range <- seq(data_tmp$year[1], tail(data_tmp$year, n = 1), 1)
  
  if (length(year_range) > 1) {
    if (method == "approx") {
      interpolated <- approx(data_tmp$year, data_tmp$mean, year_range) %>%
        as_tibble() %>%
        rename(year = x,
               mean = y) %>%
        filter(!year %in% unique(data_tmp$year)) %>%
        mutate(lake_id = lakeid_tmp,
               month_season = monthseason_tmp,
               stratum = stratum_tmp,
               taxon = taxon_tmp)
    } else if (method == "spline") {
      interpolated <- spline(data_tmp$year, data_tmp$mean,
                             xout = year_range,
                             method = "fmm") %>%
        as_tibble() %>%
        rename(year = x,
               mean = y) %>%
        filter(!year %in% unique(data_tmp$year)) %>%
        mutate(lake_id = lakeid_tmp,
               month_season = monthseason_tmp,
               stratum = stratum_tmp,
               taxon = taxon_tmp)
    }
    
    data_tmp %>%
      bind_rows(interpolated) %>%
      arrange(year) %>%
      mutate(interpolated = case_when(is.na(nmonths) ~ TRUE,
                                      TRUE ~ FALSE)) %>%
      mutate(nmonths = case_when(is.na(nmonths) ~ 0L,
                                 TRUE ~ nmonths),
             ndays = case_when(is.na(ndays) ~ 0L,
                               TRUE ~ ndays)) %>%
      return()
  }
}

# Interpolate missing data by lake, year, season, stratum (based on interpolateMonthSeasonsAllWC)
phytoplankton_split <- phytoplankton_monthseason_means %>%
  filter(month_season %in% c("Spring", "Summer", "Fall")) %>%  # Remove winter sampling dates
  group_split(lake_id, month_season, stratum, taxon)

phytoplankton_split %>%
  map(~mutate(., lakeid_monthseason_stratum_taxon = str_c(lake_id, "_", month_season, "_", stratum, "_", taxon))) %>%
  map(~pull(.,lakeid_monthseason_stratum_taxon)) %>%  # Pull out variable
  map(~as.character(.)) %>%  # Convert factor to character
  map(~unique(.)) -> names(phytoplankton_split)

phytoplankton_monthseason_interpolated <- list(data = phytoplankton_split, names = names(phytoplankton_split), method = "approx") %>%
  pmap(interpolateMonthSeasonsTaxonWC) %>%
  map_dfr(`[`, c("lake_id", "year", "month_season", "stratum", "taxon", "mean", "nmonths", "ndays", "interpolated"))

(phytoplankton_monthseason_interpolated_approx_plot <- phytoplankton_monthseason_interpolated %>%
    filter(taxon == "total_phytoplankton") %>%
    mutate(monthseason_floor = case_when(month_season == "Spring" ~ ymd(paste0(year, "-04-01")),
                                         month_season == "Summer" ~ ymd(paste0(year, "-07-01")),
                                         month_season == "Fall" ~ ymd(paste0(year, "-10-01")))) %>%
    ggplot() +
    facet_col(lake_id~factor(stratum, levels = c("Epilimnion", "Metalimnion")),
              scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = monthseason_floor, y = mean,
                   colour = as.character(interpolated),
                   shape = factor(month_season, levels = c("Winter", "Spring", "Summer", "Fall")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(labels = comma) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "dodgerblue"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_shape_manual(values = c("Spring" = 50, "Summer" = 51, "Fall" = 52)) +
    labs(y = "Biomass (\u00b5g/L)",
         colour = "Interpolated",
         shape = "Season (month)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))

# Assess month-seasons with high incidences of data interpolation
phytoplankton_monthseason_interpolated %>%
  filter(taxon == "total_phytoplankton" & interpolated) %>%
  group_by(lake_id, stratum, month_season, interpolated) %>%
  count() %>%
  ungroup() %>%
  arrange(-n)

# Summarize interpolated data by season
phytoplankton_monthseason_interpolated %>%
  group_by(lake_id, year, month_season, stratum, taxon) %>%
  summarize(mean = mean(mean),
            nseasons = n(),
            nmonths = sum(nmonths),
            ndays = sum(ndays)) %>%
  ungroup() %>%
  arrange(lake_id, year,
          factor(stratum, levels = c("Epilimnion", "Metalimnion")),
          taxon)

# Summarize interpolated data by year (based on summarizeByYearInterpolatedWC)
# Remove all fall samples for all lakes (Metalimnion fall data are sparse or absent)
phytoplankton_byyear_interpolated <- phytoplankton_monthseason_interpolated %>%
  filter(month_season != "Fall") %>%
  group_by(lake_id, year, stratum, taxon) %>%
  summarize(mean = mean(mean),
            nseasons = n(),
            nmonths = sum(nmonths),
            ndays = sum(ndays)) %>%
  ungroup() %>%
  arrange(lake_id, year,
          factor(stratum, levels = c("Epilimnion", "Metalimnion")),
          taxon)

# Plot interpolated mean annual phytoplankton biomass (including Cyanobacteria)
(phytoplankton_wcyano_byyear_interpolated_plot <- phytoplankton_byyear_interpolated %>%
    filter(taxon == "total_phytoplankton") %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    geom_point(aes(x = year_floor, y = mean,
                   colour = factor(stratum, levels = c("Epilimnion", "Metalimnion")))) +
    geom_line(aes(x = year_floor, y = mean,
                  colour = factor(stratum, levels = c("Epilimnion", "Metalimnion")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(labels = comma,
                       breaks = seq(0, max(phytoplankton_byyear_interpolated$mean), 3000)) +
    scale_colour_manual(values = palette_stratum[which(names(palette_stratum) %in% unique(phytoplankton_byyear_interpolated$stratum))]) +
    scale_fill_manual(values = palette_manipulation) +
    labs(y = "Biomass (\u00b5g/L)",
         colour = "Stratum",
         fill = "Experiment") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/phytoplankton_wcyano_byyear_interpolated_plot.pdf", phytoplankton_wcyano_byyear_interpolated_plot, width = 12, height = 10, units = "in", device = cairo_pdf)

(phytoplankton_euk_byyear_interpolated_plot <- phytoplankton_byyear_interpolated %>%
    filter(taxon %in% c("Chlorophytes", "Chrysophytes", "Cryptophytes",
                        "Diatoms", "Dinoflagellates", "Euglenophytes")) %>%
    group_by(lake_id, year, stratum) %>%
    summarize(mean = sum(mean)) %>%
    ungroup %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    geom_point(aes(x = year_floor, y = mean,
                   colour = factor(stratum, levels = c("Epilimnion", "Metalimnion")))) +
    geom_line(aes(x = year_floor, y = mean,
                  colour = factor(stratum, levels = c("Epilimnion", "Metalimnion")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(labels = comma,
                       breaks = seq(0, max(phytoplankton_byyear_interpolated$mean), 3000)) +
    scale_colour_manual(values = palette_stratum[which(names(palette_stratum) %in% unique(phytoplankton_byyear_interpolated$stratum))]) +
    scale_fill_manual(values = palette_manipulation) +
    labs(y = "Biomass (\u00b5g/L)",
         colour = "Stratum",
         fill = "Experiment") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/phytoplankton_euk_byyear_interpolated_plot.pdf", phytoplankton_euk_byyear_interpolated_plot, width = 12, height = 6, units = "in", device = cairo_pdf)

(phytoplankton_euk_byyear_interpolated_areaplot <- phytoplankton_byyear_interpolated %>%
    filter(taxon %in% c("Chlorophytes", "Chrysophytes", "Cryptophytes",
                        "Diatoms", "Dinoflagellates", "Euglenophytes")) %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id+stratum~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    scale_fill_manual(values = palette_manipulation, name = "Experiment") +
    new_scale("fill") +
    geom_area(aes(x = year_floor, y = mean,
                  fill = factor(taxon, levels = c("Chlorophytes", "Chrysophytes", "Dinoflagellates",
                                                  "Cryptophytes", "Diatoms", "Euglenophytes"))),
              stat = "identity") +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(labels = comma,
                       breaks = seq(0, max(phytoplankton_byyear_interpolated$mean), 3000)) +
    scale_fill_manual(values = palette_phytoplankton[which(names(palette_phytoplankton) %in% c("Chlorophytes", "Chrysophytes", "Dinoflagellates",
                                                                                               "Cryptophytes", "Diatoms", "Euglenophytes"))],
                      name = "Phytoplankton") +
    labs(y = "Biomass (\u00b5g/L)",
         colour = "Stratum") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/phytoplankton_euk_byyear_interpolated_areaplot.pdf", phytoplankton_euk_byyear_interpolated_areaplot, width = 12, height = 8, units = "in", device = cairo_pdf)

# Write data summary to file
# phytoplankton_byyear_interpolated %>%
#   dplyr::select(lake_id, year, stratum, taxon, mean) %>%
#   rename(biomass = mean) %>%
#   mutate(stratum = str_to_lower(stratum),
#          biomass = round(biomass, 1)) %>%
#   pivot_wider(names_from = c(taxon, stratum), values_from = biomass, names_glue = "{taxon}_{stratum}_{.value}") %>%
#   write_tsv("output/ela/ela_phytoplankton_byyear_interpolated.tsv", col_names = TRUE)
