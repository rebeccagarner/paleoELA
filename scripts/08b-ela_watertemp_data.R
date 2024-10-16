# ELA water temperature data

setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
setwd("~/Dropbox/projects/ela18s/")

# Load libraries
library(ggforce)
library(scales)

# Load environmental data
source("scripts/08-ela_env_data.R")


#### Summarize water temperature data ####
# Assess frequency of water temperature measurements
any(is.na(watertemp$temp))  # Evaluates to FALSE

# Is maximum sampling depth consistent (within lakes) across time?
(watertemp_maxdepths <- watertemp %>%
    group_by(lake_id, date) %>%
    slice_max(depth_m) %>%
    ungroup())

(watertemp_plot <- watertemp %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_line(aes(x = date, y = depth_m),
              data = watertemp_maxdepths,
              linewidth = 0.1) +
    geom_point(aes(x = date, y = depth_m, colour = temp)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_reverse(breaks = seq(0, ceiling(max(watertemp$depth_m)), 5)) +
    scale_colour_distiller(palette = "RdYlBu") +
    labs(y = "Depth (m)",
         colour = "Temperature (\u00B0C)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey85")))
#ggsave("figures/ela_watertemp.pdf", watertemp_plot, width = 18, height = 6, units = "in", device = cairo_pdf)


#### Assess water column stratification ####
# Water column is stratified when there is temperature change of >= 1 degC/m
# Define function to assign stratified water column depths to epi/meta/hypolimnion strata
assignStrata <- function(data) {
  data_tmp <- data %>%
    arrange(depth_m)
  
  lakeid_tmp <- unique(data_tmp$lake_id)
  date_tmp <- unique(data_tmp$date)
  
  #print(paste0(lakeid_tmp, "_", date_tmp))
  
  stratification <- data_tmp %>%
    mutate(temp_over_depth = (lag(temp) - temp)/(depth_m - lag(depth_m))) %>%
    mutate(temp_over_depth_over1degc = case_when(temp_over_depth >= 1 ~ TRUE,
                                                 TRUE ~ FALSE))
  
  if (any(stratification$temp_over_depth_over1degc)) {
    stratification_tmp <- "Stratified"
  } else {
    stratification_tmp <- "Mixed"
  }
  
  data_tmp <- data_tmp %>%
    mutate(stratification = stratification_tmp)
  
  if (stratification_tmp == "Stratified") {
    metalimnion_start_depth <- stratification %>%
      filter(temp_over_depth_over1degc) %>%
      slice_min(depth_m) %>%
      pull(depth_m)
    
    metalimnion_end_depth <- stratification %>%
      filter(temp_over_depth_over1degc) %>%
      slice_max(depth_m) %>%
      pull(depth_m)
    
    strata <- data_tmp %>%
      mutate(stratum = case_when(depth_m < metalimnion_start_depth ~ "Epilimnion",
                                 depth_m >= metalimnion_start_depth & depth_m <= metalimnion_end_depth ~ "Metalimnion",
                                 depth_m > metalimnion_end_depth ~ "Hypolimnion"))
  } else if (stratification_tmp == "Mixed") {
    strata <- data_tmp %>%
      mutate(stratum = "watercolumn")
  }
  
  return(strata)
}

watertemp_split <- watertemp %>%
  group_split(lake_id, date)

watertemp_split %>%
  map(~mutate(., lakeid_date = str_c(lake_id, "_", date))) %>%
  map(~pull(.,lakeid_date)) %>%  # Pull out variable
  map(~as.character(.)) %>%  # Convert factor to character
  map(~unique(.)) -> names(watertemp_split)

watertemp_strata <- list(data = watertemp_split) %>%
  pmap(assignStrata) %>%
  map_dfr(`[`, c("lake_id", "date", "depth_m", "temp", "stratification", "stratum"))

# Calculate (mean annual) number of days of stratification
(stratification_freq_plot <- watertemp_strata %>%
    distinct(lake_id, date, stratification) %>%
    ggplot() +
    geom_rect(aes(xmin = season_start_date, xmax = season_end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = factor(season, levels = names(palette_season))),
              data = formatSeasons(watertemp_strata),
              alpha = 0.1) +
    scale_fill_manual(values = palette_season, name = "Season") +
    geom_point(aes(x = date, y = lake_id, colour = stratification)) +
    scale_colour_manual(values = palette_stratification, name = "Stratification") +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_discrete(limits = rev) +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title = element_blank(),
          panel.grid = element_blank()))
#ggsave("figures/ela_stratification_frequencies.pdf", stratification_freq_plot, width = 18, height = 4, units = "in", device = cairo_pdf)

# Visualize frequency of water column thermal profiling
(watertemp_sampling_freq_plot <- watertemp %>%
    ggplot() +
    geom_rect(aes(xmin = season_start_date, xmax = season_end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = factor(season, levels = names(palette_season))),
              data = formatSeasons(watertemp),
              alpha = 0.1) +
    scale_fill_manual(values = palette_season, name = "Season") +
    geom_point(aes(x = date, y = lake_id, colour = lake_id)) +
    scale_colour_manual(values = palette_lake, guide = "none") +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_discrete(limits = rev) +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title = element_blank(),
          panel.grid = element_blank()))
#ggsave("figures/ela_watertemp_sampling_frequencies.pdf", watertemp_sampling_freq_plot, width = 18, height = 4, units = "in", device = cairo_pdf)

# Assess frequency of water temperature sampling by season
watertemp %>%
  distinct(lake_id, date) %>%
  assignSeason() %>%
  group_by(lake_id, season) %>%
  count() %>%
  arrange(lake_id, -n)

# Visualize frequency of water temperature sampling by season
(watertemp_season_freq_plot <- watertemp %>%
    distinct(lake_id, date) %>%
    assignSeason() %>%
    mutate(year_floor = floor_date(date, "year")) %>%
    group_by(lake_id, year_floor, season) %>%
    count() %>%
    arrange(lake_id, year_floor, -n) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_bar(aes(x = year_floor, y = n, fill = factor(season, levels = names(palette_season))),
             stat = "identity") +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(watertemp), 5)) +
    scale_fill_manual(values = palette_season, name = "Season") +
    labs(y = "Sampling frequency (no. days)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_watertemp_sampling_frequencies_season.pdf", watertemp_season_freq_plot, width = 14, height = 8, units = "in", device = cairo_pdf)

# Assess frequency of water temperature sampling
watertemp %>%
  calculateDataIntervals()

# Summarize intervals between sampling
watertemp %>%
  summarizeIntervals()

# Summarize data by month (surface-most temperature)
watertemp %>%
  group_by(lake_id, date) %>%
  slice_min(depth_m) %>%
  ungroup() %>%
  summarizeByMonth("temp")


#### Derive water temperature variables for various depths ####
# Calculate mean and total change in temperature for water column strata
watertemp_strata_summary <- watertemp_strata %>%
  group_by(lake_id, date, stratification, stratum) %>%
  summarize(mean_temp = mean(temp),
            min_temp = min(temp),
            max_temp = max(temp)) %>%
  ungroup() %>%
  mutate(diff_temp = max_temp - min_temp) %>%
  dplyr::select(lake_id, date, stratification, stratum, mean_temp)


#### Summarize data by time intervals ####
# For mixed water columns, assign mean water column temperature to epi/meta/hypolimnion
watertemp_strata_summary_wc <- watertemp_strata_summary %>%
  mutate(stratum = str_to_lower(stratum)) %>%
  pivot_wider(names_from = stratum, values_from = mean_temp) %>%
  mutate(epilimnion = case_when(is.na(epilimnion) & !is.na(watercolumn) ~ watercolumn,
                                TRUE ~ epilimnion),
         metalimnion = case_when(is.na(metalimnion) & !is.na(watercolumn) ~ watercolumn,
                                 TRUE ~ metalimnion),
         hypolimnion = case_when(is.na(hypolimnion) & !is.na(watercolumn) ~ watercolumn,
                                 TRUE ~ hypolimnion)) %>%
  dplyr::select(-watercolumn) %>%
  pivot_longer(!c(lake_id, date, stratification), names_to = "stratum", values_to = "mean_temp") %>%
  filter(!is.na(mean_temp)) %>%
  mutate(stratum = str_to_sentence(stratum))

# Summarize data by month and month-season
watertemp_strata_summary_wc %>%
  summarizeByMonthSeasonWC("mean_temp")

# Summarize data by year (winter sampling removed)
watertemp_strata_summary_wc %>%
  summarizeByYearWC("mean_temp")

# Interpolate missing data by lake, year, season, stratum
(watertemp_strata_summary_interpolated <- watertemp_strata_summary_wc %>%
    interpolateMonthSeasonsAllWC("mean_temp", "approx"))

(watertemp_strata_summary_interpolated_approx_plot <- watertemp_strata_summary_interpolated %>%
    mutate(monthseason_floor = case_when(month_season == "Spring" ~ ymd(paste0(year, "-04-01")),
                                         month_season == "Summer" ~ ymd(paste0(year, "-07-01")),
                                         month_season == "Fall" ~ ymd(paste0(year, "-10-01")))) %>%
    ggplot() +
    facet_col(lake_id~factor(stratum, levels = c("Epilimnion", "Metalimnion", "Hypolimnion")),
              scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = monthseason_floor, y = mean,
                   colour = as.character(interpolated),
                   shape = factor(month_season, levels = c("Winter", "Spring", "Summer", "Fall")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(watertemp_strata_summary_interpolated), 1)) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "dodgerblue"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_shape_manual(values = c("Spring" = 50, "Summer" = 51, "Fall" = 52)) +
    labs(y = "Mean temperature (\u00B0C)",
         colour = "Interpolated",
         shape = "Season (month)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))

# Summarize interpolated data by year
(watertemp_strata_byyear_interpolated <- watertemp_strata_summary_wc %>%
    summarizeByYearInterpolatedWC("mean_temp", "approx") %>%
    filter(year < 2018))

(watertemp_strata_byyear_interpolated_plot <- watertemp_strata_byyear_interpolated %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    geom_point(aes(x = year_floor, y = mean,
                   colour = factor(stratum, levels = c("Epilimnion", "Metalimnion", "Hypolimnion")))) +
    geom_line(aes(x = year_floor, y = mean,
                  colour = factor(stratum, levels = c("Epilimnion", "Metalimnion", "Hypolimnion")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(watertemp_strata_summary_wc), 1)) +
    scale_colour_manual(values = palette_stratum) +
    scale_fill_manual(values = palette_manipulation) +
    labs(y = "Mean temperature (\u00B0C)",
         colour = "Stratum",
         fill = "Experiment") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_watertemp_strata_byyear_interpolated.pdf", watertemp_strata_byyear_interpolated_plot, width = 12, height = 10, units = "in", device = cairo_pdf)

# Write data summary to file
# watertemp_strata_byyear_interpolated %>%
#   dplyr::select(lake_id, year, stratum, mean) %>%
#   rename(temp = mean) %>%
#   mutate(stratum = str_to_lower(stratum),
#          temp = round(temp, 1)) %>%
#   pivot_wider(names_from = stratum, values_from = temp, names_glue = "{stratum}_{.value}") %>%
#   write_tsv("output/ela/ela_watertemp_strata_byyear_interpolated.tsv", col_names = TRUE)
