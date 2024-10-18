# ELA Secchi disk depth data

# Load libraries
library(ggforce)
library(scales)

# Load environmental data
source("scripts/08-ela_env_data.R")


#### Format data ####
# Assess frequency of Secchi/light attenuation measurements
any(is.na(secchi$secchi_m))  # Evaluates to TRUE
any(is.na(secchi$kd))  # Evaluates to TRUE

secchi %>%
  filter(!is.na(secchi_m) & !is.na(kd)) %>%
  mutate(secchi_vs_kd = secchi_m/kd)

# Extract Secchi disk depth (SDD) sampling
sdd <- secchi %>%
  filter(!is.na(secchi_m)) %>%
  dplyr::select(lake_id, date, secchi_m)


#### Assess sampling frequency ####
# Visualize frequency of SDD sampling
(sdd_sampling_freq_plot <- sdd %>%
   ggplot() +
   geom_rect(aes(xmin = season_start_date, xmax = season_end_date,
                 ymin = -Inf, ymax = Inf,
                 fill = factor(season, levels = names(palette_season))),
             data = formatSeasons(sdd),
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
#ggsave("figures/ela_secchi_sampling_frequencies.pdf", sdd_sampling_freq_plot, width = 18, height = 4, units = "in", device = cairo_pdf)

# Visualize frequency of SDD sampling by season
(sdd_season_freq_plot <- sdd %>%
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
    scale_y_continuous(breaks = seq(0, nrow(sdd), 5)) +
    scale_fill_manual(values = palette_season, name = "Season") +
    labs(y = "Sampling frequency") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_secchi_sampling_frequencies_season.pdf", sdd_season_freq_plot, width = 14, height = 6, units = "in", device = cairo_pdf)

# Visualize frequency of SDD sampling by month
(sdd_month_freq_plot <- sdd %>%
    distinct(lake_id, date) %>%
    mutate(month = month(date, label = TRUE, abbr = FALSE)) %>%
    mutate(year_floor = floor_date(date, "year")) %>%
    group_by(lake_id, year_floor, month) %>%
    count() %>%
    arrange(lake_id, year_floor, -n) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_bar(aes(x = year_floor, y = n,
                 fill = factor(month, levels = names(palette_month))),
             stat = "identity") +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(sdd), 5)) +
    scale_fill_manual(values = palette_month) +
    labs(y = "Sampling frequency",
         fill = "Month") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_secchi_sampling_frequencies_month.pdf", sdd_month_freq_plot, width = 14, height = 6, units = "in", device = cairo_pdf)

# Assess frequency of SSD measurements
sdd %>%
  calculateDataIntervals()

# Summarize intervals between sampling
(sdd_intervals_summary <- sdd %>%
    summarizeIntervals())


#### Summarize data by time intervals ####
# Summarize data by month
sdd %>%
  summarizeByMonth("secchi_m")

(sdd_month_means <- sdd %>%
    summarizeByMonth("secchi_m") %>%
    mutate(month_floor = ymd(paste0(year, "-", str_pad(month, 2, "left", 0), "-01"))) %>%
    mutate(month = month(month, label = TRUE, abbr = FALSE)) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = month_floor, y = mean,
                   colour = factor(month, levels = names(palette_month)))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(sdd), 1)) +
    scale_colour_manual(values = palette_month) +
    labs(y = "Mean SDD (m)",
         colour = "Month") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))

# Summarize data by month and month-season
sdd %>%
  summarizeByMonthSeason("secchi_m")

# Summarize data by year
sdd %>%
  summarizeByYear("secchi_m")

(sdd_year_means_plot <- sdd %>%
    summarizeByYear("secchi_m") %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = year_floor, y = mean, colour = as.character(nseasons))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(sdd), 1)) +
    scale_colour_manual(values = c("1" = "red", "2" = "orange", "3" = "dodgerblue"),
                        guide = guide_legend(reverse = TRUE)) +
    labs(y = "Mean SDD (m)",
         colour = "No. seasons") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))

# Interpolate missing data by lake, year, season
(sdd_interpolated <- sdd %>%
    interpolateMonthSeasonsAll("secchi_m", "approx"))

(sdd_interpolated_approx_plot <- sdd_interpolated %>%
    mutate(monthseason_floor = case_when(month_season == "Spring" ~ ymd(paste0(year, "-04-01")),
                                         month_season == "Summer" ~ ymd(paste0(year, "-07-01")),
                                         month_season == "Fall" ~ ymd(paste0(year, "-10-01")))) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = monthseason_floor, y = mean,
                   colour = as.character(interpolated),
                   shape = factor(month_season, levels = c("Winter", "Spring", "Summer", "Fall")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(sdd_interpolated), 1)) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "dodgerblue"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_shape_manual(values = c("Spring" = 50, "Summer" = 51, "Fall" = 52)) +
    labs(y = "Mean SDD (m)",
         colour = "Interpolated",
         shape = "Season (month)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))

# Summarize interpolated data by year
(sdd_byyear_interpolated <- sdd %>%
    summarizeByYearInterpolated("secchi_m", "approx") %>%
    filter(year < 2018))

(sdd_byyear_interpolated_plot <- sdd_byyear_interpolated %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    geom_point(aes(x = year_floor, y = mean)) +
    geom_line(aes(x = year_floor, y = mean)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_reverse(breaks = seq(0, nrow(sdd), 1)) +
    scale_fill_manual(values = palette_manipulation) +
    labs(y = "Mean Secchi disk depth (m)",
         fill = "Experiment") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_sdd_byyear_interpolated.pdf", sdd_byyear_interpolated_plot, width = 12, height = 6, units = "in", device = cairo_pdf)

# Write data summary to file
# sdd_byyear_interpolated %>%
#   dplyr::select(lake_id, year, mean) %>%
#   rename(sdd = mean) %>%
#   mutate(sdd = round(sdd, 1)) %>%
#   write_tsv("output/ela/ela_sdd_byyear_interpolated.tsv", col_names = TRUE)
