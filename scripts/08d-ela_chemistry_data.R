# ELA chemistry data

# Load libraries
library(ggforce)
library(scales)
library(ggnewscale)

# Load environmental data
source("scripts/08-ela_env_data.R")


#### Format chemistry data ####
# Examine chemistry variables
names(chemistry)[which(!names(chemistry) %in% c("lake_id", "date", "start_depth", "stratum"))]

# Convert to long format matrix
chemistry_long <- chemistry %>%
  pivot_longer(!c(lake_id, date, start_depth, stratum),
               names_to = "var", values_to = "value") %>%
  filter(!is.na(value))


#### Format seasons ####
# Define function to assign season from date
assignSeason <- function(date) {
  month_date_numeric <- 100 * month(date) + day(date)
  season_cuts <- cut(month_date_numeric, breaks = c(0, 319, 620, 921, 1220, 1231))
  levels(season_cuts) <- c("Winter","Spring","Summer","Fall","Winter")
  return(season_cuts)
}

season_days <- tibble(season = c("Winter", "Spring", "Summer", "Fall", "Winter"),
                      season_start_day = c(1, 79, 172, 265, 355),
                      season_end_day = c(78, 171, 264, 354, 365))


##### Assess sampling frequency ####
chemistry_split <- chemistry_long %>%
  group_split(var)

chemistry_split %>%
  map(~pull(.,var)) %>%  # Pull out chemical variable
  map(~as.character(.)) %>%  # Convert factor to character
  map(~unique(.)) -> names(chemistry_split)

# Define function to visualize frequency of chemical variable sampling
plotVarFreq <- function(data, names) {
  data %>%
    distinct(lake_id, date, var) %>%
    ggplot() +
    geom_rect(aes(xmin = season_start_date, xmax = season_end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = factor(season, levels = names(palette_season))),
              data = formatSeasons(data),
              alpha = 0.1) +
    scale_fill_manual(values = palette_season, name = "Season") +
    geom_point(aes(x = date, y = lake_id, colour = lake_id)) +
    scale_colour_manual(values = palette_lake, guide = "none") +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_discrete(limits = rev) +
    ggtitle(unique(data$var)) +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title = element_blank(),
          panel.grid = element_blank())
}

# Generate individual chemical variable sampling frequency plots and save to PDF
#pdf("figures/ela_chemistry_sampling_frequencies.pdf", width = 18, height = 4)
# (allplots <- list(data = chemistry_split, names = names(chemistry_split)) %>%
#     pmap(plotVarFreq))
#dev.off()


#### Extract total nitrogen data ####
# Select nitrogen variables
nitrogen <- chemistry_long %>%
  filter(var == "tdn_ug_L" | var == "suspn_ug_L") %>%
  pivot_wider(names_from = var, values_from = value)

# Assess nitrogen variable frequencies
nitrogen %>%
  mutate(suspn_val = case_when(!is.na(suspn_ug_L) ~ TRUE,
                               TRUE ~ FALSE),
         tdn_val = case_when(!is.na(tdn_ug_L) ~ TRUE,
                             TRUE ~ FALSE)) %>%
  mutate(n_vals = case_when(suspn_val & tdn_val ~ "suspn+tdn",
                            suspn_val & !tdn_val ~ "suspn only",
                            !suspn_val & tdn_val ~ "tdn only")) %>%
  group_by(lake_id, stratum, n_vals) %>%
  count(name = "nobs") %>%
  ungroup()

# Calculate total nitrogen (mg/L) where data are available (TN = TDN + suspN)
tn <- nitrogen %>%
  mutate(tn = tdn_ug_L+ suspn_ug_L) %>%
  filter(!is.na(tn)) %>%
  dplyr::select(-c(suspn_ug_L, tdn_ug_L)) %>%
  mutate(tn = tn/1000)

# Correct suspected typo in TN data
tn <- tn %>%
  mutate(tn = case_when(lake_id == "L226N" & date == "1971-06-24" & start_depth == 0 & stratum == "epi" ~ tn + 0.2,
                        TRUE ~ tn))

# Resolve water column strata
# There are multiple observations for epi and pro for the same lake on the same date
# There are observations for epi, pro, eup, and hyp for the same lake on the same date
# Observations for eup, hyp, all, and int never have associated start_depth
# Observations for epi and pro always have associated start_depth
# Observations for epi and pro start_depth are always unique (no redundant start_depth for a lake on a specific date)
# There is only one observation for int there are other observations for the lake on that date
# There is always only one observation for each of eup, hyp, and all per lake and date
tn <- tn %>%
  filter(stratum != "int")

# Plot TN time series
(tn_plot <- tn %>%
    filter(!is.na(start_depth)) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = date, y = start_depth, colour = tn)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_reverse(breaks = seq(0, ceiling(max(tn$start_depth, na.rm = TRUE)), 5)) +
    scale_colour_distiller(palette = "RdYlBu",
                           trans = "sqrt",
                           breaks = c(1, 4, 10, 20, 30)) +
    labs(y = "Depth (m)",
         colour = "TN (mg/L)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey85")))
#ggsave("figures/ela_tn.pdf", tn_plot, width = 18, height = 6, units = "in", device = cairo_pdf)

# Filter TN observations 
tn <- filterChemObs(tn)

# Assess epilimnion observations by depth
tn %>%
  filter(stratum == "Epilimnion") %>%
  arrange(start_depth) %>%
  pivot_wider(names_from = start_depth, names_prefix = "depth", values_from = tn)

# Calculate mean TN for water column strata
tn_strata_summary <- tn %>%
  group_by(lake_id, date, stratum) %>%
  summarize(mean_tn = mean(tn),
            min_tn = min(tn),
            max_tn = max(tn)) %>%
  ungroup() %>%
  mutate(diff_tn = max_tn - min_tn) %>%
  dplyr::select(lake_id, date, stratum, mean_tn)

# Plot TN sampling frequency
(tn_freq_plot <- tn_strata_summary %>%
    mutate(season = assignSeason(date)) %>%
    mutate(year = year(date),
           year_day = yday(date)) %>%
    ggplot() +
    facet_grid(factor(stratum, levels = c("Epilimnion", "Profundal"))~
                 factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))) +
    geom_rect(aes(xmin = -Inf, xmax = Inf,
                  ymin = season_start_day, ymax = season_end_day,
                  fill = season),
              data = season_days,
              alpha = 0.1) +
    scale_fill_manual(values = palette_season, name = "Season") +
    new_scale_fill() +
    geom_tile(aes(x = year, y = year_day, fill = mean_tn)) +
    scale_y_reverse(expand = c(0,0)) +
    scale_fill_distiller(palette = "RdYlBu",
                         trans = "sqrt",
                         breaks = c(1, 4, 10, 20, 30),
                         name = "TN (mg/L)") +
    labs(y = "Day") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_tn_freq_plot.pdf", tn_freq_plot, width = 18, height = 8, units = "in", device = cairo_pdf)

# Summarize data by month and month-season
tn_strata_summary %>%
  summarizeByMonthSeason("mean_tn")

# Summarize data by year (winter sampling removed)
tn_strata_summary %>%
  summarizeByYear("mean_tn")

# Interpolate missing data by lake, year, season, stratum
(tn_strata_summary_interpolated <- tn_strata_summary %>%
    interpolateMonthSeasonsAllWC("mean_tn", "approx"))

(tn_strata_summary_interpolated_approx_plot <- tn_strata_summary_interpolated %>%
    mutate(monthseason_floor = case_when(month_season == "Spring" ~ ymd(paste0(year, "-04-01")),
                                         month_season == "Summer" ~ ymd(paste0(year, "-07-01")),
                                         month_season == "Fall" ~ ymd(paste0(year, "-10-01")))) %>%
    ggplot() +
    facet_col(lake_id~factor(stratum, levels = c("Epilimnion", "Profundal")),
              scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = monthseason_floor, y = mean,
                   colour = as.character(interpolated),
                   shape = factor(month_season, levels = c("Spring", "Summer", "Fall")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(tn_strata_summary_interpolated), 1)) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "dodgerblue"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_shape_manual(values = c("Spring" = 50, "Summer" = 51, "Fall" = 52)) +
    labs(y = "TN (mg/L)",
         colour = "Interpolated",
         shape = "Season (month)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))

# Summarize interpolated data by year
(tn_strata_byyear_interpolated <- tn_strata_summary %>%
    summarizeByYearInterpolatedWC("mean_tn", "approx") %>%
    filter(year < 2018))

(tn_strata_byyear_interpolated_plot <- tn_strata_byyear_interpolated %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    geom_point(aes(x = year_floor, y = mean,
                   colour = factor(stratum, levels = c("Epilimnion", "Profundal")))) +
    geom_line(aes(x = year_floor, y = mean,
                  colour = factor(stratum, levels = c("Epilimnion", "Profundal")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(tn_strata_summary), 1)) +
    scale_colour_manual(values = palette_stratum) +
    scale_fill_manual(values = palette_manipulation) +
    labs(y = "TN (mg/L)",
         colour = "Stratum",
         fill = "Experiment") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_tn_strata_byyear_interpolated.pdf", tn_strata_byyear_interpolated_plot, width = 12, height = 10, units = "in", device = cairo_pdf)


#### Extract total phosphorus data ####
# Select phosphorus variables
phosphorus <- chemistry_long %>%
  filter(var == "tdp_ug_L" | var == "suspp_ug_L") %>%
  pivot_wider(names_from = var, values_from = value)

# Assess phosphorus variable frequencies
phosphorus %>%
  mutate(suspp_val = case_when(!is.na(suspp_ug_L) ~ TRUE,
                               TRUE ~ FALSE),
         tdp_val = case_when(!is.na(tdp_ug_L) ~ TRUE,
                             TRUE ~ FALSE)) %>%
  mutate(n_vals = case_when(suspp_val & tdp_val ~ "suspp+tdp",
                            suspp_val & !tdp_val ~ "suspp only",
                            !suspp_val & tdp_val ~ "tdp only")) %>%
  group_by(lake_id, stratum, n_vals) %>%
  count(name = "nobs") %>%
  ungroup()

# Calculate total phosphorus (ug/L) where data are available (TP = tdp + suspp)
tp <- phosphorus %>%
  mutate(tp = tdp_ug_L+ suspp_ug_L) %>%
  filter(!is.na(tp)) %>%
  dplyr::select(-c(suspp_ug_L, tdp_ug_L))

# Resolve water column strata
# There are multiple observations for epi and pro for the same lake on the same date
# There are observations for epi, pro, eup, and hyp for the same lake on the same date
# Observations for eup, hyp, and all never have associated start_depth
# Observations for epi and pro always have associated start_depth
# Observations for epi and pro start_depth are always unique (no redundant start_depth for a lake on a specific date)
# There is always only one observation for each of eup, hyp, and all per lake and date

# Plot TP time series
(tp_plot <- tp %>%
    filter(!is.na(start_depth)) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = date, y = start_depth, colour = tp)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_reverse(breaks = seq(0, ceiling(max(tp$start_depth, na.rm = TRUE)), 5)) +
    scale_colour_distiller(palette = "RdYlBu",
                           trans = "log10",
                           breaks = c(1, 10, 100, 1000, 5000),
                           labels = comma) +
    labs(y = "Depth (m)",
         colour = "TP (\u00b5g/L)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey85")))
#ggsave("figures/ela_tp.pdf", tp_plot, width = 18, height = 6, units = "in", device = cairo_pdf)

# Filter TP observations 
tp <- filterChemObs(tp)

# Assess epilimnion observations by depth
tp %>%
  filter(stratum == "Epilimnion") %>%
  arrange(start_depth) %>%
  pivot_wider(names_from = start_depth, names_prefix = "depth", values_from = tp)

# Calculate mean TP for water column strata
tp_strata_summary <- tp %>%
  group_by(lake_id, date, stratum) %>%
  summarize(mean_tp = mean(tp),
            min_tp = min(tp),
            max_tp = max(tp)) %>%
  ungroup() %>%
  mutate(diff_tp = max_tp - min_tp) %>%
  dplyr::select(lake_id, date, stratum, mean_tp)

# Plot TP sampling frequency
(tp_freq_plot <- tp_strata_summary %>%
    mutate(season = assignSeason(date)) %>%
    mutate(year = year(date),
           year_day = yday(date)) %>%
    ggplot() +
    facet_grid(factor(stratum, levels = c("Epilimnion", "Profundal"))~
                 factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))) +
    geom_rect(aes(xmin = -Inf, xmax = Inf,
                  ymin = season_start_day, ymax = season_end_day,
                  fill = season),
              data = season_days,
              alpha = 0.1) +
    scale_fill_manual(values = palette_season, name = "Season") +
    new_scale_fill() +
    geom_tile(aes(x = year, y = year_day, fill = mean_tp)) +
    scale_y_reverse(expand = c(0,0)) +
    scale_fill_distiller(palette = "RdYlBu",
                         trans = "log10",
                         breaks = c(1, 10, 100, 1000, 5000),
                         labels = comma,
                         name = "TP (\u00b5g/L)") +
    labs(y = "Day") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_tp_freq_plot.pdf", tp_freq_plot, width = 18, height = 8, units = "in", device = cairo_pdf)

# Summarize data by month and month-season
tp_strata_summary %>%
  summarizeByMonthSeason("mean_tp")

# Summarize data by year (winter sampling removed)
tp_strata_summary %>%
  summarizeByYear("mean_tp")

# Interpolate missing data by lake, year, season, stratum
(tp_strata_summary_interpolated <- tp_strata_summary %>%
    interpolateMonthSeasonsAllWC("mean_tp", "approx"))

(tp_strata_summary_interpolated_approx_plot <- tp_strata_summary_interpolated %>%
    mutate(monthseason_floor = case_when(month_season == "Spring" ~ ymd(paste0(year, "-04-01")),
                                         month_season == "Summer" ~ ymd(paste0(year, "-07-01")),
                                         month_season == "Fall" ~ ymd(paste0(year, "-10-01")))) %>%
    ggplot() +
    facet_col(lake_id~factor(stratum, levels = c("Epilimnion", "Profundal")),
              scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = monthseason_floor, y = mean,
                   colour = as.character(interpolated),
                   shape = factor(month_season, levels = c("Spring", "Summer", "Fall")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(tp_strata_summary_interpolated), 50),
                       labels = comma) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "dodgerblue"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_shape_manual(values = c("Spring" = 50, "Summer" = 51, "Fall" = 52)) +
    labs(y = "TP (\u00b5g/L)",
         colour = "Interpolated",
         shape = "Season (month)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))

# Summarize interpolated data by year
(tp_strata_byyear_interpolated <- tp_strata_summary %>%
    summarizeByYearInterpolatedWC("mean_tp", "approx") %>%
    filter(year < 2018))

(tp_strata_byyear_interpolated_plot <- tp_strata_byyear_interpolated %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    geom_point(aes(x = year_floor, y = mean,
                   colour = factor(stratum, levels = c("Epilimnion", "Profundal")))) +
    geom_line(aes(x = year_floor, y = mean,
                  colour = factor(stratum, levels = c("Epilimnion", "Profundal")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(tp_strata_summary), 50)) +
    scale_colour_manual(values = palette_stratum) +
    scale_fill_manual(values = palette_manipulation) +
    labs(y = "TP ((\u00b5g/L)",
         colour = "Stratum",
         fill = "Experiment") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_tp_strata_byyear_interpolated.pdf", tp_strata_byyear_interpolated_plot, width = 12, height = 10, units = "in", device = cairo_pdf)


#### Extract pH data ####
# Select pH variables
ph <- chemistry_long %>%
  filter(var == "ph_situ" | var == "ph_equil") %>%
  pivot_wider(names_from = var, values_from = value)

# Assess pH variable frequencies
ph %>%
  mutate(ph_vars = case_when(!is.na(ph_situ) & !is.na(ph_equil) ~ "situ+equil",
                             !is.na(ph_situ) & is.na(ph_equil) ~ "situ only",
                             is.na(ph_situ) & !is.na(ph_equil) ~ "equil only")) %>%
  group_by(ph_vars) %>%
  count()

ph %>%
  ggplot() +
  geom_point(aes(x = ph_situ, y = ph_equil, colour = lake_id)) +
  scale_colour_manual(values = palette_lake, name = "Lake") +
  theme_classic()

# Select ph_situ variable
ph <- ph %>%
  dplyr::select(-c(ph_equil)) %>%
  filter(!is.na(ph_situ)) %>%
  rename(ph = ph_situ)

# Resolve water column strata
# There are multiple observations for epi and pro for the same lake on the same date
# There are observations for epi, pro, eup, and hyp for the same lake on the same date
# Observations for eup, hyp, and all never have associated start_depth
# Observations for epi and pro always have associated start_depth
# Observations for epi and pro start_depth are always unique (no redundant start_depth for a lake on a specific date)
# There is always only one observation for each of eup, hyp, and all per lake and date

# Plot pH time series
(ph_plot <- ph %>%
    filter(!is.na(start_depth)) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = date, y = start_depth, colour = ph)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_reverse(breaks = seq(0, ceiling(max(ph$start_depth, na.rm = TRUE)), 5)) +
    scale_colour_distiller(palette = "RdYlBu",
                           breaks = c(5, 7, 9),
                           labels = comma) +
    labs(y = "Depth (m)",
         colour = "pH") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey85")))
#ggsave("figures/ela_ph.pdf", ph_plot, width = 18, height = 6, units = "in", device = cairo_pdf)

# Filter pH observations 
ph <- filterChemObs(ph)

# Assess epilimnion observations by depth
ph %>%
  filter(stratum == "Epilimnion") %>%
  arrange(start_depth) %>%
  pivot_wider(names_from = start_depth, names_prefix = "depth", values_from = ph)

# Calculate mean pH for water column strata
ph_strata_summary <- ph %>%
  group_by(lake_id, date, stratum) %>%
  summarize(mean_ph = mean(ph),
            min_ph = min(ph),
            max_ph = max(ph)) %>%
  ungroup() %>%
  mutate(diff_ph = max_ph - min_ph) %>%
  dplyr::select(lake_id, date, stratum, mean_ph)

# Plot pH sampling frequency
(ph_freq_plot <- ph_strata_summary %>%
    mutate(season = assignSeason(date)) %>%
    mutate(year = year(date),
           year_day = yday(date)) %>%
    ggplot() +
    facet_grid(factor(stratum, levels = c("Epilimnion", "Profundal"))~
                 factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))) +
    geom_rect(aes(xmin = -Inf, xmax = Inf,
                  ymin = season_start_day, ymax = season_end_day,
                  fill = season),
              data = season_days,
              alpha = 0.1) +
    scale_fill_manual(values = palette_season, name = "Season") +
    new_scale_fill() +
    geom_tile(aes(x = year, y = year_day, fill = mean_ph)) +
    scale_y_reverse(expand = c(0,0)) +
    scale_fill_distiller(palette = "RdYlBu",
                         breaks = c(5, 7, 9),
                         name = "pH") +
    labs(y = "Day") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_ph_freq_plot.pdf", ph_freq_plot, width = 18, height = 8, units = "in", device = cairo_pdf)

# Summarize data by month and month-season
ph_strata_summary %>%
  summarizeByMonthSeason("mean_ph")

# Summarize data by year (winter sampling removed)
ph_strata_summary %>%
  summarizeByYear("mean_ph")

# Interpolate missing data by lake, year, season, stratum
(ph_strata_summary_interpolated <- ph_strata_summary %>%
    interpolateMonthSeasonsAllWC("mean_ph", "approx"))

(ph_strata_summary_interpolated_approx_plot <- ph_strata_summary_interpolated %>%
    mutate(monthseason_floor = case_when(month_season == "Spring" ~ ymd(paste0(year, "-04-01")),
                                         month_season == "Summer" ~ ymd(paste0(year, "-07-01")),
                                         month_season == "Fall" ~ ymd(paste0(year, "-10-01")))) %>%
    ggplot() +
    facet_col(lake_id~factor(stratum, levels = c("Epilimnion", "Profundal")),
              scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = monthseason_floor, y = mean,
                   colour = as.character(interpolated),
                   shape = factor(month_season, levels = c("Spring", "Summer", "Fall")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(ph_strata_summary_interpolated), 1)) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "dodgerblue"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_shape_manual(values = c("Spring" = 50, "Summer" = 51, "Fall" = 52)) +
    labs(y = "pH",
         colour = "Interpolated",
         shape = "Season (month)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))

# Summarize interpolated data by year
(ph_strata_byyear_interpolated <- ph_strata_summary %>%
    summarizeByYearInterpolatedWC("mean_ph", "approx") %>%
    filter(year < 2018))

(ph_strata_byyear_interpolated_plot <- ph_strata_byyear_interpolated %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    geom_point(aes(x = year_floor, y = mean,
                   colour = factor(stratum, levels = c("Epilimnion", "Profundal")))) +
    geom_line(aes(x = year_floor, y = mean,
                  colour = factor(stratum, levels = c("Epilimnion", "Profundal")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(ph_strata_summary), 1)) +
    scale_colour_manual(values = palette_stratum) +
    scale_fill_manual(values = palette_manipulation) +
    labs(y = "pH",
         colour = "Stratum",
         fill = "Experiment") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_ph_strata_byyear_interpolated.pdf", ph_strata_byyear_interpolated_plot, width = 12, height = 10, units = "in", device = cairo_pdf)


#### Extract oxygen data ####
# Select oxygen variable
oxygen <- chemistry_long %>%
  filter(var == "o2_mg_L") %>%
  pivot_wider(names_from = var, values_from = value) %>%
  rename(oxygen = o2_mg_L) %>%
  mutate(stratum = case_when(stratum == "epi" ~ "Epilimnion",
                             stratum == "pro" ~ "Profundal")) %>%
  arrange(lake_id, date, stratum, start_depth)

# Plot oxygen time series
(oxygen_plot <- oxygen %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = date, y = start_depth, colour = oxygen)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_reverse(breaks = seq(0, ceiling(max(oxygen$start_depth, na.rm = TRUE)), 5)) +
    scale_colour_distiller(palette = "RdYlBu",
                           direction = 2,
                           #breaks = c(5, 10, 15, 20),
                           labels = comma) +
    labs(y = "Depth (m)",
         colour = "Oxygen (mg/L)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey85")))
#ggsave("figures/ela_oxygen.pdf", oxygen_plot, width = 18, height = 6, units = "in", device = cairo_pdf)

# Assess epilimnion observations by depth
oxygen %>%
  filter(stratum == "Epilimnion") %>%
  arrange(start_depth) %>%
  pivot_wider(names_from = start_depth, names_prefix = "depth", values_from = oxygen)

# Calculate mean oxygen for water column strata
oxygen_strata_summary <- oxygen %>%
  group_by(lake_id, date, stratum) %>%
  summarize(mean_oxygen = mean(oxygen),
            min_oxygen = min(oxygen),
            max_oxygen = max(oxygen)) %>%
  ungroup() %>%
  mutate(diff_oxygen = max_oxygen - min_oxygen) %>%
  dplyr::select(lake_id, date, stratum, mean_oxygen)

# Plot oxygen sampling frequency
(oxygen_freq_plot <- oxygen_strata_summary %>%
    mutate(season = assignSeason(date)) %>%
    mutate(year = year(date),
           year_day = yday(date)) %>%
    ggplot() +
    facet_grid(factor(stratum, levels = c("Epilimnion", "Profundal"))~
                 factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))) +
    geom_rect(aes(xmin = -Inf, xmax = Inf,
                  ymin = season_start_day, ymax = season_end_day,
                  fill = season),
              data = season_days,
              alpha = 0.1) +
    scale_fill_manual(values = palette_season, name = "Season") +
    new_scale_fill() +
    geom_tile(aes(x = year, y = year_day, fill = mean_oxygen)) +
    scale_y_reverse(expand = c(0,0)) +
    scale_fill_distiller(palette = "RdYlBu",
                         direction = 2,
                         breaks = seq(0, max(oxygen_strata_summary$mean_oxygen), 5),
                         labels = comma,
                         name = "Oxygen (mg/L)") +
    labs(y = "Day") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_oxygen_freq_plot.pdf", oxygen_freq_plot, width = 18, height = 8, units = "in", device = cairo_pdf)

# Summarize data by month and month-season
oxygen_strata_summary %>%
  summarizeByMonthSeason("mean_oxygen")

# Summarize data by year (winter sampling removed)
oxygen_strata_summary %>%
  summarizeByYear("mean_oxygen")

# Interpolate missing data by lake, year, season, stratum
(oxygen_strata_summary_interpolated <- oxygen_strata_summary %>%
    interpolateMonthSeasonsAllWC("mean_oxygen", "approx"))

(oxygen_strata_summary_interpolated_approx_plot <- oxygen_strata_summary_interpolated %>%
    mutate(monthseason_floor = case_when(month_season == "Spring" ~ ymd(paste0(year, "-04-01")),
                                         month_season == "Summer" ~ ymd(paste0(year, "-07-01")),
                                         month_season == "Fall" ~ ymd(paste0(year, "-10-01")))) %>%
    ggplot() +
    facet_col(lake_id~factor(stratum, levels = c("Epilimnion", "Profundal")),
              scales = "free_y", space = "free", strip.position = "right") +
    geom_point(aes(x = monthseason_floor, y = mean,
                   colour = as.character(interpolated),
                   shape = factor(month_season, levels = c("Spring", "Summer", "Fall")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(oxygen_strata_summary_interpolated), 1)) +
    scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "dodgerblue"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_shape_manual(values = c("Spring" = 50, "Summer" = 51, "Fall" = 52)) +
    labs(y = "oxygen",
         colour = "Interpolated",
         shape = "Season (month)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))

# Summarize interpolated data by year
(oxygen_strata_byyear_interpolated <- oxygen_strata_summary %>%
    summarizeByYearInterpolatedWC("mean_oxygen", "approx") %>%
    filter(year < 2018))

(oxygen_strata_byyear_interpolated_plot <- oxygen_strata_byyear_interpolated %>%
    mutate(year_floor = ymd(paste0(year, "-01-01"))) %>%
    ggplot() +
    facet_col(lake_id~., scales = "free_y", space = "free", strip.position = "right") +
    geom_rect(aes(xmin = start_date, xmax = end_date,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations,
              alpha = 0.5) +
    geom_point(aes(x = year_floor, y = mean,
                   colour = factor(stratum, levels = c("Epilimnion", "Profundal")))) +
    geom_line(aes(x = year_floor, y = mean,
                  colour = factor(stratum, levels = c("Epilimnion", "Profundal")))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%y",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, nrow(oxygen_strata_summary), 1)) +
    scale_colour_manual(values = palette_stratum) +
    scale_fill_manual(values = palette_manipulation) +
    labs(y = "Oxygen (mg/L)",
         colour = "Stratum",
         fill = "Experiment") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank()))
#ggsave("figures/ela_oxygen_strata_byyear_interpolated.pdf", oxygen_strata_byyear_interpolated_plot, width = 12, height = 10, units = "in", device = cairo_pdf)


#### Write chemistry data to file ####
# Join chemistry data
chemistry_strata_byyear_interpolated <- tn_strata_byyear_interpolated %>%
  dplyr::select(lake_id, year, stratum, mean) %>%
  rename(tn = mean) %>%
  mutate(stratum = str_to_lower(stratum),
         tn = round(tn, 1)) %>%
  pivot_wider(names_from = stratum, values_from = tn, names_glue = "{stratum}_{.value}") %>%
  left_join(tp_strata_byyear_interpolated %>%
              dplyr::select(lake_id, year, stratum, mean) %>%
              rename(tp = mean) %>%
              mutate(stratum = str_to_lower(stratum),
                     tp = round(tp, 1)) %>%
              pivot_wider(names_from = stratum, values_from = tp, names_glue = "{stratum}_{.value}"),
            by = c("lake_id", "year")) %>%
  left_join(ph_strata_byyear_interpolated %>%
              dplyr::select(lake_id, year, stratum, mean) %>%
              rename(ph = mean) %>%
              mutate(stratum = str_to_lower(stratum),
                     ph = round(ph, 1)) %>%
              pivot_wider(names_from = stratum, values_from = ph, names_glue = "{stratum}_{.value}"),
            by = c("lake_id", "year")) %>%
  left_join(oxygen_strata_byyear_interpolated %>%
              dplyr::select(lake_id, year, stratum, mean) %>%
              rename(oxygen = mean) %>%
              mutate(stratum = str_to_lower(stratum),
                     oxygen = round(oxygen, 1)) %>%
              pivot_wider(names_from = stratum, values_from = oxygen, names_glue = "{stratum}_{.value}"),
            by = c("lake_id", "year"))

# Write data summary to file
# chemistry_strata_byyear_interpolated %>%
#   write_tsv("output/ela/ela_chemistry_strata_byyear_interpolated.tsv", col_names = TRUE)
