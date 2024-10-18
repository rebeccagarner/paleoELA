# ELA climate data

# Load libraries
library(ggforce)
library(scales)

# Load environmental data
source("scripts/08-ela_env_data.R")


#### Summarize precipitation data ####
# Assess frequency of precipitation measurements
any(is.na(precipitation$tot_precip_mm))  # Evaluates to FALSE
any(is.na(precipitation$rain_mm))  # Evaluates to FALSE
any(is.na(precipitation$snow_mm))  # Evaluates to FALSE

precipitation %>%
  calculateTimeDiffs()  # Daily measurements

# Identify months with measurements for only some days
precipitation %>%
  partMonth()

# Calculate monthly total precipitations
precipitation_month_totals <- precipitation %>%
  filter(date > "1969-06-30") %>%
  mutate(year = year(date),
         month = month(date),
         day = day(date)) %>%
  group_by(year, month) %>%
  summarize(total_precip_month = sum(tot_precip_mm)) %>%
  ungroup() %>%
  mutate(date = make_date(year, month, 15))

precipitation_month_totals %>%
  ggplot() +
  geom_rect(aes(xmin = season_start_date, xmax = season_end_date,
                ymin = -Inf, ymax = Inf,
                fill = season),
            data = formatSeasons(precipitation_month_totals),
            alpha = 0.1) +
  scale_fill_manual(values = palette_season, name = "Season") +
  geom_point(aes(x = date, y = total_precip_month)) +
  geom_line(aes(x = date, y = total_precip_month)) +
  geom_smooth(aes(x = date, y = total_precip_month)) +
  labs(y = "Total monthly precipitation (mm)") +
  scale_x_date(date_breaks = "1 year", date_labels = "%y",
               expand = c(0,0)) +
  theme_bw() %+replace%
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        panel.grid = element_blank())

# Identify years with measurements for only some months
precipitation %>%
  partYear()

# Calculate annual total precipitations
precipitation_year_totals <- precipitation %>%
  filter(date > "1969-12-31" & date < "2018-01-01") %>%
  mutate(year = year(date),
         month = month(date)) %>%
  group_by(year) %>%
  summarize(total_precip_year = sum(tot_precip_mm)) %>%
  ungroup()

(precipitation_year_totals_plot <- precipitation_year_totals %>%
    ggplot(aes(x = year, y = total_precip_year)) +
    geom_point() +
    geom_line() +
    geom_smooth() +
    labs(y = "Total annual precipitation (mm)") +
    scale_y_continuous(labels = comma) +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          panel.grid = element_blank()))
#ggsave("figures/ela_precipitation_year_totals.pdf", precipitation_year_totals_plot, width = 5, height = 4, units = "in", device = cairo_pdf)


#### Summarize air temperature data ####
# Calculate (mean annual) number of ice-free days?

# Assess frequency of air temperature measurements
any(is.na(airtemp$max_temp))  # Evaluates to TRUE
any(is.na(airtemp$min_temp))  # Evaluates to TRUE

# Assess frequency of air temperature measurements
airtemp %>%
  filter(!is.na(max_temp)) %>%
  calculateTimeDiffs()

airtemp %>%
  filter(!is.na(min_temp)) %>%
  calculateTimeDiffs()

# Identify months with measurements for only some days
airtemp %>%
  filter(!is.na(max_temp)) %>%
  partMonth()

airtemp %>%
  filter(!is.na(min_temp)) %>%
  partMonth()

# Calculate monthly mean air temperatures
airtemp_month_totals <- airtemp %>%
  filter(date > "1969-06-30") %>%
  mutate(year = year(date),
         month = month(date),
         day = day(date)) %>%
  group_by(year, month) %>%
  summarize(mean_maxairtemp_month = mean(max_temp, na.rm = TRUE),
            mean_minairtemp_month = mean(min_temp, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(date = make_date(year, month, 15))

airtemp_month_totals %>%
  pivot_longer(!c(year, month, date), names_to = "airtemp_metric", values_to = "temp") %>%
  mutate(airtemp_metric = case_when(grepl("maxairtemp", airtemp_metric) ~ "max",
                                    grepl("minairtemp", airtemp_metric) ~ "min")) %>%
  ggplot() +
  geom_rect(aes(xmin = season_start_date, xmax = season_end_date,
                ymin = -Inf, ymax = Inf,
                fill = season),
            data = formatSeasons(airtemp_month_totals),
            alpha = 0.1) +
  scale_fill_manual(values = palette_season, name = "Season") +
  geom_point(aes(x = date, y = temp, colour = airtemp_metric)) +
  geom_line(aes(x = date, y = temp, colour = airtemp_metric)) +
  geom_smooth(aes(x = date, y = temp, colour = airtemp_metric)) +
  labs(y = "Mean annual air temperature (\u00B0C)") +
  scale_x_date(date_breaks = "1 year", date_labels = "%y",
               expand = c(0,0)) +
  theme_bw() %+replace%
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        panel.grid = element_blank())

# Identify years with measurements for only some months
airtemp %>%
  partYear()

# Calculate annual air temperature extremes
airtemp_year_means <- airtemp %>%
  filter(date > "1969-12-31" & date < "2018-01-01") %>%
  mutate(year = year(date),
         month = month(date)) %>%
  group_by(year) %>%
  summarize(mean_maxairtemp_year = mean(max_temp, na.rm = TRUE),
            mean_minairtemp_year = mean(min_temp, na.rm = TRUE)) %>%
  ungroup()

(airtemp_year_means_plot <- airtemp_year_means %>%
    pivot_longer(!year, names_to = "airtemp_metric", values_to = "temp") %>%
    mutate(airtemp_metric = case_when(grepl("maxairtemp", airtemp_metric) ~ "max",
                                      grepl("minairtemp", airtemp_metric) ~ "min")) %>%
    ggplot(aes(x = year, y = temp, colour = airtemp_metric)) +
    geom_hline(yintercept = 0, size = 0.5) +
    geom_point() +
    geom_line() +
    geom_smooth() +
    labs(y = "Mean annual air temperature (\u00B0C)") +
    guides(colour = guide_legend(title = NULL,
                                 override.aes = list(fill = NA))) +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          panel.grid = element_blank()))
#ggsave("figures/ela_airtemp_year_means.pdf", airtemp_year_means_plot, width = 5, height = 5, units = "in", device = cairo_pdf)

airtemp %>%
  mutate(ydays = yday(date)) %>%
  filter(!is.na(ice)) %>%
  ggplot() +
  geom_point(aes(x = year(date), y = ydays, colour = ice)) +
  labs(y = "Calendar day",
       colour = "Ice") +
  theme_bw() %+replace%
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        panel.grid = element_blank())

airtemp %>%
  mutate(ydays = yday(date)) %>%
  filter(!is.na(first_last)) %>%
  ggplot() +
  geom_point(aes(x = year(date), y = ydays, colour = as.character(first_last))) +
  labs(y = "Calendar day",
       colour = "First/last") +
  theme_bw() %+replace%
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        panel.grid = element_blank())

ice_duration <- airtemp %>%
  mutate(year = year(date)) %>%
  filter(!is.na(ice)) %>%
  dplyr::select(year, ice, date) %>%
  pivot_wider(names_from = ice, names_prefix = "ice_", values_from = date) %>%
  dplyr::select(year, ice_off, ice_on) %>%
  mutate(ice_duration = as.numeric(difftime(ice_off, lag(ice_on))),
         ice_free = as.numeric(difftime(ice_on, ice_off)))

(ice_duration_plot <- ice_duration %>%
    ggplot(aes(x = year, y = ice_duration)) +
    geom_point() +
    geom_line() +
    geom_smooth() +
    labs(y = "Ice cover duration (days)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          panel.grid = element_blank()))
#ggsave("figures/ela_ice_duration.pdf", ice_duration_plot, width = 5, height = 4, units = "in", device = cairo_pdf)

ice_duration %>%
  ggplot(aes(x = year, y = ice_free)) +
  geom_point() +
  geom_line() +
  geom_smooth() +
  labs(y = "Ice-free period (days)") +
  theme_bw() %+replace%
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        panel.grid = element_blank())


#### Summarize data ####
(climate_byyear <- airtemp_year_means %>%
   left_join(ice_duration %>%
               dplyr::select(year, ice_duration), by = "year") %>%
   left_join(precipitation_year_totals, by = "year") %>%
   rename(maxairtemp = mean_maxairtemp_year,
          minairtemp = mean_minairtemp_year,
          precip = total_precip_year))

# climate_byyear %>%
#   mutate(maxairtemp = round(maxairtemp, 1),
#          minairtemp = round(minairtemp, 1),
#          precip = round(precip, 0)) %>%
#   write_tsv("output/ela/ela_climate_byyear.tsv", col_names = TRUE)
