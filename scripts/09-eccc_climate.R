# Environment Canada weather station climate data
# Data accessed from: https://climate.weather.gc.ca/historical_data/search_historic_data_e.html

# Load libraries
library(janitor)
library(ggforce)
library(scales)

# Load environmental data
source("scripts/08-ela_env_data.R")


### Import and format data ####
# Define function to import and format Environment Canada historical daily climate data
# Flag legend:
# - B:	More than one occurrence and estimated
# - E:	Estimated
# - M:	Missing
# - S:	More than one occurrence
# - T:	Trace
# - [empty]:	Indicates an unobserved value
# - ^:	The value displayed is based on incomplete data
importDailyClimate <- function(file) {
  read_csv(file,
           col_names = TRUE,
           col_types = cols(x = col_double(),
                            y = col_double(),
                            STATION_NAME = col_character(),
                            CLIMATE_IDENTIFIER = col_double(),
                            ID = col_character(),
                            LOCAL_DATE = col_datetime(format = ""),
                            PROVINCE_CODE = col_character(),
                            LOCAL_YEAR = col_double(),
                            LOCAL_MONTH = col_double(),
                            LOCAL_DAY = col_double(),
                            MEAN_TEMPERATURE = col_double(),
                            MEAN_TEMPERATURE_FLAG = col_character(),
                            MIN_TEMPERATURE = col_double(),
                            MIN_TEMPERATURE_FLAG = col_character(),
                            MAX_TEMPERATURE = col_double(),
                            MAX_TEMPERATURE_FLAG = col_character(),
                            TOTAL_PRECIPITATION = col_double(),
                            TOTAL_PRECIPITATION_FLAG = col_character(),
                            TOTAL_RAIN = col_double(),
                            TOTAL_RAIN_FLAG = col_character(),
                            TOTAL_SNOW = col_double(),
                            TOTAL_SNOW_FLAG = col_character(),
                            SNOW_ON_GROUND = col_double(),
                            SNOW_ON_GROUND_FLAG = col_character(),
                            DIRECTION_MAX_GUST = col_double(),
                            DIRECTION_MAX_GUST_FLAG = col_character(),
                            SPEED_MAX_GUST = col_double(),
                            SPEED_MAX_GUST_FLAG = col_character(),
                            COOLING_DEGREE_DAYS = col_double(),
                            COOLING_DEGREE_DAYS_FLAG = col_character(),
                            HEATING_DEGREE_DAYS = col_double(),
                            HEATING_DEGREE_DAYS_FLAG = col_character(),
                            MIN_REL_HUMIDITY = col_double(),
                            MIN_REL_HUMIDITY_FLAG = col_character(),
                            MAX_REL_HUMIDITY = col_double(),
                            MAX_REL_HUMIDITY_FLAG = col_character())) %>%
    clean_names() %>%
    mutate(local_date = as_date(local_date)) %>%
    rename(date = local_date,
           year = local_year,
           month = local_month,
           day = local_day) %>%
    dplyr::select(station_name,
                  date,
                  contains("temperature"),
                  contains("precipitation"), contains("rain"), contains("snow")) %>%
    dplyr::select(-contains("on_ground")) %>%
    return()
}

# Import and format Environment Canada historical climate data
kenora_1883to1939.1 <- importDailyClimate("data/eccc_historical_climate/climate-daily_kenora_6034070_records1-10000.csv")
kenora_1883to1939.2 <- importDailyClimate("data/eccc_historical_climate/climate-daily_kenora_6034070_records10001-14877.csv")
kenora_1883to1939 <- bind_rows(kenora_1883to1939.1,
                               kenora_1883to1939.2)

kenoraa_1938to2013.1 <- importDailyClimate("data/eccc_historical_climate/climate-daily_kenoraa_6034075_records1-10000.csv")
kenoraa_1938to2013.2 <- importDailyClimate("data/eccc_historical_climate/climate-daily_kenoraa_6034075_records10001-20000.csv")
kenoraa_1938to2013.3 <- importDailyClimate("data/eccc_historical_climate/climate-daily_kenoraa_6034075_records20001-27211.csv")
kenoraa_2013to2023 <- importDailyClimate("data/eccc_historical_climate/climate-daily_kenoraa_6034076_records1-3688.csv") %>%
  filter(date != "2013-01-08")  # Remove entry duplicated in a different file
kenoraa_1938to2023 <- bind_rows(kenoraa_1938to2013.1,
                                kenoraa_1938to2013.2,
                                kenoraa_1938to2013.3,
                                kenoraa_2013to2023)


#### Summarize precipitation data ####
# Combine precipitation data
precipitation <- bind_rows(kenora_1883to1939,
                           kenoraa_1938to2023) %>%
  dplyr::select(-contains("temperature")) %>%
  filter(!is.na(total_precipitation) | !is.na(total_rain) | !is.na(total_snow)) %>%
  arrange(date) %>%
  filter(date < "2018-09-01")

# Calculate total precipitations
precipitation <- precipitation %>%
  mutate(total_rain = case_when(is.na(total_rain) ~ 0,
                                TRUE ~ total_rain),
         total_snow = case_when(is.na(total_snow) ~ 0,
                                TRUE ~ total_snow)) %>%
  mutate(total_precipitation = case_when(is.na(total_precipitation) ~ total_rain + total_snow,
                                         TRUE ~ total_precipitation)) %>%
  dplyr::select(station_name, date, total_precipitation)

# Assess redundant precipitation data
precipitation_nobs <- precipitation %>%
  group_by(date) %>%
  count(name = "nobs") %>%
  ungroup() %>%
  filter(nobs > 1)

# Assess compatibility of Kenora and Kenora A stations precipitation data for overlapping dates
precipitation_kenora_vs_kenoraa <- precipitation %>%
  filter(date %in% unique(precipitation_nobs$date)) %>%
  pivot_wider(names_from = station_name, values_from = total_precipitation)

precipitation_kenora_vs_kenoraa %>%
  ggplot() +
  geom_point(aes(x = KENORA, y = `KENORA A`)) +
  theme_bw()

cor(precipitation_kenora_vs_kenoraa$KENORA, precipitation_kenora_vs_kenoraa$`KENORA A`)

# Remove redundant data
precipitation_nr <- precipitation %>%
  filter(!(date %in% unique(precipitation_nobs$date) & station_name == "KENORA")) %>%
  dplyr::select(-station_name)
rm(precipitation)

precipitation_nr %>%
  calculateTimeDiffs()  # Daily measurements

# Identify months with measurements for only some days
(precipitation_partmonth <- precipitation_nr %>%
  partMonth())

# Fill in missing precipitation data
# Calculate mean daily total precipitation for months with missing days of precipitation data
precipitation_partmonth_data <- precipitation_nr %>%
  mutate(year = year(date),
         month = month(date),
         day = day(date)) %>%
  right_join(precipitation_partmonth, by = c("year", "month"))

precipitation_partmonth_means <- precipitation_partmonth_data %>%
  group_by(year, month) %>%
  summarize(mean_total_precipitation = mean(total_precipitation)) %>%
  ungroup()

# Calculate total monthly precipitation for months with missing days of precipitation data
# by summing total precipitation for days with real data and filling in missing days with the monthly mean
precipitation_partmonth_filled <- precipitation_partmonth_data %>%
  mutate(ndays_month = case_when(month %in% c(1, 3, 5, 7, 8, 10, 12) ~ 31,
                                 month %in% c(4, 6, 9, 11) ~ 30,
                                 month == 2 ~ 28)) %>%
  mutate(ndays_missing = ndays_month - ndays) %>%
  group_by(year, month, ndays, ndays_missing) %>%
  summarize(ndays_total_precipitation = sum(total_precipitation)) %>%
  ungroup() %>%
  left_join(precipitation_partmonth_means, by = c("year", "month")) %>%
  mutate(total_precipitation = ndays_total_precipitation + ndays_missing * mean_total_precipitation) %>%
  dplyr::select(year, month, total_precipitation)

# Calculate total monthly precipitation for months without missing days
# Combine with filled-in total monthly precipitation
precipitation_month_totals <- precipitation_nr %>%
  mutate(year_month = paste(year(date), month(date), sep = "_")) %>%
  left_join(precipitation_partmonth %>%
               mutate(year_month = paste(year, month, sep = "_")) %>%
                        dplyr::select(-c(year, month)), by = c("year_month")) %>%
  filter(is.na(ndays)) %>%
  dplyr::select(-c(year_month, ndays)) %>%
  mutate(year = year(date),
         month = month(date),
         day = day(date)) %>%
  group_by(year, month) %>%
  summarize(total_precipitation = sum(total_precipitation)) %>%
  ungroup() %>%
  bind_rows(precipitation_partmonth_filled) %>%
  arrange(year, month)

# Identify years with measurements for only some months
precipitation_month_totals %>%
  group_by(year) %>%
  count(name = "nmonths") %>%
  ungroup() %>%
  filter(nmonths < 12)

# Filter precipitation data year range
precipitation_month_totals <- precipitation_month_totals %>%
  filter(year >= 1900 & year < 2018)

# Interpolate total monthly precipitation across years
precipitation_month_totals_interpolated <- interpolateMonthYears(precipitation_month_totals, "total_precipitation", "approx")

precipitation_month_totals_interpolated %>%
  #mutate(date = as_date(paste0(year, str_pad(month, 2, "left", 0), 15, sep = "-"))) %>%
  ggplot() +
  facet_wrap(~month) +
  geom_point(aes(x = year, y = total_precipitation,
                 colour = as.character(interpolated),
                 shape = as.character(month))) +
  # scale_x_date(date_breaks = "1 year", date_labels = "%y",
  #              expand = c(0,0)) +
  scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "dodgerblue"),
                      guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(values = c("1" = 49, "2" = 50, "3" = 51, "4" = 52,
                                "5" = 53, "6" = 54, "7" = 55, "8" = 56,
                                "9" = 57, "10" = 79, "11" = 78, "12" = 68)) +
  labs(y = "Total monthly precipitation (mm)",
       colour = "Interpolated",
       shape = "Month") +
  theme_bw() %+replace%
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        panel.grid = element_blank())

# Calculate annual total precipitations
precipitation_year_totals <- precipitation_month_totals_interpolated %>%
  group_by(year) %>%
  summarize(total_precip_year = sum(total_precipitation)) %>%
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
#ggsave("figures/eccc_precipitation_year_totals.pdf", precipitation_year_totals_plot, width = 5, height = 4, units = "in", device = cairo_pdf)


#### Summarize air temperature data ####
# Combine air temperature data
airtemp <- bind_rows(kenora_1883to1939,
                     kenoraa_1938to2023) %>%
  dplyr::select(station_name, date, contains("temperature")) %>%
  filter(!is.na(mean_temperature) | !is.na(min_temperature) | !is.na(max_temperature)) %>%
  arrange(date) %>%
  filter(date < "2018-09-01")

# Compare n observations for mean vs. min vs. max temperature
airtemp %>%
  select(-contains("flag")) %>%
  pivot_longer(!c(station_name, date), names_to = "var", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(var) %>%
  count(name = "nobs")

# Select mean air temperature data
airtemp <- airtemp %>%
  dplyr::select(station_name, date, mean_temperature) %>%
  filter(!is.na(mean_temperature))

# Assess redundant mean air temperature data
airtemp_nobs <- airtemp %>%
  group_by(date) %>%
  count(name = "nobs") %>%
  ungroup() %>%
  filter(nobs > 1)

# Assess compatibility of Kenora and Kenora A stations air temperature data for overlapping dates
airtemp_kenora_vs_kenoraa <- airtemp %>%
  filter(date %in% unique(airtemp_nobs$date)) %>%
  pivot_wider(names_from = station_name, values_from = mean_temperature)

airtemp_kenora_vs_kenoraa %>%
  ggplot() +
  geom_point(aes(x = KENORA, y = `KENORA A`)) +
  theme_bw()

cor(airtemp_kenora_vs_kenoraa$KENORA, airtemp_kenora_vs_kenoraa$`KENORA A`)

# Remove redundant data
airtemp_nr <- airtemp %>%
  filter(!(date %in% unique(airtemp_nobs$date) & station_name == "KENORA")) %>%
  dplyr::select(-station_name)
rm(airtemp)

airtemp_nr %>%
  calculateTimeDiffs()  # Daily measurements

# Identify months with measurements for only some days
(airtemp_partmonth <- airtemp_nr %>%
    partMonth())

# Fill in missing air temperature data
# Calculate mean daily air temperature for months with missing days of air temperature data
airtemp_partmonth_data <- airtemp_nr %>%
  mutate(year = year(date),
         month = month(date),
         day = day(date)) %>%
  right_join(airtemp_partmonth, by = c("year", "month"))

# Calculate mean air temperature for months with missing days of air temperature data
# by averaging air temperature for days with real data and filling in missing days with the monthly mean
airtemp_partmonth_filled <- airtemp_partmonth_data %>%
  mutate(ndays_month = case_when(month %in% c(1, 3, 5, 7, 8, 10, 12) ~ 31,
                                 month %in% c(4, 6, 9, 11) ~ 30,
                                 month == 2 ~ 28)) %>%
  mutate(ndays_missing = ndays_month - ndays) %>%
  group_by(year, month, ndays, ndays_missing) %>%
  summarize(ndays_sum_temperature = sum(mean_temperature),
            ndays_mean_temperature = mean(mean_temperature)) %>%
  ungroup() %>%
  mutate(mean_temperature = (ndays_sum_temperature + (ndays_mean_temperature * ndays_missing))/(ndays + ndays_missing)) %>%
  dplyr::select(year, month, mean_temperature)

# Calculate mean air temperature for months without missing days
# Combine with filled-in mean monthly air temperature
airtemp_month_means <- airtemp_nr %>%
  mutate(year_month = paste(year(date), month(date), sep = "_")) %>%
  left_join(airtemp_partmonth %>%
              mutate(year_month = paste(year, month, sep = "_")) %>%
              dplyr::select(-c(year, month)), by = c("year_month")) %>%
  filter(is.na(ndays)) %>%
  dplyr::select(-c(year_month, ndays)) %>%
  mutate(year = year(date),
         month = month(date),
         day = day(date)) %>%
  group_by(year, month) %>%
  summarize(mean_temperature = mean(mean_temperature)) %>%
  ungroup() %>%
  bind_rows(airtemp_partmonth_filled) %>%
  arrange(year, month)

# Identify years with measurements for only some months
airtemp_month_means %>%
  group_by(year) %>%
  count(name = "nmonths") %>%
  ungroup() %>%
  filter(nmonths < 12)

# Filter air temperature data year range
airtemp_month_means <- airtemp_month_means %>%
  filter(year >= 1900 & year < 2018)

# Interpolate mean monthly (mean) air temperature across years
airtemp_month_means_interpolated <- interpolateMonthYears(airtemp_month_means, "mean_temperature", "approx")

airtemp_month_means_interpolated %>%
  #mutate(date = as_date(paste0(year, str_pad(month, 2, "left", 0), 15, sep = "-"))) %>%
  ggplot() +
  facet_wrap(~month) +
  geom_point(aes(x = year, y = mean_temperature,
                 colour = as.character(interpolated),
                 shape = as.character(month))) +
  # scale_x_date(date_breaks = "1 year", date_labels = "%y",
  #              expand = c(0,0)) +
  scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "dodgerblue"),
                      guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(values = c("1" = 49, "2" = 50, "3" = 51, "4" = 52,
                                "5" = 53, "6" = 54, "7" = 55, "8" = 56,
                                "9" = 57, "10" = 79, "11" = 78, "12" = 68)) +
  labs(y = "Monthly mean air temperature (\u00B0C)",
       colour = "Interpolated",
       shape = "Month") +
  theme_bw() %+replace%
  theme(axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        panel.grid = element_blank())

# Calculate annual mean air temperature
airtemp_year_means <- airtemp_month_means_interpolated %>%
  group_by(year) %>%
  summarize(mean_temperature_year = mean(mean_temperature)) %>%
  ungroup()

(airtemp_year_means_plot <- airtemp_year_means %>%
    ggplot(aes(x = year, y = mean_temperature_year)) +
    geom_line(colour = "steelblue", alpha = 0.7) +
    geom_point(aes(colour = mean_temperature_year)) +
    geom_smooth(colour = "steelblue", alpha = 0.2) +
    labs(y = "Annual mean air temperature (\u00B0C)") +
    scale_y_continuous(labels = comma) +
    scale_colour_distiller(palette = "RdBu",
                           direction = -1,
                           limits = c(min(airtemp_year_means$mean_temperature_year), max(airtemp_year_means$mean_temperature_year)),
                           #trans = "sqrt",
                           #breaks = c(0, 1, 2),
                           name = "Temp (\u00B0C)") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          panel.grid = element_blank()))
#ggsave("figures/eccc_airtemp_year_means.pdf", airtemp_year_means_plot, width = 5, height = 4, units = "in", device = cairo_pdf)


#### Summarize data ####
(climate_byyear <- airtemp_year_means %>%
   left_join(precipitation_year_totals, by = "year") %>%
   rename(meanairtemp = mean_temperature_year,
          precip = total_precip_year))

# climate_byyear %>%
#   mutate(meanairtemp = round(meanairtemp, 1),
#          precip = round(precip, 0)) %>%
#   write_tsv("output/environmental/eccc_climate_byyear.tsv", col_names = TRUE)


#### Compare ECCC and ELA climate data ####
# Import and format ELA climate data
ela_climate <- read_tsv("output/ela/ela_climate_byyear.tsv", col_names = TRUE)

# Join ECCC and ELA climate data
climate_byyear_all <- climate_byyear %>%
  rename(eccc_precip = precip) %>%
  left_join(ela_climate %>%
              rename(ela_precip = precip), by = "year")

# Compare annual total precipitation
climate_byyear_all %>%
  ggplot() +
  geom_point(aes(x = ela_precip, eccc_precip, colour = year)) +
  geom_text(aes(x = ela_precip, eccc_precip, label = year), nudge_y = 10) +
  theme_bw()

cor(climate_byyear_all$ela_precip, climate_byyear_all$eccc_precip, use = "complete.obs")

# Compare air temperature
climate_byyear_all %>%
  ggplot() +
  geom_point(aes(x = maxairtemp, meanairtemp, colour = year)) +
  geom_text(aes(x = maxairtemp, meanairtemp, label = year), nudge_y = 0.1) +
  theme_bw()

climate_byyear_all %>%
  ggplot() +
  geom_point(aes(x = minairtemp, meanairtemp, colour = year)) +
  geom_text(aes(x = minairtemp, meanairtemp, label = year), nudge_y = 0.1) +
  theme_bw()

cor(climate_byyear_all$maxairtemp, climate_byyear_all$meanairtemp, use = "complete.obs")
cor(climate_byyear_all$minairtemp, climate_byyear_all$meanairtemp, use = "complete.obs")
