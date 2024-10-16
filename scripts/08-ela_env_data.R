# ELA environmental data

# setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
# setwd("~/Dropbox/projects/ela18s/")

# Load libraries
library(tidyverse)
library(readxl)
library(lubridate)

# Load palettes
source("scripts/00-palettes.R")


#### Import and format data ####
# Define function for filtering study lakes up until date of core sampling
filterObs <- function(data) {
  data %>%
    filter(lake_id %in% c(224, 226, 227, 373)) %>%
    filter(date < "2018-09-01") %>%
    arrange(lake_id, date) %>%
    distinct() %>%
    return()
}

# Define function to format lake IDs
formatLakeID <- function(data) {
  data_tmp <- data
  if ("lake_id" %in% names(data_tmp) & "sublocation" %in% names(data_tmp)) {
    data_tmp <- data_tmp %>%
      unite("lake_id", c(lake_id, sublocation), sep = "_") %>%
      mutate(lake_id = case_when(lake_id == "224_LA" ~ "L224",
                                 lake_id == "226_NEB" ~ "L226N",
                                 lake_id == "226_SWB" ~ "L226S",
                                 lake_id == "227_LA" ~ "L227",
                                 lake_id == "373_LA" ~ "L373"))
  } else if ("lake_id" %in% names(data_tmp) & !("sublocation" %in% names(data_tmp))) {
    data_tmp <- data_tmp %>%
      mutate(lake_id = str_c("L", lake_id))
  }
  return(data_tmp)
}

# Summarize lake maximum depths
lake_depths <- c("L226N" = 14.7,
                 "L226S" = 11.6,
                 "L227" = 10.0,
                 "L224" = 27.3,
                 "L373" = 21.2) %>%
  enframe("lake_id", "max_depth")

# Import and format experimental manipulation data
ela_manipulations <- read_tsv("data/ela/ela_manipulations.txt", col_names = TRUE, comment = "#")

ela_experiments <- bind_rows(tibble(lake_id = c(rep("L227", 2018 - 1970 + 1)),
                                    year = seq(1970, 2018, 1),
                                    experiment = c(rep("12N:P", 1974 - 1970 + 1),
                                                   rep("4N:P", 1989 - 1975 + 1),
                                                   rep("P only", 2018 - 1990 + 1))),
                             
                             tibble(lake_id = c(rep("L226N", 1980 - 1973 + 1),
                                                rep("L226S", 1980 - 1973 + 1)),
                                    year = rep(seq(1973, 1980, 1), 2),
                                    experiment = c(rep("C+N+P", 1980 - 1973 + 1),
                                                   rep("C+N", 1980 - 1973 + 1))))

# Import and format water temperature data
watertemp <- read_csv("data/ela/tempprofiles_garner.csv", col_names = TRUE) %>%
  rename(lake_id = lake,
         temp = temp_C) %>%
  filterObs() %>%
  formatLakeID()

# Import and format water chemistry data
# Filter study dates within study period (prior to core sampling)
# Filter sampling depths below lake maximum depth
chemistry <- read_csv("data/ela/all_lakes_chem_ela_sample_avgs_Garner.csv", col_names = TRUE,
                      col_types = cols(ph_equil = "d", fe_mg_L = "d")) %>%
  dplyr::select(-stn, -year, -month, -month_name, -jd) %>%
  rename(lake_id = lake,
         sublocation = subloc,
         start_depth = start_d,
         stratum = strat_layer) %>%
  filterObs() %>%
  formatLakeID() %>%
  mutate(stratum = str_to_lower(stratum)) %>%
  left_join(lake_depths, by = "lake_id") %>%
  filter(!(!is.na(start_depth) & start_depth > max_depth)) %>%
  dplyr::select(-max_depth)

# Import and format Secchi disk depth data
# kd: light attenuation coefficient
secchi <- read_csv("data/ela/light_secchi_epi_Garner.csv") %>%
  dplyr::select(-year, -station) %>%
  rename(lake_id = lake) %>%
  filterObs() %>%
  formatLakeID() %>%
  group_by(lake_id, date) %>%
  summarize(epi_depth_m = mean(epi_depth_m, na.rm = TRUE),
            secchi_m = mean(secchi_m, na.rm = TRUE),
            kd = mean(kd, na.rm = TRUE)) %>%
  ungroup()
secchi[is.na(secchi)] <- NA

# Import and format bathymetry data
bathymetry <- read_csv("data/ela/bathymetry_ELA_Garner.csv", col_names = TRUE) %>%
  rename(lake_id = lake) %>%
  filter(!is.na(lake_id)) %>%
  formatLakeID()

# Import and format air temperature data
airtemp <- read_csv("data/ela/air_temperatures_to_2021.csv", col_names = TRUE) %>%
  dplyr::select(-location, -sublocation, -station, -Year, -Month, -JD) %>%
  rename(ice = `ice-on-off`) %>%
  filter(date < "2018-09-01") %>%
  filter(!is.na(max_temp) | !is.na(min_temp)) %>%
  mutate(ice = case_when(ice == 1 ~ "off",
                         ice == 2 ~ "on")) %>%
  arrange(date)

# Import and format precipitation data
precipitation <- read_csv("data/ela/ela_precip_to_2021.csv", col_names = TRUE) %>%
  dplyr::select(-location, -sublocation, -station, -trace) %>%
  filter(date < "2018-09-01") %>%
  mutate(tot_precip_mm = rain_mm + snow_mm) %>%
  arrange(date)


#### Define functions to summarize data ####
# Define function to assess sampling intervals
calculateTimeDiffs <- function(data) {
  data %>%
    mutate(time_diff = difftime(date, lag(date), units = "days")) %>%
    filter(!is.na(time_diff)) %>%
    group_by(time_diff) %>%
    count() %>%
    return()
}

# Define function to calculate time intervals (in days) between sampling
calculateIntervals <- function(data, names) {
  data %>%
    arrange(date) %>%
    mutate(date2 = lead(date)) %>%
    filter(!is.na(date2)) %>%
    mutate(time_diff = as.numeric(difftime(date2, date, units = "days"))) %>%
    mutate(count_order = row_number()) %>%
    return()
}

# Define function to assess sampling frequency
calculateDataIntervals <- function(data) {
  data_split <- data %>%
    distinct(lake_id, date) %>%
    mutate(year = year(date)) %>%
    relocate(year, .after = lake_id) %>%
    group_split(lake_id, year)
  
  data_split %>%
    map(~mutate(., lakeid_year = str_c(lake_id, "_", year))) %>%
    map(~pull(.,lakeid_year)) %>%  # Pull out variable
    map(~as.character(.)) %>%  # Convert factor to character
    map(~unique(.)) -> names(data_split)
  
  list(data = data_split, names = names(data_split)) %>%
    pmap(calculateIntervals) %>%
    map_dfr(`[`, c("lake_id", "year", "date", "date2", "time_diff", "count_order")) %>%
    return()
}

# Define function to summarize intervals between sampling
summarizeIntervals <- function(data) {
  data %>%
    calculateDataIntervals() %>%
    group_by(lake_id, year) %>%
    summarize(start_date = min(date),
              end_date = max(date2),
              n_sampling = n(),
              min_interval = min(time_diff),
              max_interval = max(time_diff),
              mean_interval = mean(time_diff)) %>%
    ungroup() %>%
    mutate(sampling_period = as.numeric(difftime(end_date, start_date, units = "days"))) %>%
    return()
}

# Define function to identify partially sampled months (some days only)
partMonth <- function(data) {
  data %>%
    mutate(year = year(date),
           month = month(date),
           day = day(date)) %>%
    group_by(year, month) %>%
    count(name = "ndays") %>%
    ungroup() %>%
    filter(ndays < 30 & !(month == 2 & ndays >= 28)) %>%
    arrange(ndays) %>%
    return()
}

# Define function to identify partially sampled years (some months only)
partYear <- function(data) {
  data %>%
    mutate(year = year(date),
           month = month(date)) %>%
    distinct(year, month) %>%
    group_by(year) %>%
    count(name = "nmonths") %>%
    ungroup() %>%
    filter(nmonths < 12) %>%
    arrange(nmonths) %>%
    return()
}

# Define function to assign season to date
assignSeason <- function(data) {
  data %>%
    mutate(season = case_when(date >= paste0(year(date), "-03-20") & date <= paste0(year(date), "-06-20") ~ "Spring",
                              date >= paste0(year(date), "-06-21") & date <= paste0(year(date), "-09-21") ~ "Summer",
                              date >= paste0(year(date), "-09-22") & date <= paste0(year(date), "-12-20") ~ "Fall",
                              (date >= paste0(year(date), "-01-01") & date <= paste0(year(date), "-03-19")) |
                                (date >= paste0(year(date), "-12-21") & date <= paste0(year(date), "-12-31")) ~ "Winter")) %>%
    return()
}

# Define function to exclude winter sampling dates
excludeWinter <- function(data) {
  data %>%
    assignSeason() %>%
    filter(season != "Winter") %>%
    dplyr::select(-season) %>%
    return()
}

# Define function to format season data
formatSeasons <- function(data) {
  # Extract sampling start and end dates
  start_date <- min(data$date)
  end_date <- max(data$date)
  
  # Create vector of sampling years
  sampling_years <- seq(from = year(start_date), to = year(end_date), by = 1)
  
  # Format seasons
  season_dates <- tibble(year = rep(sampling_years, 5)) %>%
    arrange(year) %>%
    mutate(season = rep(c("Winter", "Spring", "Summer", "Fall", "Winter"), length(sampling_years)),
           season_start_month_day = rep(c("01-01", "03-20", "06-21", "09-22", "12-21"), length(sampling_years)),
           season_end_month_day = rep(c("03-19", "06-20", "09-21", "12-20", "12-31"), length(sampling_years))) %>%
    mutate(season_start_date = date(str_c(year, "-", season_start_month_day)),
           season_end_date = date(str_c(year, "-", season_end_month_day)))
  
  return(season_dates)
}

# Define function to summarize data by month
summarizeByMonth <- function(data, var) {
  data %>%
    mutate(year = year(date),
           month = month(date)) %>%
    group_by(lake_id, year, month) %>%
    summarize(mean = mean(eval(as.name(paste(var)))),
              sd = sd(eval(as.name(paste(var)))),
              ndays = n()) %>%
    ungroup() %>%
    arrange(lake_id, year, month) %>%
    return()
}

# Define function to summarize data by month-season
summarizeByMonthSeason <- function(data, var) {
  data %>%
    summarizeByMonth(var) %>%
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
    group_by(lake_id, year, month_season) %>%
    summarize(mean = mean(mean),
              nmonths = n(),
              ndays = sum(ndays)) %>%
    ungroup() %>%
    arrange(lake_id, year, factor(month_season, levels = names(palette_season))) %>%
    return()
}

# Define function to summarize data by year
summarizeByYear <- function(data, var) {
  data %>%
    summarizeByMonthSeason(var) %>%
    filter(month_season %in% c("Spring", "Summer", "Fall")) %>%
    group_by(lake_id, year) %>%
    summarize(mean = mean(mean),
              nseasons = n(),
              nmonths = sum(nmonths),
              ndays = sum(ndays)) %>%
    ungroup() %>%
    arrange(lake_id, year) %>%
    return()
}

# Define function to interpolate variable for missing months across years
interpolateMonthYears <- function(data, var, method = "approx") {
  data_tmp <- data %>%
    arrange(year, month)
  
  year_range <- seq(data_tmp$year[1], tail(data_tmp$year, n = 1), 1)
  
  for (i in 1:12) {
    month_tmp <- i
    data_month <- data_tmp %>%
      filter(month == month_tmp)
    
    if (method == "approx") {
      interpolated <- approx(data_month$year, pull(data_month[var]), year_range) %>%
        as_tibble() %>%
        rename(year = x,
               {{var}} := y) %>%
        filter(!year %in% unique(data_month$year)) %>%
        mutate(month = month_tmp,
               interpolated = TRUE)
    } else if (method == "spline") {
      interpolated <- spline(data_month$year, pull(data_month[var]),
                             xout = year_range,
                             method = "fmm") %>%
        as_tibble() %>%
        rename(year = x,
               {{var}} := y) %>%
        filter(!year %in% unique(data_month$year)) %>%
        mutate(month = month_tmp,
               interpolated = TRUE)
    }
    
    data_tmp <- data_tmp %>%
      bind_rows(interpolated) %>%
      mutate(interpolated = case_when(is.na(interpolated) ~ FALSE,
                                      TRUE ~ interpolated))
  }
  return(data_tmp)
}

# Define functions to interpolate variable mean for missing month-seasons
interpolateMonthSeasons <- function(data, names, method = "approx") {
  data_tmp <- data %>%
    arrange(year)
  
  lakeid_tmp <- unique(data_tmp$lake_id)
  monthseason_tmp <- unique(data_tmp$month_season)
  
  year_range <- seq(data_tmp$year[1], tail(data_tmp$year, n = 1), 1)
  
  if (method == "approx") {
    interpolated <- approx(data_tmp$year, data_tmp$mean, year_range) %>%
      as_tibble() %>%
      rename(year = x,
             mean = y) %>%
      filter(!year %in% unique(data_tmp$year)) %>%
      mutate(lake_id = lakeid_tmp,
             month_season = monthseason_tmp)
  } else if (method == "spline") {
    interpolated <- spline(data_tmp$year, data_tmp$mean,
                           xout = year_range,
                           method = "fmm") %>%
      as_tibble() %>%
      rename(year = x,
             mean = y) %>%
      filter(!year %in% unique(data_tmp$year)) %>%
      mutate(lake_id = lakeid_tmp,
             month_season = monthseason_tmp)
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

interpolateMonthSeasonsAll <- function(data, var, method = "approx") {
  data_split <- data %>%
    summarizeByMonthSeason(var) %>%
    group_split(lake_id, month_season)
  
  data_split %>%
    map(~mutate(., lakeid_monthseason = str_c(lake_id, "_", month_season))) %>%
    map(~pull(.,lakeid_monthseason)) %>%  # Pull out variable
    map(~as.character(.)) %>%  # Convert factor to character
    map(~unique(.)) -> names(data_split)
  
  list(data = data_split, names = names(data_split), method = method) %>%
    pmap(interpolateMonthSeasons) %>%
    map_dfr(`[`, c("lake_id", "year", "month_season", "mean", "nmonths", "ndays", "interpolated")) %>%
    return()
}

# Define function to summarize interpolated data by year
summarizeByYearInterpolated <- function(data, var, method = "approx") {
  data %>%
    interpolateMonthSeasonsAll(var, "approx") %>%
    group_by(lake_id, year) %>%
    summarize(mean = mean(mean),
              nseasons = n(),
              nmonths = sum(nmonths),
              ndays = sum(ndays)) %>%
    ungroup() %>%
    arrange(lake_id, year) %>%
    return()
}

# Define function to summarize data by month for water column strata
summarizeByMonthWC <- function(data, var) {
  data %>%
    mutate(year = year(date),
           month = month(date)) %>%
    group_by(lake_id, year, month, stratum) %>%
    summarize(mean = mean(eval(as.name(paste(var)))),
              sd = sd(eval(as.name(paste(var)))),
              ndays = n()) %>%
    ungroup() %>%
    arrange(lake_id, year, month,
            factor(stratum, levels = c("Epilimnion", "Metalimnion", "Hypolimnion", "Profundal"))) %>%
    return()
}

# Define function to summarize data by month-season for water column strata
summarizeByMonthSeasonWC <- function(data, var) {
  data %>%
    summarizeByMonthWC(var) %>%
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
    group_by(lake_id, year, month_season, stratum) %>%
    summarize(mean = mean(mean),
              nmonths = n(),
              ndays = sum(ndays)) %>%
    ungroup() %>%
    arrange(lake_id, year, factor(month_season, levels = names(palette_season)),
            factor(stratum, levels = c("Epilimnion", "Metalimnion", "Hypolimnion", "Profundal"))) %>%
    return()
}

# Define function to summarize data by year for water column strata
summarizeByYearWC <- function(data, var) {
  data %>%
    summarizeByMonthSeasonWC(var) %>%
    filter(month_season %in% c("Spring", "Summer", "Fall")) %>%  # Remove winter sampling dates
    group_by(lake_id, year, stratum) %>%
    summarize(mean = mean(mean),
              nseasons = n(),
              nmonths = sum(nmonths),
              ndays = sum(ndays)) %>%
    ungroup() %>%
    arrange(lake_id, year,
            factor(stratum, levels = c("Epilimnion", "Metalimnion", "Hypolimnion", "Profundal"))) %>%
    return()
}

# Define functions to interpolate variable mean for missing month-seasons
interpolateMonthSeasonsWC <- function(data, names, method = "approx") {
  data_tmp <- data %>%
    arrange(year)
  
  lakeid_tmp <- unique(data_tmp$lake_id)
  monthseason_tmp <- unique(data_tmp$month_season)
  stratum_tmp <- unique(data_tmp$stratum)
  
  #print(paste0(lakeid_tmp, "_", monthseason_tmp, "_", stratum_tmp))
  
  year_range <- seq(data_tmp$year[1], tail(data_tmp$year, n = 1), 1)
  
  if (method == "approx") {
    interpolated <- approx(data_tmp$year, data_tmp$mean, year_range) %>%
      as_tibble() %>%
      rename(year = x,
             mean = y) %>%
      filter(!year %in% unique(data_tmp$year)) %>%
      mutate(lake_id = lakeid_tmp,
             month_season = monthseason_tmp,
             stratum = stratum_tmp)
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
             stratum = stratum_tmp)
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

interpolateMonthSeasonsAllWC <- function(data, var, method = "approx") {
  data_split <- data %>%
    summarizeByMonthSeasonWC(var) %>%
    filter(month_season %in% c("Spring", "Summer", "Fall")) %>%  # Remove winter sampling dates
    group_split(lake_id, month_season, stratum)
  
  data_split %>%
    map(~mutate(., lakeid_monthseason_stratum = str_c(lake_id, "_", month_season, "_", stratum))) %>%
    map(~pull(.,lakeid_monthseason_stratum)) %>%  # Pull out variable
    map(~as.character(.)) %>%  # Convert factor to character
    map(~unique(.)) -> names(data_split)
  
  list(data = data_split, names = names(data_split), method = method) %>%
    pmap(interpolateMonthSeasonsWC) %>%
    map_dfr(`[`, c("lake_id", "year", "month_season", "stratum", "mean", "nmonths", "ndays", "interpolated")) %>%
    return()
}

# Define function to summarize interpolated data by year for water column strata
summarizeByYearInterpolatedWC <- function(data, var, method = "approx") {
  data %>%
    interpolateMonthSeasonsAllWC(var, "approx") %>%
    group_by(lake_id, year, stratum) %>%
    summarize(mean = mean(mean),
              nseasons = n(),
              nmonths = sum(nmonths),
              ndays = sum(ndays)) %>%
    ungroup() %>%
    arrange(lake_id, year,
            factor(stratum, levels = c("Epilimnion", "Metalimnion", "Hypolimnion", "Profundal"))) %>%
    return()
}

# Define function to filter water chemistry observations
filterChemObs <- function(data) {
  # Identify observations by water column stratum
  strata_obs <- data %>%
    distinct(lake_id, date, stratum) %>%
    mutate(count = TRUE) %>%
    pivot_wider(names_from = stratum, values_from = count, values_fill = FALSE)
  
  # Elimate eup and hyp observations for dates with epi/pro observations
  data_filtered <- data %>%
    left_join(strata_obs, by = c("lake_id", "date")) %>%
    filter(!(stratum == "hyp" & pro & hyp)) %>%  # Discard hyp when there is also pro
    filter(!(stratum == "eup" & epi & eup)) %>%  # Discard eup when there is also epi
    filter(!(stratum == "all" & !epi & !pro & !eup & !hyp & all))  # Discard observations when only all is available
  
  # Reassign strata based on availability of observations
  data_filtered <- data_filtered %>%
    mutate(stratum_mod = case_when(stratum == "all" & epi & !pro & !eup & !hyp & all ~ "pro",
                                   
                                   stratum == "epi" ~ "epi",
                                   stratum == "eup" & !epi & eup ~ "epi",
                                   
                                   stratum == "pro" ~ "pro",
                                   stratum == "hyp" & !pro & hyp ~ "pro")) %>%
    dplyr::select(-c(stratum, epi, pro, eup, hyp, all)) %>%
    rename(stratum = stratum_mod) %>%
    relocate(stratum, .before = start_depth) %>%
    mutate(stratum = case_when(stratum == "epi" ~ "Epilimnion",
                               stratum == "pro" ~ "Profundal")) %>%
    arrange(lake_id, date, stratum, start_depth)
  
  return(data_filtered)
}
