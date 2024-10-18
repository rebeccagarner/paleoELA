# Map lakes

# Load libraries
library(tidyverse)
library(lubridate)
library(sf)
library(maps)
library(rnaturalearth)
library(ggspatial)
library(patchwork)

# Load palettes
source("scripts/00-palettes.R")


#### Import and format data ####
# Import and format coring sites
coring <- read_delim("data/lakepulse/LakePulse_site_qc_20230419.csv", delim = ";", col_names = TRUE,
                     locale = locale(decimal_mark = ".", grouping_mark = ";", encoding = "UTF-8")) %>%
  filter(site_name == "Coring Site")

coring_l226n <- tibble(lake_id = "L226N",
                       coord_comment = coring$coord_comment[which(coring$lakepulse_id == "06-304")]) %>%
  mutate(coord_comment = str_remove(coord_comment, ".*North basin info: ")) %>%
  mutate(coord_comment = str_remove(coord_comment, ".$")) %>%
  separate(coord_comment, into = c("site_name", "sampling_date", "lat", "long", "depth"), sep = "; ") %>%
  mutate(across(everything(), ~str_remove(.x, ".*\\="))) %>%
  dplyr::select(-site_name) %>%
  mutate(sampling_date = as_date(sampling_date),
         lat = as.double(lat),
         long = as.double(long),
         depth = as.double(depth))

coring <- coring %>%
  mutate(lake_id = case_when(lakepulse_id == "06-311" ~ "L224",
                             lakepulse_id == "06-308" ~ "L373",
                             lakepulse_id == "06-304" ~ "L226S",  # L226 South coring information
                             lakepulse_id == "06-307" ~ "L227")) %>%
  relocate(lake_id, 1) %>%
  filter(!is.na(lake_id)) %>%
  dplyr::select(lake_id, sampling_date, lat, long, depth) %>%
  bind_rows(coring_l226n) %>%
  arrange(factor(lake_id, levels = names(palette_lake)))

# Import hydro features of Ontario shapefile (https://ftp.maps.canada.ca/pub/nrcan_rncan/vector/canvec/shp/Hydro/)
on_hydro <- st_read("C:/Users/Gandalf/Downloads/canvec_50K_ON_Hydro_shp/canvec_50K_ON_Hydro/waterbody_2_1.shp")
on_hydro <- st_read("~/Desktop/gis/canvec_50K_ON_Hydro_shp/canvec_50K_ON_Hydro/waterbody_2_1.shp")
st_crs(on_hydro)  # Check coordinate reference system (CRS)

# Reproject the CRS
on_hydro <- st_transform(on_hydro, crs = 4326)
# Bounding box:  xmin: -95.5 ymin: 44.99136 xmax: -74.14911 ymax: 57

on_hydro_crop <- on_hydro %>%
  st_crop(c(xmin = -93.85, xmax = -93.65, ymin = 49.6, ymax = 49.8))

# Import hydro features of Canada shapefile (https://ftp.maps.canada.ca/pub/nrcan_rncan/vector/canvec/shp/Hydro/)
canada_hydro <- st_read("C:/Users/Gandalf/Downloads/canvec_1M_CA_Hydro_shp/canvec_1M_CA_Hydro/waterbody_2.shp")
canada_hydro <- st_read("~/Desktop/gis/canvec_1M_CA_Hydro_shp/canvec_1M_CA_Hydro/waterbody_2.shp")

st_crs(canada_hydro)  # Check coordinate reference system (CRS)
canada_hydro <- st_transform(canada_hydro, crs = 3348)

# Import provinces/territories of Canada shapefile
# provinces <- st_read("C:/Users/Gandalf/Downloads/gpr_000b11a_e/gpr_000b11a_e.shp")
# st_crs(provinces)  # Check coordinate reference system (CRS)
# provinces <- st_transform(provinces, crs = 3348)

# Access and format Canadian land mass map
land <- sf::st_as_sf(maps::map("world",
                               region = c("Canada", "USA", "Greenland"),
                               plot = FALSE, fill = TRUE))
st_crs(land)  # Check CRS

# Access and format North American province/state map
canada <- ne_states(country = "canada", returnclass = "sf")
st_crs(canada)
usa <- ne_states("united states of america", returnclass = "sf")
st_crs(usa)
northamerica <- rbind(canada, usa) %>%
  st_transform(crs = 3348)


#### Format map data ####
# Define lake coordinates
lake_coords <- tibble(lake_id = c("L226", "L227", "L224", "L373"),
                      latitude = c(49.6886, 49.6879, 49.6903, 49.7438),
                      longitude = c(-93.7475, -93.6889, -93.7173, -93.7999)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Extract study lake polygons
lake_features <- st_join(lake_coords, on_hydro_crop)

lake_feature_ids <- lake_features %>%
  as_tibble() %>%
  dplyr::select(lake_id, feature_id)

lake_polygons <- on_hydro_crop %>%
  filter(feature_id %in% unique(lake_features$feature_id)) %>%
  left_join(lake_feature_ids, by = "feature_id")


#### Plot lakes map ####
# Map study lakes
(lakes_map <- ggplot() +
   geom_sf(data = land,
           fill = "#F7F7F7", colour = NA) +
   geom_sf(data = on_hydro_crop %>%
             filter(!feature_id %in% lake_features$feature_id),
           fill = "#ACE8FF", colour = NA) +
   geom_sf(data = lake_polygons,
           aes(fill = lake_id), colour = "white") +
   scale_fill_manual(values = palette_lake) +
   coord_sf(xlim = c(-93.68, -93.81), ylim = c(49.68, 49.75),
            expand = FALSE) +
   labs(fill = "Lake") +
   annotation_scale(style = "ticks") +
   annotation_north_arrow(location = "tr") +
   theme(axis.text = element_blank(),
         axis.ticks = element_blank(),
         panel.border = element_rect(colour = "black", fill = NA)))


#### Plot map of region ####
# Check map limits for Canada Lambert CRS
worldmap <- ne_countries(scale = "medium", type = "map_units", returnclass = "sf")

zoom_to <- c(-91.73, 50.75)
zoom_level <- 3.5

target_crs <- 3348

C <- 40075016.686   # ~ circumference of Earth in meters
x_span <- C / 2^zoom_level
y_span <- C / 2^(zoom_level+1)

zoom_to_xy <- st_transform(st_sfc(st_point(zoom_to), crs = 4326), crs = target_crs)

disp_window <- st_sfc(st_point(st_coordinates(zoom_to_xy - c(x_span / 2, y_span / 2))),
                      st_point(st_coordinates(zoom_to_xy + c(x_span/2, y_span/2))), crs = target_crs)

ggplot() +
  geom_sf(data = worldmap) +
  geom_sf(data = zoom_to_xy, color = 'red') +
  coord_sf(xlim = st_coordinates(disp_window)[,"X"],
           ylim = st_coordinates(disp_window)[,"Y"],
           crs = target_crs, datum = target_crs) +
  theme_bw()

# Calculate median ELA lakes coordinate
latitude_median <- median(c(49.6886, 49.6879, 49.6903, 49.7438))
longitude_median <- median(c(-93.7475, -93.6889, -93.7173, -93.7999))

ela_coords <- tibble(site = "ELA",
                     latitude = latitude_median,
                     longitude = longitude_median) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = 3348)

# Convert CRS to Canada Lambert
land_lambert <- st_transform(land, crs = 3348)

# Map region with provinces
(region_map <- ggplot() +
    geom_sf(data = northamerica,
            fill = "#F7F7F7", colour = "#d1cfcf", linewidth = 0.1) +
    geom_sf(data = canada_hydro,
            fill = "#ACE8FF", colour = NA) +
    geom_sf(data = ela_coords, colour = "black", size = 3, shape = 15) +
    coord_sf(xlim = c(5500000, 7200000), ylim = c(1200000, 2400000),
             expand = FALSE) +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          panel.grid = element_line(colour = "white"),
          panel.background = element_rect(fill = "#CFE3EB"),
          panel.border = element_rect(colour = "black", fill = NA)))

# Combine maps
(maps_all <- region_map / lakes_map)
#ggsave("figures/ela_lakes_map.pdf", maps_all, width = 5, height = 7, units = "in", device = cairo_pdf)
