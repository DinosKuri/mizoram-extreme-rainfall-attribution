# ==============================================================================
# 05_Context_map.R
# Purpose: Regional Context Map generation (Figure 2c).
# Maps elevation profiles and boundary contours for context.
# Author:
# Laldinchhuaha Khiangte
#
# Date:01-06-2026
# 2026
# ==============================================================================

library(here)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(elevatr)
library(ggplot2)
library(tidyterra)

out_dir <- here::here("outputs", "Figures")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Country boundaries
countries <- ne_countries(scale = 50, returnclass = "sf") %>%
  filter(admin %in% c("India","Bangladesh","Myanmar","China","Nepal","Bhutan"))

# Study boxes
miz_box <- st_as_sfc(st_bbox(c(xmin = 92, xmax = 93, ymin = 22, ymax = 24), crs = 4326))
cmip_box <- st_as_sfc(st_bbox(c(xmin = 91.6, xmax = 93.4, ymin = 22.1, ymax = 24.0), crs = 4326))

# Elevation data
elev_raster <- get_elev_raster(locations = st_as_sf(st_bbox(c(xmin = 80, xmax = 100,
                                                              ymin = 15, ymax = 30), crs = 4326)),
                               z = 5, clip = "locations")
elev_df <- as.data.frame(elev_raster, xy = TRUE)
colnames(elev_df) <- c("lon", "lat", "elevation")

# Plot Generation
p_context <- ggplot() +
  geom_raster(data = elev_df, aes(lon, lat, fill = elevation)) +
  scale_fill_gradientn(colors = terrain.colors(100), name = "Elev. (m)") +
  geom_sf(data = countries, fill = NA, color = "grey30", linewidth = 0.5) +
  geom_sf(data = miz_box, fill = NA, color = "#D55E00", linewidth = 1.2) +
  geom_sf(data = cmip_box, fill = NA, color = "#0072B2", linewidth = 1.0, linetype = "dashed") +
  annotate("segment", x = 88, xend = 92, y = 19, yend = 22,
           arrow = arrow(length = unit(0.3, "cm")), color = "darkblue", linewidth = 1.2) +
  annotate("text", x = 89, y = 20, label = "Monsoon flow", color = "darkblue", fontface = "italic") +
  annotate("text", x = 92.5, y = 24.2, label = "Mizoram", color = "#D55E00", fontface = "bold") +
  coord_sf(xlim = c(80, 100), ylim = c(15, 30), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude", title = "c | Regional context") +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0))

ggsave(file.path(out_dir, "Map_Regional_Context.png"), p_context, width = 7, height = 5, dpi = 600)
message("Context map generated.")
