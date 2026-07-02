# ==============================================================================
# 02_ERA5_circulation.R
# Purpose: Downloads ERA5 monthly means, computes 850 hPa moisture flux,
# performs Pettitt test, and creates composite map.
# Author:
# Laldinchhuaha Khiangte
#
# Date:01-06-2026
# 2026
# ==============================================================================

library(here)
library(ecmwfr)
library(terra)
library(trend)
library(ggplot2)
library(tidyterra)
library(patchwork)
library(dplyr)

# ------------------------------------------------------------------------------
# ENVIRONMENT & DIRECTORIES
# ------------------------------------------------------------------------------
# Security Update: Use environment variable for CDS API Key
cds_key <- Sys.getenv("CDS_API_KEY")
if (cds_key == "") {
  stop("CDS_API_KEY environment variable is not set. Please set it before running this script.")
}
wf_set_key(key = cds_key)

# Configure directories
dl_dir <- here::here("data", "ERA5")
out_dir <- here::here("outputs", "Figures")
if(!dir.exists(dl_dir)) dir.create(dl_dir, recursive = TRUE)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# DATA DOWNLOAD
# ------------------------------------------------------------------------------
file_pressure <- file.path(dl_dir, "ERA5_monsoon_pressure_levels.nc")
if (!file.exists(file_pressure)) {
  request_pressure <- list(
    product_type       = "monthly_averaged_reanalysis",
    variable           = c("u_component_of_wind","v_component_of_wind",
                           "specific_humidity","vertical_velocity"),
    pressure_level     = c("850","500"),
    year               = as.character(1950:2025),
    month              = c("06","07","08","09"),
    time               = "00:00",
    area               = c(30,80,15,100),
    format             = "netcdf",
    dataset_short_name = "reanalysis-era5-pressure-levels-monthly-means",
    target             = "ERA5_monsoon_pressure_levels.nc"
  )
  wf_request(request = request_pressure, transfer = TRUE, path = dl_dir, time_out = 3600)
}

file_surface <- file.path(dl_dir, "ERA5_monsoon_surface_level.nc")
if (!file.exists(file_surface)) {
  request_surface <- list(
    product_type       = "monthly_averaged_reanalysis",
    variable           = "mean_sea_level_pressure",
    year               = as.character(1950:2025),
    month              = c("06","07","08","09"),
    time               = "00:00",
    area               = c(30,80,15,100),
    format             = "netcdf",
    dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
    target             = "ERA5_monsoon_surface_level.nc"
  )
  wf_request(request = request_surface, transfer = TRUE, path = dl_dir, time_out = 3600)
}

# ------------------------------------------------------------------------------
# DATA PROCESSING
# ------------------------------------------------------------------------------
era_p <- rast(file_pressure)
era_s <- rast(file_surface)

u850 <- subset(era_p, grep("u_pressure_level=850", names(era_p)))
v850 <- subset(era_p, grep("v_pressure_level=850", names(era_p)))
q850 <- subset(era_p, grep("q_pressure_level=850", names(era_p)))
w500 <- subset(era_p, grep("w_pressure_level=500", names(era_p)))
mslp <- era_s

# JJAS seasonal means
years <- rep(1950:2025, each = 4)
u850_jjas <- tapp(u850, years, mean, na.rm = TRUE)
v850_jjas <- tapp(v850, years, mean, na.rm = TRUE)
q850_jjas <- tapp(q850, years, mean, na.rm = TRUE)
w500_jjas <- tapp(w500, years, mean, na.rm = TRUE)
mslp_jjas <- tapp(mslp, years, mean, na.rm = TRUE)

# Moisture flux
qu <- u850_jjas * q850_jjas
qv <- v850_jjas * q850_jjas

# ------------------------------------------------------------------------------
# STATISTICAL TESTING
# ------------------------------------------------------------------------------
inner_extent <- ext(92, 93, 22, 24)
flux_mag <- sqrt(qu^2 + qv^2)
flux_ts <- as.numeric(global(crop(flux_mag, inner_extent), "mean", na.rm = TRUE)$mean)
cat("\nPettitt test for Flux Magnitude:\n")
print(pettitt.test(flux_ts))

# Moisture convergence (central difference)
e_grad <- ext(91, 94, 21, 25)
qu_crop <- crop(qu, e_grad)
qv_crop <- crop(qv, e_grad)
filter_x <- matrix(c(-0.5,0,0.5), nrow=1)
filter_y <- matrix(c(-0.5,0,0.5), ncol=1)

mc_list <- list()
for (i in seq_len(nlyr(qu_crop))) {
  dx <- focal(qu_crop[[i]], filter_x, na.rm = TRUE)
  dy <- focal(qv_crop[[i]], filter_y, na.rm = TRUE)
  mc_list[[i]] <- -1 * (dx + dy)
}
mc <- rast(mc_list)
mc_ts <- as.numeric(global(crop(mc, inner_extent), "mean", na.rm = TRUE)$mean)
cat("\nPettitt test for Moisture Convergence:\n")
print(pettitt.test(mc_ts))

# Vertical velocity (omega) test
w_ts <- as.numeric(global(crop(w500_jjas, inner_extent), "mean", na.rm = TRUE)$mean)
cat("\nPettitt test for Vertical Velocity:\n")
print(pettitt.test(w_ts))

# ------------------------------------------------------------------------------
# COMPOSITE MAPPING
# ------------------------------------------------------------------------------
break_year <- 2003
idx_pre  <- which(1950:2025 < break_year)
idx_post <- which(1950:2025 >= break_year)

qu_pre  <- mean(qu[[idx_pre]], na.rm = TRUE)
qu_post <- mean(qu[[idx_post]], na.rm = TRUE)
qv_pre  <- mean(qv[[idx_pre]], na.rm = TRUE)
qv_post <- mean(qv[[idx_post]], na.rm = TRUE)

du <- qu_post - qu_pre
dv <- qv_post - qv_pre
d_mag <- sqrt(du^2 + dv^2)

# Visualization
df_ts <- data.frame(Year = 1950:2025, Flux = flux_ts)
mean_pre  <- mean(df_ts$Flux[df_ts$Year < break_year])
mean_post <- mean(df_ts$Flux[df_ts$Year >= break_year])
df_ts$Regime <- ifelse(df_ts$Year < break_year, mean_pre, mean_post)

plot_A <- ggplot(df_ts, aes(x = Year, y = Flux)) +
  geom_line(color = "grey50", linewidth = 0.8) +
  geom_point(color = "black", size = 1.5, alpha = 0.7) +
  geom_vline(xintercept = break_year, linetype = "dashed", color = "darkred", linewidth = 1) +
  geom_step(aes(y = Regime), color = "darkred", linewidth = 1.2) +
  annotate("text", x = 1975, y = max(df_ts$Flux), 
           label = paste("Pre-2003 Mean:", round(mean_pre, 3)), color = "darkred", fontface = "bold") +
  annotate("text", x = 2014, y = max(df_ts$Flux), 
           label = paste("Post-2003 Mean:", round(mean_post, 3)), color = "darkred", fontface = "bold") +
  theme_classic(base_size = 14) +
  labs(title = "(a) JJAS Moisture Flux Magnitude over Mizoram",
       x = "Year", y = expression("Moisture Flux Magnitude (kg m"^-1*" s"^-1*")"))

du_agg <- aggregate(du, fact = 2, fun = mean)
dv_agg <- aggregate(dv, fact = 2, fun = mean)
df_vectors <- as.data.frame(c(du_agg, dv_agg), xy = TRUE)
names(df_vectors) <- c("x","y","du","dv")
mizo_box <- data.frame(xmin = 92, xmax = 93, ymin = 22, ymax = 24)

plot_B <- ggplot() +
  geom_spatraster(data = d_mag) +
  scale_fill_viridis_c(option = "mako", name = "Anomaly\nMagnitude", alpha = 0.9) +
  geom_segment(data = df_vectors, 
               aes(x = x, y = y, xend = x + du*200, yend = y + dv*200),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), 
               color = "white", linewidth = 0.4) +
  geom_rect(data = mizo_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "red", linewidth = 1.2) +
  theme_minimal(base_size = 14) +
  labs(title = "(b) Post-2003 Moisture Flux Anomaly",
       subtitle = "Difference in 850 hPa VIMT (2003-2025 minus 1950-2002)",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey80", linetype = "dotted")) +
  coord_sf(expand = FALSE)

final_plot <- plot_A / plot_B + plot_layout(heights = c(1, 1.5))
ggsave(file.path(out_dir, "Figure_Mechanism.png"), final_plot, width = 10, height = 12, dpi = 600)

rm(era_p, era_s)
message("ERA5 circulation analysis complete.")
