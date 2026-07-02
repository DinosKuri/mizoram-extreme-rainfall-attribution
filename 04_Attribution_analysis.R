# ==============================================================================
# Script: 04_Attribution_analysis.R
#
# Project:
# Dynamic Dominance of Extreme Monsoon Rainfall Intensification
# at the Eastern Himalayan Margin
#
# Description:
# Quantifies the relative contributions of dynamic atmospheric circulation
# changes and thermodynamic moisture increases to observed extreme rainfall
# intensification using attribution diagnostics.
#
# Author:
# Laldinchhuaha Khiangte
#
# Repository:
# https://github.com/DinosKuri/mizoram-extreme-rainfall-attribution
#
# ==============================================================================

library(here)
library(terra)
library(dplyr)
library(ggplot2)
library(patchwork)
library(zyp)

# ------------------------------------------------------------------------------
# DATA LOADING
# ------------------------------------------------------------------------------
cmip_dir <- here::here("data", "CMIP6")
out_dir  <- here::here("outputs", "Figures")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

obs_file <- file.path(cmip_dir, "Observed_Rx1day_1951_2014.csv")
era_file <- file.path(cmip_dir, "ERA5_monthly_tas_JJAS.nc")

if (!file.exists(obs_file) || !file.exists(era_file)) {
  stop("Missing input files. Please ensure observed CSV and ERA5 NC are in data/CMIP6/")
}

obs <- read.csv(obs_file) %>% filter(Year >= 1951, Year <= 2014)

r_era5 <- rast(era_file)
time_era5 <- time(r_era5)
years_era5 <- as.integer(format(time_era5, "%Y"))
tas_annual <- tapp(r_era5, years_era5, mean, na.rm = TRUE)
tas_point <- terra::extract(tas_annual, cbind(92.5, 23.0))
years_tas <- as.integer(gsub("^X", "", names(tas_point)[-1]))
tas_ts <- as.numeric(tas_point[1, -1]) - 273.15

df_tas <- data.frame(Year = years_tas, Tas_JJAS = tas_ts)
scaling_df <- inner_join(obs, df_tas, by = "Year")

# ------------------------------------------------------------------------------
# PANEL A: TWO-COMPONENT DECOMPOSITION
# ------------------------------------------------------------------------------
pre_rx  <- mean(obs$Rx1day[obs$Year < 1975])
post_rx <- mean(obs$Rx1day[obs$Year >= 1975 & obs$Year <= 2014])
total_change <- post_rx - pre_rx

pre_temp  <- mean(df_tas$Tas_JJAS[df_tas$Year < 1975])
post_temp <- mean(df_tas$Tas_JJAS[df_tas$Year >= 1975 & df_tas$Year <= 2014])
delta_T <- post_temp - pre_temp

thermo <- pre_rx * 0.07 * delta_T
dynamic <- total_change - thermo

contrib_df <- data.frame(Component = c("Thermodynamic", "Dynamic"),
                         Value = c(thermo, dynamic))

p_bar <- ggplot(contrib_df, aes(x = Component, y = Value, fill = Component)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Thermodynamic" = "#0072B2", "Dynamic" = "#D55E00")) +
  labs(y = "Rx1day change (mm)", title = "a | Attribution of Rx1day increase") +
  theme_classic(base_size = 12) + theme(legend.position = "none")

ggsave(file.path(out_dir, "Figure4a_bar.png"), p_bar, width = 4, height = 5, dpi = 600)

# ------------------------------------------------------------------------------
# PANEL B: 30-YEAR TREND SCALING (Reusable Function Added)
# ------------------------------------------------------------------------------
calc_rolling_trend <- function(ts_vec, time_vec, window = 30) {
  n <- length(ts_vec)
  trends_df <- data.frame()
  for (i in 1:(n - window + 1)) {
    seg <- i:(i + window - 1)
    fit <- zyp.sen(ts_vec[seg] ~ time_vec[seg])
    trends_df <- rbind(trends_df, data.frame(Year_start = time_vec[i], Sen_slope = coef(fit)[2]))
  }
  return(trends_df)
}

trends_tas <- calc_rolling_trend(df_tas$Tas_JJAS, df_tas$Year, 30) %>% rename(Sen_slope_tas = Sen_slope)
trends_rx  <- calc_rolling_trend(obs$Rx1day, obs$Year, 30) %>% rename(Sen_slope_rx = Sen_slope)

trends_merged <- inner_join(trends_rx, trends_tas, by = "Year_start")

p_trend <- ggplot(trends_merged, aes(x = Sen_slope_tas, y = Sen_slope_rx)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", color = "#D55E00", fill = "#D55E00", alpha = 0.2) +
  labs(x = "Temperature trend (°C per 30 yr)", y = "Rx1day trend (mm per 30 yr)",
       title = "b | 30-yr trend scaling") +
  theme_classic(base_size = 12)

ggsave(file.path(out_dir, "Figure4b_trend_scaling.png"), p_trend, width = 5, height = 5, dpi = 600)

# ------------------------------------------------------------------------------
# PANEL C: INTERANNUAL TEMPERATURE SCALING
# ------------------------------------------------------------------------------
p_temp <- ggplot(scaling_df, aes(x = Tas_JJAS, y = Rx1day)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", color = "#D55E00", fill = "#D55E00", alpha = 0.2) +
  labs(x = "JJAS temperature (°C)", y = "Rx1day (mm)",
       title = "c | Interannual scaling") +
  theme_classic(base_size = 12)

ggsave(file.path(out_dir, "Figure4c_temp_scaling.png"), p_temp, width = 5, height = 5, dpi = 600)

# ------------------------------------------------------------------------------
# COMBINED FIGURE & STATISTICS
# ------------------------------------------------------------------------------
combined_attr <- p_bar / (p_trend + p_temp) + plot_annotation(tag_levels = "a")
ggsave(file.path(out_dir, "Figure4_combined.png"), combined_attr, width = 10, height = 10, dpi = 600)

# Print Summary
scaling_df$log_Rx1day <- log(scaling_df$Rx1day)
fit <- lm(log_Rx1day ~ Tas_JJAS, data = scaling_df)
scaling_rate <- (exp(coef(fit)[2]) - 1) * 100
cat(sprintf("Temperature scaling: %.1f %% per °C\n", scaling_rate))

fit_trends <- lm(Sen_slope_rx ~ Sen_slope_tas, data = trends_merged)
cat("30-yr trend scaling p-value:", summary(fit_trends)$coefficients[2,4], "\n")
cat(sprintf("Thermodynamic: %.1f mm (%.0f%%), Dynamic: %.1f mm (%.0f%%)\n",
            thermo, 100*thermo/total_change, dynamic, 100*dynamic/total_change))

message("Attribution analysis complete.")
