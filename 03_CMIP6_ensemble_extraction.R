# ==============================================================================
# Script: 03_CMIP6_ensemble_extraction.R
#
# Project:
# Dynamic Dominance of Extreme Monsoon Rainfall Intensification
# at the Eastern Himalayan Margin
#
# Description:
# Processes CMIP6 historical climate simulations, extracts precipitation
# variables, computes multi-model ensemble statistics, and prepares data
# for attribution and future analyses.
#
# Author:
# DinosKuri Laldinchhuaha Khiangte
#
# Repository:
# https://github.com/DinosKuri/mizoram-extreme-rainfall-attribution
#
# ==============================================================================

library(here)
library(terra)
library(dplyr)
library(ggplot2)
library(tidyr)

# ------------------------------------------------------------------------------
# SETUP & INPUTS
# ------------------------------------------------------------------------------
cmip_dir <- here::here("data", "CMIP6")
out_dir  <- here::here("outputs", "Figures")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

obs_file <- file.path(cmip_dir, "Observed_Rx1day_1951_2014.csv")
if (!file.exists(obs_file)) {
  stop("Missing file: Observed_Rx1day_1951_2014.csv. Ensure it is in data/CMIP6/")
}

obs <- read.csv(obs_file, stringsAsFactors = FALSE) %>%
  filter(Year >= 1951, Year <= 2014)

# ------------------------------------------------------------------------------
# DATA EXTRACTION
# ------------------------------------------------------------------------------
# 3x3 extraction box
lon_centre <- 92.5; lat_centre <- 23.0
half_step <- 0.9375
bbox_3x3 <- ext(lon_centre - half_step, lon_centre + half_step,
                lat_centre - half_step, lat_centre + half_step)

extract_3x3 <- function(nc_path) {
  r <- tryCatch(rast(nc_path), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  
  years <- as.integer(format(time(r), "%Y"))
  r_crop <- crop(r, bbox_3x3)
  if (ncell(r_crop) == 0) return(NULL)
  
  rx_ts <- global(r_crop, "mean", na.rm = TRUE)
  mod <- strsplit(basename(nc_path), "_")[[1]][3]
  data.frame(Year = years, Rx1day = rx_ts$mean, Model = mod)
}

nc_files <- list.files(cmip_dir, pattern = "rx1day.*\\.nc$", full.names = TRUE)
if(length(nc_files) == 0) {
  warning("No NetCDF files found in the CMIP6 directory. Extraction will return empty.")
}

all_rx <- bind_rows(lapply(nc_files, extract_3x3))
if(nrow(all_rx) > 0) {
  all_rx <- all_rx %>% filter(Year >= 1951, Year <= 2014)
}

# Quality filter
quality <- all_rx %>%
  group_by(Model) %>%
  summarise(n_years = n(), frac_zero = mean(Rx1day < 1, na.rm = TRUE),
            mean_rx = mean(Rx1day, na.rm = TRUE)) %>%
  filter(n_years >= 60, frac_zero < 0.10, mean_rx > 30)

good_models <- quality$Model
ensemble_clean <- all_rx %>% filter(Model %in% good_models, Rx1day > 0)
write.csv(ensemble_clean, file.path(cmip_dir, "CMIP6_historical_rx1day_ensemble_expanded.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# DETECTION FIGURES
# ------------------------------------------------------------------------------
hist_ens <- ensemble_clean %>%
  group_by(Year) %>%
  summarise(mean = mean(Rx1day), q05 = quantile(Rx1day, 0.05),
            q95 = quantile(Rx1day, 0.95), .groups = "drop")

pre75_env <- ensemble_clean %>%
  filter(Year < 1975) %>%
  summarise(q05 = quantile(Rx1day, 0.05), q95 = quantile(Rx1day, 0.95))

df_plot <- left_join(obs, hist_ens, by = "Year")

p_det <- ggplot(df_plot, aes(Year)) +
  geom_ribbon(aes(ymin = pre75_env$q05, ymax = pre75_env$q95), fill = "blue", alpha = 0.1) +
  geom_hline(yintercept = c(pre75_env$q05, pre75_env$q95), linetype = "dotted", color = "blue") +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "red", alpha = 0.15) +
  geom_line(aes(y = mean), color = "red", linewidth = 1) +
  geom_line(aes(y = Rx1day), color = "black", linewidth = 1.2) +
  geom_vline(xintercept = 1975, linetype = "dashed", color = "grey30") +
  labs(y = "Rx1day (mm)", title = "Observed vs. CMIP6 (expanded ensemble)") +
  theme_classic(base_size = 14)

ggsave(file.path(out_dir, "Detection_Figure_expanded.png"), p_det, width = 8, height = 5, dpi = 600)

# PDF shift plot
obs_shift <- obs %>% mutate(Period = ifelse(Year < 1975, "1951-1974", "1975-2014"))
p_pdf <- ggplot(obs_shift, aes(x = Rx1day, fill = Period)) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  scale_fill_manual(values = c("1951-1974" = "blue", "1975-2014" = "red")) +
  labs(title = "Distribution of observed Rx1day", x = "Rx1day (mm)") +
  theme_classic()

ggsave(file.path(out_dir, "Observed_Rx1day_PDF_shift.png"), p_pdf, width = 7, height = 5, dpi = 600)
message("CMIP6 extraction and detection figure complete.")
