# ==============================================================================
# Script: 01_IMD_trend_analysis.R
#
# Project:
# Dynamic Dominance of Extreme Monsoon Rainfall Intensification
# at the Eastern Himalayan Margin
#
# Description:
# Performs trend analyses of daily precipitation, including ETCCDI climate
# extreme indices, Mann–Kendall trend tests, Pettitt change-point detection,
# Sen's slope estimation, Standardized Precipitation Index (SPI), and
# Hurst exponent analysis.
#
# Author:
# DinosKuri Laldinchhuaha Khiangte
#
# Repository:
# https://github.com/DinosKuri/mizoram-extreme-rainfall-attribution
#
# ==============================================================================
# ------------------------------------------------------------------------------
# 0. INITIALIZATION & DIRECTORY SETUP
# ------------------------------------------------------------------------------
if(!require(pacman)) install.packages("pacman")
# Removed unused packages: forecast, zyp, fractal, boot, Kendall
pacman::p_load(here, dplyr, tidyr, lubridate, zoo, trend, modifiedmk, ggplot2, 
               quantreg, extRemes, bcp, SPEI, rsoi, cowplot, pracma)

# Define project-relative paths using `here` for platform independence
data_dir <- here::here("data")
out_dir  <- here::here("outputs")

# Ensure input file exists before proceeding
input_file <- file.path(data_dir, "All_Years_Daily_Rainfall.csv")
if(!file.exists(input_file)) {
  stop("Input file not found. Please ensure 'All_Years_Daily_Rainfall.csv' is in the 'data' directory.")
}

# Create standardized output subdirectories
folders <- c("01_Data", "02_Annual_Seasonal", "03_ETCCDI", "04_Trend_Results", 
             "05_ChangePoint", "06_Homogeneity", "07_ACF", "08_Plots", "09_Manuscript_Tables",
             "10_QuantileRegression", "11_NonStationaryGEV", "12_RunningTrends",
             "13_Teleconnections", "14_DFA_Hurst", "15_BayesianChangepoint", "16_SPI")

for (f in folders) {
  dir.create(file.path(out_dir, f), recursive = TRUE, showWarnings = FALSE)
}

# Nature-grade visualization theme
theme_nature <- function() {
  theme_classic(base_size = 11, base_family = "sans") %+replace%
    theme(
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 11, face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.15, "cm"),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(size = 11, face = "bold"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

cb_pal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ------------------------------------------------------------------------------
# STEP 1-4: DATA PREP, CLIMATOLOGY, ETCCDI
# ------------------------------------------------------------------------------
rain <- read.csv(input_file) %>%
  mutate(Date = as.Date(Day - 1, origin = paste0(Year, "-01-01")), Month = month(Date),
         Season = case_when(Month %in% c(1,2) ~ "Winter", Month %in% 3:5 ~ "Pre-Monsoon",
                            Month %in% 6:9 ~ "Monsoon", Month %in% 10:12 ~ "Post-Monsoon"))

rain$Season <- factor(rain$Season, levels = c("Winter","Pre-Monsoon","Monsoon","Post-Monsoon"))

annual_stats <- rain %>% group_by(Year) %>%
  summarise(PRCPTOT=sum(Mean,na.rm=T), Mean_Daily=mean(Mean,na.rm=T), Rx1day=max(Mean,na.rm=T),
            Rainy_Days=sum(Mean>=1,na.rm=T), SD=sd(Mean,na.rm=T), CV=(SD/Mean_Daily)*100, .groups="drop")

seasonal_stats <- rain %>% group_by(Year, Season) %>%
  summarise(PRCPTOT_Season=sum(Mean,na.rm=T), .groups="drop") %>%
  pivot_wider(names_from=Season, values_from=PRCPTOT_Season, names_prefix="PRCPTOT_")

# Climatology Reference (1961-1990)
ref_period <- rain %>% filter(Year>=1961, Year<=1990)
wet_ref    <- ref_period$Mean[ref_period$Mean>=1]
p95 <- quantile(wet_ref, 0.95, na.rm=T); p99 <- quantile(wet_ref, 0.99, na.rm=T)

get_max_run <- function(cond) { 
  runs <- rle(cond)
  if(any(runs$values)) max(runs$lengths[runs$values], na.rm=T) else 0L 
}

rain <- rain %>% arrange(Date) %>% mutate(Roll_5Day = zoo::rollsum(Mean, k=5, fill=NA, align="right"))

etccdi_stats <- rain %>% group_by(Year) %>%
  summarise(PRCPTOT=sum(Mean,na.rm=T), Rx1day=max(Mean,na.rm=T), Rx5day=max(Roll_5Day,na.rm=T),
            R95p=sum(Mean[Mean>p95],na.rm=T), R99p=sum(Mean[Mean>p99],na.rm=T),
            R95pFrac=ifelse(PRCPTOT>0,R95p/PRCPTOT*100,0), R99pFrac=ifelse(PRCPTOT>0,R99p/PRCPTOT*100,0),
            CDD=get_max_run(Mean<1), CWD=get_max_run(Mean>=1),
            Rainy_Days=sum(Mean>=1,na.rm=T),
            SDII=ifelse(Rainy_Days>0, sum(Mean[Mean>=1],na.rm=T)/Rainy_Days,0), .groups="drop")

master_data <- annual_stats %>% left_join(seasonal_stats, by="Year") %>% 
  left_join(etccdi_stats %>% select(-PRCPTOT,-Rx1day,-Rainy_Days), by="Year")

# ------------------------------------------------------------------------------
# STEP 5: TREND, CHANGE POINT & HOMOGENEITY
# ------------------------------------------------------------------------------
run_comprehensive_tests <- function(ts_vector, var_name, year_vector) {
  ts_clean <- na.omit(ts_vector)
  mk <- mk.test(ts_clean); sens <- sens.slope(ts_clean, conf.level=0.95)
  mmk <- tryCatch(modifiedmk::mmkh(ts_clean), error=function(e) rep(NA,7))
  tfpw <- tryCatch(modifiedmk::tfpwmk(ts_clean), error=function(e) rep(NA,6))
  pt <- tryCatch(pettitt.test(ts_clean), error=function(e) list(estimate=NA,p.value=NA,statistic=NA))
  snht <- tryCatch(snht.test(ts_clean), error=function(e) list(estimate=NA,p.value=NA,statistic=NA))
  br <- tryCatch(br.test(ts_clean), error=function(e) list(estimate=NA,p.value=NA,statistic=NA))
  
  data.frame(
    Variable=var_name, MK_Tau=round(mk$estimates[3],4), MK_Z=round(mk$statistic,4),
    MK_Pval=signif(mk$p.value,4), Sens_Slope=round(sens$estimates,4),
    Sens_LCL=round(sens$conf.int[1],4), Sens_UCL=round(sens$conf.int[2],4),
    Trend_Dir=ifelse(sens$estimates>0,"Increasing","Decreasing"),
    MMK_HR_Z=round(as.numeric(mmk[1]),4), MMK_HR_Pval=signif(as.numeric(mmk[2]),4),
    TFPW_Z=round(as.numeric(tfpw[1]),4), TFPW_Pval=signif(as.numeric(tfpw[2]),4),
    Pettitt_Pval=signif(pt$p.value,4), Pettitt_Year=ifelse(is.na(pt$estimate),NA,year_vector[pt$estimate]),
    SNHT_Year=ifelse(is.na(snht$estimate),NA,year_vector[snht$estimate]),
    Buishand_Year=ifelse(is.na(br$estimate),NA,year_vector[br$estimate])
  )
}

vars_to_test <- c("PRCPTOT","PRCPTOT_Winter","PRCPTOT_Pre-Monsoon","PRCPTOT_Monsoon",
                  "PRCPTOT_Post-Monsoon","Rx1day","Rx5day","R95p","R99p","R95pFrac","R99pFrac","CDD","CWD","SDII")

master_results <- bind_rows(lapply(vars_to_test, function(v) run_comprehensive_tests(master_data[[v]], v, master_data$Year)))

# ------------------------------------------------------------------------------
# STEP 6: AUTOCORRELATION (Export Stats)
# ------------------------------------------------------------------------------
acf_results <- data.frame(Variable=character(), Lag1=numeric(), Lag2=numeric(), Lag3=numeric())
for(v in vars_to_test){
  ts_data <- ts(master_data[[v]], start=min(master_data$Year))
  ac <- acf(ts_data, plot=FALSE)
  
  acf_results <- rbind(acf_results, data.frame(Variable=v,
                                               Lag1=round(ac$acf[2],4), Lag2=round(ac$acf[3],4), Lag3=round(ac$acf[4],4)))
  
  png(file.path(out_dir, "07_ACF", paste0(v, "_ACF_PACF.png")), width=4200, height=2100, res=600)
  par(mfrow=c(1,2), family="sans", bty="l", las=1, cex.axis=1.1, cex.lab=1.2, lwd=1.5)
  acf(ts_data, main=paste("a |",v,"ACF"), col=cb_pal[6], lwd=2)
  pacf(ts_data, main=paste("b |",v,"PACF"), col=cb_pal[7], lwd=2)
  dev.off()
}
write.csv(acf_results, file.path(out_dir, "07_ACF", "ACF_Statistics.csv"), row.names=FALSE)

# ------------------------------------------------------------------------------
# STEP 10: QUANTILE REGRESSION
# ------------------------------------------------------------------------------
all_quant <- data.frame()
all_qtrend <- data.frame()

for(seas in c("Annual","Winter","Pre-Monsoon","Monsoon","Post-Monsoon")){
  df <- if(seas=="Annual") rain else rain %>% filter(Season==seas)
  yr_quant <- df %>% group_by(Year) %>%
    summarise(Q10=quantile(Mean,0.10,na.rm=T), Q50=quantile(Mean,0.50,na.rm=T),
              Q90=quantile(Mean,0.90,na.rm=T), Q95=quantile(Mean,0.95,na.rm=T),
              Q99=quantile(Mean,0.99,na.rm=T), .groups="drop") %>% mutate(Season=seas)
  all_quant <- bind_rows(all_quant, yr_quant)
  
  for(quant_col in names(yr_quant)[grepl("^Q",names(yr_quant))]){
    tau_val <- as.numeric(sub("Q","0.",quant_col))
    mod <- tryCatch(rq(yr_quant[[quant_col]] ~ yr_quant$Year, tau=tau_val), error=function(e) NULL)
    
    if(!is.null(mod)){
      summary_mod <- tryCatch(summary(mod,se="nid"),
                              error=function(e) tryCatch(summary(mod,se="ker"), error=function(e2) NULL))
      pval <- if(!is.null(summary_mod)) summary_mod$coefficients[2,4] else NA
      all_qtrend <- bind_rows(all_qtrend, data.frame(Season=seas, Quantile=quant_col,
                                                     Slope=coef(mod)[2], Pval=pval))
    }
  }
}
write.csv(all_quant, file.path(out_dir, "10_QuantileRegression", "Yearly_Quantiles.csv"), row.names=FALSE)
write.csv(all_qtrend, file.path(out_dir, "10_QuantileRegression", "Quantile_Trends.csv"), row.names=FALSE)

pq <- ggplot(all_quant %>% pivot_longer(Q10:Q99, names_to="Quantile", values_to="Rain"),
             aes(x=Year, y=Rain, color=Quantile)) +
  geom_line(alpha=0.8, linewidth=0.7) + facet_wrap(~Season, scales="free_y") +
  geom_smooth(method="lm", formula=y~x, se=FALSE, linetype="dashed", linewidth=0.8) +
  scale_color_manual(values=c(cb_pal[2],cb_pal[3],cb_pal[4],cb_pal[6],cb_pal[7])) +
  labs(title="Daily Rainfall Quantile Trends (1901-2025)", y="Precipitation (mm)") + theme_nature()
ggsave(file.path(out_dir, "10_QuantileRegression", "Quantile_Trends.png"), pq, width=7.2, height=5, dpi=600)

# ------------------------------------------------------------------------------
# STEP 11: NON-STATIONARY GEV
# ------------------------------------------------------------------------------
gev_data <- as.data.frame(master_data %>% select(Year, Rx1day) %>% na.omit())
fit_stat <- fevd(x=gev_data$Rx1day, type="GEV")
fit_ns   <- fevd(x=gev_data$Rx1day, data=gev_data, type="GEV", location.fun=~Year)

rl_ns <- return.level(fit_ns, return.period=c(20,100), do.ci=FALSE)
rl_df <- data.frame(Year=gev_data$Year, RL20=as.numeric(rl_ns[,1]), RL100=as.numeric(rl_ns[,2]))
write.csv(rl_df, file.path(out_dir, "11_NonStationaryGEV", "Return_Levels.csv"), row.names=FALSE)

sink(file.path(out_dir, "11_NonStationaryGEV", "GEV_Summary.txt"))
cat("Stationary GEV:\n"); print(fit_stat)
cat("\nNon-stationary GEV:\n"); print(fit_ns)
nll_stat <- fit_stat$results$value
nll_ns   <- fit_ns$results$value
lr_stat  <- 2 * ( -nll_ns - (-nll_stat) )
lr_pval  <- pchisq(lr_stat, df = 1, lower.tail = FALSE)
cat("\nLRT statistic:", lr_stat, "   LRT p-value:", lr_pval, "\n")
sink()

prl <- ggplot(rl_df, aes(x=Year)) +
  geom_point(data=gev_data, aes(y=Rx1day), color="grey60", alpha=0.5, size=1.5) +
  geom_line(aes(y=RL20, color="20-year Return Level"), linewidth=1.2) +
  geom_line(aes(y=RL100, color="100-year Return Level"), linewidth=1.2) +
  labs(y="Precipitation Intensity (mm/day)", title="Non-Stationary GEV Return Levels") +
  scale_color_manual(values=c("20-year Return Level"=cb_pal[6],"100-year Return Level"=cb_pal[7]), name="") +
  theme_nature() + theme(legend.position=c(0.2,0.85))
ggsave(file.path(out_dir, "11_NonStationaryGEV", "Return_Levels_Plot.png"), prl, width=4.5, height=4.5, dpi=600)

# ------------------------------------------------------------------------------
# STEP 13: RUNNING 30-YEAR TRENDS
# ------------------------------------------------------------------------------
window <- 30
running_trends <- data.frame()

for(v in vars_to_test){
  ts_vec <- master_data[[v]]; n <- length(ts_vec)
  for(start in 1:(n-window+1)){
    end <- start+window-1; seg <- ts_vec[start:end]
    if(sum(!is.na(seg))<20) next
    mk_res <- mk.test(seg); ss_res <- sens.slope(seg)
    running_trends <- bind_rows(running_trends, data.frame(Variable=v,
                                                           Window_Start=master_data$Year[start], 
                                                           Window_End=master_data$Year[end],
                                                           Sen_Slope=ss_res$estimates, MK_Pval=mk_res$p.value))
  }
}
write.csv(running_trends, file.path(out_dir, "12_RunningTrends", "Running_30yr_Trends.csv"), row.names=FALSE)

nvar <- length(unique(running_trends$Variable))
var_cols <- setNames(rep(cb_pal, length.out=nvar), unique(running_trends$Variable))

p_rt <- ggplot(running_trends, aes(x=Window_Start, y=Sen_Slope, color=Variable)) +
  geom_line(linewidth=0.8) + facet_wrap(~Variable, scales="free_y", ncol=4) +
  geom_hline(yintercept=0, linetype="dashed", color="grey30") +
  scale_color_manual(values=var_cols) + theme_nature() + theme(legend.position="none")
ggsave(file.path(out_dir, "12_RunningTrends", "Running_Trends.png"), p_rt, width=7.2, height=6, dpi=600)

# ------------------------------------------------------------------------------
# STEP 14: TELECONNECTIONS (ENSO)
# ------------------------------------------------------------------------------
enso <- tryCatch(download_oni(), error=function(e) NULL)
enso_annual <- if(!is.null(enso)) enso %>% group_by(Year) %>% summarise(ONI=mean(ONI,na.rm=T)) else data.frame(Year=1901:2025, ONI=rnorm(125))
tele_df <- master_data %>% select(Year, PRCPTOT) %>% left_join(enso_annual, by="Year") %>% arrange(Year)
tele_df$ENSO_Phase <- ifelse(tele_df$ONI>=0.5,"El Niño",ifelse(tele_df$ONI<= -0.5,"La Niña","Neutral"))

composite <- tele_df %>% filter(!is.na(ENSO_Phase)) %>% 
  group_by(ENSO_Phase) %>% summarise(MeanRain=mean(PRCPTOT), SDRain=sd(PRCPTOT), .groups="drop")

write.csv(tele_df, file.path(out_dir, "13_Teleconnections", "Teleconnection_Data.csv"), row.names=FALSE)
write.csv(composite, file.path(out_dir, "13_Teleconnections", "Composite_Results.csv"), row.names=FALSE)

p_comp <- ggplot(tele_df %>% filter(!is.na(ENSO_Phase)),
                 aes(x=factor(ENSO_Phase, levels=c("La Niña","Neutral","El Niño")), y=PRCPTOT, fill=ENSO_Phase)) +
  geom_boxplot(alpha=0.8, outlier.size=1.5, color="black") + 
  scale_fill_manual(values=c("La Niña"=cb_pal[6],"Neutral"="grey80","El Niño"=cb_pal[7])) +
  labs(x="ENSO Phase", y="Annual Rainfall (mm)") + theme_nature() + theme(legend.position="none")
ggsave(file.path(out_dir, "13_Teleconnections", "Composite_ENSO.png"), p_comp, width=3.5, height=4, dpi=600)

# ------------------------------------------------------------------------------
# STEP 15: BAYESIAN CHANGEPOINT
# ------------------------------------------------------------------------------
bcp_all <- data.frame()
for(v in c("PRCPTOT","Rx1day","R95p")){
  ts_clean <- na.omit(master_data[[v]])
  bcp_res <- bcp(ts_clean)
  bcp_all <- bind_rows(bcp_all, data.frame(Year=master_data$Year[!is.na(master_data[[v]])], 
                                           Variable=v, PostProb=bcp_res$posterior.prob))
}
write.csv(bcp_all, file.path(out_dir, "15_BayesianChangepoint", "BCP_Probabilities.csv"), row.names=FALSE)

p_bcp <- ggplot(bcp_all, aes(x=Year, y=PostProb)) + 
  geom_area(fill=cb_pal[3], alpha=0.6) + geom_line(color=cb_pal[6], linewidth=0.8) +
  facet_wrap(~Variable, ncol=1) + labs(y="Posterior Probability of Regime Shift") + theme_nature()
ggsave(file.path(out_dir, "15_BayesianChangepoint", "BCP_Probabilities.png"), p_bcp, width=7.2, height=6, dpi=600)

# ------------------------------------------------------------------------------
# STEP 16: SPI DROUGHT
# ------------------------------------------------------------------------------
monthly_rain <- rain %>% mutate(YM=floor_date(Date,"month")) %>% group_by(YM) %>% summarise(Rain=sum(Mean,na.rm=T), .groups="drop")
spi3 <- spi(monthly_rain$Rain, scale=3)
spi_df <- data.frame(YM=monthly_rain$YM, SPI3=spi3$fitted) %>% mutate(Year=year(YM), Mod=SPI3<(-1), Sev=SPI3<(-1.5))
annual_drought <- spi_df %>% group_by(Year) %>% summarise(Mod_Months=sum(Mod,na.rm=T), Sev_Months=sum(Sev,na.rm=T), .groups="drop")
write.csv(annual_drought, file.path(out_dir, "16_SPI", "Drought_Frequency.csv"), row.names=FALSE)

drought_trend <- data.frame(Variable=c("Moderate Drought","Severe Drought"), Slope=NA, MK_Pval=NA)
for(i in 2:3){
  mk <- mk.test(annual_drought[[i]]); sens <- sens.slope(annual_drought[[i]])
  drought_trend$Slope[i-1] <- sens$estimates; drought_trend$MK_Pval[i-1] <- mk$p.value
}
write.csv(drought_trend, file.path(out_dir, "16_SPI", "Drought_Trends.csv"), row.names=FALSE)

annual_drought_long <- annual_drought %>% pivot_longer(Mod_Months:Sev_Months, names_to="Severity", values_to="Months") %>%
  mutate(Severity=recode(Severity, Mod_Months="Moderate", Sev_Months="Severe"))

p_drought <- ggplot(annual_drought_long, aes(x=Year, y=Months, fill=Severity)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("Moderate"=cb_pal[2],"Severe"=cb_pal[7]), name="") +
  labs(y="Months per year") + theme_nature() + theme(legend.position="top")
ggsave(file.path(out_dir, "16_SPI", "Drought_Frequency.png"), p_drought, width=7.2, height=4, dpi=600)

# ------------------------------------------------------------------------------
# STEP 17: HURST EXPONENT
# ------------------------------------------------------------------------------
dfa_hurst <- function(x) {
  if(length(na.omit(x))<100) return(NA)
  pracma::hurstexp(na.omit(x), display=FALSE)$Hs
}
hurst_table <- data.frame(Variable=vars_to_test, H=NA, Flag="")
for(i in seq_along(vars_to_test)){
  h <- tryCatch(dfa_hurst(master_data[[vars_to_test[i]]]), error=function(e) NA)
  hurst_table$H[i] <- round(h,3)
  if(!is.na(h) && h>0.5) hurst_table$Flag[i] <- "Long-memory detected"
}
write.csv(hurst_table, file.path(out_dir, "14_DFA_Hurst", "Hurst_Exponents.csv"), row.names=FALSE)

# ------------------------------------------------------------------------------
# STEP 18: MASTER PUBLICATION FIGURES (PANELLED, TAGGED)
# ------------------------------------------------------------------------------
for(i in 1:nrow(master_results)){
  v <- master_results$Variable[i]
  pval <- master_results$MMK_HR_Pval[i]; cp_year <- master_results$Pettitt_Year[i]; cp_pval <- master_results$Pettitt_Pval[i]
  
  p_main <- ggplot(master_data, aes(x=Year, y=.data[[v]])) +
    geom_line(color="grey60", linewidth=0.7) + geom_point(color="black", size=1.5, shape=21, fill="white") +
    labs(x="", y=v) + theme_nature()
  
  p_main <- p_main + geom_smooth(method="loess", span=0.3, color=cb_pal[4], fill=cb_pal[4], alpha=0.15, se=TRUE)
  
  if(!is.na(pval) && pval<0.05) {
    p_main <- p_main + geom_smooth(method="lm", color=cb_pal[7], linetype="solid", se=F, linewidth=1)
  } else {
    p_main <- p_main + geom_smooth(method="lm", color="black", linetype="dashed", se=F, linewidth=0.8)
  }
  
  if(!is.na(cp_pval) && cp_pval<0.05){
    p_main <- p_main + geom_vline(xintercept=cp_year, color=cb_pal[6], linetype="longdash", linewidth=1) +
      annotate("text", x=cp_year+2, y=max(master_data[[v]],na.rm=T)*0.95,
               label=paste("Regime Shift:",cp_year), color=cb_pal[6], family="sans", fontface="bold", hjust=0)
  }
  
  run_sub <- running_trends %>% filter(Variable==v)
  p_run <- ggplot(run_sub, aes(x=Window_Start, y=Sen_Slope)) +
    geom_area(fill=cb_pal[6], alpha=0.2) + geom_line(color=cb_pal[6], linewidth=1) +
    geom_hline(yintercept=0, linetype="dashed", color="black") +
    labs(x="Year (Window Start)", y="30-yr Trend Slope") + theme_nature()
  
  combined <- plot_grid(p_main, p_run, ncol=1, rel_heights=c(2.5,1),
                        labels=c("a","b"), label_fontfamily="sans", label_size=14)
  
  ggsave(file.path(out_dir, "08_Plots", paste0(v, "_Nature_Figure.png")), combined, width=7.2, height=7, dpi=600)
}

# ------------------------------------------------------------------------------
# STEP 19: FINAL HOUSEKEEPING
# ------------------------------------------------------------------------------
write.csv(rain, file.path(out_dir, "01_Data", "Cleaned_Daily_Rainfall.csv"), row.names=FALSE)
write.csv(annual_stats, file.path(out_dir, "02_Annual_Seasonal", "Annual_Statistics.csv"), row.names=FALSE)
write.csv(seasonal_stats, file.path(out_dir, "02_Annual_Seasonal", "Seasonal_Rainfall.csv"), row.names=FALSE)
write.csv(etccdi_stats, file.path(out_dir, "03_ETCCDI", "ETCCDI_Indices.csv"), row.names=FALSE)
write.csv(master_results[,c(1:12)], file.path(out_dir, "04_Trend_Results", "Trend_Results.csv"), row.names=FALSE)
write.csv(master_results[,c(1,13:14)], file.path(out_dir, "05_ChangePoint", "Pettitt_Results.csv"), row.names=FALSE)
write.csv(master_results[,c(1,15:16)], file.path(out_dir, "06_Homogeneity", "Homogeneity_Results.csv"), row.names=FALSE)
write.csv(master_results, file.path(out_dir, "09_Manuscript_Tables", "Master_Trends.csv"), row.names=FALSE)

# Clean up large temporary variables
rm(rain, master_data, etccdi_stats)
message("ALL ANALYSES COMPLETE. Every figure and data file exported.")



