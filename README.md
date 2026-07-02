# Mizoram Extreme Rainfall Attribution

This repository contains the R scripts used in the study:

*Dynamic Dominance of Extreme Monsoon Rainfall Intensification at the Eastern Himalayan Margin*

## Repository Structure

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── scripts/
│   └── README.md
├── 01_IMD_trend_analysis.R
├── 02_ERA5_circulation.R
├── 03_CMIP6_ensemble_extraction.R
├── 04_Attribution_analysis.R
└── 05_Context_map.R
```

## Script Description

| Script | Purpose |
|---------|---------|
| `01_IMD_trend_analysis.R` | Trend analysis, ETCCDI indices, Mann–Kendall test, Pettitt test, Sen's slope, SPI, and Hurst exponent estimation |
| `02_ERA5_circulation.R` | ERA5 atmospheric circulation and moisture-flux analysis |
| `03_CMIP6_ensemble_extraction.R` | CMIP6 model processing and ensemble analysis |
| `04_Attribution_analysis.R` | Dynamic and thermodynamic attribution analyses |
| `05_Context_map.R` | Generation of the study area and contextual maps |

## Software Requirements

- R ≥ 4.4

The required R packages are automatically loaded within each script.

## Main R Packages

The analyses rely primarily on the following R packages:

- terra
- sf
- raster
- ggplot2
- dplyr
- modifiedmk
- trend
- quantreg
- extRemes
- SPEI
- lubridate

## Data Sources

This repository does not redistribute proprietary datasets.

The analyses require the following datasets:

- India Meteorological Department (IMD) daily gridded rainfall
- ERA5 reanalysis (Copernicus Climate Data Store)
- CMIP6 historical climate simulations (Earth System Grid Federation)

Please obtain these datasets directly from their official providers before reproducing the analyses.

## Citation

If you use these scripts, please cite both the associated publication and this GitHub repository.

## License

This project is distributed under the MIT License.

## Contact

For questions regarding the scripts, methodology, or reproducibility of the analyses, please open a GitHub Issue or contact the corresponding author.
