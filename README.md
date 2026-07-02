[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21132748.svg)](https://doi.org/10.5281/zenodo.21132748)

# Mizoram Extreme Rainfall Attribution

This repository contains the **reproducible R analysis scripts** accompanying the manuscript:

**_Dynamic Dominance of Extreme Monsoon Rainfall Intensification at the Eastern Himalayan Margin_**

The repository provides the software required to reproduce the principal analyses presented in the manuscript, including precipitation trend analysis, climate extreme indices, atmospheric circulation diagnostics, climate model processing, attribution analyses, and study-area mapping. The analyses are based on observations from the India Meteorological Department (IMD), ERA5 reanalysis, and CMIP6 climate model simulations.

---

## Repository Status

**Current release:** v1.0.0

**Software DOI:** https://doi.org/10.5281/zenodo.21132748

---

## Repository Structure

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── scripts/
│   └── README.md
├── data/
│   └── README.md
├── docs/
│   └── README.md
├── 01_IMD_trend_analysis.R
├── 02_ERA5_circulation.R
├── 03_CMIP6_ensemble_extraction.R
├── 04_Attribution_analysis.R
└── 05_Context_map.R
```

---

## Script Description

| Script | Description |
|---------|-------------|
| `01_IMD_trend_analysis.R` | Computes precipitation climatology, ETCCDI climate extreme indices, Mann–Kendall trend tests, Pettitt change-point detection, Sen's slope estimation, Standardized Precipitation Index (SPI), and Hurst exponent analyses using IMD daily rainfall observations. |
| `02_ERA5_circulation.R` | Processes ERA5 reanalysis data to investigate atmospheric circulation, geopotential height, wind fields, vertically integrated moisture transport, moisture-flux convergence, and circulation anomalies associated with extreme rainfall events. |
| `03_CMIP6_ensemble_extraction.R` | Processes CMIP6 historical climate simulations, extracts precipitation variables, performs quality control, and generates multi-model ensemble datasets for climate analysis. |
| `04_Attribution_analysis.R` | Quantifies the relative contributions of dynamic atmospheric circulation changes and thermodynamic moisture increases to observed extreme rainfall intensification. |
| `05_Context_map.R` | Produces publication-quality maps of the study region and associated geographic context used throughout the manuscript. |

---

## Software Requirements

- **R ≥ 4.4**

The required R packages are loaded within the individual scripts.

---

## Main R Packages

The analyses primarily rely on the following R packages:

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

---

## Reproducibility

This repository has been developed to facilitate computational reproducibility of the analyses presented in the accompanying manuscript.

To reproduce the analyses, users should:

1. Obtain the required datasets from their original data providers.
2. Update local file paths where necessary.
3. Install the required R packages.
4. Execute the scripts in the order appropriate for the intended analysis.

---

## Data Sources

This repository **does not redistribute proprietary datasets**.

The analyses are based on the following datasets:

- India Meteorological Department (IMD) daily gridded rainfall
- ERA5 reanalysis (Copernicus Climate Data Store)
- CMIP6 historical climate simulations (Earth System Grid Federation)

Please obtain these datasets directly from their respective providers before reproducing the analyses.

---

## Citation

If you use this software, please cite both:

1. The associated journal publication.
2. The archived software record available through Zenodo.

**Software DOI**

https://doi.org/10.5281/zenodo.21132748

---

## License

This project is distributed under the **MIT License**.

---

## Contact

For questions regarding the scripts, computational workflow, or reproducibility of the analyses, please open a GitHub Issue or contact the corresponding author.

---

## Acknowledgements

The analyses presented in this repository make use of datasets provided by:

- India Meteorological Department (IMD)
- Copernicus Climate Change Service (ERA5)
- Earth System Grid Federation (CMIP6)

The authors gratefully acknowledge these data providers for making the datasets available to the scientific community.
