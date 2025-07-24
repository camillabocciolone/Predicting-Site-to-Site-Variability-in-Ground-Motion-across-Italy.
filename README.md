# Predicting Site-to-Site Variability in Ground Motion across Italy

## Overview

This project was developed for the course *Applied Statistics* at Politecnico di Milano (A.Y. 2024–2025), in collaboration with INGV researchers. The objective was to model the site-to-site variability term (δS2S) of Ground Motion Models (GMMs) using geological, topographic, and seismic information collected from over 900 stations in Italy.

We explored the explanatory power of various geophysical proxies and statistical models to build spatial predictions of δS2S.

## Project Outcome

The project results are summarized in the final poster:

📎 `Final Poster group 11.pdf`

No written report was required for this project.

## Repository Structure

This repository contains all code files developed in R for the project, divided by topic:

- `analysis on lito classes`: MANOVA and ANOVA to test the effect of lithological categories
- `exploratory analysis`: initial data exploration, summary statistics and visualizations
- `functional data analysis`: FPCA applied to δS2S over the spectrum of periods
- `geostatistics`: kriging and spatial interpolation techniques (ordinary and universal)
- `linear models`: forward selection and multiple regression on covariates
- `regional analysis`: subdivision of Italy into macro-regions and local δS2S behavior

> ⚠️ **Note**: The dataset used in this study is **private** and cannot be published due to confidentiality restrictions.

## Methodology Summary

- **Multivariate Analysis (MANOVA/ANOVA)** to investigate the influence of categorical covariates
- **Linear Models** with forward selection for covariate importance
- **Geostatistics** (ordinary and universal kriging) for spatial modeling
- **Functional Principal Component Analysis (FPCA)** for spectral decomposition

## Tools & Languages

- **R** with packages: `gstat`, `geoR`, `fda`, `MASS`, `ggplot2`, `sp`, etc.
- **RStudio** for development

## Credits

- **Camilla Bocciolone**, Mattéo Bevilacqua, Martina Maria Calandro, Mattia De Bartolomeis, Pınar Nur Özkaplan  
  *Group 11 - Applied Statistics, Politecnico di Milano*  
- **Supervisor**: Prof. Alessandra Menafoglio  
- **INGV Advisors**: Dr. Sara Sgobba, Dr. Giulio Brunelli, Dr. Giovanni Lanzano


