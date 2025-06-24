# ParametricSurvFit <img src="man/figures/ParametricSurvFit_sticker.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/gosukehommaEX/ParametricSurvFit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gosukehommaEX/ParametricSurvFit/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/ParametricSurvFit)](https://CRAN.R-project.org/package=ParametricSurvFit)
<!-- badges: end -->

## Overview

`ParametricSurvFit` is a comprehensive R package designed for fitting multiple parametric survival distributions to clinical trial data. It provides streamlined functions for analyzing ADaM datasets commonly used in pharmaceutical research, with support for stratified analysis, covariance matrix extraction, formatted output tables, and advanced visualization capabilities including both Kaplan-Meier curves and Q-Q plots.

## Features

- **Multiple Distribution Support**: Fit exponential, Weibull, log-normal, log-logistic, Gompertz, generalized gamma, and gamma distributions
- **ADaM Dataset Integration**: Direct support for ADSL and ADTTE datasets
- **Stratified Analysis**: Flexible stratification by any categorical variable
- **Covariance Matrix Extraction**: Extract variance-covariance matrices from fitted models
- **Formatted Output**: Professional tables with kableExtra formatting
- **Parameter Extraction**: Comprehensive extraction of distribution parameters with confidence intervals
- **Population Filtering**: Support for ITT, Safety, and Randomized populations
- **🆕 Kaplan-Meier Visualization**: Create KM curves with parametric distribution overlays
- **🆕 Q-Q Plot Analysis**: S(t)-based Q-Q plots with native flexsurv confidence intervals
- **🆕 Professional Plotting**: Risk tables, customizable colors, and parameter annotations

## Installation

You can install the development version of ParametricSurvFit from GitHub:

```r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install ParametricSurvFit
devtools::install_github("gosukehommaEX/ParametricSurvFit")
```

## Dependencies

The package requires the following R packages:
- `dplyr` (>= 1.0.0)
- `purrr` (>= 0.3.0)
- `tidyr` (>= 1.1.0)
- `flexsurv` (>= 2.0)
- `survival` (>= 3.2.0)
- `haven` (>= 2.4.0)
- `kableExtra` (>= 1.3.0)
- `rlang` (>= 0.4.0)
- `stats` (>= 3.5.0)
- `🆕 ggplot2` (>= 3.3.0)
- `🆕 survminer` (>= 0.4.0)
- `🆕 grid` (>= 3.5.0)
- `🆕 scales` (>= 1.1.0)

## Quick Start

```r
library(ParametricSurvFit)

# 1. Create survival dataset from ADaM files
surv_data <- DataParametricSurv(
  adsl_path = "path/to/adsl.sas7bdat",
  adtte_path = "path/to/adtte.sas7bdat",
  population = "ITTFL",
  variable = "OS",
  stratify_by = "SEX"
)

# 2. Fit multiple parametric distributions
results <- FitSurvMods(
  dataset = surv_data,
  distributions = c("exp", "weibull", "lnorm", "llogis"),
  table_caption = "Parametric Models for Overall Survival"
)

# 3. Extract covariance matrix for a specific distribution
cov_matrix <- ExtractCovarianceMatrix(
  dataset = surv_data,
  distribution = "weibull",
  stratify_reference_level = "M",  # Use Male as reference
  create_kable = TRUE
)

# 4. 🆕 Create Kaplan-Meier plots with parametric overlays
km_plots <- ParametricSurvKM(
  dataset = surv_data,
  distribution = "weibull",
  figure_caption = "Overall Survival Analysis",
  stratify_name = "SEX",
  control_arm = "Placebo",
  time_scale = "Months"
)

# 5. 🆕 Create S(t)-based Q-Q plots with native flexsurv confidence intervals
qq_plots <- PlotSurvivalQQ(
  dataset = surv_data,
  distribution = "weibull",
  time_scale = "Months",
  confidence_level = 0.95
)

# 6. View results
print(results)
print(cov_matrix)
print(km_plots)  # Display all stratified plots
print(qq_plots)  # Display all Q-Q plots
```

## Main Functions

### `DataParametricSurv()`
Creates a standardized dataset for parametric survival analysis from ADaM datasets.

### `FitSurvMods()`
Fits multiple parametric survival distributions and creates formatted output tables.

### `ExtractParams()`
Extracts parameters, standard errors, and confidence intervals from fitted models.

### `ExtractCovarianceMatrix()`
Extracts variance-covariance matrices from fitted parametric survival models.

### `ParametricSurvKM()` 🆕
Creates Kaplan-Meier survival curves with overlaid parametric distribution fits.

**Key Parameters:**
- `dataset`: Output from `DataParametricSurv()`
- `distribution`: Single distribution name to overlay (e.g., "weibull", "exp")
- `control_arm`: Specify which ARM should be treated as control (red color)
- `time_scale`: Label for time axis ("Months", "Years", "Weeks", "Days")
- `time_max`: Maximum time for x-axis
- `return_plots`: Return plot objects vs. display directly

**Features:**
- Separate panels for each stratification level
- Color-coded treatment arms (control=red, treatment=blue)
- Risk tables with proper alignment
- Parameter estimates displayed on plots
- Customizable time scales and axes
- Professional survminer-based styling

### `PlotSurvivalQQ()` 🆕
Creates S(t)-based Q-Q plots with log(-log(S)) transformation using native flexsurv confidence intervals.

**Key Parameters:**
- `dataset`: Output from `DataParametricSurv()`
- `distribution`: Distribution to fit and evaluate
- `confidence_level`: Confidence level for intervals (default: 0.95)
- `stratify_name`: Original stratification variable name
- `time_scale`: Time axis label

**Features:**
- Native flexsurv confidence intervals (no complex transformations)
- Separate plots for each ARM*STRATIFY combination
- Dynamic text box sizing for parameter display
- Log(-log(S)) transformation for linearity assessment
- Events plotted as circles, censored observations as vertical lines
- Data-driven axis scaling

**Advantages over SAS LIFEREG approach:**
- **Statistical consistency**: Uses S(t) directly without F(t) transformations
- **Implementation simplicity**: Leverages flexsurv's native capabilities
- **Numerical stability**: Avoids parameter conversion errors
- **Regulatory compliance**: Statistically valid and well-documented methodology

## Supported Distributions

| Distribution | Parameters | Description |
|--------------|------------|-------------|
| `exp` | Rate, Intercept | Exponential distribution |
| `weibull` | Shape, Scale, Intercept | Weibull distribution |
| `lnorm` | Intercept (meanlog), Scale (sdlog) | Log-normal distribution |
| `llogis` | Intercept, Scale | Log-logistic distribution |
| `gompertz` | Intercept, Rate, Shape | Gompertz distribution |
| `gengamma` | Intercept (mu), Scale (sigma), Shape (Q) | Generalized gamma distribution |
| `gamma` | Intercept, Shape, Rate | Gamma distribution |

## Example Workflows

### Basic Analysis with Both Visualization Types
```r
# Simple overall survival analysis without stratification
surv_data <- DataParametricSurv(
  adsl_path = "adsl.sas7bdat",
  adtte_path = "adtte.sas7bdat",
  population = "ITTFL",
  variable = "OS",
  stratify_by = NULL  # No stratification
)

# Fit models and create both types of plots
results <- FitSurvMods(surv_data, distributions = c("exp", "weibull"))
km_plots <- ParametricSurvKM(surv_data, distribution = "weibull")
qq_plots <- PlotSurvivalQQ(surv_data, distribution = "weibull")
```

### Stratified Analysis with Custom Visualization
```r
# Analysis stratified by sex
surv_data_sex <- DataParametricSurv(
  adsl_path = "adsl.sas7bdat",
  adtte_path = "adtte.sas7bdat",
  population = "ITTFL",
  variable = "OS",
  stratify_by = "SEX"
)

# Create customized plots
km_plots_custom <- ParametricSurvKM(
  dataset = surv_data_sex,
  distribution = "weibull",
  figure_caption = "Overall Survival by Sex",
  stratify_name = "SEX",
  control_arm = "Placebo",
  time_scale = "Years",
  time_max = 5
)

qq_plots_custom <- PlotSurvivalQQ(
  dataset = surv_data_sex,
  distribution = "weibull",
  stratify_name = "SEX",
  confidence_level = 0.90,
  time_scale = "Years"
)

# Access individual plots for each sex
male_km_plot <- km_plots_custom[["SEX_M"]]
female_qq_plot <- qq_plots_custom[["SEX_F"]]

# Save plots
ggsave("OS_Male_km_weibull.png", male_km_plot$plot, width = 12, height = 8)
ggsave("OS_Female_qq_weibull.png", female_qq_plot, width = 10, height = 8)
```

### Complete Analysis Pipeline
```r
# 1. Prepare data with region stratification
surv_data_region <- DataParametricSurv(
  adsl_path = "adsl.sas7bdat",
  adtte_path = "adtte.sas7bdat",
  population = "SAFFL",
  variable = "PFS",
  stratify_by = "REGION"
)

# 2. Fit all distributions
results_all <- FitSurvMods(
  dataset = surv_data_region,
  table_caption = "PFS Analysis by Region"
)

# 3. Extract covariance for best-fitting distribution
cov_weibull <- ExtractCovarianceMatrix(
  dataset = surv_data_region,
  distribution = "weibull",
  stratify_reference_level = "US",
  create_kable = TRUE
)

# 4. Create visualizations
km_plots_weibull <- ParametricSurvKM(
  dataset = surv_data_region,
  distribution = "weibull",
  figure_caption = "Progression-Free Survival by Region",
  stratify_name = "REGION",
  control_arm = "Control",
  time_scale = "Months"
)

qq_plots_weibull <- PlotSurvivalQQ(
  dataset = surv_data_region,
  distribution = "weibull",
  stratify_name = "REGION",
  time_scale = "Months"
)

# 5. Display results
print(results_all)
print(cov_weibull)
print(km_plots_weibull)
print(qq_plots_weibull)
```

## Covariance Matrix Features

The `ExtractCovarianceMatrix()` function provides:
- **Distribution-specific handling**: Adapts to each distribution's parameter structure
- **Shape parameter control**: Optional inclusion of shape parameters
- **Reference level specification**: Control stratification encoding (e.g., "M" vs "F" for SEX)
- **Formatted output**: Professional kable tables with merged cells
- **Raw matrix option**: For advanced statistical analysis

### Distributions with Shape Parameters
- **Gamma**: `include_shape = TRUE` includes the shape parameter
- **Generalized Gamma**: `include_shape = TRUE` includes the Q parameter
- **Gompertz**: `include_shape = TRUE` includes the shape parameter

### Common Stratification Reference Levels
- **SEX**: Use `"M"` (Male) or `"F"` (Female)
- **REGION**: Use specific region codes like `"US"`, `"EU"`, etc.
- **AGE_GROUP**: Use `"<65"`, `">=65"`, etc.

## Visualization Features

### KM Plots (`ParametricSurvKM`)
- **Professional KM curves** using survminer package
- **Parametric overlays** as dashed lines matching treatment colors
- **Risk tables** properly aligned with x-axis
- **Parameter annotations** showing estimates on each plot
- **Color coding**: Control arm (red), Treatment arm (blue)
- **Flexible time scales**: Months, Years, Weeks, Days
- **Stratified panels**: Separate plots for each stratification level
- **Exportable plots**: High-quality ggplot2 objects ready for publication

### Q-Q Plots (`PlotSurvivalQQ`) 🆕
- **S(t)-based analysis**: Uses survival function directly (not cumulative failure)
- **Native flexsurv confidence intervals**: No complex parameter transformations
- **Log(-log(S)) transformation**: Optimal for assessing distributional fit
- **Dynamic text boxes**: Parameter information with auto-sizing
- **ARM*STRATIFY combinations**: Separate plots for each subgroup
- **Professional styling**: Publication-ready appearance
- **Event/censoring indicators**: Clear visualization of data structure

## Statistical Approach

### Traditional SAS LIFEREG vs. flexsurv Native Approach

**SAS LIFEREG Approach:**
- Models F(t) = 1 - S(t) directly
- Requires complex parameter transformations
- Delta method applied to F(t)

**Our flexsurv Native Approach:**
- Models S(t) directly (more natural for survival analysis)
- Uses flexsurv's optimized confidence intervals
- Statistically equivalent but computationally more stable
- Better suited for regulatory submissions requiring transparency

Both approaches are statistically valid, but the flexsurv native approach offers:
- **Implementation simplicity**
- **Numerical stability**
- **Statistical consistency**
- **Regulatory transparency**

## Output Format

The package provides professionally formatted tables and plots with:
- Grouped headers by distribution
- Parameter estimates with standard errors
- 95% confidence intervals
- Stratified results when applicable
- Interactive hover effects
- Exportable formats (HTML, PDF, PNG, etc.)
- Publication-ready visualizations

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this package in your research, please cite:

```
Homma, G. (2025). ParametricSurvFit: Parametric Survival Model Fitting and Analysis. 
R package version 0.1.0. https://github.com/gosukehommaEX/ParametricSurvFit
```

## Acknowledgments

- Built using the excellent [flexsurv](https://github.com/chjackson/flexsurv) package
- Table formatting powered by [kableExtra](https://github.com/haozhu233/kableExtra)
- Data manipulation using the [tidyverse](https://www.tidyverse.org/) ecosystem
- Survival plotting using [survminer](https://github.com/kassambara/survminer) package
- Visualization with [ggplot2](https://ggplot2.tidyverse.org/)
- Statistical graphics with [grid](https://stat.ethz.ch/R-manual/R-devel/library/grid/html/00Index.html) package
